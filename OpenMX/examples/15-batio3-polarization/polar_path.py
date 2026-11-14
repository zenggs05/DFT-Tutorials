#!/usr/bin/env python3
"""Spontaneous polarization from polB outputs along a lambda path.

Each positional argument is "lambda:file" (file = polB output). The script
reads the total polarization component --comp (default z), the cell volume
and the lattice vectors printed by polB, removes branch jumps by adding
integer multiples of the polarization quantum  Q = f * e * R_comp / Omega
so that P(lambda) is continuous, and prints Delta P = P(1) - P(0).

f = 1 is the safe default (any odd ionic valence or a spin-polarized case);
for spin-degenerate electrons and only even valence charges, the values can
differ only by multiples of 2*e*R/Omega (use --factor 2 to see that lattice).

Example:
    python3 polar_path.py 0.00:polB_bto_lam0.00.out 0.25:polB_bto_lam0.25.out \
        0.50:polB_bto_lam0.50.out 0.75:polB_bto_lam0.75.out 1.00:polB_bto_lam1.00.out
"""
import argparse
import re
import sys

UC_CM2_PER_E_BOHR2 = 5721.4766   # 1 e/bohr^2 in muC/cm^2
ROW = re.compile(r"^\s*P([xyz])\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")
TV = re.compile(r"^\s*tv([123])\s*=\s*(\S+)\s+(\S+)\s+(\S+)")


def parse(path):
    P, tv, vol = {}, {}, None
    with open(path, errors="replace") as fh:
        for line in fh:
            m = re.search(r"Cell_Volume\s*=\s*([0-9.Ee+-]+)", line)
            if m:
                vol = float(m.group(1))                     # Bohr^3
            m = TV.match(line)
            if m:
                tv[int(m.group(1))] = [float(m.group(i)) for i in (2, 3, 4)]   # Bohr
            m = ROW.match(line)
            if m and m.group(1) not in P:
                try:
                    P[m.group(1)] = [float(m.group(i)) for i in (2, 3, 4, 5)]
                except ValueError:
                    pass
    if len(P) != 3 or vol is None or len(tv) != 3:
        sys.exit(f"{path}: incomplete polB output (need Px/Py/Pz, Cell_Volume, tv1-3)")
    return P, tv, vol


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("points", nargs="+", help="lambda:file pairs")
    ap.add_argument("--comp", default="z", choices="xyz")
    ap.add_argument("--R", type=int, default=3, choices=(1, 2, 3),
                    help="lattice vector that defines the quantum (default 3 = c)")
    ap.add_argument("--factor", type=float, default=1.0,
                    help="quantum multiplier f (see the description above)")
    args = ap.parse_args()

    k = "xyz".index(args.comp)
    data = []
    for item in args.points:
        if ":" not in item:
            sys.exit(f"bad argument '{item}', expected lambda:file")
        lam, path = item.split(":", 1)
        P, tv, vol = parse(path)
        data.append((float(lam), P[args.comp], tv, vol, path))
    data.sort(key=lambda t: t[0])

    tv, vol = data[0][2], data[0][3]
    q = args.factor * abs(tv[args.R][k]) / vol * UC_CM2_PER_E_BOHR2
    if q == 0:
        sys.exit(f"lattice vector {args.R} has no {args.comp} component; choose --R")
    print(f"Polarization quantum along {args.comp} (f={args.factor:g}, R=tv{args.R}): "
          f"{q:.3f} muC/cm^2")
    print(f"{'lambda':>7} {'core':>12} {'electron':>12} {'total(raw)':>12} "
          f"{'total(unwrapped)':>17}  file")
    prev = None
    unwrapped = []
    for lam, p, _, _, path in data:
        tot = p[3]
        if prev is not None:
            tot -= q * round((tot - prev) / q)
        unwrapped.append(tot)
        prev = tot
        print(f"{lam:7.3f} {p[1]:12.4f} {p[2]:12.4f} {p[3]:12.4f} {tot:17.4f}  {path}")

    dP = unwrapped[-1] - unwrapped[0]
    print()
    print(f"Delta P = P(lambda={data[-1][0]:g}) - P(lambda={data[0][0]:g}) = {dP:.3f} muC/cm^2")
    print(f"(defined modulo {q:.3f} muC/cm^2; the branch that is continuous along "
          "the path is the physical one)")
    if len(data) >= 3:
        steps = [b - a for a, b in zip(unwrapped, unwrapped[1:])]
        if max(abs(s) for s in steps) > 0.3 * q:
            print("WARNING: a step along the path exceeds 0.3 quantum; add more lambda points.")


if __name__ == "__main__":
    main()
