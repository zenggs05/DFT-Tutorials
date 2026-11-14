#!/usr/bin/env python3
"""Dipole moments under an external field and finite-difference polarizability.

Each positional argument is "F:file.out", F = field value in GV/m that was
given as the a-component of scf.Electric.Field, file = OpenMX .out file.
The script reads the "Dipole moment (Debye)" block (Total row), fits
mu_alpha(F) with a straight line (or a central difference when only +F/-F
are given) and prints alpha_{alpha,a} in atomic units and in Ang^3.

Sign check: the polarizability of a molecule must be positive. OpenMX
defines the sign of scf.Electric.Field "with respect to the field applied to
electrons"; if alpha_xx comes out negative, a positive keyword value
corresponds to a physical field along -x, and the script tells you so.

Example:
    python3 dipole_alpha.py 0:nben_f0.out -2:nben_m2.out 2:nben_p2.out
"""
import argparse
import sys

DEBYE_PER_E_BOHR = 2.541746      # 1 e*bohr in Debye
AU_FIELD_PER_GV_M = 0.1944688e-2 # 1 GV/m in Hartree/(e*bohr), as in OpenMX Input_std.c
ANG3_PER_AU = 0.1481847          # 1 bohr^3 in Ang^3


def read_dipole(path):
    with open(path, errors="replace") as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if "Dipole moment (Debye)" in line:
            for l2 in lines[i + 1:i + 15]:
                f = l2.split()
                if f and f[0] == "Total" and len(f) >= 4:
                    return [float(x) for x in f[1:4]]
    sys.exit(f"{path}: 'Dipole moment (Debye)' block with a 'Total' row not found")


def slope(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def parse_args_with_negative(ap):
    """argparse treats "-0.05:file" as an unknown option; collect such items too."""
    args, extra = ap.parse_known_args()
    bad = [e for e in extra if ":" not in e]
    if bad:
        ap.error("unrecognized arguments: " + " ".join(bad))
    args.points = list(args.points) + extra
    if not args.points:
        ap.error("no value:file arguments given")
    return args


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("points", nargs="*", help="F:file.out pairs (F in GV/m)")
    ap.add_argument("--fmax", type=float, default=3.0,
                    help="use only |F| <= fmax GV/m in the linear fit (default 3)")
    args = parse_args_with_negative(ap)

    pts = []
    for item in args.points:
        if ":" not in item:
            sys.exit(f"bad argument '{item}', expected F:file")
        f, path = item.split(":", 1)
        pts.append((float(f), read_dipole(path), path))
    pts.sort(key=lambda t: t[0])

    print(f"{'F (GV/m)':>9} {'Dx':>12} {'Dy':>12} {'Dz':>12}   (Debye, Total)  file")
    for f, d, path in pts:
        print(f"{f:9.3f} " + " ".join(f"{v:12.6f}" for v in d) + f"   {path}")

    fit = [p for p in pts if abs(p[0]) <= args.fmax]
    if len({p[0] for p in fit}) < 2:
        sys.exit(f"need at least two different fields with |F| <= {args.fmax} GV/m")
    xs = [p[0] * AU_FIELD_PER_GV_M for p in fit]
    print()
    print(f"Linear fit over {len(fit)} points with |F| <= {args.fmax} GV/m:")
    for k, c in enumerate("xyz"):
        ys = [p[1][k] / DEBYE_PER_E_BOHR for p in fit]
        a = slope(xs, ys)
        print(f"  alpha_{c}a = {a:10.3f} bohr^3 = {a * ANG3_PER_AU:9.3f} Ang^3")
        if c == "x" and a < 0:
            print("  -> alpha_xa < 0: a positive scf.Electric.Field value acts as a physical")
            print("     field along -x here; report |alpha| and keep this sign convention in mind.")

    big = [p for p in pts if abs(p[0]) > args.fmax]
    if big and fit:
        print()
        print("Points outside the fit window (non-linear response / field too strong?):")
        for f, d, _ in big:
            print(f"  F = {f:7.2f} GV/m  Dx = {d[0]:10.5f} Debye")


if __name__ == "__main__":
    main()
