#!/usr/bin/env python3
"""Born effective charge from polB outputs (OpenMX, chapter 15).

Each positional argument is "u:file", where u is the displacement (Ang) of
one atom along Cartesian direction --dir and file is the text output of
    polB System.Name.scfout < polB.in > file

Z*_{alpha,beta} = dD_alpha / du_beta   (D = cell dipole, in e*Ang)
                = (V/e) dP_alpha / du_beta

The script fits D_alpha(u) linearly for alpha = x, y, z (total, core and
electron parts separately), cross-checks with the polarization block and
warns when the points are not on a straight line (possible branch jump).

Example:
    python3 born_charge.py --atom Na --dir x \
        -0.05:polB_NaCl_dx_m0.05.out 0.00:polB_NaCl_dx_p0.00.out \
        0.05:polB_NaCl_dx_p0.05.out
"""
import argparse
import re
import sys

DEBYE_PER_EANG = 4.80320471      # 1 e*Ang in Debye
BOHR = 0.529177210903            # Ang
UC_CM2_PER_E_BOHR2 = 5721.4766   # 1 e/bohr^2 in muC/cm^2

ROW = re.compile(r"^\s*([DP])([xyz])\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")


def parse_polb(path):
    """Return dict with 'D' and 'P' -> {comp: (bg, core, el, tot)}, and volume."""
    data = {"D": {}, "P": {}, "vol": None}
    with open(path, errors="replace") as fh:
        for line in fh:
            m = re.search(r"Cell_Volume\s*=\s*([0-9.Ee+-]+)", line)
            if m:
                data["vol"] = float(m.group(1))       # Bohr^3
            m = ROW.match(line)
            if m:
                kind, comp = m.group(1), m.group(2)
                try:
                    vals = tuple(float(m.group(i)) for i in range(3, 7))
                except ValueError:
                    continue
                data[kind].setdefault(comp, vals)      # keep first occurrence
    for kind in ("D", "P"):
        if len(data[kind]) != 3:
            sys.exit(f"{path}: could not find the {kind}x/{kind}y/{kind}z rows "
                     "(did polB finish? check the file)")
    return data


def slope(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        sys.exit("all displacements are identical")
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
    ap.add_argument("points", nargs="*", help="u:file pairs (u in Ang)")
    ap.add_argument("--dir", default="x", choices="xyz",
                    help="Cartesian direction of the displacement (beta)")
    ap.add_argument("--atom", default="?", help="label used in the printout")
    ap.add_argument("--tol", type=float, default=0.05,
                    help="relative tolerance for the linearity check")
    args = parse_args_with_negative(ap)

    pts = []
    for item in args.points:
        if ":" not in item:
            sys.exit(f"bad argument '{item}', expected u:file")
        u, path = item.split(":", 1)
        pts.append((float(u), parse_polb(path), path))
    if len(pts) < 2:
        sys.exit("need at least two displacements")
    pts.sort(key=lambda t: t[0])
    us = [p[0] for p in pts]

    print(f"Born effective charge of {args.atom}, displacement along {args.dir}")
    print(f"{'u (Ang)':>9} {'Dx_tot':>14} {'Dy_tot':>14} {'Dz_tot':>14}   (Debye)  file")
    for u, d, path in pts:
        print(f"{u:9.4f} " + " ".join(f"{d['D'][c][3]:14.6f}" for c in "xyz") + f"   {path}")

    print()
    print(f"{'alpha':>5} {'Z*_tot':>10} {'core':>10} {'electron':>10} {'Z*(from P)':>11}")
    warn = []
    for c in "xyz":
        res = []
        for col in (3, 1, 2):          # total, core, electron
            ys = [p[1]["D"][c][col] / DEBYE_PER_EANG for p in pts]
            res.append(slope(us, ys))
        # cross-check with polarization block: Z = V dP / (e du)
        vol = pts[0][1]["vol"]
        zp = float("nan")
        if vol:
            ys = [p[1]["P"][c][3] / UC_CM2_PER_E_BOHR2 * vol for p in pts]   # e*bohr
            zp = slope([u / BOHR for u in us], ys)
        print(f"{c + args.dir:>5} {res[0]:10.4f} {res[1]:10.4f} {res[2]:10.4f} {zp:11.4f}")
        # linearity check with one-sided slopes
        if len(pts) >= 3:
            ys = [p[1]["D"][c][3] / DEBYE_PER_EANG for p in pts]
            s1 = (ys[1] - ys[0]) / (us[1] - us[0])
            s2 = (ys[-1] - ys[-2]) / (us[-1] - us[-2])
            scale = max(abs(s1), abs(s2), 1e-3)
            if abs(s1 - s2) / scale > args.tol and scale > 0.05:
                warn.append(f"Z*_{c}{args.dir}: one-sided slopes {s1:.3f} vs {s2:.3f} differ; "
                            "check for a polarization-quantum jump or too large displacement")
    print()
    print("Note: Z*_tot = core + electron; the core part equals the valence charge of the")
    print("      pseudopotential (semicore electrons included), the electron part screens it.")
    for w in warn:
        print("WARNING:", w)


if __name__ == "__main__":
    main()
