#!/usr/bin/env python3
"""Check the smearing entropy term T*S of VASP runs (ISMEAR/SIGMA test).

VASP writes the electronic entropy contribution of the whole cell as
    entropy T*S    EENTRO =        -0.00123456
in OUTCAR. The usual criterion is |T*S| / NIONS < 1 meV/atom.

Usage
  python vasp_tsigma_check.py OUTCAR
  python vasp_tsigma_check.py sigma_*/OUTCAR --thresh 1.0
"""
import argparse
import sys


def parse(path):
    eentro = nions = e0 = sigma = ismear = None
    with open(path, errors="replace") as fh:
        for line in fh:
            if "EENTRO" in line:
                eentro = float(line.split()[-1])
            elif "NIONS" in line and nions is None:
                nions = int(line.split()[-1])
            elif "energy  without entropy" in line:
                e0 = float(line.split()[-1])
            elif (ismear is None and line.strip().startswith("ISMEAR")
                  and "SIGMA" in line):
                # "   ISMEAR =     1;   SIGMA  =   0.20  broadening in eV ..."
                parts = line.replace(";", " ").split()
                try:
                    ismear = int(parts[parts.index("ISMEAR") + 2])
                    sigma = float(parts[parts.index("SIGMA") + 2])
                except (ValueError, IndexError):
                    pass
    return eentro, nions, e0, ismear, sigma


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outcars", nargs="+")
    p.add_argument("--thresh", type=float, default=1.0, help="meV/atom (default 1.0)")
    a = p.parse_args()
    print(f"{'file':40s} {'ISMEAR':>6s} {'SIGMA':>6s} {'T*S/atom(meV)':>14s} {'E0(eV)':>14s}  ok?")
    bad = 0
    for f in a.outcars:
        try:
            ts, n, e0, ism, sig = parse(f)
        except OSError as err:
            print(f"{f}: {err}", file=sys.stderr)
            bad += 1
            continue
        if ts is None or not n:
            print(f"{f:40s}  no EENTRO/NIONS found (run not started or ISMEAR=-5 without smearing?)")
            continue
        per = 1000.0 * ts / n
        ok = abs(per) < a.thresh
        bad += 0 if ok else 1
        print(f"{f:40s} {str(ism):>6s} {str(sig):>6s} {per:14.3f} {e0 if e0 is not None else float('nan'):14.6f}  "
              f"{'yes' if ok else 'NO -> lower SIGMA (and densify k-points)'}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
