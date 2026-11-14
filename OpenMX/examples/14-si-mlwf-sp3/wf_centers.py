#!/usr/bin/env python3
"""Summarize an OpenMX MLWF run from its standard output (the .std file).

OpenMX tags the optimization log with markers (see manual 47.2.3):
    ---> DISE   disentanglement loop: Iter, Omega_I, Delta_I
    ---> SPRD   spread components: Omega_I, Omega_D, Omega_OD, Tot_Omega (Ang^2)
    ---> CENT   "WF n ( x, y, z) | spread" lines, centres in Ang
This script prints the last DISE/SPRD values, whether "CONVERGENCE ACHIEVED"
appeared, and the final centre and spread of every Wannier function.

Optionally, --points gives reference positions (Ang), e.g. bond midpoints, and
--cell the three lattice vectors (Ang, 9 numbers); the distance from each centre
to the nearest periodic image of a reference point is then printed.

Usage:
    python3 wf_centers.py si_mlwf_sp3.std
    python3 wf_centers.py si_mlwf_bond4.std \
        --cell 0 2.715 2.715  2.715 0 2.715  2.715 2.715 0 \
        --points "0.67875,0.67875,0.67875;0.67875,-0.67875,-0.67875;-0.67875,0.67875,-0.67875;-0.67875,-0.67875,0.67875"
Standard library only (numpy not needed).
"""
import argparse
import itertools
import math
import re
import sys

NUM = r"[-+]?\d+\.?\d*(?:[EeDd][-+]?\d+)?"
RE_CENT = re.compile(r"WF\s+(\d+)\s*\(\s*(%s)\s*,\s*(%s)\s*,\s*(%s)\s*\)\s*\|\s*(%s)" % ((NUM,) * 5))


def parse(fname):
    dise, sprd, cent, nconv = None, None, {}, 0
    with open(fname, errors="replace") as f:
        for line in f:
            if "CONVERGENCE ACHIEVED" in line and "SPRD" in line:
                nconv += 1
            if "---> DISE" in line:
                nums = re.findall(NUM, line.split("--->")[0])
                if len(nums) >= 3:
                    dise = (int(float(nums[0])), float(nums[1]), float(nums[2]))
            elif "---> SPRD" in line:
                nums = re.findall(NUM, line.split("--->")[0])
                if len(nums) >= 5:
                    sprd = [float(v) for v in nums[-4:]]
            elif "CENT" in line:
                m = RE_CENT.search(line)
                if m:
                    cent[int(m.group(1))] = [float(m.group(i)) for i in range(2, 6)]
    return dise, sprd, cent, nconv


def nearest(c, points, cell):
    best = (1e9, None)
    shifts = [(0, 0, 0)] if cell is None else itertools.product((-1, 0, 1), repeat=3)
    for n1, n2, n3 in shifts:
        t = [0.0, 0.0, 0.0]
        if cell is not None:
            t = [n1 * cell[0][i] + n2 * cell[1][i] + n3 * cell[2][i] for i in range(3)]
        for ip, p in enumerate(points):
            d = math.dist(c, [p[i] + t[i] for i in range(3)])
            if d < best[0]:
                best = (d, ip + 1)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("std", help="standard output of the OpenMX MLWF run")
    ap.add_argument("--points", help='reference points "x,y,z;x,y,z;..." in Ang')
    ap.add_argument("--cell", type=float, nargs=9, help="lattice vectors a1 a2 a3 in Ang")
    a = ap.parse_args()

    dise, sprd, cent, nconv = parse(a.std)
    if not cent:
        sys.exit("no '---> CENT' lines found: is Wannier.Func.Calc on and level.of.stdout >= 1?")
    if dise:
        print("last DISE step %d: Omega_I = %.6f Ang^2 (Delta_I = %.2e)" % dise)
    else:
        print("no DISE lines (isolated bands or disentanglement skipped)")
    if sprd:
        print("last SPRD: Omega_I = %.5f  Omega_D = %.5f  Omega_OD = %.5f  Total = %.5f Ang^2" % tuple(sprd))
    print("'CONVERGENCE ACHIEVED' (SPRD) found: %s" % ("yes" if nconv else "NO - check Max.Steps / criteria"))

    points = None
    cell = None
    if a.points:
        points = [[float(v) for v in p.split(",")] for p in a.points.split(";") if p.strip()]
    if a.cell:
        cell = [a.cell[0:3], a.cell[3:6], a.cell[6:9]]
    print("\n WF        x (Ang)     y (Ang)     z (Ang)   spread (Ang^2)" +
          ("   d_nearest_ref (Ang)" if points else ""))
    total = 0.0
    for n in sorted(cent):
        x, y, z, s = cent[n]
        total += s
        extra = ""
        if points:
            d, ip = nearest((x, y, z), points, cell)
            extra = "   %8.4f (point %d)" % (d, ip)
        print("%4d  %11.5f %11.5f %11.5f %12.5f%s" % (n, x, y, z, s, extra))
    print("sum of spreads = %.5f Ang^2 (%d functions)" % (total, len(cent)))


if __name__ == "__main__":
    main()
