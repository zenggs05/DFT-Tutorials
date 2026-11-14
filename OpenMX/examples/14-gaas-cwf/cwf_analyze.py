#!/usr/bin/env python3
"""Summarize an OpenMX 4.0 CWF calculation.

Reads
  * System.Name.out : the "Populations evaluated by CWF" section
                      (DM_func = minimized distance measure, CWF populations)
  * System.Name.CWF.Dis_vs_H (CWF.Dis.vs.H on): one line per hopping
        collinear     : spin GA GB l1 l2 l3 i j distance(Ang) t(eV)
        non-collinear : GA GB l1 l2 l3 i j distance(Ang) t(eV)   (no spin column)
and prints DM_func, the populations per atom, the on-site energies of the
CWFs and the largest |t| in distance shells (to judge how short-ranged the
tight-binding model is).

Usage:
    python3 cwf_analyze.py gaas_cwf_ao                 # System.Name (reads .out and .CWF.Dis_vs_H)
    python3 cwf_analyze.py gaas_cwf_ao --rmax 10 --plot gaas_cwf_ao_t.png
numpy is not needed; matplotlib only for --plot.
"""
import argparse
import os
import re
import sys

NUM = r"[-+]?\d+\.?\d*(?:[EeDd][-+]?\d+)?"


def parse_out(fname):
    res = dict(dm=None, atoms=[], sums=[])
    if not os.path.exists(fname):
        print("(no %s)" % fname)
        return res
    inside = False
    with open(fname, errors="replace") as f:
        for line in f:
            if "Populations evaluated by CWF" in line and "Sum" not in line:
                inside, res["atoms"] = True, []
                continue
            if not inside:
                continue
            if "decomposed populations" in line.lower():
                inside = False
                continue
            if line.strip().startswith("DM_func"):
                nums = re.findall(NUM, line)
                res["dm"] = (float(nums[0]), float(nums[1]))
            elif "=" in line and ("up" in line or "total" in line):
                res["sums"].append(line.rstrip())
            else:
                m = re.match(r"\s*(\d+)\s+([A-Za-z]\S*)\s+(%s)\s+(%s)\s+(%s)\s+(%s)\s*$" % ((NUM,) * 4), line)
                if m:
                    res["atoms"].append((int(m.group(1)), m.group(2),
                                         *[float(m.group(i)) for i in range(3, 7)]))
    return res


def parse_dis_vs_h(fname):
    onsite, hops = [], []
    with open(fname, errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) == 10:
                spin, p = int(p[0]), p[1:]
            elif len(p) == 9:
                spin = 0
            else:
                continue
            ga, gb, l1, l2, l3, i, j = (int(v) for v in p[:7])
            d, t = float(p[7]), float(p[8])
            if ga == gb and (l1, l2, l3) == (0, 0, 0) and i == j:
                onsite.append((spin, ga, i, t))
            else:
                hops.append((d, abs(t)))
    return onsite, hops


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("name", help="System.Name of the CWF run")
    ap.add_argument("--rmax", type=float, default=8.0, help="largest distance to tabulate (Ang)")
    ap.add_argument("--dr", type=float, default=0.05, help="shell width (Ang)")
    ap.add_argument("--plot", help="write log|t| vs distance (png/pdf)")
    a = ap.parse_args()

    out = parse_out(a.name + ".out")
    if out["dm"]:
        print("DM_func = %.6f   (%.6f per CWF and spin)" % out["dm"])
        print("  rule of thumb (manual 47.1.1): ~0.16 is small, > 0.7 means the guiding")
        print("  functions do not match the symmetry of the selected Bloch states")
    else:
        print("no DM_func line in %s.out (is CWF.Calc on?)" % a.name)
    if out["atoms"]:
        print("\nCWF populations (from .out):")
        print("  atom        up         down        sum        diff")
        for n, el, up, dn, s, df in out["atoms"]:
            print("  %3d %-4s %10.5f %10.5f %10.5f %10.5f" % (n, el, up, dn, s, df))
        for line in out["sums"]:
            print(line)

    fdh = a.name + ".CWF.Dis_vs_H"
    if not os.path.exists(fdh):
        sys.exit("\n%s not found: set CWF.Dis.vs.H on" % fdh)
    onsite, hops = parse_dis_vs_h(fdh)
    print("\non-site energies of the CWFs (eV, as written in %s):" % fdh)
    for spin, ga, i, t in onsite:
        print("  spin %d  atom %3d  CWF %2d : %9.4f" % (spin, ga, i, t))

    shells = {}
    for d, t in hops:
        if d <= a.rmax:
            key = round(d / a.dr) * a.dr
            shells[key] = max(shells.get(key, 0.0), t)
    print("\nlargest |t| per distance shell (eV):")
    print("   d (Ang)     max|t|")
    for d in sorted(shells):
        print("  %7.3f  %10.5f" % (d, shells[d]))

    if a.plot:
        import math
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        xs = [d for d, t in hops if t > 1e-8]
        ys = [math.log10(t) for d, t in hops if t > 1e-8]
        fig, ax = plt.subplots(figsize=(4.5, 3.5))
        ax.plot(xs, ys, ".", ms=2)
        ax.set_xlabel("distance between CWF centres (Ang)")
        ax.set_ylabel("log10 |t| (eV)")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=200)
        print("figure written to", a.plot)


if __name__ == "__main__":
    main()
