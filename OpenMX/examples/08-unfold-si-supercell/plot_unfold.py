#!/usr/bin/env python3
"""Plot OpenMX band-unfolding results (System.Name.unfold_tot*, unfold_orb*).

Files written by OpenMX when Unfolding.Electronic.Band = on (manual Sec. 55):
    <prefix>.unfold_totup / _totdn (collinear) or <prefix>.unfold_tot (non-collinear)
        col 1: distance along the path (Bohr^-1), col 2: E - E_F (eV), col 3: weight
    <prefix>.unfold_orbup / _orbdn / _orb : col 3.. orbital-resolved weights,
        column = (orbital serial number in <prefix>.out) + 2
    <prefix>.unfold_plotexample : gnuplot script; its "set xtics(...)" line gives
        the positions of the high-symmetry points, which this script reuses.

Usage:
    python3 plot_unfold.py si8_unfold
    python3 plot_unfold.py si7ge1_unfold --band ../08-band-cdia-si/si.BANDDAT1 --out unfold.png
    python3 plot_unfold.py si7ge1_unfold --orb-cols 3,4,5,6 --out unfold_ge.png
Requires numpy and matplotlib.
"""
import argparse
import os
import re

import numpy as np


def read_xtics(fname):
    if not os.path.exists(fname):
        return []
    with open(fname) as f:
        for line in f:
            if line.strip().lower().startswith("set xtics"):
                return [(float(x), lab) for lab, x in re.findall(r"'([^']*)'\s+([-+0-9.eE]+)", line)]
    return []


def read_banddat(fname):
    """Read bandgnu13 BANDDAT: blocks separated by blank lines.
    The first blocks are vertical lines at the segment boundaries (constant x); skip them."""
    blocks, cur = [], []
    with open(fname) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                cur.append((float(p[0]), float(p[1])))
            elif cur:
                blocks.append(np.array(cur))
                cur = []
    if cur:
        blocks.append(np.array(cur))
    return [b for b in blocks if np.ptp(b[:, 0]) > 1e-8]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="System.Name of the unfolding run")
    ap.add_argument("--band", default="", help="optional BANDDAT1 of the primitive cell to overlay")
    ap.add_argument("--orb-cols", default="", help="sum these 1-based columns of unfold_orb* instead of the total weight")
    ap.add_argument("--scale", type=float, default=20.0, help="marker size per unit weight")
    ap.add_argument("--emin", type=float, default=-13.0)
    ap.add_argument("--emax", type=float, default=6.0)
    ap.add_argument("--out", default="", help="output image (default <prefix>_unfold.png)")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    kind = "orb" if args.orb_cols else "tot"
    files = [(args.prefix + ".unfold_%sup" % kind, "tab:red", "up"),
             (args.prefix + ".unfold_%sdn" % kind, "tab:blue", "down"),
             (args.prefix + ".unfold_%s" % kind, "tab:red", "nc")]
    files = [f for f in files if os.path.exists(f[0])]
    if not files:
        raise SystemExit("no %s.unfold_%s* file found" % (args.prefix, kind))

    fig, ax = plt.subplots(figsize=(6, 4.5))
    if args.band:
        for b in read_banddat(args.band):
            ax.plot(b[:, 0], b[:, 1], color="0.6", lw=0.8)
    for fname, color, tag in files:
        d = np.loadtxt(fname, ndmin=2)
        if args.orb_cols:
            cols = [int(c) - 1 for c in args.orb_cols.split(",")]
            w = d[:, cols].sum(axis=1)
        else:
            w = d[:, 2]
        sel = (d[:, 1] >= args.emin) & (d[:, 1] <= args.emax)
        ax.scatter(d[sel, 0], d[sel, 1], s=args.scale * np.clip(w[sel], 0, None),
                   color=color, alpha=0.6, lw=0, label=tag)
        print("%s: %d points, total weight range %.3f - %.3f" % (fname, len(d), w.min(), w.max()))
    ticks = read_xtics(args.prefix + ".unfold_plotexample")
    for x, _ in ticks:
        ax.axvline(x, color="k", lw=0.5)
    if ticks:
        ax.set_xticks([x for x, _ in ticks])
        ax.set_xticklabels([r"$\Gamma$" if lab.upper() == "G" else lab for _, lab in ticks])
        ax.set_xlim(ticks[0][0], ticks[-1][0])
    ax.axhline(0.0, color="gray", ls="--", lw=0.8)
    ax.set_ylim(args.emin, args.emax)
    ax.set_ylabel("E - E$_F$ (eV)")
    if len(files) > 1:
        ax.legend(loc="lower right")
    out = args.out or args.prefix + "_unfold.png"
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print("figure written to", out)


if __name__ == "__main__":
    main()
