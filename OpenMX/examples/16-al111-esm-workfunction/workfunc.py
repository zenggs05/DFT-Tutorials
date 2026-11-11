#!/usr/bin/env python3
"""Work function of a slab from an OpenMX planar-averaged potential.

Input 1: the 1D file written by gcube2oned (columns: index, position / Ang,
         bc-plane averaged potential / Hartree), e.g.
             ./gcube2oned al111_esm_on1_2e.v0.cube 1 > 1d_pot.txt
Input 2: the OpenMX *.out file of the same run (chemical potential, Hartree).

The vacuum level is taken as the mean potential inside a window near each
cell boundary (ESM on1: vacuum on both sides), and

    Phi = (V_vac - mu) * 27.211386  [eV]

is printed for the left (x ~ 0) and the right (x ~ |a|) vacuum separately,
so that asymmetric slabs (adsorbates, polar surfaces) are handled too.

Usage:
    python workfunc.py 1d_pot.txt al111_esm_on1_2e.out
    python workfunc.py 1d_pot.txt al111_esm_on1_2e.out --left 3 8 --right 42 47 --plot wf.png

Only numpy is required; matplotlib is needed for --plot.
"""
import argparse
import re
import sys

import numpy as np

HARTREE2EV = 27.211386


def read_mu(outfile):
    """Return the last 'Chemical potential (Hartree)' value found in an .out file."""
    pat = re.compile(r"^\s*Chemical potential \(Hartree\)\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)")
    mu = None
    with open(outfile, errors="replace") as f:
        for line in f:
            m = pat.match(line)
            if m:
                mu = float(m.group(1))
    if mu is None:
        sys.exit(f"error: 'Chemical potential (Hartree)' not found in {outfile}")
    return mu


def plateau(x, v, lo, hi, label):
    sel = (x >= lo) & (x <= hi)
    if sel.sum() < 3:
        sys.exit(f"error: fewer than 3 grid points in the {label} window [{lo}, {hi}] Ang")
    return v[sel].mean(), v[sel].std(), sel.sum()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pot1d", help="output of gcube2oned (3 columns)")
    ap.add_argument("outfile", help="OpenMX System.Name.out")
    ap.add_argument("--left", nargs=2, type=float, metavar=("X1", "X2"),
                    help="vacuum window near x=0 in Ang (default: 5%%-15%% of the cell)")
    ap.add_argument("--right", nargs=2, type=float, metavar=("X1", "X2"),
                    help="vacuum window near x=|a| in Ang (default: 85%%-95%% of the cell)")
    ap.add_argument("--plot", metavar="PNG", help="write a figure of the averaged potential")
    args = ap.parse_args()

    data = np.loadtxt(args.pot1d)
    if data.ndim != 2 or data.shape[1] < 3:
        sys.exit("error: expected at least 3 columns (index, x/Ang, V/Hartree)")
    # column 1 is not a reliable index for axis 2/3 in gcube2oned; use columns 2 and 3 only
    x, v = data[:, 1], data[:, 2]
    length = x[-1] + (x[1] - x[0])
    left = args.left or (0.05 * length, 0.15 * length)
    right = args.right or (0.85 * length, 0.95 * length)

    mu = read_mu(args.outfile)
    vl, sl, nl = plateau(x, v, *left, "left")
    vr, sr, nr = plateau(x, v, *right, "right")

    print(f"# cell length along the averaging axis : {length:10.4f} Ang")
    print(f"# chemical potential mu               : {mu: .8f} Hartree = {mu*HARTREE2EV: .4f} eV")
    for name, vv, ss, nn, win in (("left ", vl, sl, nl, left), ("right", vr, sr, nr, right)):
        print(f"{name} vacuum [{win[0]:6.2f},{win[1]:6.2f}] Ang: V = {vv: .6f} Ha "
              f"(std {ss:.1e}, {nn} pts)  ->  Phi = {(vv - mu)*HARTREE2EV:7.3f} eV")
    if max(sl, sr) * HARTREE2EV > 0.01:
        print("warning: potential is not flat (> 10 meV) inside a window; "
              "enlarge the vacuum or move the window away from the slab")
    dv = (vr - vl) * HARTREE2EV
    print(f"vacuum-level difference V_right - V_left = {dv:7.3f} eV "
          "(non-zero for asymmetric/polar slabs)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 3.5))
        ax.plot(x, (v - mu) * HARTREE2EV, lw=1.2, color="#1f4e79", label=r"$\bar V(x)-\mu$")
        ax.axhline(0.0, color="0.5", lw=0.8, ls="--")
        for win in (left, right):
            ax.axvspan(*win, color="#f2c14e", alpha=0.3)
        ax.set_xlabel("x along a-axis (Ang)")
        ax.set_ylabel("energy relative to mu (eV)")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200)
        print(f"figure written to {args.plot}")


if __name__ == "__main__":
    main()
