#!/usr/bin/env python3
"""1D free-energy profile from a CP2K metadynamics HILLS file.

CP2K writes <PROJECT>-HILLS.metadynLog with one line per Gaussian:
    time[fs]  s_1..s_n  width_1..width_n  height[hartree]
(CV values and widths in CP2K internal units, e.g. bohr for distances).

Free energy estimate:  F(s) = - sum_j h_j * exp(-0.5 * ((s - s_j)/w_j)^2)  + const
For well-tempered runs CP2K stores the heights already multiplied by
(T + DELTA_T)/DELTA_T (see src/metadynamics.F), so the same formula applies to
standard and well-tempered metadynamics.  Walls are not included in HILLS.

Examples:
    python3 fes_from_hills.py nacl_wtmetad_dft-HILLS.metadynLog --bohr2ang
    python3 fes_from_hills.py HILLS.metadynLog --bohr2ang --grid 2.0 6.0 401 --blocks 4 -o fes.dat
For 2D/3D profiles use CP2K's "graph" tool or "plumed sum_hills".
Needs numpy (matplotlib only for --png).
"""
import argparse
import sys

import numpy as np

BOHR_ANG = 0.529177210903
HARTREE = {"kJ/mol": 2625.4996394799, "kcal/mol": 627.5094740631, "eV": 27.211386245988}


def read_hills(path):
    rows = []
    with open(path, errors="replace") as fh:
        for line in fh:
            p = line.split()
            if not p or p[0].startswith("#"):
                continue
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                continue
    if not rows:
        sys.exit(f"no data in {path}")
    a = np.array(rows)
    ncv = (a.shape[1] - 2) // 2
    if a.shape[1] != 2 * ncv + 2:
        sys.exit(f"unexpected number of columns ({a.shape[1]}) in {path}")
    return a, ncv


def fes(grid, centers, widths, heights, period=None):
    f = np.zeros_like(grid)
    for c, w, h in zip(centers, widths, heights):
        d = grid - c
        if period:
            d = d - period * np.round(d / period)
        f -= h * np.exp(-0.5 * (d / w) ** 2)
    return f


def extrema(x, y):
    out = []
    for i in range(1, len(y) - 1):
        if y[i] < y[i - 1] and y[i] <= y[i + 1]:
            out.append(("min", x[i], y[i]))
        elif y[i] > y[i - 1] and y[i] >= y[i + 1]:
            out.append(("max", x[i], y[i]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("hills", help="<PROJECT>-HILLS.metadynLog")
    ap.add_argument("--cv", type=int, default=1, help="which CV to project on if the file has several (sum over others is NOT done)")
    ap.add_argument("--bohr2ang", action="store_true", help="CV is a distance: convert bohr -> Angstrom")
    ap.add_argument("--periodic", type=float, default=None, help="period of the CV in its (converted) unit, e.g. 6.2831853")
    ap.add_argument("--grid", nargs=3, type=float, metavar=("MIN", "MAX", "N"), help="grid (default: range of centers +/- 3 widths, 400 points)")
    ap.add_argument("--unit", choices=list(HARTREE), default="kJ/mol")
    ap.add_argument("--blocks", type=int, default=1, help="also write FES after 1/B, 2/B, ... of the hills (convergence check)")
    ap.add_argument("-o", "--out", default="fes.dat")
    ap.add_argument("--png", help="plot to this file (needs matplotlib)")
    args = ap.parse_args()

    a, ncv = read_hills(args.hills)
    if ncv != 1 and args.cv > ncv:
        sys.exit(f"file has {ncv} CVs")
    k = args.cv
    time = a[:, 0]
    s = a[:, k].copy()
    w = a[:, ncv + k].copy()
    h = a[:, -1] * HARTREE[args.unit]
    if ncv > 1:
        print(f"# warning: {ncv} CVs in file; projecting naively on CV {k} (other CVs ignored)", file=sys.stderr)
    if args.bohr2ang:
        s *= BOHR_ANG
        w *= BOHR_ANG
    if args.grid:
        x = np.linspace(args.grid[0], args.grid[1], int(args.grid[2]))
    else:
        x = np.linspace(s.min() - 3 * w.max(), s.max() + 3 * w.max(), 400)

    nb = max(1, args.blocks)
    cols, labels = [], []
    for b in range(1, nb + 1):
        n = int(round(len(h) * b / nb))
        f = fes(x, s[:n], w[:n], h[:n], args.periodic)
        cols.append(f - f.min())
        labels.append(f"F_{n}hills")
    data = np.column_stack([x] + cols)
    np.savetxt(args.out, data, fmt="%14.6f",
               header=f"CV({'A' if args.bohr2ang else 'internal'}) " + " ".join(labels) + f"   [{args.unit}]")
    print(f"# {len(h)} hills, t = {time[0]:.1f} ... {time[-1]:.1f} fs, written {args.out}")
    print(f"# last Gaussian height: {h[-1]:.4f} {args.unit} (first: {h[0]:.4f}) -- should decay in WT-MetaD")
    for kind, xi, yi in extrema(x, cols[-1]):
        print(f"  local {kind}: CV = {xi:8.3f}   F = {yi:8.3f} {args.unit}")

    if args.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.8, 3.4))
        for c, lab in zip(cols, labels):
            ax.plot(x, c, label=lab)
        ax.set_xlabel("CV (A)" if args.bohr2ang else "CV (internal units)")
        ax.set_ylabel(f"F ({args.unit})")
        if nb > 1:
            ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(args.png, dpi=200)
        print(f"# plot written to {args.png}")


if __name__ == "__main__":
    main()
