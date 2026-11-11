#!/usr/bin/env python3
"""plot_conv.py - analyse and plot a convergence scan written by run_encut.sh / run_kpts.sh

The input is a whitespace separated table whose first line starts with '#'
and lists the column names, e.g.

    # ENCUT natoms E0 E0_per_atom TOTEN EENTRO pressure_kB maxF_eVA nkpts time_s

For every row the deviation from a reference row is computed:
  * energy per atom          (meV/atom)   threshold --thr-e    (default 1.0)
  * maximum force            (meV/A)      threshold --thr-f    (default 5.0)
  * pressure / stress        (kB)         threshold --thr-p    (default 1.0)
The reference is the most converged row: largest x (ENCUT, N) by default,
or smallest x with --ref min (KSPACING scans).

A value of x is called converged when it AND every more accurate value lie
within the thresholds. The figure is written to <input>.png.

Usage:
    python plot_conv.py encut_POSCAR.dat
    python plot_conv.py encut_POSCAR_displaced.dat --thr-f 2
    python plot_conv.py kpts_kspacing.dat --ref min
Requires: python >= 3.9, numpy, matplotlib
"""
import argparse
import math
import sys

import numpy as np


def read_table(path):
    with open(path, encoding="ascii") as fh:
        header = fh.readline()
        if not header.startswith("#"):
            sys.exit(f"{path}: first line must be a '#' header with column names")
        names = header[1:].split()
        rows = []
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            vals = []
            for tok in line.split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    vals.append(math.nan)
            rows.append(vals)
    data = np.array(rows, dtype=float)
    if data.ndim != 2 or data.shape[1] != len(names):
        sys.exit(f"{path}: inconsistent number of columns")
    return names, data


def converged_from(x, dev, thr, ref):
    """Smallest-effort x such that it and all more accurate points are within thr."""
    order = np.argsort(x) if ref == "max" else np.argsort(-x)   # least -> most accurate
    ok = np.abs(dev[order]) <= thr
    first = None
    for i in range(len(order)):
        if np.all(ok[i:]) and not np.any(np.isnan(dev[order][i:])):
            first = x[order][i]
            break
    return first


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("table", help="data file written by run_encut.sh or run_kpts.sh")
    ap.add_argument("--ref", choices=["max", "min"], default="max",
                    help="reference row: largest x (default) or smallest x (KSPACING)")
    ap.add_argument("--thr-e", type=float, default=1.0, help="energy threshold in meV/atom")
    ap.add_argument("--thr-f", type=float, default=5.0, help="force threshold in meV/A")
    ap.add_argument("--thr-p", type=float, default=1.0, help="pressure threshold in kB")
    ap.add_argument("--no-plot", action="store_true", help="print the table only")
    args = ap.parse_args()

    names, d = read_table(args.table)
    col = {n: i for i, n in enumerate(names)}
    x = d[:, 0]
    good = ~np.isnan(d[:, col["E0_per_atom"]])
    if good.sum() < 2:
        sys.exit("need at least two finished calculations")
    d, x = d[good], x[good]
    iref = int(np.argmax(x) if args.ref == "max" else np.argmin(x))

    de = 1000.0 * (d[:, col["E0_per_atom"]] - d[iref, col["E0_per_atom"]])   # meV/atom
    df = 1000.0 * (d[:, col["maxF_eVA"]] - d[iref, col["maxF_eVA"]])         # meV/A
    dp = d[:, col["pressure_kB"]] - d[iref, col["pressure_kB"]]               # kB

    print(f"reference: {names[0]} = {x[iref]:g}")
    print(f"{names[0]:>10s} {'nkpts':>6s} {'E0/atom (eV)':>14s} {'dE (meV/at)':>12s} "
          f"{'P (kB)':>9s} {'dP (kB)':>8s} {'Fmax (eV/A)':>12s} {'dF (meV/A)':>11s}")
    order = np.argsort(x) if args.ref == "max" else np.argsort(-x)
    for i in order:
        print(f"{x[i]:10g} {d[i, col['nkpts']]:6.0f} {d[i, col['E0_per_atom']]:14.6f} {de[i]:12.3f} "
              f"{d[i, col['pressure_kB']]:9.2f} {dp[i]:8.2f} {d[i, col['maxF_eVA']]:12.5f} {df[i]:11.2f}")

    for label, dev, thr, unit in (("energy", de, args.thr_e, "meV/atom"),
                                  ("pressure", dp, args.thr_p, "kB"),
                                  ("max. force", df, args.thr_f, "meV/A")):
        xc = converged_from(x, dev, thr, args.ref)
        msg = f"{xc:g}" if xc is not None else "not reached in this scan"
        print(f"{label:>10s} converged within {thr:g} {unit}: {names[0]} = {msg}")

    if args.no_plot:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(5.5, 8.0), sharex=True)
    panels = ((de, args.thr_e, "dE (meV/atom)"),
              (dp, args.thr_p, "dP (kB)"),
              (df, args.thr_f, "dFmax (meV/A)"))
    for ax, (dev, thr, ylabel) in zip(axes, panels):
        ax.axhspan(-thr, thr, color="0.85", zorder=0, label=f"+/- {thr:g}")
        ax.plot(x[order], dev[order], "o-", color="C0", zorder=2)
        ax.axhline(0.0, color="0.4", lw=0.8)
        ax.set_ylabel(ylabel)
        ax.legend(loc="best", frameon=False)
    axes[-1].set_xlabel(names[0])
    if args.ref == "min":
        axes[-1].invert_xaxis()
    fig.suptitle(f"Convergence: {args.table}")
    fig.tight_layout()
    png = args.table.rsplit(".", 1)[0] + ".png"
    fig.savefig(png, dpi=150)
    print(f"figure written to {png}")


if __name__ == "__main__":
    main()
