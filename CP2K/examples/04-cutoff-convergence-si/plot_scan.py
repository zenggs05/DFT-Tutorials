#!/usr/bin/env python3
"""Plot a CUTOFF or REL_CUTOFF scan written by collect_scan.sh.

Example:
    python plot_scan.py cutoff_scan.csv --threshold 0.03 -o cutoff_scan.png
    python plot_scan.py relcutoff_scan.csv -o relcutoff_scan.png

The y axis is |E - E_ref| in meV/atom (log scale), where E_ref is the energy
at the highest scanned value (last row of the CSV). Requires numpy, matplotlib.
"""
import argparse
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read_csv(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{path}: no data rows")
    xkey = list(rows[0].keys())[0]
    x = np.array([float(r[xkey]) for r in rows])
    de = np.array([float(r["dE_meV_per_atom"]) for r in rows])
    conv = [r["scf_converged"] for r in rows]
    return xkey, x, de, conv


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="cutoff_scan.csv or relcutoff_scan.csv")
    ap.add_argument("--threshold", type=float, default=0.03,
                    help="convergence threshold in meV/atom (default 0.03 ~ 1e-6 Ha/atom)")
    ap.add_argument("-o", "--output", default=None, help="output image (default: <csv>.png)")
    args = ap.parse_args()

    xkey, x, de, conv = read_csv(args.csv)
    # drop the reference point (exactly zero) from the log plot
    mask = np.abs(de) > 0
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.semilogy(x[mask], np.abs(de[mask]), "o-", color="#1f5aa6", label="|E - E_ref|")
    bad = [i for i, c in enumerate(conv) if c != "yes" and mask[i]]
    if bad:
        ax.semilogy(x[bad], np.abs(de[bad]), "x", color="#c0392b", ms=9, label="SCF not converged")
    ax.axhline(args.threshold, ls="--", color="0.4", lw=1,
               label=f"threshold {args.threshold:g} meV/atom")
    ax.set_xlabel(xkey.replace("_Ry", "") + " (Ry)")
    ax.set_ylabel("|dE| (meV/atom)")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = args.output or args.csv.rsplit(".", 1)[0] + ".png"
    fig.savefig(out, dpi=200)
    print(f"saved {out}")

    # first value from which all following points stay below the threshold
    ok = np.abs(de) < args.threshold
    for i in range(len(x)):
        if ok[i:].all():
            print(f"converged (< {args.threshold:g} meV/atom) from {xkey} = {x[i]:g} Ry")
            break


if __name__ == "__main__":
    main()
