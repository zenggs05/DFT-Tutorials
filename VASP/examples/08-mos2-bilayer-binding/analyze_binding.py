#!/usr/bin/env python3
"""Interlayer binding energy of a bilayer from a rigid distance scan.

    E_b(d) = [E_bilayer(d) - 2 E_mono] / A

Prints meV/A^2 and meV per formula unit (MoS2), fits a cubic spline/polynomial
around the minimum to estimate d_eq, and saves binding_curve.png.

Usage: python analyze_binding.py --mono mono/OUTCAR --dirs d_* --a 3.19 [--nfu 1]
"""
import argparse
import math
import os
import re
import sys

import numpy as np


def energy(outcar):
    e = None
    with open(outcar) as fh:
        for line in fh:
            m = re.search(r"energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)", line)
            if m:
                e = float(m.group(1))
    if e is None:
        sys.exit(f"no energy in {outcar}")
    return e


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mono", required=True)
    p.add_argument("--dirs", nargs="+", required=True)
    p.add_argument("--a", type=float, required=True, help="hexagonal lattice constant (A)")
    p.add_argument("--nfu", type=int, default=1, help="formula units per layer in the cell")
    p.add_argument("--no-plot", action="store_true")
    o = p.parse_args()

    area = math.sqrt(3.0) / 2.0 * o.a ** 2
    e_mono = energy(o.mono)
    rows = []
    for d in o.dirs:
        m = re.search(r"d_(\d+(?:\.\d+)?)$", os.path.normpath(d))
        out = os.path.join(d, "OUTCAR")
        if not m or not os.path.exists(out):
            continue
        eb = energy(out) - 2.0 * e_mono                     # eV per cell
        rows.append((float(m.group(1)), eb * 1000.0 / area, eb * 1000.0 / o.nfu))
    if not rows:
        sys.exit("no finished d_* directories")
    rows.sort()
    print(f"{'d(Mo-Mo) A':>11s} {'E_b meV/A^2':>13s} {'E_b meV/f.u.':>13s}")
    for d, ea, ef in rows:
        print(f"{d:11.2f} {ea:13.2f} {ef:13.1f}")

    d, ea, _ = map(np.array, zip(*rows))
    i = int(np.argmin(ea))
    lo, hi = max(i - 2, 0), min(i + 3, len(d))
    if hi - lo >= 4:
        c = np.polyfit(d[lo:hi], ea[lo:hi], 3)
        dd = np.linspace(d[lo], d[hi - 1], 400)
        fit = np.polyval(c, dd)
        k = int(np.argmin(fit))
        print(f"\ncubic fit around the minimum: d_eq ~ {dd[k]:.3f} A, "
              f"E_b ~ {fit[k]:.2f} meV/A^2")
    else:
        print("\nnot enough points around the minimum for a fit")

    if not o.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.5, 3.4))
        ax.plot(d, ea, "o-", ms=4)
        ax.axhline(0.0, color="grey", lw=0.8)
        ax.set_xlabel("Mo-Mo interlayer distance d (A)")
        ax.set_ylabel("E_b (meV/A$^2$)")
        fig.tight_layout()
        fig.savefig("binding_curve.png", dpi=200)
        print("saved binding_curve.png")


if __name__ == "__main__":
    main()
