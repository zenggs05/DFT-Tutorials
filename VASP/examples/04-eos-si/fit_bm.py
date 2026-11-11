#!/usr/bin/env python3
"""Fit a third-order Birch-Murnaghan equation of state to E(V) data.

The 3rd-order BM energy is exactly a cubic polynomial in x = V**(-2/3):
    E(V) = c0 + c1*x + c2*x**2 + c3*x**3
so a linear least-squares fit is sufficient (no initial guess needed).
From the polynomial we get V0, E0, B0 = V*E''(V) and B0' = -1 - V*E'''/E''.

Input file: two columns, volume (A^3 per cell) and energy (eV per cell).
Lines starting with '#' are ignored.

Usage:
    python fit_bm.py eos.dat
    python fit_bm.py eos.dat --lattice fcc --plot
    python fit_bm.py eos.dat --natoms 2

Dependencies: numpy (required), matplotlib (only for --plot).
"""
import argparse
import sys

import numpy as np

EV_A3_TO_GPA = 160.21766208

# primitive-cell volume -> conventional cubic lattice constant
LATTICE_FACTOR = {"sc": 1.0, "bcc": 2.0, "fcc": 4.0}


def fit_bm3(volumes, energies):
    """Return (coeffs, V0, E0, B0[eV/A^3], B0p) for the 3rd-order BM EOS."""
    x = volumes ** (-2.0 / 3.0)
    coeffs = np.polyfit(x, energies, 3)          # highest power first
    poly = np.poly1d(coeffs)
    dpoly = poly.deriv(1)
    roots = dpoly.r
    # physical root: real, positive x, positive curvature, inside (or near) the data
    xmin, xmax = x.min(), x.max()
    cands = []
    for r in roots:
        if abs(r.imag) > 1e-12 or r.real <= 0:
            continue
        r = r.real
        if poly.deriv(2)(r) <= 0:
            continue
        cands.append(r)
    if not cands:
        raise RuntimeError("no energy minimum found in the fitted curve; extend the volume range")
    xc = 0.5 * (xmin + xmax)
    x0 = min(cands, key=lambda r: abs(r - xc))
    v0 = x0 ** (-1.5)
    e0 = poly(x0)

    # derivatives with respect to V, using E(V) = sum_n a_n V^(-2n/3)
    a = coeffs[::-1]                              # a[0] + a[1] x + a[2] x^2 + a[3] x^3

    def dE(v, order):
        total = 0.0
        for n, an in enumerate(a):
            p = -2.0 * n / 3.0
            fac = 1.0
            for k in range(order):
                fac *= (p - k)
            total += an * fac * v ** (p - order)
        return total

    e2 = dE(v0, 2)
    e3 = dE(v0, 3)
    b0 = v0 * e2
    b0p = -1.0 - v0 * e3 / e2
    return coeffs, v0, e0, b0, b0p


def main():
    ap = argparse.ArgumentParser(description="3rd-order Birch-Murnaghan fit of E(V) data")
    ap.add_argument("datafile", help="two columns: V (A^3/cell)  E (eV/cell)")
    ap.add_argument("--lattice", choices=sorted(LATTICE_FACTOR), default=None,
                    help="if V is the volume of a primitive sc/bcc/fcc cell, also print the cubic lattice constant")
    ap.add_argument("--natoms", type=int, default=None, help="atoms per cell, to print per-atom values")
    ap.add_argument("--plot", action="store_true", help="save eos_fit.png (needs matplotlib)")
    args = ap.parse_args()

    data = np.loadtxt(args.datafile, comments="#", ndmin=2)
    if data.shape[0] < 5:
        sys.exit("need at least 5 (V, E) points for a meaningful 4-parameter fit")
    order = np.argsort(data[:, 0])
    V, E = data[order, 0], data[order, 1]

    coeffs, v0, e0, b0, b0p = fit_bm3(V, E)
    resid = E - np.poly1d(coeffs)(V ** (-2.0 / 3.0))
    rms = np.sqrt(np.mean(resid ** 2))

    print("Birch-Murnaghan (3rd order) fit of %d points" % len(V))
    print("  V0  = %12.5f A^3/cell" % v0)
    print("  E0  = %12.6f eV/cell" % e0)
    print("  B0  = %12.2f GPa" % (b0 * EV_A3_TO_GPA))
    print("  B0' = %12.3f" % b0p)
    print("  RMS residual = %.2e eV" % rms)
    if args.natoms:
        print("  V0/atom = %.5f A^3, E0/atom = %.6f eV" % (v0 / args.natoms, e0 / args.natoms))
    if args.lattice:
        a0 = (LATTICE_FACTOR[args.lattice] * v0) ** (1.0 / 3.0)
        print("  cubic lattice constant a0 (%s) = %.5f A" % (args.lattice, a0))

    # sanity checks
    if not (V.min() < v0 < V.max()):
        print("WARNING: V0 lies outside the scanned range -> shift the reference volume and rescan")
    else:
        frac = (v0 - V.min()) / (V.max() - V.min())
        if frac < 0.2 or frac > 0.8:
            print("WARNING: V0 is close to the edge of the scanned range (%.0f%%); consider recentering" % (100 * frac))
    if rms > 1e-3:
        print("WARNING: RMS residual > 1 meV: check SCF convergence / magnetic state of outlier points")
    print("\n  V(A^3)        E(eV)        E-fit(meV)")
    for v, e, r in zip(V, E, resid):
        print("  %10.4f  %14.6f  %10.3f" % (v, e, 1000.0 * r))

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            sys.exit("matplotlib not available; rerun without --plot")
        vv = np.linspace(V.min() * 0.98, V.max() * 1.02, 300)
        ee = np.poly1d(coeffs)(vv ** (-2.0 / 3.0))
        fig, ax = plt.subplots(figsize=(4.5, 3.4))
        ax.plot(vv, ee, "-", lw=1.5, label="BM3 fit")
        ax.plot(V, E, "o", ms=5, label="VASP")
        ax.axvline(v0, ls="--", lw=0.8, color="gray")
        ax.set_xlabel(r"Volume ($\AA^3$/cell)")
        ax.set_ylabel("Energy (eV/cell)")
        ax.set_title(r"$V_0$=%.3f $\AA^3$, $B_0$=%.1f GPa" % (v0, b0 * EV_A3_TO_GPA), fontsize=9)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig("eos_fit.png", dpi=200)
        print("\nsaved eos_fit.png")


if __name__ == "__main__":
    main()
