#!/usr/bin/env python3
"""Transmission spectrum, Landauer conductance and Seebeck coefficient from OpenMX NEGF output.

OpenMX Step 3 writes one file System.Name.tran<i2>_<i3> per (k2, k3) point of NEGF.tran.Kgrid.
Columns (see the file header): iw, Re E (Hartree), Im E (Hartree), Re E (eV), Im E (eV),
Re T_up, Im T_up, Re T_down, Im T_down.  E is measured from the chemical potential of the
LEFT lead.  This script averages all given k files with equal weights (uniform grid),
sums the two spin channels, and evaluates at temperature T (zero-bias linear response):

    L_n = int T(E) (E - mu)^n (-df/dE) dE
    G   = (e^2/h) L_0            (T summed over spins; G0 = e^2/h per spin channel)
    S   = -(1/(e T)) L_1 / L_0    (Seebeck coefficient, V/K)

Usage:
    python plot_tran.py negf-chain.tran0_0
    python plot_tran.py negf-graphene.tran*_* --temp 300 --mu 0.0 --plot tran.png --csv tran_avg.csv
Requires numpy; matplotlib only for --plot.
"""
import argparse
import glob

import numpy as np

KB_EV = 8.617333262e-5       # eV/K


def load(fn):
    d = np.loadtxt(fn, comments="#")
    return d[:, 3], d[:, 5], d[:, 7]          # E (eV), T_up, T_down


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="System.Name.tran*_* files (shell globs allowed)")
    ap.add_argument("--temp", type=float, default=300.0, help="temperature in K for G and S")
    ap.add_argument("--mu", type=float, nargs="*", default=[0.0],
                    help="chemical potential(s) in eV relative to mu(left lead) for G and S")
    ap.add_argument("--plot", metavar="PNG")
    ap.add_argument("--csv", metavar="CSV", help="write E, T_up, T_down, T_total (k-averaged)")
    args = ap.parse_args()

    files = sorted({f for p in args.files for f in (glob.glob(p) or [p])})
    E, tu, td = load(files[0])
    tu, td = tu.copy(), td.copy()
    for fn in files[1:]:
        e2, a, b = load(fn)
        if len(e2) != len(E) or not np.allclose(e2, E):
            raise SystemExit(f"energy grid of {fn} differs from {files[0]}")
        tu += a
        td += b
    tu /= len(files)
    td /= len(files)
    tt = tu + td
    print(f"# {len(files)} k-point file(s) averaged; E grid {E[0]:.3f} .. {E[-1]:.3f} eV, {len(E)} points")

    kt = KB_EV * args.temp
    print(f"# T = {args.temp:.1f} K")
    print("#   mu-mu_L(eV)   T(mu)    G (e^2/h)     S (uV/K)")
    for mu in args.mu:
        x = (E - mu) / kt
        mdf = np.where(np.abs(x) < 40, 1.0 / (4.0 * kt * np.cosh(np.clip(x, -40, 40) / 2.0) ** 2), 0.0)
        l0 = np.trapezoid(tt * mdf, E) if hasattr(np, "trapezoid") else np.trapz(tt * mdf, E)
        l1 = np.trapezoid(tt * (E - mu) * mdf, E) if hasattr(np, "trapezoid") else np.trapz(tt * (E - mu) * mdf, E)
        t_at = np.interp(mu, E, tt)
        seebeck = -(l1 / l0) / args.temp * 1e6 if l0 > 1e-12 else float("nan")
        print(f"  {mu:12.4f} {t_at:8.4f} {l0:12.5f} {seebeck:12.3f}")

    if args.csv:
        np.savetxt(args.csv, np.column_stack([E, tu, td, tt]), fmt="%.6e",
                   header="E_eV T_up T_down T_total", comments="# ")
        print(f"k-averaged transmission written to {args.csv}")
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.plot(E, tu, lw=1.2, color="#1f4e79", label="spin up")
        if np.max(np.abs(td - tu)) > 1e-6:
            ax.plot(E, td, lw=1.2, ls="--", color="#c0504d", label="spin down")
        ax.axvline(0.0, color="0.5", lw=0.8, ls=":")
        ax.set_xlabel(r"$E-\mu_L$ (eV)")
        ax.set_ylabel("transmission per spin")
        ax.set_xlim(E[0], E[-1])
        ax.set_ylim(bottom=0)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200)
        print(f"figure written to {args.plot}")


if __name__ == "__main__":
    main()
