#!/usr/bin/env python3
"""Fit 2D (in-plane) elastic constants from the energy-strain runs created by
run_strain.sh.

    python fit_2d_elastic.py                 # reads ./POSCAR and strain_*/*/OUTCAR
    python fit_2d_elastic.py --deg 3 --plot

Energy per unit area of the undeformed cell A0 (eps in Voigt form,
engineering shear gamma = eps6):

    dE/A0 = 1/2 C11 e1^2 + 1/2 C22 e2^2 + C12 e1 e2 + 1/2 C66 e6^2

    xx      -> dE/A0 = 1/2 C11 d^2
    yy      -> dE/A0 = 1/2 C22 d^2
    biaxial -> dE/A0 = 1/2 (C11 + C22 + 2 C12) d^2
    shear   -> dE/A0 = 1/2 C66 d^2

A polynomial E(d) = c0 + c1 d + c2 d^2 (+ c3 d^3) is fitted for each mode;
C = 2 c2 / A0 in eV/A^2, converted to N/m (1 eV/A^2 = 16.0218 N/m).
A large linear coefficient c1 means the reference cell was not relaxed.
Only numpy is required (matplotlib for --plot).
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

EV_A2_TO_N_M = 16.02176634


def area_from_poscar(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
    lat = lat * scale if scale > 0 else lat
    return float(np.linalg.norm(np.cross(lat[0], lat[1])))


def final_energy(outcar):
    """Last 'energy(sigma->0)' value in OUTCAR (eV)."""
    val = None
    with open(outcar, errors="replace") as fh:
        for line in fh:
            if "energy(sigma->0)" in line:
                val = float(line.split("=")[-1])
    return val


def collect(mode):
    data = []
    for d in sorted(glob.glob(f"strain_{mode}/*")):
        out = os.path.join(d, "OUTCAR")
        if not os.path.isfile(out):
            continue
        e = final_energy(out)
        if e is None:
            print(f"warning: no energy in {out}", file=sys.stderr)
            continue
        with open(out, errors="replace") as fh:
            if "reached required accuracy" not in fh.read():
                print(f"warning: ionic relaxation not converged in {d}", file=sys.stderr)
        data.append((float(os.path.basename(d)), e))
    return np.array(sorted(data))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--poscar", default="POSCAR", help="undeformed cell (for A0)")
    ap.add_argument("--deg", type=int, default=2, choices=(2, 3),
                    help="polynomial degree of the E(d) fit")
    ap.add_argument("--plot", action="store_true", help="write energy_strain.png")
    args = ap.parse_args()

    a0 = area_from_poscar(args.poscar)
    print(f"A0 = {a0:.4f} A^2 (from {args.poscar})\n")
    print(f"{'mode':>8} {'npts':>5} {'c1 (eV)':>10} {'2c2/A0 (N/m)':>14} {'max|res| (meV)':>15}")
    k = {}
    fits = {}
    for mode in ("xx", "yy", "biaxial", "shear"):
        arr = collect(mode)
        if len(arr) < args.deg + 2:
            print(f"{mode:>8}  not enough points ({len(arr)})")
            continue
        d, e = arr[:, 0], arr[:, 1]
        coef = np.polyfit(d, e, args.deg)          # highest power first
        c2, c1 = coef[-3], coef[-2]
        res = e - np.polyval(coef, d)
        k[mode] = 2 * c2 / a0 * EV_A2_TO_N_M
        fits[mode] = (d, e, coef)
        print(f"{mode:>8} {len(d):5d} {c1:10.4f} {k[mode]:14.2f} {1000 * np.abs(res).max():15.3f}")

    if not {"xx", "yy", "biaxial"} <= set(k):
        sys.exit("\nneed at least the xx, yy and biaxial series")
    c11, c22 = k["xx"], k["yy"]
    c12 = (k["biaxial"] - c11 - c22) / 2
    print(f"\nC11 = {c11:8.2f} N/m\nC22 = {c22:8.2f} N/m\nC12 = {c12:8.2f} N/m")
    if "shear" in k:
        c66 = k["shear"]
        print(f"C66 = {c66:8.2f} N/m   (hexagonal check: (C11-C12)/2 = {(c11 - c12) / 2:.2f})")
    else:
        c66 = None

    det = c11 * c22 - c12 ** 2
    print(f"\nY_x  = (C11 C22 - C12^2)/C22 = {det / c22:8.2f} N/m")
    print(f"Y_y  = (C11 C22 - C12^2)/C11 = {det / c11:8.2f} N/m")
    print(f"nu_xy = C12/C22 = {c12 / c22:.3f}    nu_yx = C12/C11 = {c12 / c11:.3f}")
    ok = c11 > 0 and c22 > 0 and det > 0 and (c66 is None or c66 > 0)
    print("2D Born criteria (C11>0, C22>0, C11*C22>C12^2, C66>0):", "satisfied" if ok else "VIOLATED")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 4))
        for mode, (d, e, coef) in fits.items():
            x = np.linspace(d.min(), d.max(), 100)
            e0 = np.polyval(coef, 0.0)
            ax.plot(100 * d, 1000 * (e - e0), "o", label=mode)
            ax.plot(100 * x, 1000 * (np.polyval(coef, x) - e0), "-", color=ax.lines[-1].get_color())
        ax.set_xlabel("strain (%)")
        ax.set_ylabel("E - E0 (meV/cell)")
        ax.legend()
        fig.tight_layout()
        fig.savefig("energy_strain.png", dpi=200)
        print("\nwrote energy_strain.png")


if __name__ == "__main__":
    main()
