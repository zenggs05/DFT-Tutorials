#!/usr/bin/env python3
"""Planar average of a (difference) density cube along the third cell axis.

Usage:
    python cdd_planar.py co_cu_dAB.cube [--png cdd_z.png] [--zref 13.1]

Output (stdout), three columns:
    z (Ang)   dlambda(z) (e/Ang)   dQ(z) (e, accumulated from the bottom of the cell)

dlambda(z) is the density integrated over the plane spanned by the first two
grid vectors (line density); dQ(z) is the number of electrons gained (+) or
lost (-) below height z. The cube file is assumed to use Bohr units and the
third grid vector to be parallel to the z axis (slab normal along c).
"""
import argparse

import numpy as np

BOHR = 0.529177210903  # Angstrom


def read_cube(fname):
    with open(fname) as f:
        f.readline()
        f.readline()
        tok = f.readline().split()
        natom = abs(int(tok[0]))
        origin = np.array(tok[1:4], float)
        ngrid, vec = [], []
        for _ in range(3):
            tok = f.readline().split()
            ngrid.append(int(tok[0]))
            vec.append([float(x) for x in tok[1:4]])
        for _ in range(natom):
            f.readline()
        data = np.array(f.read().split(), float).reshape(ngrid)
    return origin, np.array(vec), data


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cube")
    p.add_argument("--png", default=None, help="save a plot (needs matplotlib)")
    p.add_argument("--zref", type=float, default=None,
                   help="print dQ at this height (Ang), e.g. midway between surface and adsorbate")
    a = p.parse_args()

    origin, vec, rho = read_cube(a.cube)
    dv = abs(np.linalg.det(vec))                       # Bohr^3
    area = np.linalg.norm(np.cross(vec[0], vec[1]))    # Bohr^2 per in-plane voxel
    dz = dv / area                                     # Bohr, spacing along the normal
    n3 = rho.shape[2]
    z = (origin[2] + np.arange(n3) * vec[2][2]) * BOHR  # Ang
    lam = rho.sum(axis=(0, 1)) * area                  # e/Bohr
    dq = np.cumsum(lam) * dz                           # e
    lam_ang = lam / BOHR                               # e/Ang

    print("# z(Ang)   dlambda(e/Ang)   dQ(e)")
    for zi, li, qi in zip(z, lam_ang, dq):
        print(f"{zi:10.4f} {li:15.6e} {qi:12.6f}")
    print(f"# total integral over the cell = {dq[-1]: .6f} e (should be ~0 for rho(AB)-rho(A)-rho(B))")
    if a.zref is not None:
        k = int(np.argmin(abs(z - a.zref)))
        print(f"# dQ at z = {z[k]:.3f} Ang : {dq[k]: .6f} e")

    if a.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=(6, 3.5))
        ax1.plot(z, lam_ang, color="tab:blue", label=r"$\Delta\lambda(z)$")
        ax1.axhline(0.0, color="gray", lw=0.6)
        ax1.set_xlabel("z (Angstrom)")
        ax1.set_ylabel(r"$\Delta\lambda$ (e/Angstrom)", color="tab:blue")
        ax2 = ax1.twinx()
        ax2.plot(z, dq, color="tab:red", ls="--", label=r"$\Delta Q(z)$")
        ax2.set_ylabel(r"$\Delta Q$ (e)", color="tab:red")
        fig.tight_layout()
        fig.savefig(a.png, dpi=200)


if __name__ == "__main__":
    main()
