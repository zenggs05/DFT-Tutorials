#!/usr/bin/env python3
"""Dirac-point position of graphene from an OpenMX *.Band file.

The .Band file (written when Band.dispersion = on) contains, in Hartree:
    line 1 : number_of_bands  SpinP_switch  chemical_potential
    line 2 : 9 numbers, reciprocal lattice vectors
    line 3 : number of k-paths, followed by one line per path
    then, for every k point (and spin): "n k1 k2 k3" and a line with n eigenvalues.

The script finds the k point closest to K (default k1 k2 k3 = 0 1/3 1/3), takes
the two eigenvalues closest to the chemical potential there, and reports their
mean E_D relative to mu in eV.  For electron doping E_D - mu < 0.

Usage:
    python dirac_point.py graphene_q0.Band graphene_qm002_bg.Band graphene_qm002_on3.Band
    python dirac_point.py x.Band --k 0 0.3333333 0.3333333
"""
import argparse

import numpy as np

HARTREE2EV = 27.211386


def read_band(path):
    with open(path) as f:
        tok = f.read().split("\n")
    head = tok[0].split()
    nband, spin, mu = int(head[0]), int(head[1]), float(head[2])
    npath = int(tok[2].split()[0])
    i = 3 + npath
    kpts, eigs = [], []
    while i + 1 < len(tok):
        kl = tok[i].split()
        if len(kl) < 4:
            i += 1
            continue
        kpts.append([float(v) for v in kl[1:4]])
        eigs.append([float(v) for v in tok[i + 1].split()[:nband]])
        i += 2
    nspin = spin + 1 if spin in (0, 1) else 1
    return mu, np.array(kpts)[::nspin], np.array(eigs)[::nspin]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("band", nargs="+", help="System.Name.Band file(s)")
    ap.add_argument("--k", nargs=3, type=float, default=[0.0, 1 / 3, 1 / 3], help="K point (k1 k2 k3)")
    args = ap.parse_args()
    kt = np.array(args.k)
    print(f"{'file':40s} {'mu (eV)':>10s} {'E_D-mu (eV)':>12s} {'gap at K (eV)':>14s}")
    for path in args.band:
        mu, k, e = read_band(path)
        d = np.linalg.norm(k - kt, axis=1)
        ik = int(np.argmin(d))
        ek = np.sort((e[ik] - mu) * HARTREE2EV)
        two = ek[np.argsort(np.abs(ek))[:2]]
        print(f"{path:40s} {mu*HARTREE2EV:10.4f} {two.mean():12.4f} {abs(two[0]-two[1]):14.4f}"
              + ("" if d[ik] < 1e-3 else f"   (nearest k at distance {d[ik]:.3g})"))


if __name__ == "__main__":
    main()
