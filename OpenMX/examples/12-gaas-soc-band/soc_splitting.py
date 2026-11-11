#!/usr/bin/env python3
"""Read an OpenMX System.Name.Band file and print the levels near the chemical
potential at labelled k points, grouped by degeneracy, together with the
splitting between the two highest valence groups and the two lowest conduction
groups. Useful to quantify spin-orbit splittings (e.g. Gamma_8 - Gamma_7 in GaAs,
the valence-band splitting at K in monolayer MoS2) and band gaps.

Usage:
    python soc_splitting.py gaas_soc.Band [--at G L] [--nshow 6] [--tol 0.005]

Format of *.Band (as read by bandgnu13.c): header "nband nspin ChemP(Hartree)",
reciprocal vectors (9 numbers), number of paths, one line per path
"n k1 k2 k3 k1' k2' k3' label1 label2", then for every k point and spin:
"n kx ky kz" (fractional) followed by n eigenvalues in Hartree (not shifted).
"""
import argparse

import numpy as np

HARTREE = 27.211386245988  # eV


def read_band(fname):
    tok = open(fname).read().split()
    p = 0

    def nxt():
        nonlocal p
        p += 1
        return tok[p - 1]

    nband, mspin, chemp = int(nxt()), int(nxt()), float(nxt())
    rtv = np.array([float(nxt()) for _ in range(9)]).reshape(3, 3)
    npath = int(nxt())
    paths = []
    for _ in range(npath):
        n = int(nxt())
        k0 = np.array([float(nxt()) for _ in range(3)])
        k1 = np.array([float(nxt()) for _ in range(3)])
        paths.append((n, k0, k1, nxt(), nxt()))
    kpts, eig = [], [[] for _ in range(mspin + 1)]
    for n, _, _, _, _ in paths:
        for _ in range(n):
            for s in range(mspin + 1):
                n1 = int(nxt())
                k = np.array([float(nxt()) for _ in range(3)])
                e = np.array([float(nxt()) for _ in range(n1)])
                if s == 0:
                    kpts.append(k)
                eig[s].append(e)
    nmin = min(len(e) for es in eig for e in es)
    eig = np.array([[e[:nmin] for e in es] for es in eig])   # (spin, k, band)
    return chemp, rtv, paths, np.array(kpts), (eig - chemp) * HARTREE


def groups(levels, tol):
    out = []
    for e in np.sort(levels):
        if out and abs(e - out[-1][-1]) < tol:
            out[-1].append(e)
        else:
            out.append([e])
    return [(float(np.mean(g)), len(g)) for g in out]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("band")
    ap.add_argument("--at", nargs="*", default=None, help="k-point labels (default: all labels)")
    ap.add_argument("--nshow", type=int, default=6, help="number of level groups shown below/above E_F")
    ap.add_argument("--tol", type=float, default=0.005, help="degeneracy tolerance in eV")
    a = ap.parse_args()

    chemp, rtv, paths, kpts, e = read_band(a.band)
    print(f"file: {a.band}   spin blocks: {e.shape[0]}   k points: {e.shape[1]}   bands: {e.shape[2]}")
    print(f"chemical potential = {chemp:.6f} Hartree; energies below are relative to it (eV)")

    labels = {}
    for _, k0, k1, l0, l1 in paths:
        labels.setdefault(l0, k0)
        labels.setdefault(l1, k1)
    want = a.at if a.at else list(labels)

    for s in range(e.shape[0]):
        occ = e[s][e[s] < 0.0]
        emp = e[s][e[s] >= 0.0]
        if occ.size and emp.size:
            print(f"\n[spin block {s}] VBM = {occ.max():.4f} eV, CBM = {emp.min():.4f} eV, "
                  f"gap (along the path) = {emp.min() - occ.max():.4f} eV")
        for lab in want:
            if lab not in labels:
                print(f"label {lab} not found; available: {' '.join(labels)}")
                continue
            d = np.linalg.norm(kpts - labels[lab], axis=1)
            ik = int(np.argmin(d))
            if d[ik] > 1e-4:
                print(f"label {lab}: no k point within 1e-4 of {labels[lab]}")
                continue
            g = groups(e[s, ik], a.tol)
            val = [x for x in g if x[0] < 0.0][-a.nshow:]
            con = [x for x in g if x[0] >= 0.0][:a.nshow]
            print(f"\n  k = {lab} {np.round(labels[lab], 4)}   (energy eV, degeneracy)")
            for en, deg in reversed(con):
                print(f"     {en:10.4f}  x{deg}")
            print("     ---------- E_F ----------")
            for en, deg in reversed(val):
                print(f"     {en:10.4f}  x{deg}")
            if len(val) >= 2:
                print(f"     splitting of the two highest valence groups : {val[-1][0] - val[-2][0]:.4f} eV")
            if len(con) >= 2:
                print(f"     splitting of the two lowest conduction groups: {con[1][0] - con[0][0]:.4f} eV")


if __name__ == "__main__":
    main()
