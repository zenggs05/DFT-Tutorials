#!/usr/bin/env python3
"""Surface energy of a symmetric slab series (e.g. Cu(100), N = 3, 5, 7 ...).

    gamma(N) = [E_slab(N) - N * eps_bulk] / (2 A)

eps_bulk is taken (1) from the oriented bulk cell obulk/ (energy per atom,
same in-plane k mesh) and (2) from a Fiorentini-Methfessel linear fit
E_slab(N) = 2 A gamma + N eps over the thickest slabs. Both relaxed (slab_N)
and unrelaxed (slab_N_unrel) slabs are reported if present.

Usage: python surface_energy.py [--root .] [--nfit 3]
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

EV_A2_TO_J_M2 = 16.021766


def energy(outcar):
    e = None
    if not os.path.exists(outcar):
        return None
    with open(outcar) as fh:
        for line in fh:
            m = re.search(r"energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)", line)
            if m:
                e = float(m.group(1))
    return e


def natoms_and_area(poscar):
    lines = open(poscar).read().splitlines()
    s = float(lines[1].split()[0])
    cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * abs(s)
    counts_line = lines[6] if not lines[5].split()[0].isdigit() else lines[5]
    nat = sum(int(x) for x in counts_line.split())
    area = np.linalg.norm(np.cross(cell[0], cell[1]))
    return nat, area


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default=".")
    p.add_argument("--nfit", type=int, default=3, help="number of thickest slabs in the linear fit")
    a = p.parse_args()

    ob_e = energy(os.path.join(a.root, "obulk", "OUTCAR"))
    if ob_e is None:
        sys.exit("obulk/OUTCAR missing or unfinished")
    ob_n, _ = natoms_and_area(os.path.join(a.root, "obulk", "POSCAR"))
    eps_bulk = ob_e / ob_n
    print(f"oriented bulk: {ob_n} atoms, eps_bulk = {eps_bulk:.6f} eV/atom")

    for kind, pattern in (("relaxed", r"slab_(\d+)$"), ("unrelaxed", r"slab_(\d+)_unrel$")):
        rows = []
        for d in glob.glob(os.path.join(a.root, "slab_*")):
            m = re.search(pattern, os.path.normpath(d))
            if not m:
                continue
            e = energy(os.path.join(d, "OUTCAR"))
            if e is None:
                continue
            n, area = natoms_and_area(os.path.join(d, "POSCAR"))
            rows.append((n, e, area))
        if not rows:
            continue
        rows.sort()
        print(f"\n== {kind} slabs ==")
        print(f"{'N':>4s} {'E_slab (eV)':>14s} {'gamma (eV/A^2)':>15s} {'gamma (J/m^2)':>14s}")
        for n, e, area in rows:
            g = (e - n * eps_bulk) / (2.0 * area)
            print(f"{n:4d} {e:14.6f} {g:15.5f} {g * EV_A2_TO_J_M2:14.4f}")
        if len(rows) >= 2:
            sub = rows[-a.nfit:] if len(rows) >= a.nfit else rows
            ns = np.array([r[0] for r in sub], float)
            es = np.array([r[1] for r in sub])
            slope, icpt = np.polyfit(ns, es, 1)
            area = sub[-1][2]
            g_fm = icpt / (2.0 * area)
            print(f"linear fit over N = {', '.join(str(int(x)) for x in ns)}: "
                  f"eps_fit = {slope:.6f} eV/atom (bulk ref {eps_bulk:.6f}, "
                  f"diff {1000 * (slope - eps_bulk):+.2f} meV)")
            print(f"Fiorentini-Methfessel gamma = {g_fm:.5f} eV/A^2 = "
                  f"{g_fm * EV_A2_TO_J_M2:.4f} J/m^2")


if __name__ == "__main__":
    main()
