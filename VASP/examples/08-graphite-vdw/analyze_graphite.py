#!/usr/bin/env python3
"""Summarise the graphite vdW comparison produced by run_compare.sh.

For each method directory <m>/ it reads
    <m>/relax/CONTCAR        -> a, c
    <m>/graphite/OUTCAR      -> E(graphite, 4 atoms)
    <m>/graphene/OUTCAR      -> E(graphene, 2 atoms)
and prints the interlayer binding energy
    E_b = E(graphene)/2 - E(graphite)/4          (meV/atom, > 0 = bound)
also per in-plane area (meV/A^2).

Usage: python analyze_graphite.py PBE D3BJ optB86b
"""
import os
import re
import sys

import numpy as np


def lattice(path):
    lines = open(path).read().splitlines()
    s = float(lines[1].split()[0])
    cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)])
    if s < 0:
        s = (-s / abs(np.linalg.det(cell))) ** (1 / 3)
    return cell * s


def energy(outcar):
    e = None
    with open(outcar) as fh:
        for line in fh:
            m = re.search(r"energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)", line)
            if m:
                e = float(m.group(1))
    return e


def main(methods):
    print(f"{'method':10s} {'a (A)':>8s} {'c (A)':>8s} {'E_b (meV/atom)':>15s} {'E_b (meV/A^2)':>14s}")
    for m in methods:
        try:
            cell = lattice(os.path.join(m, "relax", "CONTCAR"))
            eg = energy(os.path.join(m, "graphite", "OUTCAR"))
            e1 = energy(os.path.join(m, "graphene", "OUTCAR"))
        except OSError as err:
            print(f"{m:10s} missing files ({err.filename})")
            continue
        a = np.linalg.norm(cell[0])
        c = np.linalg.norm(cell[2])
        if eg is None or e1 is None:
            print(f"{m:10s} {a:8.4f} {c:8.4f}   (static runs not finished)")
            continue
        eb_atom = (e1 / 2.0 - eg / 4.0) * 1000.0
        area_per_atom = abs(np.cross(cell[0], cell[1])[2]) / 2.0   # 2 atoms per layer cell
        print(f"{m:10s} {a:8.4f} {c:8.4f} {eb_atom:15.1f} {eb_atom / area_per_atom:14.2f}")
    print("\nreference: exp. a = 2.464 A, c = 6.711 A; "
          "RPA binding ~48 meV/atom (Lebegue et al., PRL 105, 196401, 2010)")


if __name__ == "__main__":
    main(sys.argv[1:] or ["PBE", "D3BJ", "optB86b"])
