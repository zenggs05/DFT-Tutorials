#!/usr/bin/env python3
"""Build a graphene monolayer with the SAME in-plane lattice as a (relaxed) graphite
CONTCAR/POSCAR, with a large vacuum along c. Used as the reference for the
interlayer binding energy.

Usage: python make_graphene.py CONTCAR [--c 20.0] > POSCAR_graphene
"""
import argparse

import numpy as np


def read_lattice(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    scale = float(lines[1].split()[0])
    cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)])
    if scale < 0:                       # negative scale = target volume
        scale = (-scale / abs(np.linalg.det(cell))) ** (1.0 / 3.0)
    return cell * scale


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("structure", help="graphite CONTCAR or POSCAR")
    p.add_argument("--c", type=float, default=20.0, help="cell length along z (A)")
    a = p.parse_args()

    cell = read_lattice(a.structure)
    a1, a2 = cell[0], cell[1]
    print("graphene from %s, a = %.5f A, c = %.2f A" % (a.structure, np.linalg.norm(a1), a.c))
    print("1.0")
    for v in (a1, a2, (0.0, 0.0, a.c)):
        print("  %14.10f %14.10f %14.10f" % tuple(v))
    print("C")
    print("2")
    print("Direct")
    print("  0.0000000000  0.0000000000  0.5000000000")
    print("  0.3333333333  0.6666666667  0.5000000000")


if __name__ == "__main__":
    main()
