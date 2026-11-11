#!/usr/bin/env python3
"""Integrate one or more Gaussian cube files written by OpenMX.

Usage:
    python cube_integrate.py h2o.tden.cube [h2o.dden.cube ...]

For a *.tden.cube the integral is the number of valence electrons,
for *.dden.cube it should be close to zero, for *.sden.cube it is the
total spin moment (muB). Units in the cube file: Bohr and e/Bohr^3.
"""
import sys

import numpy as np


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
        atoms = [f.readline().split() for _ in range(natom)]
        data = np.array(f.read().split(), float)
    if data.size != ngrid[0] * ngrid[1] * ngrid[2]:
        raise ValueError(f"{fname}: expected {np.prod(ngrid)} values, got {data.size}")
    return origin, np.array(vec), atoms, data.reshape(ngrid)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    for fname in sys.argv[1:]:
        origin, vec, atoms, rho = read_cube(fname)
        dv = abs(np.linalg.det(vec))
        print(f"{fname}")
        print(f"  atoms = {len(atoms)}, grid = {rho.shape}, voxel volume = {dv:.6e} Bohr^3")
        print(f"  integral      = {rho.sum() * dv: .6f}")
        print(f"  min / max     = {rho.min(): .4e} / {rho.max(): .4e}")


if __name__ == "__main__":
    main()
