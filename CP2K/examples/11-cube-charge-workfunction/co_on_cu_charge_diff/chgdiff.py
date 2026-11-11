#!/usr/bin/env python3
"""Charge density difference  drho = rho(AB) - rho(A) - rho(B)  from CP2K cube files.

Writes
  diff.cube        3D difference density (same header/grid as the AB cube)
  diff_planar.dat  z [Angstrom], planar-averaged drho(z) [e/Angstrom],
                   charge displacement dQ(z) = integral_0^z drho(z') dz' [e]

Usage:  python chgdiff.py AB.cube A.cube B.cube
Needs numpy. The three cubes must come from runs with identical CELL, CUTOFF,
REL_CUTOFF and STRIDE, otherwise the grids differ and the script stops.
"""
import sys

import numpy as np

BOHR = 0.529177210903  # Angstrom


def read_cube(fname):
    with open(fname) as f:
        lines = f.readlines()
    natoms = abs(int(lines[2].split()[0]))
    grid = [lines[3 + i].split() for i in range(3)]
    n = [int(g[0]) for g in grid]
    vec = np.array([[float(x) for x in g[1:4]] for g in grid])   # bohr, voxel vectors
    header = lines[:6 + natoms]
    values = np.array(" ".join(lines[6 + natoms:]).split(), dtype=float)
    if values.size != n[0] * n[1] * n[2]:
        raise SystemExit("%s: unexpected number of values (multiple orbitals per cube?)" % fname)
    return header, n, vec, values.reshape(n)                      # x slowest, z fastest


def write_cube(fname, header, data):
    with open(fname, "w") as f:
        f.write(" charge density difference written by chgdiff.py\n")
        f.write(" rho(AB) - rho(A) - rho(B), e/bohr^3\n")
        f.writelines(header[2:])
        flat = data.reshape(data.shape[0] * data.shape[1], data.shape[2])
        for row in flat:                                          # 6 values per line per z-column
            for i in range(0, row.size, 6):
                f.write("".join("%13.5E" % v for v in row[i:i + 6]) + "\n")


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    hab, nab, vab, rab = read_cube(sys.argv[1])
    _, na, va, ra = read_cube(sys.argv[2])
    _, nb, vb, rb = read_cube(sys.argv[3])
    if not (nab == na == nb and np.allclose(vab, va) and np.allclose(vab, vb)):
        raise SystemExit("grids differ: use identical CELL/CUTOFF/REL_CUTOFF/STRIDE")
    drho = rab - ra - rb                                          # e/bohr^3
    write_cube("diff.cube", hab, drho)

    dvol = abs(np.linalg.det(vab))                                # bohr^3 per voxel
    print("integral of drho (should be ~0): %.4e e" % (drho.sum() * dvol))
    dz = np.linalg.norm(vab[2]) * BOHR                            # Angstrom
    area = np.linalg.norm(np.cross(vab[0] * nab[0], vab[1] * nab[1])) * BOHR ** 2
    # planar average: electrons per Angstrom along z
    line = drho.sum(axis=(0, 1)) * dvol / dz
    z = np.arange(nab[2]) * dz
    dq = np.cumsum(line) * dz
    np.savetxt("diff_planar.dat", np.column_stack([z, line, dq]),
               header="z[Angstrom]  drho_planar[e/Angstrom]  dQ(z)[e]   (cell area %.4f A^2)" % area,
               fmt="%12.6f")
    print("wrote diff.cube and diff_planar.dat")


if __name__ == "__main__":
    main()
