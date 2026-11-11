#!/usr/bin/env python3
"""Generate Cu(100) inputs for a surface-energy convergence study.

From the relaxed bulk (CONTCAR of the fcc primitive cell) it builds
  obulk/            oriented bulk: 1x1 (100) surface cell, 2 layers, periodic
  slab_N/           N-layer slab, central layer fixed (relaxed run)
  slab_N_unrel/     same slab, static (unrelaxed surface energy)
for every N in --layers, each with POSCAR and a KPOINTS whose in-plane mesh is
identical for bulk and slabs.

Usage: python make_slabs.py bulk/CONTCAR --layers 3 5 7 9 11 --vacuum 15 --kin 16 --kz 12
Requires ASE.
"""
import argparse
import os

import numpy as np
from ase.build import fcc100
from ase.constraints import FixAtoms
from ase.io import write


def lattice_constant(contcar):
    """fcc primitive cell: V = a^3 / 4."""
    lines = open(contcar).read().splitlines()
    s = float(lines[1].split()[0])
    cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)])
    vol = abs(np.linalg.det(cell)) * (s ** 3 if s > 0 else 1.0)
    if s < 0:
        vol = -s
    return (4.0 * vol) ** (1.0 / 3.0)


def kpoints(path, kin, kz, comment):
    with open(path, "w", newline="\n") as fh:
        fh.write(f"{comment}\n0\nGamma\n{kin} {kin} {kz}\n0 0 0\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("contcar", help="relaxed fcc primitive cell (bulk/CONTCAR)")
    p.add_argument("--layers", type=int, nargs="+", default=[3, 5, 7, 9, 11])
    p.add_argument("--vacuum", type=float, default=15.0, help="total vacuum (A)")
    p.add_argument("--kin", type=int, default=16, help="in-plane k mesh")
    p.add_argument("--kz", type=int, default=12, help="k mesh along z for the oriented bulk")
    p.add_argument("--element", default="Cu")
    a = p.parse_args()

    alat = lattice_constant(a.contcar)
    print(f"lattice constant from {a.contcar}: a = {alat:.5f} A")

    os.makedirs("obulk", exist_ok=True)
    ob = fcc100(a.element, size=(1, 1, 2), a=alat, periodic=True)
    write("obulk/POSCAR", ob, format="vasp", direct=True)
    kpoints("obulk/KPOINTS", a.kin, a.kz, "oriented bulk: same in-plane mesh as the slabs")
    print(f"obulk: {len(ob)} atoms, c = {ob.cell[2, 2]:.4f} A")

    for n in a.layers:
        slab = fcc100(a.element, size=(1, 1, n), a=alat, vacuum=a.vacuum / 2.0)
        z = slab.positions[:, 2]
        layers = np.unique(np.round(z, 4))
        mid = layers[len(layers) // 2]
        unrel = slab.copy()
        slab.set_constraint(FixAtoms(indices=np.where(np.abs(z - mid) < 1e-3)[0]))
        for d, at in ((f"slab_{n}", slab), (f"slab_{n}_unrel", unrel)):
            os.makedirs(d, exist_ok=True)
            write(f"{d}/POSCAR", at, format="vasp", direct=True)
            kpoints(f"{d}/KPOINTS", a.kin, 1, f"{n}-layer slab, 1 k-point along z")
        print(f"slab_{n}: {len(slab)} atoms, c = {slab.cell[2, 2]:.3f} A, "
              f"central layer z = {mid:.3f} A fixed")


if __name__ == "__main__":
    main()
