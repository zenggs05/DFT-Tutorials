#!/usr/bin/env python3
"""Build a p(3x3) Cu(111) slab with 4 layers, 15 A vacuum and the two
bottom layers fixed (Selective dynamics F F F) -> POSCAR_Cu111_slab.

  python build_cu111_slab.py                 # defaults below
  python build_cu111_slab.py --a 3.635 --layers 5 --fix 2 --vacuum 20

The lattice constant should be YOUR bulk value optimized with the same
functional/ENCUT/POTCAR as the slab calculation (PBE gives about 3.63 A
for Cu; experiment at room temperature is about 3.615 A). Using a value
from another functional strains the slab.
"""
import argparse

import numpy as np
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io import write


def build(a, size, layers, vacuum, nfix):
    # ASE adds `vacuum` on EACH side of the slab, so the total vacuum gap
    # is 2*vacuum. We therefore pass half of the requested total.
    slab = fcc111("Cu", size=(size, size, layers), a=a,
                  vacuum=vacuum / 2.0, orthogonal=False)
    # ASE tags layers 1 (top) ... layers (bottom)
    tags = slab.get_tags()
    fixed = tags > layers - nfix
    slab.set_constraint(FixAtoms(mask=fixed))
    return slab, fixed


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a", type=float, default=3.63,
                   help="bulk fcc lattice constant in A (default 3.63)")
    p.add_argument("--size", type=int, default=3, help="in-plane repeat")
    p.add_argument("--layers", type=int, default=4)
    p.add_argument("--vacuum", type=float, default=15.0,
                   help="TOTAL vacuum thickness in A")
    p.add_argument("--fix", type=int, default=2,
                   help="number of bottom layers to fix")
    p.add_argument("-o", "--output", default="POSCAR_Cu111_slab")
    args = p.parse_args()

    slab, fixed = build(args.a, args.size, args.layers, args.vacuum,
                        args.fix)
    z = slab.positions[:, 2]
    thickness = z.max() - z.min()
    gap = slab.cell[2, 2] - thickness
    d111 = args.a / np.sqrt(3.0)
    with open(args.output, "w", newline="\n") as f:
        write(f, slab, format="vasp", direct=True)

    print("Cu(111) p(%dx%d), %d layers, %d atoms (%d fixed)"
          % (args.size, args.size, args.layers, len(slab), fixed.sum()))
    print("a(bulk) = %.4f A, d(111) = %.4f A" % (args.a, d111))
    print("in-plane |a| = |b| = %.4f A, gamma = %.1f deg"
          % (slab.cell.lengths()[0], slab.cell.angles()[2]))
    print("c = %.3f A, slab thickness = %.3f A, vacuum gap = %.3f A"
          % (slab.cell[2, 2], thickness, gap))
    print("written:", args.output)


if __name__ == "__main__":
    main()
