#!/usr/bin/env python3
"""Rutile TiO2(110) slab with pymatgen SlabGenerator.

Steps
  1. build bulk rutile from its space group (P4_2/mnm, No. 136):
     Ti at 2a (0,0,0), O at 4f (u,u,0)
  2. enumerate all (110) terminations and report stoichiometry,
     polarity and symmetry of each candidate
  3. keep a stoichiometric, non-polar, symmetric slab, optionally
     repeat it in-plane, fix the bottom atoms and write POSCAR

  python build_slab_TiO2_110.py                       # 1x1, ~12 A thick
  python build_slab_TiO2_110.py --thick 16 --repeat 2 1 --fix-trilayers 2

Default lattice: a = 4.594 A, c = 2.959 A, u = 0.3048 (room-temperature
neutron diffraction, Burdett et al., JACS 109, 3639 (1987)). For
production runs replace them by your own optimized bulk values.
"""
import argparse

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp import Poscar


def rutile(a, c, u):
    return Structure.from_spacegroup("P4_2/mnm", Lattice.tetragonal(a, c),
                                     ["Ti", "O"], [[0, 0, 0], [u, u, 0]])


def describe(slab):
    comp = slab.composition.reduced_formula
    return ("%-6s %3d atoms  symmetric=%-5s polar=%-5s shift=%.4f"
            % (comp, len(slab), slab.is_symmetric(), slab.is_polar(),
               slab.shift))


def trilayer_index(slab, cation="Ti", tol=0.3):
    """Assign every atom to the nearest cation (Ti2O2) plane.

    Rutile (110) is a stack of O / Ti2O2 / O "trilayers"; counting
    trilayers is the usual way to quote the slab thickness.
    """
    zc = sorted(s.coords[2] for s in slab if s.specie.symbol == cation)
    planes = []
    for z in zc:
        if not planes or z - planes[-1] > tol:
            planes.append(z)
    idx = [int(np.argmin([abs(s.coords[2] - zp) for zp in planes]))
           for s in slab]
    return idx, len(planes)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a", type=float, default=4.594)
    p.add_argument("--c", type=float, default=2.959)
    p.add_argument("--u", type=float, default=0.3048)
    p.add_argument("--thick", type=float, default=12.0,
                   help="minimum slab thickness in A")
    p.add_argument("--vacuum", type=float, default=15.0,
                   help="minimum vacuum thickness in A")
    p.add_argument("--repeat", type=int, nargs=2, default=[1, 1],
                   help="in-plane supercell, e.g. 2 1")
    p.add_argument("--fix-trilayers", type=int, default=2,
                   help="number of bottom O-Ti2O2-O trilayers to fix")
    p.add_argument("-o", "--output", default="POSCAR_TiO2_110")
    args = p.parse_args()

    bulk = rutile(args.a, args.c, args.u)
    print("Bulk:", bulk.composition.reduced_formula,
          bulk.get_space_group_info(), "%d atoms" % len(bulk))

    gen = SlabGenerator(bulk, (1, 1, 0), min_slab_size=args.thick,
                        min_vacuum_size=args.vacuum, center_slab=True,
                        lll_reduce=True)
    # Every distinct "shift" is a different termination. Optional filter:
    #   gen.get_slabs(bonds={("Ti", "O"): 2.2}, max_broken_bonds=2)
    # discards cuts that break too many Ti-O bonds.
    slabs = gen.get_slabs()
    print("\nCandidate (110) terminations:")
    for i, s in enumerate(slabs):
        print("  [%d] %s" % (i, describe(s)))

    good = [s for s in slabs
            if s.composition.reduced_formula == "TiO2"
            and s.is_symmetric() and not s.is_polar()]
    if not good:
        raise SystemExit("No stoichiometric non-polar symmetric slab found;"
                         " try another --thick value.")
    slab = good[0].get_orthogonal_c_slab()
    slab.make_supercell([args.repeat[0], args.repeat[1], 1])

    slab = slab.get_sorted_structure()  # group elements for POSCAR
    layer, nlayers = trilayer_index(slab)
    sd = [[False] * 3 if li < args.fix_trilayers else [True] * 3
          for li in layer]
    Poscar(slab, selective_dynamics=sd,
           comment="TiO2(110) slab %dx%d" % tuple(args.repeat)).write_file(
        args.output)

    z = np.array([s.coords[2] for s in slab])
    a, b, c = slab.lattice.abc
    print("\nChosen slab: %s" % describe(slab))
    print("trilayers = %d, thickness = %.2f A" % (nlayers, z.max() - z.min()))
    print("surface cell %.3f x %.3f A, c = %.3f A, gamma = %.1f deg"
          % (a, b, c, slab.lattice.gamma))
    print("vacuum gap = %.2f A, fixed atoms = %d / %d (bottom %d trilayers)"
          % (c - (z.max() - z.min()), sum(1 for f in sd if not f[0]),
             len(slab), args.fix_trilayers))
    print("written:", args.output)


if __name__ == "__main__":
    main()
