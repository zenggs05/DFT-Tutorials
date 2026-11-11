#!/usr/bin/env python3
"""Symmetry analysis with spglib (and pymatgen if installed).

  python check_symmetry.py Si.cif
  python check_symmetry.py CONTCAR --rattle 0.01     # simulate a relaxed,
                                                       # slightly noisy cell
What it shows
  1. space group found for a range of symprec values (tolerance in Angstrom)
  2. primitive and conventional standard cells (spglib)
  3. the same via pymatgen SpacegroupAnalyzer (optional)
  4. symmetry-inequivalent sites (Wyckoff letters) -> useful for defects
"""
import argparse

import numpy as np
import spglib
from ase.io import read


def spg_cell(atoms):
    return (atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("infile")
    p.add_argument("--rattle", type=float, default=0.0,
                   help="random displacement (A) added to every atom")
    p.add_argument("--seed", type=int, default=7)
    args = p.parse_args()

    atoms = read(args.infile)
    if args.rattle > 0:
        atoms.rattle(stdev=args.rattle, seed=args.seed)
        print("Atoms rattled with stdev = %.3f A" % args.rattle)
    cell = spg_cell(atoms)
    print("spglib version:", spglib.__version__)
    print("Input: %s, %d atoms" % (atoms.get_chemical_formula(), len(atoms)))

    print("\n[1] symprec scan")
    print("%10s  %-14s %4s" % ("symprec", "space group", "No."))
    for symprec in (1e-5, 1e-4, 1e-3, 1e-2, 1e-1):
        ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
        if ds is None:
            print("%10.0e  %s" % (symprec, "failed"))
            continue
        print("%10.0e  %-14s %4d" % (symprec, ds.international, ds.number))

    symprec = 1e-3 if args.rattle == 0 else max(10 * args.rattle, 1e-3)
    ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
    print("\n[2] spglib with symprec = %g: %s (No. %d), %d operations"
          % (symprec, ds.international, ds.number, len(ds.rotations)))
    prim = spglib.standardize_cell(cell, to_primitive=True, symprec=symprec)
    conv = spglib.standardize_cell(cell, to_primitive=False, symprec=symprec)
    for label, (lat, pos, num) in (("primitive", prim),
                                   ("conventional", conv)):
        lengths = np.linalg.norm(lat, axis=1)
        print("  %-12s %3d atoms, |a|,|b|,|c| = %s A"
              % (label, len(num), np.round(lengths, 4).tolist()))

    print("\n[3] Wyckoff positions / equivalent atoms")
    eq = ds.equivalent_atoms
    for rep in sorted(set(eq)):
        n = int(np.sum(eq == rep))
        print("  atom %3d (%s): Wyckoff %s, multiplicity %d"
              % (rep, atoms[rep].symbol, ds.wyckoffs[rep], n))

    try:
        from pymatgen.io.ase import AseAtomsAdaptor
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    except ImportError:
        print("\n[4] pymatgen not installed - skipped")
        return
    st = AseAtomsAdaptor.get_structure(atoms)
    sga = SpacegroupAnalyzer(st, symprec=symprec)
    ps = sga.get_primitive_standard_structure()
    cs = sga.get_conventional_standard_structure()
    print("\n[4] pymatgen SpacegroupAnalyzer: %s (%d)"
          % (sga.get_space_group_symbol(), sga.get_space_group_number()))
    print("  primitive standard   : %d atoms, abc = %s"
          % (len(ps), np.round(ps.lattice.abc, 4).tolist()))
    print("  conventional standard: %d atoms, abc = %s"
          % (len(cs), np.round(cs.lattice.abc, 4).tolist()))
    sym = sga.get_symmetrized_structure()
    print("  inequivalent site groups:",
          [len(g) for g in sym.equivalent_indices])


if __name__ == "__main__":
    main()
