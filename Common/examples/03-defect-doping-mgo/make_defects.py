#!/usr/bin/env python3
"""Point defects and doping in a MgO supercell (ASE + spglib, pymatgen
optional for de-duplicating dopant pairs).

  python make_defects.py                 # 2x2x2 conventional cell, 64 atoms
  python make_defects.py --n 3 --a 4.25  # 216 atoms, PBE-like a

Writes (POSCAR format, files *.vasp):
  MgO_perfect.vasp     perfect supercell (reference for formation energies)
  MgO_VO.vasp          oxygen vacancy V_O
  MgO_AlMg.vasp        Al on a Mg site (substitution, n-type dopant)
  MgO_Li_i.vasp        Li on a tetrahedral interstitial site
  MgO_2AlMg_<k>.vasp   symmetry-distinct Al-Al pair configurations

Default a = 4.212 A is close to the room-temperature experimental value
of MgO; for production use the lattice constant optimized with your own
functional so that the supercell is not strained.
"""
import argparse
import itertools

import numpy as np
import spglib
from ase.build import bulk
from ase.geometry import get_distances


def lf_write(fn, atoms):
    from ase.io import write
    with open(fn, "w", newline="\n") as f:
        write(f, atoms, format="vasp", direct=True, sort=True)


def inequivalent(atoms, symprec=1e-3):
    cell = (atoms.cell[:], atoms.get_scaled_positions(), atoms.numbers)
    ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
    reps = sorted(set(ds.equivalent_atoms))
    return ds, {r: int(np.sum(ds.equivalent_atoms == r)) for r in reps}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a", type=float, default=4.212)
    p.add_argument("--n", type=int, default=2,
                   help="repeat of the 8-atom cubic cell (n x n x n)")
    p.add_argument("--max-pairs", type=int, default=6)
    args = p.parse_args()

    unit = bulk("MgO", "rocksalt", a=args.a, cubic=True)      # 8 atoms
    sc = unit.repeat((args.n, args.n, args.n))
    L = sc.cell.lengths()[0]
    print("Perfect supercell: %s, %d atoms, L = %.3f A"
          % (sc.get_chemical_formula(), len(sc), L))
    print("Defect-defect periodic image distance = %.3f A "
          "(aim for >= 10 A in production)" % L)

    ds, groups = inequivalent(sc)
    print("Space group %s (%d); inequivalent sites:"
          % (ds.international, ds.number))
    for rep, mult in groups.items():
        print("   atom %3d  %-2s  Wyckoff %s  multiplicity %d"
              % (rep, sc[rep].symbol, ds.wyckoffs[rep], mult))
    lf_write("MgO_perfect.vasp", sc)

    # 1) vacancy: remove the first O (all O are equivalent here)
    i_o = [a.index for a in sc if a.symbol == "O"][0]
    vac = sc.copy()
    del vac[i_o]
    lf_write("MgO_VO.vasp", vac)
    print("V_O   : removed O #%d -> %s" % (i_o, vac.get_chemical_formula()))

    # 2) substitution: Al on the first Mg site
    mg_idx = [a.index for a in sc if a.symbol == "Mg"]
    sub = sc.copy()
    sub[mg_idx[0]].symbol = "Al"
    lf_write("MgO_AlMg.vasp", sub)
    print("Al_Mg : Mg #%d -> Al -> %s" % (mg_idx[0],
                                          sub.get_chemical_formula()))

    # 3) interstitial: Li at the tetrahedral site (1/4,1/4,1/4) of the
    #    first cubic cell (fractional coordinates of the supercell / n)
    inter = sc.copy()
    pos = np.array([0.25, 0.25, 0.25]) / args.n @ sc.cell[:]
    inter.append("Li")
    inter.positions[-1] = pos
    d = get_distances(pos, inter.positions[:-1], cell=inter.cell,
                      pbc=True)[1].min()
    lf_write("MgO_Li_i.vasp", inter)
    print("Li_i  : tetrahedral site, nearest neighbour %.3f A" % d)

    # 4) two Al_Mg dopants: list pairs by Al-Al distance and remove
    #    symmetry-equivalent configurations
    first = mg_idx[0]
    dists = get_distances(sc.positions[first], sc.positions[mg_idx],
                          cell=sc.cell, pbc=True)[1][0]
    shells = {}
    for j, dj in zip(mg_idx, dists):
        if j == first:
            continue
        shells.setdefault(round(float(dj), 3), []).append(j)
    print("Al-Al pair candidates (fix Al#1 on Mg #%d):" % first)
    try:
        from pymatgen.analysis.structure_matcher import StructureMatcher
        from pymatgen.io.ase import AseAtomsAdaptor
        matcher = StructureMatcher()
    except ImportError:
        matcher = None
    unique = []
    for k, (dist, js) in enumerate(sorted(shells.items())):
        if k >= args.max_pairs:
            break
        cands = []
        for j in js:
            s2 = sc.copy()
            s2[first].symbol = "Al"
            s2[j].symbol = "Al"
            cands.append(s2)
        if matcher is not None:
            structs = [AseAtomsAdaptor.get_structure(c) for c in cands]
            groups = matcher.group_structures(structs)
            reps = [AseAtomsAdaptor.get_atoms(g[0]) for g in groups]
        else:
            reps = cands[:1]
        for r in reps:
            unique.append((dist, r))
            fn = "MgO_2AlMg_%d.vasp" % len(unique)
            lf_write(fn, r)
            print("   d(Al-Al) = %6.3f A  (%2d equivalent j)  -> %s"
                  % (dist, len(js), fn))
    if matcher is None:
        print("   (pymatgen missing: one structure per distance shell)")


if __name__ == "__main__":
    main()
