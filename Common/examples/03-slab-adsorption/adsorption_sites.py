#!/usr/bin/env python3
"""Adsorption sites and initial adsorbate geometries on Cu(111).

Part A (pymatgen): AdsorbateSiteFinder finds the symmetry-distinct
  top / bridge / hollow sites of any slab read from POSCAR and builds one
  CO structure per site.
Part B (ASE): add_adsorbate() with the named sites of ase.build.fcc111
  ('ontop', 'bridge', 'fcc', 'hcp') - quick and transparent.

  python adsorption_sites.py                      # uses POSCAR_Cu111_slab
  python adsorption_sites.py --slab CONTCAR_clean_slab

Output: ads_pmg_<site>_<i>.vasp and ads_ase_<name>.vasp (POSCAR format).
The heights used here are only STARTING GUESSES (about 1.2 A for O in a
hollow, 1.85 A for C of CO on top); every structure must be relaxed.
"""
import argparse

from ase import Atoms
from ase.build import add_adsorbate, fcc111
from ase.constraints import FixAtoms
from ase.io import write


def part_pymatgen(slab_file):
    from pymatgen.analysis.adsorption import AdsorbateSiteFinder
    from pymatgen.core import Molecule, Structure

    slab = Structure.from_file(slab_file)
    # selective_dynamics=False: keep the F F F / T T T flags read from the
    # POSCAR (True would re-fix ALL non-surface atoms) and give the
    # adsorbate atoms T T T.
    asf = AdsorbateSiteFinder(slab, selective_dynamics=False)
    sites = asf.find_adsorption_sites(distance=1.85,
                                      positions=("ontop", "bridge", "hollow"))
    print("[A] pymatgen AdsorbateSiteFinder on %s (%d atoms)"
          % (slab_file, len(slab)))
    for kind in ("ontop", "bridge", "hollow"):
        print("    %-7s %d site(s)" % (kind, len(sites[kind])))
        for xyz in sites[kind]:
            print("            x,y,z = %7.3f %7.3f %7.3f" % tuple(xyz))
    # CO molecule, C at the origin and O above it (C-O = 1.15 A)
    co = Molecule(["C", "O"], [[0, 0, 0], [0, 0, 1.15]])
    n = 0
    for kind in ("ontop", "bridge", "hollow"):
        for i, xyz in enumerate(sites[kind]):
            st = asf.add_adsorbate(co, xyz, translate=False, reorient=True)
            fn = "ads_pmg_CO_%s_%d.vasp" % (kind, i)
            st.to(filename=fn, fmt="poscar")
            n += 1
    print("    wrote %d CO structures (ads_pmg_CO_*.vasp)" % n)


def part_ase(a, size, layers, nfix, vacuum):
    print("\n[B] ASE add_adsorbate on a freshly built Cu(111) slab")
    jobs = [("O_fcc", "O", 1.2, "fcc"),
            ("O_hcp", "O", 1.2, "hcp"),
            ("CO_top", "CO", 1.85, "ontop")]
    for name, ads, height, site in jobs:
        slab = fcc111("Cu", size=(size, size, layers), a=a,
                      vacuum=vacuum / 2.0)
        tags = slab.get_tags()
        slab.set_constraint(FixAtoms(mask=tags > layers - nfix))
        if ads == "CO":
            mol = Atoms("CO", positions=[(0, 0, 0), (0, 0, 1.15)])
            add_adsorbate(slab, mol, height, site, mol_index=0)
        else:
            add_adsorbate(slab, ads, height, site)
        # The cell is NOT changed: using the same cell as the clean slab
        # keeps E_ads = E(slab+ads) - E(slab) - E(ads) consistent.
        fn = "ads_ase_%s.vasp" % name
        with open(fn, "w", newline="\n") as f:
            write(f, slab, format="vasp", direct=True)
        print("    %-7s site=%-6s height=%.2f A -> %s (%d atoms)"
              % (name, site, height, fn, len(slab)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--slab", default="POSCAR_Cu111_slab")
    p.add_argument("--a", type=float, default=3.63)
    p.add_argument("--size", type=int, default=3)
    p.add_argument("--layers", type=int, default=4)
    p.add_argument("--fix", type=int, default=2)
    p.add_argument("--vacuum", type=float, default=15.0)
    p.add_argument("--skip-pymatgen", action="store_true")
    args = p.parse_args()

    if not args.skip_pymatgen:
        try:
            part_pymatgen(args.slab)
        except ImportError:
            print("pymatgen not installed - part A skipped")
    part_ase(args.a, args.size, args.layers, args.fix, args.vacuum)


if __name__ == "__main__":
    main()
