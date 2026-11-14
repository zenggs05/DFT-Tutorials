#!/usr/bin/env python3
"""Query several crystal databases at once through OPTIMADE (pymatgen).

OPTIMADE is a common REST API implemented by Materials Project, OQMD,
Alexandria, COD, NOMAD and many more providers, so one filter works for all
of them.  No API key is needed for public providers.

Usage:
  python pmg_optimade_query.py --elements Ti O --nelements 2
  python pmg_optimade_query.py --providers mp oqmd --formula TiO2 --max 5
  python pmg_optimade_query.py --filter 'elements HAS ALL "Ga","N" AND nsites<=8'

Structures are written as CIF files into optimade_cif/<provider>/<id>.cif.
"""
import argparse
import os
import re

from pymatgen.ext.optimade import OptimadeRester


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--providers", nargs="+", default=["mp", "oqmd"],
                    help="pymatgen aliases (mp, oqmd, alexandria, cod, nmd, ...)")
    ap.add_argument("--elements", nargs="+")
    ap.add_argument("--nelements", type=int)
    ap.add_argument("--formula", help="reduced formula in Hill order, e.g. O2Ti for TiO2")
    ap.add_argument("--filter", help="raw OPTIMADE filter string (overrides the options above)")
    ap.add_argument("--max", type=int, default=10, help="max structures written per provider")
    ap.add_argument("--outdir", default="optimade_cif")
    a = ap.parse_args()
    if not (a.filter or a.elements or a.nelements or a.formula):
        ap.error("give --elements/--nelements/--formula or --filter "
                 "(an empty query would try to download whole databases)")

    with OptimadeRester(a.providers) as opt:
        if a.filter:
            results = opt.get_structures_with_filter(a.filter)
        else:
            hill = None
            if a.formula:
                # OPTIMADE uses Hill order (C, H, then alphabetical): TiO2 -> O2Ti
                from pymatgen.core import Composition
                hill = Composition(a.formula).reduced_composition.hill_formula.replace(" ", "")
            results = opt.get_structures(elements=a.elements, nelements=a.nelements,
                                         chemical_formula_hill=hill)

    for provider, structs in results.items():
        print(f"[{provider}] {len(structs)} structures")
        d = os.path.join(a.outdir, provider)
        os.makedirs(d, exist_ok=True)
        for i, (sid, s) in enumerate(sorted(structs.items())):
            if i >= a.max:
                break
            safe = re.sub(r"[^A-Za-z0-9_.-]", "_", str(sid))
            if hasattr(s, "lattice"):           # periodic -> CIF
                s.to(filename=os.path.join(d, safe + ".cif"))
                print(f"  {sid:<24} {s.composition.reduced_formula:<10} "
                      f"nsites={len(s):<4} V={s.volume:9.2f} A^3")
            else:                               # molecules -> xyz
                s.to(filename=os.path.join(d, safe + ".xyz"))
    print("CIF files written under", a.outdir)


if __name__ == "__main__":
    main()
