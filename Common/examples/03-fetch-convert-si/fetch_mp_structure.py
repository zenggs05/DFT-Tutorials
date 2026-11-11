#!/usr/bin/env python3
"""Download a crystal structure from the Materials Project (new API, 2025+).

The legacy REST API (materialsproject.org/rest/v2) was retired in 2025-09;
old 16-character keys no longer work. Get a new key from
https://next-gen.materialsproject.org/dashboard (it is 32 characters long).

Two client back-ends are supported:
  * mp-api    : pip install mp-api      ; key from $MP_API_KEY
  * pymatgen  : pymatgen.ext.matproj    ; key from $PMG_MAPI_KEY
                (or PMG_MAPI_KEY in ~/.config/.pmgrc.yaml)

Examples
  export MP_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  python fetch_mp_structure.py mp-149                 # Si, primitive cell
  python fetch_mp_structure.py mp-149 --conventional  # Si, 8-atom cubic cell
  python fetch_mp_structure.py --formula TiO2         # list candidates first
  python fetch_mp_structure.py mp-2657 --backend pymatgen

Outputs: <prefix>.cif and POSCAR_<prefix> (prefix defaults to the mp-id).
"""
import argparse
import os
import sys

FIELDS = ["material_id", "formula_pretty", "symmetry", "nsites",
          "energy_above_hull", "band_gap", "is_stable"]


def get_key(backend):
    """Return the API key or exit with a helpful message."""
    names = ["MP_API_KEY", "PMG_MAPI_KEY"]
    if backend == "pymatgen":
        names = ["PMG_MAPI_KEY", "MP_API_KEY"]
    for name in names:
        key = os.environ.get(name, "").strip()
        if key:
            break
    else:
        key = ""
    if not key and backend == "pymatgen":
        # pymatgen can also read the key from ~/.config/.pmgrc.yaml
        try:
            from pymatgen.core import SETTINGS
            key = SETTINGS.get("PMG_MAPI_KEY", "") or ""
        except Exception:
            key = ""
    if not key:
        sys.exit("No API key found. Set MP_API_KEY (mp-api) or PMG_MAPI_KEY "
                 "(pymatgen). Keys are shown on the MP dashboard.")
    if len(key) != 32:
        print("WARNING: key length is %d, new-API keys have 32 characters. "
              "Legacy keys stopped working in 2025-09." % len(key))
    return key


def field(doc, name):
    """Read a field from an mp-api document object or a plain dict."""
    if isinstance(doc, dict):
        return doc.get(name)
    return getattr(doc, name, None)


def open_rester(backend, key):
    if backend == "mp-api":
        try:
            from mp_api.client import MPRester
        except ImportError:
            sys.exit("mp-api is not installed: pip install mp-api "
                     "(or use --backend pymatgen)")
        return MPRester(key)
    from pymatgen.ext.matproj import MPRester
    return MPRester(key)


def list_formula(mpr, backend, formula):
    """Print all MP entries with this reduced formula, most stable first."""
    if backend == "mp-api":
        docs = mpr.materials.summary.search(formula=formula, fields=FIELDS)
    else:
        docs = mpr.materials.summary.search(formula=formula, _fields=FIELDS)
    rows = []
    for d in docs:
        sym = field(d, "symmetry")
        if isinstance(sym, dict):
            sg = sym.get("symbol")
        else:
            sg = getattr(sym, "symbol", "?")
        ehull = field(d, "energy_above_hull")
        rows.append((ehull if ehull is not None else 99.0,
                     str(field(d, "material_id")), str(sg),
                     field(d, "nsites"), field(d, "band_gap")))
    rows.sort()
    print("%-14s %-12s %6s %12s %9s" % ("material_id", "spacegroup",
                                        "nsites", "E_hull(eV)", "gap(eV)"))
    for ehull, mpid, sg, nsites, gap in rows:
        print("%-14s %-12s %6s %12.3f %9.3f" % (mpid, sg, nsites, ehull,
                                                 gap or 0.0))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mpid", nargs="?", default=None,
                   help="Materials Project id, e.g. mp-149")
    p.add_argument("--formula", help="list all entries of a formula and exit")
    p.add_argument("--backend", choices=["mp-api", "pymatgen"],
                   default="mp-api")
    p.add_argument("--conventional", action="store_true",
                   help="return the conventional standard cell")
    p.add_argument("--prefix", help="output file prefix (default: mpid)")
    args = p.parse_args()

    if not args.mpid and not args.formula:
        p.error("give an mp-id or --formula")

    key = get_key(args.backend)
    with open_rester(args.backend, key) as mpr:
        if args.formula:
            list_formula(mpr, args.backend, args.formula)
            return
        # same method name and keyword in both clients
        st = mpr.get_structure_by_material_id(
            args.mpid, conventional_unit_cell=args.conventional)

    prefix = args.prefix or args.mpid
    st.to(filename=prefix + ".cif")
    st.to(filename="POSCAR_" + prefix, fmt="poscar")
    print("Formula      :", st.composition.reduced_formula)
    print("Sites        :", len(st))
    print("Lattice abc  : %.4f %.4f %.4f" % st.lattice.abc)
    print("Angles       : %.2f %.2f %.2f" % st.lattice.angles)
    print("Space group  : %s (%d)" % st.get_space_group_info(symprec=0.01))
    print("Written      : %s.cif  POSCAR_%s" % (prefix, prefix))
    print("Note: MP structures are relaxed with PBE (or r2SCAN); PBE "
          "lattice constants are typically ~1% larger than experiment.")


if __name__ == "__main__":
    main()
