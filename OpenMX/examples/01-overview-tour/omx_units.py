#!/usr/bin/env python3
"""Unit conversions frequently needed with OpenMX input/output.

DFT-Tutorials, OpenMX chapter 1 (01-overview.md).  Standard library only, Python >= 3.9.
Constants: CODATA 2018.

OpenMX conventions (OpenMX 4.0 manual): energies in *.out are in Hartree,
forces in Hartree/Bohr, scf.energycutoff in Ry, dipole moments in Debye,
scf.Electric.Field in GV/m, Dos.Erange and Hubbard U in eV.

Examples
--------
  python omx_units.py 1 Ha eV                 # 27.211386 eV
  python omx_units.py -- -8.033540953061 Ha eV
  python omx_units.py 3e-4 Ha/Bohr eV/Ang     # default MD.Opt.criterion
  python omx_units.py 0.5 V/Ang GV/m          # electric field
  python omx_units.py 1.85 D e*Ang            # dipole moment
  python omx_units.py --table                 # print a cheat sheet
  python omx_units.py --list                  # list supported units
"""
import argparse
import sys

HA_EV = 27.211386245988          # Hartree in eV
BOHR_ANG = 0.529177210903        # Bohr in Angstrom
HA_KCALMOL = 627.5094740631      # Hartree in kcal/mol
HA_KJMOL = 2625.4996394799       # Hartree in kJ/mol
HA_K = 315775.02480407           # Hartree / k_B in Kelvin
HA_CM1 = 219474.6313632          # Hartree in cm^-1
EBOHR_DEBYE = 2.541746473        # e*Bohr in Debye
EFIELD_AU_VM = 5.14220674763e11  # atomic unit of electric field in V/m

# every unit is expressed in a base unit of its dimension
UNITS = {
    "energy": {          # base: eV
        "eV": 1.0, "meV": 1e-3, "Ha": HA_EV, "Hartree": HA_EV, "mHa": 1e-3 * HA_EV,
        "Ry": HA_EV / 2, "Rydberg": HA_EV / 2, "kcal/mol": HA_EV / HA_KCALMOL,
        "kJ/mol": HA_EV / HA_KJMOL, "K": HA_EV / HA_K, "cm-1": HA_EV / HA_CM1,
    },
    "length": {          # base: Angstrom
        "Ang": 1.0, "A": 1.0, "nm": 10.0, "pm": 0.01, "Bohr": BOHR_ANG, "AU": BOHR_ANG,
    },
    "force": {           # base: eV/Ang
        "eV/Ang": 1.0, "eV/A": 1.0, "Ha/Bohr": HA_EV / BOHR_ANG, "Ry/Bohr": HA_EV / 2 / BOHR_ANG,
        "meV/Ang": 1e-3, "nN": 1.0 / 1.602176634,
    },
    "field": {           # base: V/m
        "V/m": 1.0, "GV/m": 1e9, "V/Ang": 1e10, "V/nm": 1e9, "MV/cm": 1e8, "au_field": EFIELD_AU_VM,
    },
    "dipole": {          # base: Debye
        "D": 1.0, "Debye": 1.0, "e*Bohr": EBOHR_DEBYE, "e*Ang": EBOHR_DEBYE / BOHR_ANG,
    },
}


def find(unit):
    for dim, table in UNITS.items():
        if unit in table:
            return dim, table[unit]
    raise SystemExit(f"unknown unit '{unit}'. Use --list to see the supported units.")


def convert(value, src, dst):
    d1, f1 = find(src)
    d2, f2 = find(dst)
    if d1 != d2:
        raise SystemExit(f"cannot convert {d1} ({src}) to {d2} ({dst})")
    return value * f1 / f2


def table():
    rows = [
        (1, "Ha", "eV"), (1, "Ha", "Ry"), (1, "Ry", "eV"), (1, "meV", "Ha"),
        (1, "Ha", "kcal/mol"), (1, "Ha", "kJ/mol"), (300, "K", "meV"),
        (1, "Bohr", "Ang"), (1, "Ang", "Bohr"),
        (1, "Ha/Bohr", "eV/Ang"), (3e-4, "Ha/Bohr", "eV/Ang"), (0.01, "eV/Ang", "Ha/Bohr"),
        (1, "V/Ang", "GV/m"), (1, "au_field", "GV/m"),
        (1, "e*Bohr", "D"), (1, "D", "e*Ang"),
    ]
    for v, a, b in rows:
        print(f"{v:>10g} {a:<9s} = {convert(v, a, b):<16.8g} {b}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("value", nargs="?", type=float, help="number to convert")
    ap.add_argument("src", nargs="?", help="source unit, e.g. Ha")
    ap.add_argument("dst", nargs="?", help="target unit, e.g. eV")
    ap.add_argument("--table", action="store_true", help="print a cheat sheet of common conversions")
    ap.add_argument("--list", action="store_true", help="list supported units")
    args = ap.parse_args()
    if args.list:
        for dim, tab in UNITS.items():
            print(f"{dim:7s}: {', '.join(tab)}")
        return
    if args.table or args.value is None:
        table()
        return
    if args.src is None or args.dst is None:
        ap.error("give: value source_unit target_unit")
    print(f"{args.value:g} {args.src} = {convert(args.value, args.src, args.dst):.10g} {args.dst}")


if __name__ == "__main__":
    sys.exit(main())
