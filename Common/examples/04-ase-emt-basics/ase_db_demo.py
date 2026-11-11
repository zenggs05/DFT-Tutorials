#!/usr/bin/env python3
"""Store and query results with ase.db (SQLite), then export to pandas/CSV.

EMT is used to compute fcc lattice constants of 6 metals so that the demo
runs anywhere.  Replace the calculator to archive real DFT results.

Usage:
    python ase_db_demo.py            # creates metals.db
    ase db metals.db                 # list rows on the command line
    ase db metals.db 'a0>3.9' -c id,metal,a0,B_GPa
    ase gui metals.db@metal=Cu       # open the selected row(s) in the GUI
"""
import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.db import connect
from ase.eos import EquationOfState
from ase.units import kJ

METALS = {"Al": 4.05, "Ni": 3.52, "Cu": 3.61, "Pd": 3.89, "Ag": 4.09, "Pt": 3.92}


def fit_metal(symbol, a_guess):
    atoms0 = bulk(symbol, "fcc", a=a_guess)
    vols, ens = [], []
    for s in np.linspace(0.97, 1.03, 7):
        at = atoms0.copy()
        at.set_cell(atoms0.cell * s, scale_atoms=True)
        at.calc = EMT()
        vols.append(at.get_volume())
        ens.append(at.get_potential_energy())
    v0, e0, B = EquationOfState(vols, ens).fit()
    a0 = (4 * v0) ** (1 / 3)
    best = bulk(symbol, "fcc", a=a0)
    best.calc = EMT()
    best.get_potential_energy()       # attach energy/forces so db stores them
    return best, a0, B / kJ * 1e24, vols, ens


def main():
    db = connect("metals.db")         # .db -> SQLite; .json also works
    for sym, a in METALS.items():
        if db.count(metal=sym) > 0:  # skip if already stored (idempotent re-runs)
            continue
        atoms, a0, b_gpa, vols, ens = fit_metal(sym, a)
        # key_value_pairs: searchable scalars; data: arbitrary arrays (not searchable).
        # Reserved names (id, energy, formula, natoms, ...) cannot be used as keys.
        db.write(atoms, metal=sym, a0=round(a0, 4), B_GPa=round(b_gpa, 1),
                 method="EMT", project="demo",
                 data={"volumes": vols, "energies": ens})

    print(f"{'id':>3} {'formula':>7} {'a0(A)':>8} {'B(GPa)':>8} {'E(eV)':>9}")
    for row in db.select(project="demo", sort="a0"):
        print(f"{row.id:>3} {row.formula:>7} {row.a0:8.4f} {row.B_GPa:8.1f} {row.energy:9.4f}")

    # Query with conditions; row.toatoms() gives back an Atoms object
    big = [r.metal for r in db.select("a0>3.9")]
    print("a0 > 3.9 A :", big)
    row = db.get(metal="Cu")
    print("Cu EOS points stored in row.data:", len(row.data["volumes"]))

    try:
        import pandas as pd
        df = pd.DataFrame([{"formula": r.formula, "a0": r.a0, "B_GPa": r.B_GPa,
                            "energy": r.energy} for r in db.select(project="demo")])
        df.to_csv("metals.csv", index=False)
        print("exported metals.csv with", len(df), "rows")
    except ImportError:
        print("pandas not installed: skip CSV export")


if __name__ == "__main__":
    main()
