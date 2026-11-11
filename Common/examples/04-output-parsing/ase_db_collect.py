#!/usr/bin/env python3
"""Walk a project tree, read every finished VASP / Gaussian run with ASE and
archive final structures + energies + metadata in an ase.db (SQLite).
Then export a CSV (pandas) for plotting and an extxyz for sharing.

Usage:
  python ase_db_collect.py ROOT --db project.db
  ase db project.db code=vasp -c id,formula,energy,max_force,dir
  python ase_db_collect.py ROOT --db project.db --csv summary.csv --extxyz all.xyz
"""
import argparse
import os

from ase.db import connect
from ase.io import read, write


def candidates(root):
    for d, _dirs, files in os.walk(root):
        if "OUTCAR" in files:
            yield os.path.join(d, "OUTCAR"), "vasp-out"
        for f in files:
            if f.endswith((".log", ".out")):
                p = os.path.join(d, f)
                with open(p, errors="ignore") as fh:
                    head = fh.read(5000)
                if "Entering Gaussian System" in head:
                    yield p, "gaussian-out"


def finished(path, fmt):
    with open(path, errors="ignore") as f:
        txt = f.read()
    if fmt == "vasp-out":
        return "General timing and accounting" in txt
    return "Normal termination" in txt[-3000:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--db", default="project.db")
    ap.add_argument("--csv")
    ap.add_argument("--extxyz")
    a = ap.parse_args()

    db = connect(a.db)
    n_new = 0
    for path, fmt in candidates(a.root):
        rel = os.path.relpath(os.path.dirname(path), a.root)
        if db.count(dir=rel, code=fmt.split("-")[0]) > 0:
            continue                          # already archived
        if not finished(path, fmt):
            print("skip (not finished):", path)
            continue
        atoms = read(path, format=fmt, index=-1)
        fmax = float(abs(atoms.get_forces()).max()) if "forces" in atoms.calc.results else -1.0
        # NOTE: names such as energy, fmax, formula, natoms are reserved by ase.db
        db.write(atoms, dir=rel, code=fmt.split("-")[0], max_force=round(fmax, 4))
        n_new += 1
    print(f"added {n_new} rows; database now has {db.count()} rows")

    if a.csv:
        import pandas as pd
        rows = [{"id": r.id, "dir": r.dir, "code": r.code, "formula": r.formula,
                 "energy_eV": r.get("energy"), "energy_per_atom": r.get("energy", 0) / r.natoms,
                 "max_force": r.max_force} for r in db.select()]
        pd.DataFrame(rows).to_csv(a.csv, index=False)
        print("written", a.csv)
    if a.extxyz:
        write(a.extxyz, [r.toatoms(add_additional_information=True) for r in db.select()],
              format="extxyz")
        print("written", a.extxyz)


if __name__ == "__main__":
    main()
