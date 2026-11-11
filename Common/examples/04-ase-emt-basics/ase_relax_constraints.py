#!/usr/bin/env python3
"""Relax a Cu adatom on a Cu(100) slab with FixAtoms; compare BFGS and FIRE.

Runs with the EMT potential (no DFT code required).  The same pattern
(build -> constrain -> attach calculator -> optimizer.run(fmax)) is used
with Vasp/CP2K/OpenMX calculators; only the calculator line changes.

Usage:
    python ase_relax_constraints.py
Outputs:
    relax_bfgs.traj, relax_fire.traj, relax_*.log, CONTCAR_bfgs (POSCAR format)
"""
from ase.build import add_adsorbate, fcc100
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import write
from ase.optimize import BFGS, FIRE


def build_slab():
    # 3x3 surface cell, 4 layers, 7.5 A vacuum on each side (15 A total)
    slab = fcc100("Cu", size=(3, 3, 4), vacuum=7.5)
    add_adsorbate(slab, "Cu", height=1.6, position="hollow")
    # Fix the two bottom layers (tag 3 and 4; adsorbate has tag 0)
    mask = [atom.tag > 2 for atom in slab]
    slab.set_constraint(FixAtoms(mask=mask))
    return slab


def run(opt_class, name, fmax=0.02):
    atoms = build_slab()
    atoms.calc = EMT()
    opt = opt_class(atoms, trajectory=f"relax_{name}.traj", logfile=f"relax_{name}.log")
    converged = opt.run(fmax=fmax, steps=500)
    e = atoms.get_potential_energy()
    nsteps = opt.get_number_of_steps()
    print(f"{name:5s} converged={converged}  steps={nsteps:3d}  E={e:.4f} eV")
    return atoms


if __name__ == "__main__":
    final = run(BFGS, "bfgs")
    run(FIRE, "fire")
    # Selective dynamics flags (F F F) are written for the fixed atoms
    write("CONTCAR_bfgs", final, format="vasp", direct=True)
    print("fixed atoms:", len(final.constraints[0].index), "of", len(final))
