#!/usr/bin/env python3
"""Single-point energy and forces of an isolated H2O molecule with CP2K via ASE.

The ASE CP2K calculator talks to a running "cp2k.psmp -s" (shell mode)
process through pipes, so a sequence of energy/force calls (optimisation,
NEB, MD in ASE) has no start-up cost.  Anything the simple keywords do
not cover goes into the 'inp' template (a normal CP2K input fragment).

Requirements:
  export CP2K_DATA_DIR=/path/to/cp2k/data         # BASIS_MOLOPT, GTH_POTENTIALS
  export ASE_CP2K_COMMAND="mpirun -np 4 cp2k.psmp -s"
Usage:
  python ase_cp2k_water_sp.py
  python ase_cp2k_water_sp.py --relax             # also relax with ASE BFGS
"""
import argparse

from ase.build import molecule
from ase.calculators.cp2k import CP2K
from ase.optimize import BFGS

# Extra input merged into the generated one.  Works with CP2K 2023.x - 2026.x.
INP = """
&FORCE_EVAL
  &DFT
    &SCF
      EPS_SCF 1.0E-7
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF
        MAX_SCF 10
        EPS_SCF 1.0E-7
      &END OUTER_SCF
    &END SCF
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &MGRID
      REL_CUTOFF 50
    &END MGRID
  &END DFT
&END FORCE_EVAL
"""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--relax", action="store_true")
    args = p.parse_args()

    h2o = molecule("H2O")
    h2o.center(vacuum=5.0)        # ~10 A box
    h2o.pbc = False               # -> ASE adds POISSON PERIODIC NONE + PSOLVER MT

    calc = CP2K(
        label="h2o",              # PROJECT name -> h2o-*.restart etc.
        xc="PBE",
        basis_set="DZVP-MOLOPT-SR-GTH",
        basis_set_file="BASIS_MOLOPT",
        pseudo_potential="GTH-PBE",
        potential_file="GTH_POTENTIALS",
        cutoff=400 * 13.605693,   # ASE expects eV: 400 Ry
        max_scf=30,               # inner OT loop; OUTER_SCF in INP handles the rest
        stress_tensor=False,      # no stress for an isolated molecule
        inp=INP,
    )
    h2o.calc = calc
    with calc:                    # closes the cp2k shell at the end
        e = h2o.get_potential_energy()
        f = h2o.get_forces()
        print(f"E = {e:.6f} eV")
        print("forces (eV/A):")
        print(f.round(4))
        if args.relax:
            BFGS(h2o, trajectory="h2o_cp2k.traj").run(fmax=0.02)
            print(f"relaxed E = {h2o.get_potential_energy():.6f} eV")
            print("O-H distance (A):", round(h2o.get_distance(0, 1), 4))


if __name__ == "__main__":
    main()
