#!/usr/bin/env python3
"""Optimise an isolated H2O molecule with ASE driving CP2K through its shell mode.

Chapter 17 of the CP2K tutorial (ASE CP2K calculator).

How the CP2K process is started (first match wins, see ASE's cp2k.py):
  1. the --command argument of this script
  2. the environment variable ASE_CP2K_COMMAND
  3. the ASE default: "cp2k.psmp -s" (ASE >= 3.29) or "cp2k_shell" (ASE <= 3.28)

ASE <= 3.28 asserts that the command string contains "cp2k_shell". With a modern
CMake build of CP2K use the wrapper ase_cp2k_shell_wrapper.sh from this directory
(its file name contains "cp2k_shell") or upgrade ASE.

Everything that is not covered by the calculator keywords (UZH basis sets, OT
settings, ...) is passed in through the input template INP below. Keywords set by
the calculator must NOT appear in the template as well, otherwise ASE raises
"Keyword ... already present".

Requirements: python >= 3.9 (ASE >= 3.27 needs >= 3.10), ase, numpy, a CP2K binary,
CP2K_DATA_DIR pointing to the CP2K data directory.
"""
import argparse
import os
import sys

import ase
from ase.build import molecule
from ase.calculators.cp2k import CP2K
from ase.io import write
from ase.optimize import BFGS
from ase.units import Rydberg

# Input template: only what the calculator cannot express itself.
# Do not put PROJECT, PRINT_LEVEL, CUTOFF, MAX_SCF, XC, CELL or COORD here.
INP = """
&FORCE_EVAL
  &DFT
    &MGRID
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      EPS_SCF 1.0E-7
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-7
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
  &END DFT
  &SUBSYS
    &KIND H
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
"""


def ase_version_tuple():
    parts = []
    for p in ase.__version__.split(".")[:2]:
        digits = "".join(ch for ch in p if ch.isdigit())
        parts.append(int(digits) if digits else 0)
    return tuple(parts)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--command", default=None,
                    help="command that starts CP2K in shell mode, "
                         "e.g. 'mpirun -np 4 cp2k.psmp -s'")
    ap.add_argument("--fmax", type=float, default=0.01,
                    help="force convergence for BFGS in eV/A (default 0.01)")
    ap.add_argument("--cutoff", type=float, default=500.0,
                    help="CUTOFF of the finest grid in Ry (default 500)")
    ap.add_argument("--vacuum", type=float, default=5.0,
                    help="vacuum around the molecule in A (default 5.0)")
    args = ap.parse_args()

    command = args.command or os.environ.get("ASE_CP2K_COMMAND")
    if ase_version_tuple() < (3, 29):
        if command is None or "cp2k_shell" not in command:
            sys.exit(f"ASE {ase.__version__} requires 'cp2k_shell' in the command. "
                     "Set ASE_CP2K_COMMAND to the wrapper ase_cp2k_shell_wrapper.sh "
                     "or upgrade to ASE >= 3.29.")
    if "CP2K_DATA_DIR" not in os.environ:
        print("WARNING: CP2K_DATA_DIR is not set; BASIS_MOLOPT_UZH may not be found.",
              file=sys.stderr)

    atoms = molecule("H2O")
    atoms.center(vacuum=args.vacuum)      # sets a box; pbc stays False
    atoms.pbc = False                     # -> ASE adds PERIODIC NONE + MT Poisson solver

    calc = CP2K(
        label="h2o_ase",                  # becomes PROJECT; files h2o_ase.inp / h2o_ase.out
        command=command,                  # None -> ASE_CP2K_COMMAND or ASE default
        inp=INP,
        xc="PBE",
        cutoff=args.cutoff * Rydberg,     # ASE expects eV and writes CUTOFF [eV]
        max_scf=30,                       # inner SCF steps (OT), outer loop set in INP
        basis_set_file="BASIS_MOLOPT_UZH",
        potential_file="POTENTIAL_UZH",
        basis_set=None,                   # None: KIND basis comes from the template
        pseudo_potential=None,            # None: KIND potential comes from the template
        stress_tensor=False,              # no stress for an isolated molecule
        charge=0,
        uks=False,
        print_level="MEDIUM",
    )
    print(f"ASE {ase.__version__}; CP2K command: {calc.command}")

    with calc:                            # closes the CP2K shell at the end
        atoms.calc = calc
        e0 = atoms.get_potential_energy()
        print(f"initial energy: {e0:.6f} eV")
        opt = BFGS(atoms, trajectory="h2o_opt.traj", logfile="h2o_opt.log")
        opt.run(fmax=args.fmax, steps=100)
        e1 = atoms.get_potential_energy()

    d_oh = atoms.get_distance(0, 1)
    a_hoh = atoms.get_angle(1, 0, 2)
    print(f"final energy  : {e1:.6f} eV  (relaxation {e1 - e0:.4f} eV)")
    print(f"r(O-H) = {d_oh:.4f} A, angle(H-O-H) = {a_hoh:.2f} deg")
    write("h2o_opt.xyz", atoms)
    print("final structure written to h2o_opt.xyz, BFGS log in h2o_opt.log")


if __name__ == "__main__":
    main()
