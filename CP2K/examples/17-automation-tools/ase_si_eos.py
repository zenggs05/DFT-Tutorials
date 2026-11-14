#!/usr/bin/env python3
"""Equation of state (E-V curve) of bulk Si with ASE driving CP2K (k-points via template).

Chapter 17 of the CP2K tutorial (ASE CP2K calculator, periodic system).

For every lattice constant a new CP2K calculator (and a new CP2K shell process) is
created. This is deliberate: when only the cell of an existing calculator changes,
ASE sends SET_CELL to the running CP2K shell and the existing force environment is
reused (as in a cell optimization, the grids are set up for the first cell). A fresh
force environment per volume gives the usual "fixed CUTOFF" E-V curve and makes every
point independent of the previous one.

Everything the calculator cannot express (k-points, diagonalization, mixing, UZH
basis/potential names) is in the template INP. Keywords that ASE writes itself
(PROJECT, PRINT_LEVEL, CUTOFF, MAX_SCF, XC, CELL, COORD, ...) must not be in INP.

CP2K 2026.2 and earlier: k-points require diagonalization (no OT).

Usage:
    export CP2K_DATA_DIR=/path/to/cp2k/data
    export ASE_CP2K_COMMAND="cp2k.psmp -s"          # ASE >= 3.29
    python ase_si_eos.py                            # 7 points, +-3 % in a
    python ase_si_eos.py --kpts 8 --cutoff 600 --npoints 9 --range 0.04

Requirements: python >= 3.10, ase, numpy, (matplotlib for the plot), a CP2K binary.
"""
import argparse
import os
import sys

import numpy as np

import ase
from ase.build import bulk
from ase.calculators.cp2k import CP2K
from ase.eos import EquationOfState
from ase.units import Rydberg, kJ

INP_TEMPLATE = """
&FORCE_EVAL
  &DFT
    &MGRID
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK {k} {k} {k}
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-8
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
        NBROYDEN 8
      &END MIXING
    &END SCF
  &END DFT
  &SUBSYS
    &KIND Si
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE-q4
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
                    help="command that starts CP2K in shell mode "
                         "(default: ASE_CP2K_COMMAND or the ASE default)")
    ap.add_argument("--a0", type=float, default=5.431, help="central lattice constant in A")
    ap.add_argument("--range", type=float, default=0.03,
                    help="relative variation of a, e.g. 0.03 = +-3 %% (default)")
    ap.add_argument("--npoints", type=int, default=7, help="number of E-V points (default 7)")
    ap.add_argument("--kpts", type=int, default=6, help="Monkhorst-Pack grid n x n x n (default 6)")
    ap.add_argument("--cutoff", type=float, default=500.0, help="CUTOFF in Ry (default 500)")
    args = ap.parse_args()

    command = args.command or os.environ.get("ASE_CP2K_COMMAND")
    if ase_version_tuple() < (3, 29) and (command is None or "cp2k_shell" not in command):
        sys.exit(f"ASE {ase.__version__} requires 'cp2k_shell' in the command; use "
                 "ase_cp2k_shell_wrapper.sh or upgrade to ASE >= 3.29.")
    if "CP2K_DATA_DIR" not in os.environ:
        print("WARNING: CP2K_DATA_DIR is not set; BASIS_MOLOPT_UZH may not be found.",
              file=sys.stderr)

    inp = INP_TEMPLATE.format(k=args.kpts)
    scales = np.linspace(1.0 - args.range, 1.0 + args.range, args.npoints)
    volumes, energies = [], []

    with open("si_eos.dat", "w") as fout:
        fout.write("# a[A]  V/atom[A^3]  E/atom[eV]\n")
        for i, s in enumerate(scales):
            a = args.a0 * s
            atoms = bulk("Si", "diamond", a=a)        # 2-atom primitive cell, pbc=True
            calc = CP2K(
                label=f"si_eos_{i:02d}",
                command=command,
                inp=inp,
                xc="PBE",
                cutoff=args.cutoff * Rydberg,          # ASE expects eV
                max_scf=100,
                basis_set_file="BASIS_MOLOPT_UZH",
                potential_file="POTENTIAL_UZH",
                basis_set=None,                        # taken from &KIND in the template
                pseudo_potential=None,
                stress_tensor=False,                   # E-V fit needs energies only
                print_level="MEDIUM",
            )
            with calc:
                atoms.calc = calc
                e = atoms.get_potential_energy()
            n = len(atoms)
            volumes.append(atoms.get_volume() / n)
            energies.append(e / n)
            fout.write(f"{a:10.5f} {volumes[-1]:12.5f} {energies[-1]:16.8f}\n")
            fout.flush()
            print(f"a = {a:.4f} A  V/atom = {volumes[-1]:.4f} A^3  E/atom = {energies[-1]:.6f} eV")

    eos = EquationOfState(volumes, energies, eos="birchmurnaghan")
    v0, e0, b = eos.fit()
    a_eq = (4.0 * 2 * v0) ** (1.0 / 3.0)          # V_prim = a^3 / 4 for fcc, 2 atoms per cell
    print(f"V0 = {v0:.4f} A^3/atom, a0 = {a_eq:.4f} A, B = {b / kJ * 1.0e24:.1f} GPa")
    try:
        eos.plot("si_eos.png")
        print("plot written to si_eos.png, data in si_eos.dat")
    except Exception as exc:                       # matplotlib missing or no display
        print(f"no plot ({exc}); data in si_eos.dat")


if __name__ == "__main__":
    main()
