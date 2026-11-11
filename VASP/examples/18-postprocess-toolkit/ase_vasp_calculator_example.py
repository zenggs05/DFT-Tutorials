#!/usr/bin/env python3
"""Drive VASP from ASE: E-V scan of Si and an IDPP-interpolated NEB path (input only).

Environment needed for real runs:
    export ASE_VASP_COMMAND="mpirun -np 8 vasp_std"
    export VASP_PP_PATH=/path/to/potentials        # contains potpaw_PBE/ (ASE uses it for xc='PBE')

Examples
--------
    python ase_vasp_calculator_example.py --dry-run      # only write INCAR/KPOINTS/POSCAR/POTCAR, VASP is not run
    python ase_vasp_calculator_example.py                # E-V scan with VASP, fit with ASE EquationOfState
    python ase_vasp_calculator_example.py --neb-demo     # write neb/00..06 POSCARs (IDPP interpolation, no VASP)

Requires ASE >= 3.23 (NEB lives in ase.mep); tested with the API of ASE 3.29.
"""
import argparse
import os

import numpy as np
from ase.build import bulk
from ase.calculators.vasp import Vasp
from ase.eos import EquationOfState
from ase.io import write


def make_calc(directory):
    # keyword names are the lower-case INCAR tags; ASE writes INCAR, KPOINTS, POSCAR, POTCAR
    return Vasp(directory=directory,
                xc="PBE",                 # selects potpaw_PBE in $VASP_PP_PATH
                setups={"Si": ""},        # standard Si POTCAR (use e.g. {"Ga": "_d"} for Ga_d)
                prec="Accurate",
                encut=400,
                ediff=1e-6,
                ismear=0,
                sigma=0.05,
                kpts=(8, 8, 8),
                gamma=True,               # Gamma-centred mesh
                lwave=False,
                lcharg=False)


def ev_scan(dry_run):
    a0 = 5.431
    volumes, energies = [], []
    for scale in np.linspace(0.97, 1.03, 7):
        a = a0 * scale
        atoms = bulk("Si", "diamond", a=a)
        d = "ev_%.4f" % a
        calc = make_calc(d)
        atoms.calc = calc
        if dry_run:
            # ASE needs the POTCAR library even to write inputs (it builds POTCAR and the INCAR together)
            if not os.environ.get("VASP_PP_PATH"):
                raise SystemExit("set VASP_PP_PATH (directory containing potpaw_PBE/) first")
            os.makedirs(d, exist_ok=True)
            calc.write_input(atoms)               # INCAR, KPOINTS, POSCAR, POTCAR; VASP is not run
            print("wrote inputs in", d)
            continue
        e = atoms.get_potential_energy()          # runs VASP in folder d
        volumes.append(atoms.get_volume())
        energies.append(e)
        print("a = %.4f A  V = %.3f A^3  E = %.6f eV" % (a, volumes[-1], e))
    if dry_run:
        return
    eos = EquationOfState(volumes, energies, eos="birchmurnaghan")
    v0, e0, b = eos.fit()
    a_fit = (4 * v0) ** (1 / 3)                  # fcc primitive cell holds 1/4 of the cubic cell
    print("V0 = %.3f A^3, a0 = %.4f A, B0 = %.1f GPa" % (v0, a_fit, b * 160.21766208))
    eos.plot("si_eos.png")


def neb_demo():
    """Linear vs IDPP interpolation for a vacancy hop in a small Si cell (geometry only)."""
    from ase.mep import NEB

    sc = bulk("Si", "diamond", a=5.431, cubic=True)       # 8 atoms
    initial = sc.copy()
    del initial[0]                                          # vacancy on site 0
    final = sc.copy()
    final.positions[1] = sc.positions[0]                    # nearest neighbour (site 1) hops into the vacancy
    del final[0]                                            # same atom order as 'initial'
    images = [initial] + [initial.copy() for _ in range(5)] + [final]
    neb = NEB(images, method="improvedtangent")
    neb.interpolate("idpp")                                 # image-dependent pair potential
    for i, im in enumerate(images):
        d = "neb/%02d" % i
        os.makedirs(d, exist_ok=True)
        write(os.path.join(d, "POSCAR"), im, format="vasp", direct=True)
    print("wrote neb/00 ... neb/%02d (IMAGES = %d for VASP/VTST)" % (len(images) - 1, len(images) - 2))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dry-run", action="store_true", help="write input files only")
    p.add_argument("--neb-demo", action="store_true", help="write IDPP-interpolated NEB images")
    args = p.parse_args()
    if args.neb_demo:
        neb_demo()
    else:
        ev_scan(args.dry_run)


if __name__ == "__main__":
    main()
