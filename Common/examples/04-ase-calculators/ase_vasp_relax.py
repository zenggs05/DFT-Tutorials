#!/usr/bin/env python3
"""Relax bulk Si with VASP driven by ASE, in two different ways.

  --mode vasp : VASP's own optimizer (IBRION=2, ISIF=3).  ASE only writes
                INCAR/KPOINTS/POSCAR/POTCAR, runs VASP once and reads CONTCAR.
  --mode ase  : ASE optimizer (BFGS + FrechetCellFilter); VASP does single
                points (NSW=0) and returns energy/forces/stress each step.

Requirements (see README):
  export VASP_PP_PATH=/path/to/pp           # contains potpaw_PBE/ (or potpaw_PBE.64/)
  export ASE_VASP_COMMAND="mpirun -np 16 vasp_std"
Usage:
  python ase_vasp_relax.py --mode vasp
  python ase_vasp_relax.py --mode ase --fmax 0.01
  python ase_vasp_relax.py --dry-run          # only write the inputs to ./si_relax
"""
import argparse
import os

from ase.build import bulk
from ase.calculators.vasp import Vasp
from ase.filters import FrechetCellFilter
from ase.io import write
from ase.optimize import BFGS

COMMON = dict(
    xc="PBE",              # selects potpaw_PBE (+ GGA=PE)
    setups="recommended",  # VASP-recommended POTCAR variants (Si -> Si)
    encut=520,             # ENMAX(Si) is ~245 eV; >=1.3x and higher for stress (Pulay)
    prec="Accurate",
    kpts=(8, 8, 8),
    gamma=True,            # Gamma-centred mesh
    ismear=0, sigma=0.05,  # Gaussian smearing for a semiconductor
    ediff=1e-7,
    lreal=False,
    lwave=False, lcharg=False,
    ncore=4,               # adjust to your node (see Common chapter 2)
)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["vasp", "ase"], default="vasp")
    p.add_argument("--fmax", type=float, default=0.01, help="eV/A, ASE mode only")
    p.add_argument("--dry-run", action="store_true", help="write inputs, do not run")
    args = p.parse_args()

    si = bulk("Si", "diamond", a=5.40)       # deliberately off the PBE minimum
    workdir = f"si_relax_{args.mode}"

    if args.mode == "vasp":
        calc = Vasp(directory=workdir, ibrion=2, isif=3, nsw=60, ediffg=-0.005,
                    **COMMON)
    else:
        calc = Vasp(directory=workdir, ibrion=-1, nsw=0, isif=2, **COMMON)
    si.calc = calc

    if args.dry_run:
        os.makedirs(workdir, exist_ok=True)
        calc.write_input(si)
        print("inputs written to", workdir, ":", sorted(os.listdir(workdir)))
        return

    if args.mode == "vasp":
        e = si.get_potential_energy()      # runs VASP; atoms updated from CONTCAR
        outcar = os.path.join(workdir, "OUTCAR")
        with open(outcar) as f:
            done = "reached required accuracy" in f.read()
        print("VASP ionic relaxation converged:", done)
    else:
        ecf = FrechetCellFilter(si)        # relax positions + cell (ASE >= 3.23)
        opt = BFGS(ecf, trajectory="si_ase_relax.traj", logfile="si_ase_relax.log")
        opt.run(fmax=args.fmax, steps=50)
        e = si.get_potential_energy()

    a = si.cell.cellpar()[0] * 2 ** 0.5     # primitive fcc vector -> cubic a
    print(f"E = {e:.6f} eV   a(cubic) = {a:.4f} A")
    write("CONTCAR_ase_" + args.mode, si, format="vasp", direct=True)


if __name__ == "__main__":
    main()
