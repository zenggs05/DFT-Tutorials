#!/usr/bin/env python3
"""Single-point SCF of diamond Si with OpenMX through ASE.

ASE's OpenMX interface is useful to generate the structure blocks and to
read results back; for anything advanced (DFT+U, NEGF, Wannier ...) write the
.dat file by hand (see OpenMX chapters) and only use ASE for parsing.

Unit conventions of the ASE OpenMX calculator (important!):
  energy_cutoff : eV   (written as scf.energycutoff in Ry)
  convergence   : eV   (written as scf.criterion in Hartree)

Requirements:
  export OPENMX_DFT_DATA_PATH=/path/to/openmx4.0/DFT_DATA19
  export ASE_OPENMX_COMMAND=openmx                 # needed by ASE >= 3.23
Usage:
  python ase_openmx_si.py --dry-run               # write si.dat only
  python ase_openmx_si.py --np 8                  # run: mpirun -np 8 openmx si.dat -nt 1
"""
import argparse
import os

from ase.build import bulk
from ase.calculators.openmx import OpenMX
from ase.units import Ha, Ry


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--np", type=int, default=4, help="MPI processes")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    si = bulk("Si", "diamond", a=5.43)
    os.makedirs("omx", exist_ok=True)
    calc = OpenMX(
        label="omx/si",
        xc="PBE",                         # -> scf.XcType GGA-PBE
        energy_cutoff=220 * Ry,           # 220 Ry, given in eV
        convergence=1.0e-7 * Ha,          # 1e-7 Hartree, given in eV
        kpts=(8, 8, 8),                   # scf.Kgrid
        maxiter=100,
        eigensolver="Band",
        dft_data_year="19",               # DFT_DATA19 (OpenMX 3.9 / 4.0)
        # Basis "Standard" level for Si in DFT_DATA19 (see OpenMX chapter 4)
        definition_of_atomic_species=[["Si", "Si7.0-s2p2d1", "Si_PBE19"]],
        mpi={"processes": args.np, "threads": 1},
        # ASE would otherwise write an absolute "<cwd>/<label>" path here;
        # "./" keeps all outputs next to si.dat (ASE runs OpenMX inside omx/)
        system_currentdirectory="./",
    )
    si.calc = calc

    if args.dry_run:
        calc.write_input(si)
        print(open("omx/si.dat").read())
        return

    e = si.get_potential_energy()
    print(f"E = {e:.6f} eV")
    with open("omx/si.out") as f:
        txt = f.read()
    print("Utot line:", [l.strip() for l in txt.splitlines() if "Utot." in l][-1:])


if __name__ == "__main__":
    main()
