#!/usr/bin/env python3
"""Pre-optimise a structure with a machine-learned potential, then hand it to DFT.

Typical use: a rough POSCAR/CIF/xyz (built by hand, from a database or by
substitution) is relaxed with MACE-MP (or DeePMD / NEP) in seconds to minutes;
the result is written as the starting geometry of the real DFT relaxation.
The MLIP energy is NOT a DFT energy - never report or mix it with DFT numbers.

Examples
  python mace_ase_preopt.py POSCAR_guess                       # positions only, MACE-MP medium, CPU
  python mace_ase_preopt.py MoS2.cif --cell --fmax 0.05 -o POSCAR_preopt
  python mace_ase_preopt.py slab.vasp --calc mace-mp --d3 --device cuda
  python mace_ase_preopt.py mol.xyz --calc mace-off -o mol_preopt.xyz
  python mace_ase_preopt.py POSCAR --calc nep --model nep.txt
  python mace_ase_preopt.py POSCAR --calc emt                    # toy test without any MLIP package

Constraints present in the input (e.g. Selective dynamics F F F in a POSCAR)
are kept, so fixed slab layers stay fixed.
"""
import argparse
import sys

import numpy as np
from ase.io import read, write
from ase.optimize import BFGS, FIRE

from mlip_calc import KINDS, make_calculator


def cell_filter(atoms):
    try:
        from ase.filters import FrechetCellFilter      # ASE >= 3.23
        return FrechetCellFilter(atoms)
    except ImportError:
        from ase.constraints import ExpCellFilter      # older ASE
        return ExpCellFilter(atoms)


def max_force(atoms):
    f = atoms.get_forces()
    return float(np.sqrt((f ** 2).sum(axis=1)).max()) if len(f) else 0.0


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("structure", help="any file ASE can read (POSCAR, CIF, xyz, gjf, ...)")
    p.add_argument("--calc", default="mace-mp", choices=KINDS)
    p.add_argument("--model", default=None, help="model name or file (mace-mp: small/medium/large or a path)")
    p.add_argument("--device", default="cpu", help="cpu, cuda or mps")
    p.add_argument("--dtype", default="float64", choices=["float32", "float64"])
    p.add_argument("--d3", action="store_true", help="add DFT-D3(BJ) to MACE-MP (needs torch-dftd)")
    p.add_argument("--cell", action="store_true", help="also relax the cell (bulk / 2D in-plane: check vacuum!)")
    p.add_argument("--fmax", type=float, default=0.05, help="eV/Angstrom (default 0.05)")
    p.add_argument("--steps", type=int, default=500)
    p.add_argument("--opt", default="FIRE", choices=["FIRE", "BFGS"])
    p.add_argument("-o", "--output", default="POSCAR_preopt")
    p.add_argument("--traj", default="preopt.traj")
    a = p.parse_args()

    atoms = read(a.structure)
    atoms.calc = make_calculator(a.calc, a.model, a.device, a.dtype, a.d3)
    e0, f0, v0 = atoms.get_potential_energy(), max_force(atoms), atoms.get_volume() if atoms.pbc.any() else 0.0
    print(f"start : E = {e0:.4f} eV  ({e0 / len(atoms):.4f} eV/atom)  max|F| = {f0:.3f} eV/A")

    target = cell_filter(atoms) if a.cell else atoms
    Opt = FIRE if a.opt == "FIRE" else BFGS
    opt = Opt(target, trajectory=a.traj, logfile="-")
    converged = opt.run(fmax=a.fmax, steps=a.steps)

    e1, f1 = atoms.get_potential_energy(), max_force(atoms)
    print(f"end   : E = {e1:.4f} eV  ({e1 / len(atoms):.4f} eV/atom)  max|F| = {f1:.3f} eV/A")
    if a.cell and atoms.pbc.any():
        print(f"volume: {v0:.3f} -> {atoms.get_volume():.3f} A^3 ({100 * (atoms.get_volume() / v0 - 1):+.2f} %)")
    if not converged:
        print(f"WARNING: not converged in {a.steps} steps - inspect {a.traj} (ase gui {a.traj})", file=sys.stderr)

    # Short interatomic distances usually mean the MLIP left its training domain
    if len(atoms) > 1:
        d = atoms.get_all_distances(mic=bool(atoms.pbc.any()))
        dmin = d[np.triu_indices(len(atoms), 1)].min()
        print(f"shortest interatomic distance: {dmin:.3f} A")
        if dmin < 0.7:
            print("WARNING: atoms collapsed - do not use this geometry", file=sys.stderr)

    atoms.calc = None
    if a.output.upper().startswith("POSCAR") or a.output.endswith(".vasp"):
        write(a.output, atoms, format="vasp", direct=True, sort=False)
    else:
        write(a.output, atoms)
    print(f"written: {a.output}  -> next step: DFT relaxation starting from this file")


if __name__ == "__main__":
    main()
