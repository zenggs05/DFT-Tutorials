#!/usr/bin/env python3
"""Equation of state (E-V) of fcc Cu with ASE.

Default calculator is EMT (an empirical potential shipped with ASE), so the
script runs anywhere in a few seconds.  Use --calc vasp to repeat the same
workflow with VASP (needs VASP_PP_PATH and ASE_VASP_COMMAND, see README).

Usage:
    python ase_eos_cu.py                 # EMT, prints a0 and B, writes cu_eos.png
    python ase_eos_cu.py --calc vasp     # same scan with VASP (one dir per volume)
"""
import argparse

import numpy as np
from ase.build import bulk
from ase.eos import EquationOfState
from ase.io import write
from ase.units import kJ


def make_calc(name, tag):
    """Return a fresh calculator for one volume point."""
    if name == "emt":
        from ase.calculators.emt import EMT
        return EMT()
    if name == "vasp":
        from ase.calculators.vasp import Vasp
        # One directory per strain so that outputs do not overwrite each other.
        return Vasp(directory=f"eos_{tag}", xc="PBE", setups="recommended",
                    encut=500, kpts=(14, 14, 14), gamma=True,
                    ismear=1, sigma=0.1, ediff=1e-6,
                    lwave=False, lcharg=False)
    raise ValueError(f"unknown calculator: {name}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--calc", default="emt", choices=["emt", "vasp"])
    parser.add_argument("--a0", type=float, default=3.61, help="initial lattice constant (A)")
    parser.add_argument("--eos", default="birchmurnaghan",
                        help="EOS form: sj, birchmurnaghan, murnaghan, vinet ...")
    args = parser.parse_args()

    cu0 = bulk("Cu", "fcc", a=args.a0)          # 1-atom primitive cell
    cell0 = cu0.get_cell()
    volumes, energies, frames = [], [], []
    # Scale the lattice by -3% ... +3% (volume range about +-9%).
    for i, s in enumerate(np.linspace(0.97, 1.03, 7)):
        atoms = cu0.copy()
        atoms.set_cell(cell0 * s, scale_atoms=True)
        atoms.calc = make_calc(args.calc, f"{i:02d}")
        e = atoms.get_potential_energy()        # eV per cell
        volumes.append(atoms.get_volume())      # A^3 per cell
        energies.append(e)
        frames.append(atoms)
        print(f"scale={s:.3f}  V={volumes[-1]:8.4f} A^3  E={e:12.6f} eV")

    write("cu_eos.traj", frames)               # view with: ase gui cu_eos.traj
    eos = EquationOfState(volumes, energies, eos=args.eos)
    v0, e0, B = eos.fit()                      # B in eV/A^3
    a0 = (4.0 * v0) ** (1.0 / 3.0)             # fcc: V_conv = 4 * V_prim
    print(f"\nEOS = {args.eos}")
    print(f"V0  = {v0:.4f} A^3/atom   a0 = {a0:.4f} A")
    print(f"E0  = {e0:.6f} eV          B  = {B / kJ * 1.0e24:.1f} GPa")
    eos.plot("cu_eos.png")
    print("figure written: cu_eos.png")


if __name__ == "__main__":
    main()
