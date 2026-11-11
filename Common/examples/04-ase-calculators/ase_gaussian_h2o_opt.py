#!/usr/bin/env python3
"""Optimise H2O with Gaussian 16 through ASE, then run a frequency job.

Gaussian's own optimizer is used (GaussianOptimizer -> 'opt' keyword), which
is much more efficient than driving Gaussian with an ASE optimizer.

Requirements:
  export ASE_GAUSSIAN_COMMAND="g16 < PREFIX.com > PREFIX.log"
  export GAUSS_SCRDIR=/local/scratch/$USER        # fast local disk
Usage:
  python ase_gaussian_h2o_opt.py              # opt + freq
  python ase_gaussian_h2o_opt.py --dry-run    # only write h2o_opt.com
"""
import argparse
import os

from ase.build import molecule
from ase.calculators.gaussian import Gaussian, GaussianOptimizer
from ase.io import read, write

LINK0 = dict(mem="4GB", nprocshared=8)       # %mem and %nprocshared
LEVEL = dict(method="B3LYP", basis="def2SVP",
             EmpiricalDispersion="GD3BJ",     # free-form route keyword=value
             scf="tight", charge=0, mult=1)


def normal_termination(logfile):
    with open(logfile) as f:
        tail = f.read()[-2000:]
    return "Normal termination" in tail


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()
    os.makedirs("g16", exist_ok=True)

    h2o = molecule("H2O")
    h2o.positions[1, 2] += 0.10             # distort a little so opt has work to do

    if args.dry_run:
        # Write the input that the optimizer would produce, without running g16
        write("g16/h2o_opt.com", h2o, format="gaussian-in", chk="h2o.chk",
              opt="tight", **LINK0, **LEVEL)
        print(open("g16/h2o_opt.com").read())
        return

    calc_opt = Gaussian(label="g16/h2o_opt", chk="h2o.chk", **LINK0, **LEVEL)
    opt = GaussianOptimizer(h2o, calc_opt)
    ok = opt.run(fmax="tight", steps=100)   # -> opt=(tight,maxcycle=100)
    print("opt finished:", ok, " normal termination:",
          normal_termination("g16/h2o_opt.log"))

    # Frequency job on the optimised geometry (same level of theory!)
    calc_freq = Gaussian(label="g16/h2o_freq", chk="h2o.chk", freq="", **LINK0, **LEVEL)
    h2o.calc = calc_freq
    e = h2o.get_potential_energy()          # eV (ASE converts from Hartree)
    print(f"E(SCF) = {e:.6f} eV")

    last = read("g16/h2o_freq.log", format="gaussian-out")
    print("final geometry (A):")
    print(last.get_positions().round(4))
    # Frequencies are not parsed by ASE; use cclib (see 04-output-parsing)


if __name__ == "__main__":
    main()
