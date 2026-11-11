#!/usr/bin/env python3
"""Run VASP under custodian: automatic detection and correction of common errors.

Two modes:
  single  : one VASP run in the current directory with error handlers
  double  : AFLOW-style double relaxation (relax1 -> relax2), handy for ISIF = 3

Every correction is logged in custodian.json (which error, which INCAR/KPOINTS change).
Requires: pip install custodian pymatgen

Examples (inside a job script, VASP inputs in the current directory)
--------
    python custodian_run.py --cmd "mpirun -np 32 vasp_std"
    python custodian_run.py --cmd "srun vasp_std" --mode double --max-errors 5
"""
import argparse
import shlex

from custodian.custodian import Custodian
from custodian.vasp.handlers import (
    MeshSymmetryErrorHandler,
    NonConvergingErrorHandler,
    PotimErrorHandler,
    UnconvergedErrorHandler,
    VaspErrorHandler,
)
from custodian.vasp.jobs import VaspJob
from custodian.vasp.validators import VasprunXMLValidator


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cmd", required=True, help='VASP command, e.g. "mpirun -np 32 vasp_std"')
    p.add_argument("--mode", choices=["single", "double"], default="single")
    p.add_argument("--max-errors", type=int, default=5)
    args = p.parse_args()
    vasp_cmd = shlex.split(args.cmd)

    handlers = [
        VaspErrorHandler(),            # parses vasp.out for ZBRENT, BRMIX, EDDDAV, Sub-Space-Matrix, ...
        MeshSymmetryErrorHandler(),
        UnconvergedErrorHandler(),     # electronic / ionic steps not converged
        NonConvergingErrorHandler(),
        PotimErrorHandler(),
    ]
    validators = [VasprunXMLValidator()]

    if args.mode == "double":
        jobs = VaspJob.double_relaxation_run(vasp_cmd)
    else:
        jobs = [VaspJob(vasp_cmd)]     # stdout goes to vasp.out, which the handlers read

    c = Custodian(handlers, jobs, validators=validators, max_errors=args.max_errors)
    c.run()


if __name__ == "__main__":
    main()
