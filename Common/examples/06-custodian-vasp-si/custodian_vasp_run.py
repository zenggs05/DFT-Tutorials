#!/usr/bin/env python
"""Run a VASP relax -> relax -> static sequence under custodian.

Every step is watched by custodian error handlers.  When a known VASP error
(ZBRENT, EDDDAV, BRMIX, SCF not converged, ...) is detected, custodian
backs up the files to error.N.tar.gz, edits INCAR/POSCAR and restarts VASP.

Usage (inside a Slurm allocation, NOT on the login node):
    export VASP_CMD="srun vasp_std"            # or "mpirun -np 64 vasp_std"
    python custodian_vasp_run.py               # uses INCAR KPOINTS POSCAR POTCAR in cwd

Tested against custodian 2026.10.21 (Python >= 3.11).
"""
import os
import shlex
import sys

from custodian.custodian import Custodian
from custodian.vasp.handlers import (
    FrozenJobErrorHandler,
    MeshSymmetryErrorHandler,
    NonConvergingErrorHandler,
    PositiveEnergyErrorHandler,
    PotimErrorHandler,
    StdErrHandler,
    UnconvergedErrorHandler,
    VaspErrorHandler,
    WalltimeHandler,
)
from custodian.vasp.jobs import VaspJob
from custodian.vasp.validators import VaspFilesValidator, VasprunXMLValidator


def main() -> int:
    for name in ("INCAR", "KPOINTS", "POSCAR", "POTCAR"):
        if not os.path.isfile(name):
            print(f"[error] {name} not found in {os.getcwd()}", file=sys.stderr)
            return 1

    # custodian needs the command as a list of arguments, not one string
    vasp_cmd = shlex.split(os.environ.get("VASP_CMD", "srun vasp_std"))

    handlers = [
        VaspErrorHandler(),             # monitor: scans vasp.out for ~50 known messages
        MeshSymmetryErrorHandler(),     # k-mesh does not respect the symmetry
        UnconvergedErrorHandler(),      # after the run: electronic/ionic not converged
        NonConvergingErrorHandler(),    # monitor: last 10 ionic steps all hit NELM
        PotimErrorHandler(),            # monitor: energy jumps -> reduce POTIM
        PositiveEnergyErrorHandler(),   # monitor: positive total energy
        FrozenJobErrorHandler(timeout=7200),  # vasp.out not updated for 2 h
        StdErrHandler(),                # std_err.txt: MPI / memory problems
    ]
    # Writes a STOPCAR shortly before the Slurm time limit (needs SBATCH_TIMELIMIT
    # or an explicit wall_time in seconds); otherwise the handler does nothing.
    if "SBATCH_TIMELIMIT" in os.environ or "CUSTODIAN_WALLTIME" in os.environ:
        wall = os.environ.get("CUSTODIAN_WALLTIME")
        handlers.append(WalltimeHandler(wall_time=int(wall) if wall else None, buffer_time=600))

    copy_contcar = {"file": "CONTCAR", "action": {"_file_copy": {"dest": "POSCAR"}}}
    jobs = [
        # 1st relaxation: outputs are copied to *.relax1, originals backed up as *.orig
        VaspJob(vasp_cmd, final=False, suffix=".relax1"),
        # 2nd relaxation starts from CONTCAR: removes the Pulay-stress bias of run 1
        VaspJob(vasp_cmd, final=False, suffix=".relax2", backup=False,
                settings_override=[copy_contcar]),
        # static run on the relaxed cell: accurate total energy / charge density
        VaspJob(vasp_cmd, final=True, suffix=".static", backup=False,
                settings_override=[
                    copy_contcar,
                    {"dict": "INCAR",
                     "action": {"_set": {"NSW": 0, "IBRION": -1, "ISIF": 2, "ALGO": "Normal"}}},
                ]),
    ]

    validators = [VasprunXMLValidator(), VaspFilesValidator()]

    c = Custodian(
        handlers,
        jobs,
        validators=validators,
        max_errors=10,          # total corrections allowed for the whole sequence
        max_errors_per_job=5,   # corrections allowed for one VaspJob
        polling_time_step=10,   # check every 10 s whether VASP is still running
        monitor_freq=30,        # run monitor handlers every 30 x 10 s = 5 min
    )
    c.run()                     # raises CustodianError subclasses on failure
    print("custodian finished; see custodian.json for the list of corrections")
    return 0


if __name__ == "__main__":
    sys.exit(main())
