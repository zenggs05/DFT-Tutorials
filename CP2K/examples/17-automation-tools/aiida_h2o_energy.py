#!/usr/bin/env python3
"""Single-point energy of H2O through AiiDA + aiida-cp2k (Cp2kBaseWorkChain).

Chapter 17 of the CP2K tutorial (AiiDA introduction).

Prerequisites (once):
    pip install aiida-cp2k            # pulls aiida-core >= 2.7 (Python >= 3.10)
    verdi presto                      # quick profile (or: verdi profile setup ...)
    verdi computer setup / verdi computer configure core.ssh <label>   # remote cluster
    verdi code create core.code.installed       # interactive; executable = cp2k.psmp,
                                                # default plugin = cp2k
    verdi code list                   # e.g. cp2k@mycluster

Run:
    python aiida_h2o_energy.py cp2k@mycluster              # run in the foreground
    python aiida_h2o_energy.py cp2k@mycluster --submit     # submit to the daemon
    verdi process list -a ; verdi process show <PK>
    verdi calcjob outputcat <PK of the CalcJob> | tail     # look at aiida.out

Notes
  * aiida-cp2k writes PROJECT aiida, the cell vectors A/B/C and an XYZ coordinate
    file from the StructureData. It does NOT translate pbc into PERIODIC: for an
    isolated molecule PERIODIC NONE (CELL and POISSON) must be in the parameters.
  * By default the basis/potential files are taken from the data directory of the
    CP2K installation on the remote computer. With --with-data-files the local
    files from $CP2K_DATA_DIR are uploaded as SinglefileData (stored in provenance).
  * Cp2kBaseWorkChain restarts failed calculations for some known errors
    (e.g. walltime exceeded); see the aiida-cp2k documentation for details.
"""
import argparse
import os
import sys

from ase.build import molecule

from aiida import load_profile
from aiida.engine import run_get_node, submit
from aiida.orm import Dict, SinglefileData, StructureData, load_code
from aiida.plugins import WorkflowFactory

PARAMETERS = {
    "GLOBAL": {
        "RUN_TYPE": "ENERGY_FORCE",
        "PRINT_LEVEL": "MEDIUM",
    },
    "FORCE_EVAL": {
        "METHOD": "QUICKSTEP",
        "DFT": {
            "BASIS_SET_FILE_NAME": "BASIS_MOLOPT_UZH",
            "POTENTIAL_FILE_NAME": "POTENTIAL_UZH",
            "MGRID": {"CUTOFF": 500, "REL_CUTOFF": 60},
            "QS": {"EPS_DEFAULT": 1.0e-12},
            "POISSON": {"PERIODIC": "NONE", "POISSON_SOLVER": "MT"},
            "SCF": {
                "SCF_GUESS": "ATOMIC",
                "EPS_SCF": 1.0e-7,
                "MAX_SCF": 30,
                "OT": {"_": "ON", "MINIMIZER": "DIIS", "PRECONDITIONER": "FULL_ALL"},
                "OUTER_SCF": {"_": "ON", "EPS_SCF": 1.0e-7, "MAX_SCF": 10},
            },
            "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},
        },
        "SUBSYS": {
            "CELL": {"PERIODIC": "NONE"},
            "KIND": [
                {"_": "O", "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q6", "POTENTIAL": "GTH-PBE-q6"},
                {"_": "H", "BASIS_SET": "DZVP-MOLOPT-PBE-GTH-q1", "POTENTIAL": "GTH-PBE-q1"},
            ],
        },
        "PRINT": {"FORCES": {"_": "ON"}},
    },
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("code", help="AiiDA code label, e.g. cp2k@mycluster")
    ap.add_argument("--submit", action="store_true",
                    help="submit to the AiiDA daemon instead of running in the foreground")
    ap.add_argument("--nprocs", type=int, default=4, help="MPI processes (default 4)")
    ap.add_argument("--walltime", type=int, default=1800, help="max wallclock in s (default 1800)")
    ap.add_argument("--with-data-files", action="store_true",
                    help="upload BASIS_MOLOPT_UZH/POTENTIAL_UZH from $CP2K_DATA_DIR")
    args = ap.parse_args()

    load_profile()
    code = load_code(args.code)

    atoms = molecule("H2O")
    atoms.center(vacuum=5.0)
    atoms.pbc = False
    structure = StructureData(ase=atoms)

    Cp2kBaseWorkChain = WorkflowFactory("cp2k.base")
    builder = Cp2kBaseWorkChain.get_builder()
    builder.cp2k.code = code
    builder.cp2k.structure = structure
    builder.cp2k.parameters = Dict(PARAMETERS)
    if args.with_data_files:
        data_dir = os.environ.get("CP2K_DATA_DIR")
        if not data_dir:
            sys.exit("--with-data-files needs CP2K_DATA_DIR")
        builder.cp2k.file = {
            "basis": SinglefileData(file=os.path.join(data_dir, "BASIS_MOLOPT_UZH")),
            "pseudo": SinglefileData(file=os.path.join(data_dir, "POTENTIAL_UZH")),
        }
    builder.cp2k.metadata.options.resources = {
        "num_machines": 1,
        "num_mpiprocs_per_machine": args.nprocs,
    }
    builder.cp2k.metadata.options.max_wallclock_seconds = args.walltime

    if args.submit:
        node = submit(builder)
        print(f"submitted Cp2kBaseWorkChain <{node.pk}>; follow with 'verdi process list -a'")
        return

    results, node = run_get_node(builder)
    print(f"Cp2kBaseWorkChain <{node.pk}> finished, exit status {node.exit_status}")
    out = results.get("output_parameters")
    if out is not None:
        d = out.get_dict()
        print(f"energy = {d.get('energy')} {d.get('energy_units', '')}")


if __name__ == "__main__":
    main()
