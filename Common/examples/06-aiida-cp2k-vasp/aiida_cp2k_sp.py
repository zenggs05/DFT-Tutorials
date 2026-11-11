#!/usr/bin/env python
"""AiiDA + aiida-cp2k: single-point PBE energy of one H2O molecule.

The CP2K input is written as a nested Python dict; aiida-cp2k turns it into
aiida.inp, adds the coordinates (aiida.coords.xyz) and the cell, submits the job
through the configured Computer/scheduler and parses aiida.out.

Prerequisites
  pip install aiida-core aiida-cp2k          (aiida-core 2.9.x, aiida-cp2k 2.x)
  verdi presto                                # profile + localhost computer
  verdi code create core.code.installed --config cp2k_code.yaml
Usage
  verdi run aiida_cp2k_sp.py                  # blocking run (no daemon needed)
  python aiida_cp2k_sp.py --submit            # submit to the daemon instead
"""
import argparse

import ase.build
from aiida import load_profile, orm
from aiida.engine import run_get_node, submit
from aiida.plugins import CalculationFactory

CODE_LABEL = "cp2k@localhost"   # label@computer as created with "verdi code create"


def build_parameters() -> orm.Dict:
    """CP2K input as nested sections. Keys are section/keyword names."""
    return orm.Dict(
        {
            "GLOBAL": {"RUN_TYPE": "ENERGY", "PRINT_LEVEL": "MEDIUM"},
            "FORCE_EVAL": {
                "METHOD": "Quickstep",
                "DFT": {
                    # File names are looked up in the CP2K data directory
                    "BASIS_SET_FILE_NAME": "BASIS_MOLOPT",
                    "POTENTIAL_FILE_NAME": "GTH_POTENTIALS",
                    "QS": {"EPS_DEFAULT": 1.0e-12},
                    "MGRID": {"CUTOFF": 400, "REL_CUTOFF": 50},
                    "POISSON": {"PERIODIC": "NONE", "POISSON_SOLVER": "MT"},
                    "SCF": {"SCF_GUESS": "ATOMIC", "EPS_SCF": 1.0e-6, "MAX_SCF": 50,
                            "OT": {"MINIMIZER": "DIIS", "PRECONDITIONER": "FULL_ALL"}},
                    "XC": {"XC_FUNCTIONAL": {"_": "PBE"}},   # "_" = section parameter
                },
                "SUBSYS": {
                    "CELL": {"PERIODIC": "NONE"},   # A/B/C are added from the structure
                    "KIND": [
                        {"_": "O", "BASIS_SET": "DZVP-MOLOPT-SR-GTH", "POTENTIAL": "GTH-PBE-q6"},
                        {"_": "H", "BASIS_SET": "DZVP-MOLOPT-SR-GTH", "POTENTIAL": "GTH-PBE-q1"},
                    ],
                },
            },
        }
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--submit", action="store_true", help="submit to the daemon")
    p.add_argument("--code", default=CODE_LABEL)
    args = p.parse_args()

    load_profile()                                   # default profile (verdi presto)

    atoms = ase.build.molecule("H2O")
    atoms.center(vacuum=5.0)                         # ~10 A box, non-periodic
    structure = orm.StructureData(ase=atoms)

    Cp2kCalculation = CalculationFactory("cp2k")
    builder = Cp2kCalculation.get_builder()
    builder.code = orm.load_code(args.code)
    builder.structure = structure
    builder.parameters = build_parameters()
    builder.metadata.label = "H2O PBE single point"
    builder.metadata.options.resources = {"num_machines": 1, "num_mpiprocs_per_machine": 4}
    builder.metadata.options.max_wallclock_seconds = 30 * 60

    if args.submit:
        node = submit(builder)
        print(f"submitted CalcJob pk={node.pk}; follow with: verdi process list -a")
        return

    results, node = run_get_node(builder)
    print(f"CalcJob pk={node.pk}  exit_status={node.exit_status}")
    out = results["output_parameters"].get_dict()
    energy_ha = out.get("energy")
    if energy_ha is not None:
        print(f"E(PBE) = {energy_ha:.8f} Ha = {energy_ha * 27.211386:.5f} eV")
    print("provenance graph: verdi node graph generate", node.pk)


if __name__ == "__main__":
    main()
