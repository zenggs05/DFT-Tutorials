#!/usr/bin/env python
"""AiiDA + aiida-gaussian: B3LYP/6-31G(d) optimisation of H2O with GaussianBaseWorkChain.

GaussianBaseWorkChain wraps GaussianCalculation in a restart loop.  Its process
handlers react to known failures, e.g. an SCF convergence failure is retried
with SCF=YQC and then SCF=XQC; the optional ASyTop handler rounds coordinates.
Every attempt is a separate CalcJob node in the provenance graph.

Prerequisites
  pip install aiida-core aiida-gaussian      (aiida-gaussian 2.x parses with cclib)
  verdi presto
  verdi code create core.code.installed --config gaussian_code.yaml
Usage
  verdi run aiida_gaussian_base.py           # blocking run
  python aiida_gaussian_base.py --submit     # through the daemon
"""
import argparse

import ase.build
from aiida import load_profile, orm
from aiida.engine import run_get_node, submit
from aiida.plugins import WorkflowFactory

NPROC = 8          # cores used by Gaussian (%nprocshared)
MEM_MB = 4000      # %mem in MB


def build_parameters() -> orm.Dict:
    """Gaussian input as a dict; aiida-gaussian writes aiida.inp from it."""
    return orm.Dict(
        {
            "link0_parameters": {
                "%chk": "aiida.chk",
                "%mem": f"{MEM_MB}MB",
                "%nprocshared": NPROC,
            },
            "functional": "B3LYP",
            "basis_set": "6-31G(d)",
            "charge": 0,
            "multiplicity": 1,
            "route_parameters": {
                "opt": None,                                # keyword without value
                "scf": {"conver": 8, "maxcycle": 128},      # -> scf=(conver=8,maxcycle=128)
            },
        }
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--code", default="g16@localhost")
    p.add_argument("--submit", action="store_true")
    args = p.parse_args()

    load_profile()
    GaussianBaseWorkChain = WorkflowFactory("gaussian.base")

    builder = GaussianBaseWorkChain.get_builder()
    builder.gaussian.code = orm.load_code(args.code)
    builder.gaussian.structure = orm.StructureData(ase=ase.build.molecule("H2O"))
    builder.gaussian.parameters = build_parameters()
    builder.gaussian.metadata.options.resources = {"num_machines": 1, "tot_num_mpiprocs": NPROC}
    # ask the scheduler for 25 % more memory than %mem
    builder.gaussian.metadata.options.max_memory_kb = int(1.25 * MEM_MB) * 1024
    builder.gaussian.metadata.options.max_wallclock_seconds = 2 * 3600
    builder.metadata.label = "H2O B3LYP/6-31G(d) opt"

    if args.submit:
        node = submit(builder)
        print(f"submitted WorkChain pk={node.pk}; verdi process report {node.pk}")
        return

    results, node = run_get_node(builder)
    print(f"WorkChain pk={node.pk} exit_status={node.exit_status}")
    out = results["output_parameters"].get_dict()
    scf = out.get("scfenergies")
    if scf:
        print(f"final SCF energy = {scf[-1]:.6f} eV (cclib units)")


if __name__ == "__main__":
    main()
