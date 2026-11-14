#!/usr/bin/env python
"""AiiDA + aiida-vasp (5.x): relax bulk Si with VaspRelaxWorkChain from a protocol.

The "balanced" protocol of aiida-vasp sets PREC/ENCUT/EDIFF/k-spacing (and GGA = PS,
i.e. PBEsol, which is overridden to PBE below), restarts
the relaxation until the volume change is below 1 %, and finishes with a static
run.  Every intermediate VASP run is a node in the provenance graph.

Prerequisites
  pip install aiida-core aiida-vasp
  verdi presto  (or a PostgreSQL profile for production)
  aiida-vasp potcar uploadfamily --path /path/to/potpaw_PBE.54.tgz \
        --name PBE.54 --description "PBE.54 from the VASP portal"
  verdi code create core.code.installed --config vasp_code.yaml
Usage
  verdi run aiida_vasp_relax.py               # blocking
  python aiida_vasp_relax.py --submit         # through the daemon
"""
import argparse

from aiida import load_profile, orm
from aiida.engine import run_get_node, submit
from aiida.plugins import WorkflowFactory
from ase.build import bulk


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--code", default="vasp_std@localhost")
    p.add_argument("--submit", action="store_true")
    args = p.parse_args()

    load_profile()
    VaspRelaxWorkChain = WorkflowFactory("vasp.relax")

    structure = orm.StructureData(ase=bulk("Si", "diamond", a=5.43))

    builder = VaspRelaxWorkChain.get_builder_from_protocol(
        code=orm.load_code(args.code),
        structure=structure,
        protocol="balanced",
        # Same nesting as the input namespaces; merged on top of the protocol
        # NOTE: the "balanced" protocol of aiida-vasp 5.1 uses GGA = PS (PBEsol);
        # "gga": "pe" switches back to PBE, consistent with the PBE.54 POTCARs.
        overrides={"vasp": {"parameters": {"incar": {"encut": 400, "ispin": 1, "ncore": 4,
                                                     "gga": "pe"}}}},
        relax_settings={"force_cutoff": 0.01, "steps": 40},
        options={
            "resources": {"num_machines": 1, "num_mpiprocs_per_machine": 32},
            "max_wallclock_seconds": 3600,
            "queue_name": "cpu",           # Slurm partition
        },
    )
    builder.metadata.label = "Si relax (balanced)"

    if args.submit:
        node = submit(builder)
        print(f"submitted WorkChain pk={node.pk}; verdi process status {node.pk}")
        return

    results, node = run_get_node(builder)
    print(f"WorkChain pk={node.pk} exit_status={node.exit_status}")
    if "relax" in results and "structure" in results["relax"]:
        cell = results["relax"]["structure"].get_ase().cell.cellpar()
        print(f"relaxed |a1| = {cell[0]:.4f} A  (primitive cell vector)")
    if "misc" in results:
        print("total_energies:", results["misc"].get_dict().get("total_energies"))


if __name__ == "__main__":
    main()
