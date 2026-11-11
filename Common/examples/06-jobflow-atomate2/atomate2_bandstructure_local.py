#!/usr/bin/env python
"""atomate2: relax + static + uniform/line-mode band structure of Si, run locally.

What happens (RelaxBandStructureMaker):
    DoubleRelaxMaker (2 x "relax") -> "static" -> "non-scf uniform" + "non-scf line"
Every VASP step is run through custodian with atomate2's default handlers.

Prerequisites
  * Python >= 3.12, pip install atomate2   (0.1.x at the time of writing)
  * ~/.atomate2.yaml with VASP_CMD (see atomate2.yaml.example)
  * POTCARs configured for pymatgen:  pmg config --add PMG_VASP_PSP_DIR /path/to/psp
  * run inside a Slurm allocation (VASP is started by custodian)

Usage:
    python atomate2_bandstructure_local.py            # creates job_* folders in cwd
    python atomate2_bandstructure_local.py --dry-run  # only print the flow graph
"""
import argparse

from jobflow import run_locally
from pymatgen.core import Lattice, Structure

from atomate2.vasp.flows.core import RelaxBandStructureMaker
from atomate2.vasp.powerups import add_metadata_to_flow, update_user_incar_settings


def si_structure() -> Structure:
    """Primitive diamond Si (a = 5.43 A); replace by Structure.from_file('POSCAR')."""
    a_prim = 5.43 / 2**0.5                      # fcc primitive vector length
    lattice = Lattice.from_parameters(a_prim, a_prim, a_prim, 60, 60, 60)
    return Structure(lattice, ["Si", "Si"], [[0, 0, 0], [0.25, 0.25, 0.25]])


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dry-run", action="store_true", help="print the jobs and exit")
    p.add_argument("--ncore", type=int, default=4, help="NCORE written to every INCAR")
    args = p.parse_args()

    structure = si_structure()
    flow = RelaxBandStructureMaker().make(structure)

    # Powerups modify all VASP makers inside the flow (returns a modified copy)
    flow = update_user_incar_settings(flow, {"NCORE": args.ncore, "LWAVE": False})
    flow = add_metadata_to_flow(flow, {"project": "dft-tutorial", "system": "Si"})

    if args.dry_run:
        for job, parents in flow.iterflow():
            print(f"{job.name:28s} <- {len(parents)} parent job(s)")
        return

    # MemoryStore is used unless ~/.jobflow.yaml defines a MongoDB JOB_STORE
    responses = run_locally(flow, create_folders=True, ensure_success=True)

    for job, _ in flow.iterflow():
        resp = responses.get(job.uuid, {}).get(1)
        if resp is None:
            continue
        doc = resp.output                           # emmet TaskDoc
        out = getattr(doc, "output", None)
        gap = getattr(out, "bandgap", None)
        energy = getattr(out, "energy", None)
        print(f"{job.name:28s} E = {energy}  gap = {gap}  dir = {getattr(doc, 'dir_name', '')}")


if __name__ == "__main__":
    main()
