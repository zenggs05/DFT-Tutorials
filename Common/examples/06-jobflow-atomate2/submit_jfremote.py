#!/usr/bin/env python
"""Submit an atomate2 flow to the jobflow-remote queue (it does NOT run it here).

The runner daemon ("jf runner start") picks the jobs up, uploads inputs to the
worker, submits them to Slurm, downloads outputs and stores the results.

Usage:
    python submit_jfremote.py POSCAR              # relax + band structure
    python submit_jfremote.py POSCAR --phonon     # phonopy workflow instead
    jf job list                                   # follow the progress
    jf flow list
"""
import argparse

from jobflow_remote import submit_flow
from pymatgen.core import Structure

from atomate2.vasp.flows.core import RelaxBandStructureMaker
from atomate2.vasp.powerups import update_user_incar_settings


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("structure", help="POSCAR / CIF file")
    p.add_argument("--phonon", action="store_true", help="submit PhononMaker instead")
    p.add_argument("--worker", default="hpc_cpu", help="worker name in the project file")
    args = p.parse_args()

    structure = Structure.from_file(args.structure)
    if args.phonon:
        # needs: pip install "atomate2[phonons]"  (phonopy + seekpath)
        from atomate2.vasp.flows.phonons import PhononMaker

        # min_length: supercell edge >= 15 A; displacement in Angstrom
        flow = PhononMaker(min_length=15.0, displacement=0.01).make(structure)
    else:
        flow = RelaxBandStructureMaker().make(structure)
    flow = update_user_incar_settings(flow, {"NCORE": 8})

    ids = submit_flow(flow, worker=args.worker, exec_config="vasp6")
    print(f"submitted {len(ids)} job(s): {ids}")


if __name__ == "__main__":
    main()
