#!/usr/bin/env python3
"""Minimal atomate2 workflow: double relaxation + uniform and line-mode band structures.

The flow is executed locally with jobflow (one job after the other, each in its own folder).
Configure once:
  * ~/.atomate2.yaml   ->  VASP_CMD: mpirun -np 32 vasp_std   (or env ATOMATE2_VASP_CMD)
  * pymatgen POTCAR   ->  pmg config -p <unzipped_potcars> <pmg_potcars>
                          pmg config --add PMG_VASP_PSP_DIR <pmg_potcars>
Run this script inside a batch job (it calls VASP several times).
For many structures or remote execution use jobflow-remote or FireWorks instead (see chapter 6 of Common).

Requires: pip install atomate2 (pulls in jobflow, pymatgen, custodian)
"""
import argparse

from jobflow import run_locally
from pymatgen.core import Structure

from atomate2.vasp.flows.core import RelaxBandStructureMaker


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("structure", nargs="?", default="POSCAR", help="POSCAR / CIF of the material")
    args = p.parse_args()

    structure = Structure.from_file(args.structure)
    flow = RelaxBandStructureMaker().make(structure)
    # create_folders=True: every job runs in its own job_<uuid> directory
    responses = run_locally(flow, create_folders=True)
    print("finished %d jobs" % len(responses))


if __name__ == "__main__":
    main()
