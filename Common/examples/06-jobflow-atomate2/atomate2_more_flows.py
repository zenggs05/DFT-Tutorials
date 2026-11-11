#!/usr/bin/env python
"""Build (and optionally run) atomate2 PhononMaker / VaspLobsterMaker flows.

Without --run the script only prints the job graph, so it is a cheap way to see
what a flow will do before spending CPU hours.

    pip install "atomate2[phonons,lobster]"     # phonopy, seekpath, lobsterpy, ijson
    python atomate2_more_flows.py phonon         # Si phonons: list the jobs
    python atomate2_more_flows.py lobster        # GaAs COHP: list the jobs
    python atomate2_more_flows.py phonon --run   # run locally (inside a Slurm job)
"""
import argparse

from jobflow import run_locally
from pymatgen.core import Lattice, Structure


def fcc_prim(a: float, species: list) -> Structure:
    """Diamond / zinc-blende primitive cell with conventional lattice constant a."""
    b = a / 2**0.5
    lat = Lattice.from_parameters(b, b, b, 60, 60, 60)
    return Structure(lat, species, [[0, 0, 0], [0.25, 0.25, 0.25]])


def phonon_flow():
    from atomate2.vasp.flows.phonons import PhononMaker

    si = fcc_prim(5.43, ["Si", "Si"])
    maker = PhononMaker(
        min_length=15.0,          # supercell edges >= 15 A (phonopy supercell search)
        displacement=0.01,        # finite displacement in A
        born_maker=None,          # Si is non-polar: skip the DFPT Born-charge job
        store_force_constants=False,
    )
    return maker.make(si)


def lobster_flow():
    from atomate2.vasp.flows.lobster import VaspLobsterMaker

    gaas = fcc_prim(5.65, ["Ga", "As"])
    # relax -> static with ISYM = 0 and enough bands -> LOBSTER for several basis sets
    return VaspLobsterMaker().make(gaas)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("which", choices=["phonon", "lobster"])
    p.add_argument("--run", action="store_true", help="run with run_locally()")
    args = p.parse_args()

    flow = phonon_flow() if args.which == "phonon" else lobster_flow()
    for job, parents in flow.iterflow():
        print(f"{job.name:40s} parents={len(parents)}")
    if args.run:
        run_locally(flow, create_folders=True, ensure_success=True)


if __name__ == "__main__":
    main()
