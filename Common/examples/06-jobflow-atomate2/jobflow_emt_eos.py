#!/usr/bin/env python
"""A first jobflow workflow that runs anywhere (no DFT code needed).

Flow:  N x energy_at(a)  -->  fit_eos([...])      (fan-in, like a convergence scan)
The "calculator" is ASE's EMT potential, so the whole flow takes seconds.
Replace energy_at() by a VASP/CP2K job and the structure of the flow stays the same.

Usage:
    pip install jobflow ase
    python jobflow_emt_eos.py              # results kept in memory (MemoryStore)
    python jobflow_emt_eos.py --folders    # one job_<time>-<id> folder per job
"""
import argparse

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.eos import EquationOfState
from ase.units import kJ
from jobflow import Flow, job, run_locally


@job
def energy_at(a: float, element: str = "Cu") -> dict:
    """Single-point energy of an fcc primitive cell with lattice constant a (Angstrom)."""
    atoms = bulk(element, "fcc", a=a)
    atoms.calc = EMT()
    return {"a": a, "volume": atoms.get_volume(), "energy": atoms.get_potential_energy()}


@job
def fit_eos(points: list) -> dict:
    """Birch-Murnaghan fit of E(V); 'points' arrives as a list of resolved dicts."""
    vols = np.array([p["volume"] for p in points])
    ens = np.array([p["energy"] for p in points])
    eos = EquationOfState(vols, ens, eos="birchmurnaghan")
    v0, e0, b = eos.fit()
    a0 = (4.0 * v0) ** (1.0 / 3.0)          # primitive fcc: V = a^3 / 4
    return {"a0": a0, "v0": v0, "e0": e0, "B_GPa": b / kJ * 1.0e24}


def build_flow(element: str = "Cu", a_guess: float = 3.6) -> Flow:
    scan = [energy_at(a, element) for a in np.linspace(0.95, 1.05, 7) * a_guess]
    for j in scan:
        j.name = f"emt a={j.function_args[0]:.3f}"
    fit = fit_eos([j.output for j in scan])  # OutputReferences, resolved at run time
    fit.name = "eos fit"
    return Flow([*scan, fit], output=fit.output, name=f"{element} EMT EOS")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--folders", action="store_true", help="create one folder per job")
    args = p.parse_args()

    flow = build_flow()
    responses = run_locally(flow, create_folders=args.folders, ensure_success=True)

    fit_uuid = flow.jobs[-1].uuid
    res = responses[fit_uuid][1].output        # [uuid][job index] -> Response
    print(f"a0 = {res['a0']:.4f} A   B = {res['B_GPa']:.1f} GPa   (EMT, not DFT)")


if __name__ == "__main__":
    main()
