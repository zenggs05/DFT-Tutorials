#!/usr/bin/env python
"""quacc recipes in 20 lines: relax a Cu(111) slab and its bulk with EMT.

quacc recipes are plain functions decorated with @job/@flow.  With no workflow
engine configured (the default, QUACC_WORKFLOW_ENGINE unset) they run eagerly,
so this script needs no database and no scheduler.  Swap the import to
quacc.recipes.vasp.core (or .gaussian.core) and the same call runs VASP/Gaussian.

Usage:
    pip install quacc              # quacc 1.x, Python >= 3.11
    python quacc_emt_demo.py
"""
from ase.build import bulk, fcc111

from quacc.recipes.emt.core import relax_job, static_job


def main() -> None:
    cu = bulk("Cu", "fcc", a=3.6)
    bulk_res = relax_job(cu, relax_cell=True)        # BFGS + cell filter
    e_bulk = bulk_res["results"]["energy"] / len(bulk_res["atoms"])
    print(f"bulk: a = {bulk_res['atoms'].cell.cellpar()[0] * 2**0.5:.4f} A,"
          f" E/atom = {e_bulk:.4f} eV  (dir: {bulk_res['dir_name']})")

    a_relaxed = bulk_res["atoms"].cell.cellpar()[0] * 2**0.5   # fcc primitive -> cubic a
    slab = fcc111("Cu", size=(2, 2, 4), vacuum=8.0, a=a_relaxed)  # consistent with e_bulk
    slab_res = relax_job(slab, opt_params={"fmax": 0.02})
    single = static_job(slab_res["atoms"])
    n = len(single["atoms"])
    area = slab_res["atoms"].cell.area(2)
    gamma = (single["results"]["energy"] - n * e_bulk) / (2 * area)
    print(f"Cu(111) surface energy (EMT, illustrative only): {gamma:.4f} eV/A^2")


if __name__ == "__main__":
    main()
