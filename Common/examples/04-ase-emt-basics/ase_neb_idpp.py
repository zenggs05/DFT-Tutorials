#!/usr/bin/env python3
"""CI-NEB for a Cu adatom hopping between two hollow sites on Cu(100).

Uses EMT so it finishes in well under a minute.  Shows the full ASE NEB
workflow: relax end points -> IDPP interpolation -> NEB (explicit method)
-> climbing image -> barrier with NEBTools.

ASE >= 3.28 changed the default NEB method from 'aseneb' to
'improvedtangent'.  We always pass method= explicitly so that results do
not depend on the ASE version.

Usage:
    python ase_neb_idpp.py [--nimages 5] [--method improvedtangent]
"""
import argparse

from ase.build import add_adsorbate, fcc100
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.mep import NEB, NEBTools
from ase.optimize import BFGS, FIRE


def endpoint(offset):
    """Adatom in a hollow site; offset=(dx, dy) in A from the default hollow."""
    slab = fcc100("Cu", size=(3, 3, 3), vacuum=7.5)
    add_adsorbate(slab, "Cu", height=1.6, position="hollow", offset=offset)
    slab.set_constraint(FixAtoms(mask=[a.tag > 1 for a in slab]))  # only top layer + adatom move
    slab.calc = EMT()
    BFGS(slab, logfile=None).run(fmax=0.01)
    return slab


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--nimages", type=int, default=5, help="number of intermediate images")
    p.add_argument("--method", default="improvedtangent",
                   choices=["improvedtangent", "aseneb", "eb", "spline", "string"])
    args = p.parse_args()

    initial = endpoint((0, 0))
    final = endpoint((1, 0))   # offset in units of the surface lattice vector -> next hollow
    print(f"E(initial) = {initial.get_potential_energy():.4f} eV")
    print(f"E(final)   = {final.get_potential_energy():.4f} eV")

    images = [initial]
    for _ in range(args.nimages):
        img = initial.copy()
        img.calc = EMT()
        images.append(img)
    images.append(final)

    neb = NEB(images, k=0.1, climb=False, method=args.method)
    neb.interpolate(method="idpp")            # IDPP initial path (better than linear)
    write("neb_initial_path.traj", images)

    # Stage 1: plain NEB to get a reasonable path; stage 2: switch on the climbing image
    FIRE(neb, trajectory="neb.traj", logfile="neb.log").run(fmax=0.10)
    neb.climb = True
    FIRE(neb, trajectory="neb.traj", logfile="neb.log", append_trajectory=True).run(fmax=0.03)

    nt = NEBTools(images)
    barrier, de = nt.get_barrier()
    print(f"method = {args.method}")
    print(f"forward barrier = {barrier:.3f} eV   reaction energy = {de:.3f} eV")
    fig = nt.plot_band()
    fig.savefig("neb_barrier.png", dpi=200)
    write("neb_final_path.traj", images)
    print("written: neb_final_path.traj (ase gui neb_final_path.traj), neb_barrier.png")


if __name__ == "__main__":
    main()
