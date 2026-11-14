#!/usr/bin/env python3
"""Build a cheap NEB initial path with a machine-learned potential.

Workflow: DFT-relaxed initial and final states -> IDPP interpolation ->
CI-NEB with MACE-MP / DeePMD / NEP -> images written as 00/POSCAR ... NN/POSCAR
for a VASP+VTST CI-NEB (or converted to CP2K / OpenMX by yourself).
The MLIP barrier is only a guess; the DFT NEB decides.

Examples
  python mlip_neb_guess.py IS/CONTCAR FS/CONTCAR --nimages 5
  python mlip_neb_guess.py is.vasp fs.vasp --nimages 7 --calc mace-mp --device cuda --outdir neb_init
  python mlip_neb_guess.py is.vasp fs.vasp --no-neb            # IDPP only, no MLIP needed

Atom order must be identical in both end points (same as for nebmake.pl).
"""
import argparse
import os
import sys

import numpy as np
from ase.io import read, write
from ase.optimize import FIRE

from mlip_calc import KINDS, make_calculator

try:
    from ase.mep import NEB            # ASE >= 3.23
except ImportError:
    from ase.neb import NEB            # older ASE


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("initial")
    p.add_argument("final")
    p.add_argument("--nimages", type=int, default=5, help="intermediate images (IMAGES in VASP)")
    p.add_argument("--calc", default="mace-mp", choices=KINDS)
    p.add_argument("--model", default=None)
    p.add_argument("--device", default="cpu")
    p.add_argument("--dtype", default="float64", choices=["float32", "float64"])
    p.add_argument("--method", default="improvedtangent",
                   help="NEB tangent method; ASE >= 3.28 default is improvedtangent, older ASE used aseneb")
    p.add_argument("--fmax", type=float, default=0.1, help="eV/Angstrom, loose on purpose")
    p.add_argument("--steps", type=int, default=300)
    p.add_argument("--no-neb", action="store_true", help="only write the IDPP path")
    p.add_argument("--outdir", default="neb_init")
    a = p.parse_args()

    ini, fin = read(a.initial), read(a.final)
    if ini.get_chemical_symbols() != fin.get_chemical_symbols():
        sys.exit("ERROR: atom order/species differ between initial and final structure")
    images = [ini] + [ini.copy() for _ in range(a.nimages)] + [fin]

    neb = NEB(images, climb=not a.no_neb, method=a.method, allow_shared_calculator=True)
    neb.interpolate(method="idpp", mic=True)
    print(f"IDPP path with {a.nimages} intermediate images built")

    if not a.no_neb:
        calc = make_calculator(a.calc, a.model, a.device, a.dtype)
        for im in images:
            im.calc = calc
        opt = FIRE(neb, trajectory="neb_mlip.traj", logfile="-")
        ok = opt.run(fmax=a.fmax, steps=a.steps)
        energies = np.array([im.get_potential_energy() for im in images])
        rel = energies - energies[0]
        print("image  E - E_IS (eV, MLIP)")
        for i, e in enumerate(rel):
            print(f"{i:5d}  {e:10.4f}")
        print(f"MLIP forward barrier ~ {rel.max():.3f} eV at image {int(rel.argmax())} (guess only)")
        if not ok:
            print("WARNING: MLIP NEB not converged; the path may still be a useful start", file=sys.stderr)

    os.makedirs(a.outdir, exist_ok=True)
    for i, im in enumerate(images):
        d = os.path.join(a.outdir, f"{i:02d}")
        os.makedirs(d, exist_ok=True)
        im.calc = None
        write(os.path.join(d, "POSCAR"), im, format="vasp", direct=True, sort=False)
    print(f"written {a.outdir}/00 ... {a.outdir}/{len(images) - 1:02d}; copy the DFT OUTCARs of the end points into 00 and the last dir")


if __name__ == "__main__":
    main()
