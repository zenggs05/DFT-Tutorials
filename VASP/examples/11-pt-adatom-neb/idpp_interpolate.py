#!/usr/bin/env python3
"""Build the NEB image directories 00 ... NN for VASP (linear or IDPP path).

    python idpp_interpolate.py POSCAR_IS POSCAR_FS -n 3            # IDPP (default)
    python idpp_interpolate.py IS/CONTCAR FS/CONTCAR -n 5 --method linear
    python idpp_interpolate.py POSCAR_IS POSCAR_FS -n 3 --check-only

What it does
  1. reads the relaxed initial (IS) and final (FS) structures,
  2. checks that both have the same cell, the same element order and the same
     selective-dynamics flags, and reports the atoms that move most
     (a jump > --warn Angstrom usually means the atom order is not consistent),
  3. interpolates N intermediate images with the minimum-image convention
     (atoms that cross a periodic boundary are handled correctly),
  4. writes 00/POSCAR ... NN/POSCAR (NN = N+1) and prints the distance between
     neighbouring images (like dist.pl of the VTST scripts).

The image-dependent pair potential (IDPP) method of Smidstrup et al.,
J. Chem. Phys. 140, 214106 (2014) usually gives a much better starting path
than linear interpolation when several atoms move.

Requires ASE (>= 3.23 uses ase.mep; older versions ase.neb are also accepted).
"""
import argparse
import os
import sys

import numpy as np

try:
    from ase.io import read, write
    from ase.constraints import FixAtoms, FixScaled
except ImportError:
    sys.exit("ERROR: this script needs ASE (pip install ase)")
try:                                   # ASE >= 3.23
    from ase.mep import NEB
except ImportError:                    # older ASE
    from ase.neb import NEB


def fixed_mask(atoms):
    """Boolean array: True where an atom has any frozen Cartesian component."""
    m = np.zeros(len(atoms), dtype=bool)
    for c in atoms.constraints:
        if isinstance(c, FixAtoms):
            m[c.index] = True
        elif isinstance(c, FixScaled):
            m[c.index if np.ndim(c.index) else [c.index]] = True
    return m


def mic_displacement(a, b):
    """Displacement b - a for every atom with the minimum-image convention."""
    d = b.get_positions() - a.get_positions()
    frac = np.linalg.solve(a.cell.T, d.T).T
    pbc = a.get_pbc()
    frac[:, pbc] -= np.round(frac[:, pbc])
    return frac @ a.cell.array


def check_endpoints(ini, fin, warn):
    ok = True
    if len(ini) != len(fin):
        print("ERROR: different number of atoms (%d vs %d)" % (len(ini), len(fin)))
        return False
    if ini.get_chemical_symbols() != fin.get_chemical_symbols():
        print("ERROR: element order differs between IS and FS")
        ok = False
    if not np.allclose(ini.cell.array, fin.cell.array, atol=1e-4):
        print("WARNING: IS and FS cells differ - NEB in VASP keeps the cell fixed")
    if not np.array_equal(fixed_mask(ini), fixed_mask(fin)):
        print("WARNING: selective-dynamics flags differ between IS and FS")
    d = mic_displacement(ini, fin)
    dn = np.linalg.norm(d, axis=1)
    order = np.argsort(dn)[::-1]
    print("total IS->FS distance (sqrt of sum of squares): %.3f A" % np.sqrt((dn ** 2).sum()))
    print("atoms that move most (1-based index, element, |d| in A):")
    for k in order[:min(5, len(order))]:
        print("   %4d %-3s %7.3f" % (k + 1, ini[k].symbol, dn[k]))
    if dn.max() > warn:
        print("WARNING: an atom moves more than %.1f A - check that the atom "
              "order of IS and FS is identical" % warn)
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("initial")
    ap.add_argument("final")
    ap.add_argument("-n", "--nimages", type=int, default=3,
                    help="number of intermediate images (= IMAGES in INCAR)")
    ap.add_argument("--method", choices=["idpp", "linear"], default="idpp")
    ap.add_argument("--warn", type=float, default=3.0,
                    help="warn if one atom moves more than this (A)")
    ap.add_argument("--check-only", action="store_true",
                    help="only compare IS and FS, write nothing")
    args = ap.parse_args()

    ini = read(args.initial, format="vasp")
    fin = read(args.final, format="vasp")
    if not check_endpoints(ini, fin, args.warn):
        sys.exit(1)
    if args.check_only:
        return

    # put FS atoms on the periodic image closest to IS (avoids crossing the cell)
    fin = fin.copy()
    fin.set_positions(ini.get_positions() + mic_displacement(ini, fin))

    images = [ini] + [ini.copy() for _ in range(args.nimages)] + [fin]
    try:                               # silence the default-method warning of new ASE
        neb = NEB(images, method="improvedtangent")
    except TypeError:
        neb = NEB(images)
    method = "linear" if args.method == "linear" else "idpp"
    try:      # interpolate positions freely; IS/FS already satisfy the constraints
        neb.interpolate(method=method, mic=True, apply_constraint=False)
    except TypeError:                  # very old ASE without apply_constraint
        neb.interpolate(method=method, mic=True)

    nn = args.nimages + 1
    for i, img in enumerate(images):
        d = "%02d" % i
        os.makedirs(d, exist_ok=True)
        img.set_constraint(ini.constraints)
        write(os.path.join(d, "POSCAR"), img, format="vasp", direct=True)
    print("wrote %02d ... %02d (%s interpolation, IMAGES = %d)"
          % (0, nn, args.method, args.nimages))
    print("distance between neighbouring images (A):")
    for i in range(nn):
        dn = np.linalg.norm(mic_displacement(images[i], images[i + 1]), axis=1)
        print("   %02d -> %02d : %6.3f   (largest single-atom move %.3f)"
              % (i, i + 1, np.sqrt((dn ** 2).sum()), dn.max()))
    # shortest interatomic distance per image: catches overlapping atoms
    for i, img in enumerate(images[1:-1], 1):
        dm = img.get_all_distances(mic=True)
        np.fill_diagonal(dm, 1e9)
        if dm.min() < 1.0:
            print("WARNING: image %02d has atoms closer than 1.0 A (%.2f A)" % (i, dm.min()))
    print("next: copy the OUTCAR of the relaxed IS into 00/ and of the FS into %02d/" % nn)


if __name__ == "__main__":
    main()
