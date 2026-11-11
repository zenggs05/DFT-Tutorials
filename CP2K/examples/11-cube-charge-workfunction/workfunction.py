#!/usr/bin/env python3
"""Work function from a CP2K V_HARTREE cube:  Phi = V_vac - E_F.

  * planar (xy) average of the Hartree potential along z, hartree -> eV
  * E_F is read from the last "Fermi energy:" line of the CP2K output (hartree)
  * V_vac is taken as the mean potential in a window around the lowest-density
    region of the vacuum; with a surface dipole correction the two sides of an
    asymmetric slab differ, both plateaus are reported.

Usage:  python workfunction.py pt111_wf-v_hartree-1_0.cube pt111_workfunction.out
Writes  vh_planar.dat  (z [Angstrom], V(z) [eV], V(z) - E_F [eV]).
Needs numpy only. Assumes the c axis is parallel to z (slab normal along z).
"""
import re
import sys

import numpy as np

HA2EV = 27.211386245988
BOHR = 0.529177210903


def read_cube(fname):
    with open(fname) as f:
        lines = f.readlines()
    natoms = abs(int(lines[2].split()[0]))
    grid = [lines[3 + i].split() for i in range(3)]
    n = [int(g[0]) for g in grid]
    vec = np.array([[float(x) for x in g[1:4]] for g in grid])
    atoms_z = np.array([float(lines[6 + i].split()[4]) for i in range(natoms)]) * BOHR
    data = np.array(" ".join(lines[6 + natoms:]).split(), dtype=float).reshape(n)
    return n, vec, atoms_z, data


def fermi_energy(outfile):
    ef = None
    with open(outfile) as f:
        for line in f:
            m = re.search(r"Fermi energy:\s+(-?\d+\.\d+)", line)
            if m:
                ef = float(m.group(1))
    if ef is None:
        raise SystemExit("no 'Fermi energy:' line found (smearing switched off?)")
    return ef * HA2EV


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    n, vec, atoms_z, v = read_cube(sys.argv[1])
    ef = fermi_energy(sys.argv[2])
    dz = np.linalg.norm(vec[2]) * BOHR
    cz = dz * n[2]
    z = np.arange(n[2]) * dz
    vz = v.mean(axis=(0, 1)) * HA2EV          # CP2K sign: electron potential energy
    np.savetxt("vh_planar.dat", np.column_stack([z, vz, vz - ef]),
               header="z[Angstrom]  V_H(z)[eV]  V_H(z)-E_F[eV]", fmt="%12.6f")

    # vacuum region = largest gap between atomic z positions (periodic in z)
    za = np.sort(np.mod(atoms_z, cz))
    gaps = np.diff(np.append(za, za[0] + cz))
    i = int(np.argmax(gaps))
    lo, hi = za[i], za[i] + gaps[i]
    print("E_F = %.4f eV, vacuum from z = %.2f to %.2f Angstrom (%.1f A thick)"
          % (ef, lo, hi % cz, gaps[i]))
    # two windows of 2 Angstrom, 4 Angstrom away from the outermost atoms
    for label, (a, b) in (("upper side", (lo + 4.0, lo + 6.0)), ("lower side", (hi - 6.0, hi - 4.0))):
        sel = ((z - a) % cz) <= (b - a)
        vvac = vz[sel].mean()
        print("%-10s  V_vac = %8.4f eV   Phi = V_vac - E_F = %.3f eV   (spread %.3f eV)"
              % (label, vvac, vvac - ef, vz[sel].max() - vz[sel].min()))
    print("wrote vh_planar.dat; check the plateau by plotting columns 1 and 3")


if __name__ == "__main__":
    main()
