#!/usr/bin/env python3
"""Turn a molecule (.xyz) into starting inputs for Gaussian, CP2K, VASP
and OpenMX, and check that charge and multiplicity are consistent.

  python xyz_to_inputs.py h2o.xyz
  python xyz_to_inputs.py paracetamol.xyz --task "opt freq" --nproc 16
  python xyz_to_inputs.py nh4.xyz --charge 1 --padding 8

charge/mult are taken from the xyz comment line if it contains
"charge=<q> mult=<m>" (as written by smiles_to_3d.py), else from the
command line (default 0 / 1).

Outputs (prefix = xyz file name):
  <prefix>.gjf              Gaussian input (Link0 + route + geometry)
  <prefix>_cp2k_subsys.inc  CP2K &CELL/&COORD/&TOPOLOGY for an isolated
                            molecule (PERIODIC NONE)
  POSCAR_<prefix>_box       molecule centred in a cubic box (VASP,
                            use a Gamma-only KPOINTS and vasp_gam)
  <prefix>_openmx_atoms.dat OpenMX Atoms.* blocks, same box, Angstrom
"""
import argparse
import os
import re
import sys

import numpy as np
from ase.data import atomic_numbers
from ase.io import read, write

VALENCE_PBE19 = {"H": 1.0, "C": 4.0, "N": 5.0, "O": 6.0, "S": 6.0,
                 "Si": 4.0}


def parse_comment(fn):
    with open(fn) as f:
        f.readline()
        comment = f.readline()
    q = re.search(r"charge=(-?\d+)", comment)
    m = re.search(r"mult=(\d+)", comment)
    return (int(q.group(1)) if q else None, int(m.group(1)) if m else None)


def check_spin(atoms, charge, mult):
    nelec = sum(atomic_numbers[s] for s in atoms.get_chemical_symbols())
    nelec -= charge
    ok = (nelec % 2 == 0 and mult % 2 == 1) or \
         (nelec % 2 == 1 and mult % 2 == 0)
    return nelec, ok


def write_gjf(fn, atoms, args, charge, mult, title):
    base = os.path.splitext(os.path.basename(fn))[0]
    with open(fn, "w", newline="\n") as f:
        f.write("%%nprocshared=%d\n%%mem=%s\n%%chk=%s.chk\n"
                % (args.nproc, args.mem, base))
        f.write("#p %s %s/%s %s\n\n" % (args.task, args.method, args.basis,
                                       args.extra))
        f.write("%s\n\n%d %d\n" % (title, charge, mult))
        for s, (x, y, z) in zip(atoms.get_chemical_symbols(),
                                atoms.positions):
            f.write("%-2s %13.6f %13.6f %13.6f\n" % (s, x, y, z))
        f.write("\n")      # a blank line terminates the molecule spec


def boxed(atoms, padding):
    """Cubic box: largest extent + 2*padding, molecule at the centre."""
    mol = atoms.copy()
    extent = np.ptp(mol.positions, axis=0).max()
    L = extent + 2.0 * padding
    mol.set_cell([L, L, L])
    mol.center()
    mol.pbc = True
    return mol, extent, L


def write_cp2k(fn, mol, extent):
    with open(fn, "w", newline="\n") as f:
        f.write("# CP2K &SUBSYS content for an isolated molecule\n")
        f.write("# Pair it with &DFT/&POISSON: PERIODIC NONE and PSOLVER MT\n"
                "# (Martyna-Tuckerman needs a cell roughly twice the size of\n"
                "# the molecule; here extent = %.2f A) or PSOLVER WAVELET.\n"
                % extent)
        f.write("&CELL\n  ABC [angstrom] %.4f %.4f %.4f\n  PERIODIC NONE\n"
                "&END CELL\n" % tuple(mol.cell.lengths()))
        f.write("&COORD\n")
        for s, (x, y, z) in zip(mol.get_chemical_symbols(), mol.positions):
            f.write("  %-2s %13.6f %13.6f %13.6f\n" % (s, x, y, z))
        f.write("&END COORD\n")
        f.write("&TOPOLOGY\n  &CENTER_COORDINATES\n  &END CENTER_COORDINATES\n"
                "&END TOPOLOGY\n")


def write_openmx(fn, mol):
    syms = mol.get_chemical_symbols()
    with open(fn, "w", newline="\n") as f:
        f.write("# OpenMX blocks; species names must match "
                "<Definition.of.Atomic.Species\n")
        f.write("Atoms.Number  %d\n" % len(mol))
        f.write("Atoms.SpeciesAndCoordinates.Unit  Ang\n")
        f.write("<Atoms.SpeciesAndCoordinates\n")
        for i, (s, (x, y, z)) in enumerate(zip(syms, mol.positions), 1):
            v = VALENCE_PBE19.get(s)
            occ = ("%5.2f  %5.2f" % (v / 2, v / 2)) if v else "  ?.?    ?.?"
            f.write("  %3d  %-2s %12.6f %12.6f %12.6f  %s\n"
                    % (i, s, x, y, z, occ))
        f.write("Atoms.SpeciesAndCoordinates>\n")
        f.write("Atoms.UnitVectors.Unit  Ang\n<Atoms.UnitVectors\n")
        for v in mol.cell:
            f.write("  %12.6f %12.6f %12.6f\n" % tuple(v))
        f.write("Atoms.UnitVectors>\n")
        f.write("# isolated molecule: scf.Kgrid 1 1 1 (Gamma only)\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("xyz")
    p.add_argument("--charge", type=int, default=None)
    p.add_argument("--mult", type=int, default=None)
    p.add_argument("--method", default="B3LYP")
    p.add_argument("--basis", default="def2SVP")
    p.add_argument("--task", default="opt freq")
    p.add_argument("--extra", default="EmpiricalDispersion=GD3BJ")
    p.add_argument("--nproc", type=int, default=8)
    p.add_argument("--mem", default="8GB")
    p.add_argument("--padding", type=float, default=7.5,
                   help="vacuum on each side of the box (A)")
    args = p.parse_args()

    atoms = read(args.xyz)
    q0, m0 = parse_comment(args.xyz)
    charge = args.charge if args.charge is not None else (q0 or 0)
    mult = args.mult if args.mult is not None else (m0 or 1)
    nelec, ok = check_spin(atoms, charge, mult)
    print("%s: %s, %d atoms, charge %d, multiplicity %d, %d electrons"
          % (args.xyz, atoms.get_chemical_formula(), len(atoms), charge,
             mult, nelec))
    if not ok:
        sys.exit("ERROR: multiplicity %d is impossible with %d electrons "
                 "(Gaussian would stop with 'The combination of "
                 "multiplicity ... and ... electrons is impossible')."
                 % (mult, nelec))

    prefix = os.path.splitext(os.path.basename(args.xyz))[0]
    write_gjf(prefix + ".gjf", atoms, args, charge, mult,
              "%s generated by xyz_to_inputs.py" % prefix)
    mol, extent, L = boxed(atoms, args.padding)
    write_cp2k(prefix + "_cp2k_subsys.inc", mol, extent)
    with open("POSCAR_%s_box" % prefix, "w", newline="\n") as f:
        write(f, mol, format="vasp", direct=False, sort=True)
    write_openmx(prefix + "_openmx_atoms.dat", mol)
    print("molecule extent %.2f A -> cubic box %.2f A (vacuum >= %.1f A)"
          % (extent, L, 2 * args.padding))
    print("written: %s.gjf %s_cp2k_subsys.inc POSCAR_%s_box "
          "%s_openmx_atoms.dat" % (prefix, prefix, prefix, prefix))


if __name__ == "__main__":
    main()
