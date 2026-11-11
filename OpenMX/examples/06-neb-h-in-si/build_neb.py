#!/usr/bin/env python3
"""Build an OpenMX NEB input from two relaxed end-point files (*.dat#).

The coordinate block of the initial-state file replaces <Atoms.SpeciesAndCoordinates
of the NEB template, and the coordinate block of the final-state file replaces
<NEB.Atoms.SpeciesAndCoordinates. Checks performed:
  * both end points use Ang or AU (FRAC breaks NEB, see OpenMX manual 51.8)
  * same number of atoms and same species order
  * minimum-image continuity: if an atom of the final state is more than half a
    lattice vector away from its initial position (wrapped through the boundary),
    it is shifted back by a lattice vector and a warning is printed

Usage (defaults shown):
  python build_neb.py --ini H_Si8_T1_opt.dat# --fin H_Si8_T2_opt.dat# \
                      --template H_diffusion_NEB.dat --out H_diffusion_NEB_relaxed.dat
Only numpy is required.
"""
import argparse
import re
import sys

import numpy as np

BOHR = 0.529177210903  # Angstrom


def find_block(text, name):
    pat = re.compile(r"^[ \t]*<%s[ \t]*\n(.*?)^[ \t]*%s>[ \t]*$" % (re.escape(name), re.escape(name)),
                     re.M | re.S | re.I)
    m = pat.search(text)
    if not m:
        sys.exit("ERROR: block <%s ... %s> not found" % (name, name))
    return m


def keyword(text, name):
    m = re.search(r"^[ \t]*%s[ \t]+(\S+)" % re.escape(name), text, re.M | re.I)
    return m.group(1) if m else None


def parse_atoms(block):
    atoms = []
    for line in block.splitlines():
        tok = line.split("#", 1)[0].split()
        if len(tok) < 5:
            continue
        atoms.append({"idx": tok[0], "sp": tok[1],
                      "xyz": np.array([float(t) for t in tok[2:5]]), "rest": tok[5:]})
    return atoms


def format_atoms(atoms):
    out = []
    for a in atoms:
        out.append("  %-3s %-4s %12.6f %12.6f %12.6f   %s" % (a["idx"], a["sp"], a["xyz"][0], a["xyz"][1],
                                                           a["xyz"][2], "  ".join(a["rest"])))
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ini", default="H_Si8_T1_opt.dat#", help="relaxed initial state (*.dat#)")
    ap.add_argument("--fin", default="H_Si8_T2_opt.dat#", help="relaxed final state (*.dat#)")
    ap.add_argument("--template", default="H_diffusion_NEB.dat")
    ap.add_argument("--out", default="H_diffusion_NEB_relaxed.dat")
    args = ap.parse_args()

    texts = {}
    for key in ("ini", "fin", "template"):
        with open(getattr(args, key)) as fh:
            texts[key] = fh.read()

    units = {k: (keyword(texts[k], "Atoms.SpeciesAndCoordinates.Unit") or "").lower() for k in texts}
    for k in ("ini", "fin", "template"):
        if units[k] not in ("ang", "au"):
            sys.exit("ERROR: %s uses coordinate unit '%s'; NEB needs Ang or AU" % (k, units[k]))
    if not (units["ini"] == units["fin"] == units["template"]):
        sys.exit("ERROR: coordinate units differ: %s" % units)

    ini = parse_atoms(find_block(texts["ini"], "Atoms.SpeciesAndCoordinates").group(1))
    fin = parse_atoms(find_block(texts["fin"], "Atoms.SpeciesAndCoordinates").group(1))
    if len(ini) != len(fin) or [a["sp"] for a in ini] != [a["sp"] for a in fin]:
        sys.exit("ERROR: the two end points differ in number or order of atoms")

    cell = np.array([[float(x) for x in l.split()[:3]]
                     for l in find_block(texts["template"], "Atoms.UnitVectors").group(1).splitlines()
                     if len(l.split()) >= 3])
    cell_unit = (keyword(texts["template"], "Atoms.UnitVectors.Unit") or "ang").lower()
    if cell_unit == "ang" and units["ini"] == "au":
        cell = cell / BOHR
    elif cell_unit == "au" and units["ini"] == "ang":
        cell = cell * BOHR
    inv = np.linalg.inv(cell)

    print("# atom  species  |displacement| between end points (%s)" % units["ini"])
    for a, b in zip(ini, fin):
        frac = (b["xyz"] - a["xyz"]) @ inv
        shift = np.round(frac)
        if np.any(shift != 0):
            b["xyz"] = b["xyz"] - shift @ cell
            print("WARNING: atom %s was wrapped through the cell boundary; shifted by %s lattice vectors"
                  % (b["idx"], (-shift).astype(int)), file=sys.stderr)
        print("  %-5s %-7s %.4f" % (a["idx"], a["sp"], np.linalg.norm(b["xyz"] - a["xyz"])))

    out = texts["template"]
    m = find_block(out, "NEB.Atoms.SpeciesAndCoordinates")
    out = out[:m.start(1)] + format_atoms(fin) + out[m.end(1):]
    m = find_block(out, "Atoms.SpeciesAndCoordinates")
    out = out[:m.start(1)] + format_atoms(ini) + out[m.end(1):]
    out = out.replace("The coordinates below are the UNRELAXED ideal sites.",
                      "End-point coordinates were taken from %s and %s by build_neb.py." % (args.ini, args.fin))
    with open(args.out, "w", newline="\n") as fh:
        fh.write(out)
    print("written: %s" % args.out)


if __name__ == "__main__":
    main()
