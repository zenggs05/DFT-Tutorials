#!/usr/bin/env python3
"""Write Gaussian 16 inputs with pymatgen.io.gaussian and read a log file back.

Usage:
  python pmg_gaussian_input.py                      # writes benzene_opt_freq.gjf
  python pmg_gaussian_input.py --xyz mol.xyz --charge 0 --mult 1
  python pmg_gaussian_input.py --parse job.log      # summary of a finished job
"""
import argparse

from pymatgen.core import Molecule
from pymatgen.io.gaussian import GaussianInput, GaussianOutput


def benzene():
    """Planar benzene, C-C 1.39 A, C-H 1.09 A (D6h), built from scratch."""
    import math
    species, coords = [], []
    for i in range(6):
        t = math.radians(60 * i)
        species.append("C")
        coords.append([1.39 * math.cos(t), 1.39 * math.sin(t), 0.0])
    for i in range(6):
        t = math.radians(60 * i)
        species.append("H")
        coords.append([2.48 * math.cos(t), 2.48 * math.sin(t), 0.0])
    return Molecule(species, coords)


def write_input(mol, charge, mult, fname):
    gi = GaussianInput(
        mol, charge=charge, spin_multiplicity=mult,
        title="benzene opt+freq B3LYP-D3(BJ)/def2-SVP",
        functional="B3LYP", basis_set="def2SVP",
        # route keywords: "" for a bare keyword, a string for keyword=value
        # (G16 already uses the UltraFine integration grid by default)
        route_parameters={"opt": "", "freq": "", "EmpiricalDispersion": "GD3BJ",
                          "scf": "tight"},
        link0_parameters={"%mem": "8GB", "%nprocshared": "16", "%chk": "benzene.chk"},
        dieze_tag="#P",
    )
    text = gi.to_str(cart_coords=True)
    if not text.endswith("\n\n"):          # Gaussian needs a terminating blank line
        text = text.rstrip("\n") + "\n\n"
    with open(fname, "w", newline="\n") as f:
        f.write(text)
    print(text)


def parse(logfile):
    go = GaussianOutput(logfile)
    print("normal termination :", go.properly_terminated)
    print("SCF energy (Ha)    :", go.final_energy)
    print("stationary point   :", getattr(go, "stationary_type", "n/a"))
    if go.frequencies:
        freqs = [m["frequency"] for m in go.frequencies[-1]]
        nimag = sum(1 for f in freqs if f < 0)
        print(f"frequencies (cm-1) : {freqs[:6]} ... ({len(freqs)} modes, {nimag} imaginary)")
    print("final geometry:")
    print(go.final_structure.to(fmt="xyz"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", help="input molecule (xyz); default: built-in benzene")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--out", default="benzene_opt_freq.gjf")
    ap.add_argument("--parse", metavar="LOG", help="parse a Gaussian .log/.out instead")
    a = ap.parse_args()
    if a.parse:
        parse(a.parse)
        return
    mol = Molecule.from_file(a.xyz) if a.xyz else benzene()
    write_input(mol, a.charge, a.mult, a.out)


if __name__ == "__main__":
    main()
