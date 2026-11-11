#!/usr/bin/env python3
"""Generate and edit CP2K inputs with pymatgen.io.cp2k; parse a CP2K output.

pymatgen's CP2K sets (StaticSet, RelaxSet, CellOptSet, Hybrid*) produce a
complete input from a Structure, but their defaults are "high-throughput"
choices (UKS on, &DFT/&PRINT/&PDOS block, OT+IRAC ...).  Always inspect and
edit the result: here we switch off UKS for closed-shell water, drop the
print block (CP2K 2026.2 restructured &DFT/&PRINT/&PDOS) and make the
molecule non-periodic.

Usage:
  python pmg_cp2k_input.py                       # writes h2o_static.inp
  python pmg_cp2k_input.py --parse h2o.out       # summary of a finished run
"""
import argparse
import warnings

from pymatgen.core import Molecule
from pymatgen.io.cp2k.inputs import Cp2kInput
from pymatgen.io.cp2k.sets import StaticSet


def build_input():
    mol = Molecule(["O", "H", "H"], [[0.0, 0.0, 0.1193],
                                     [0.0, 0.7632, -0.4770],
                                     [0.0, -0.7632, -0.4770]])
    box = mol.get_boxed_structure(10.0, 10.0, 10.0)
    bp = {  # explicit names -> no need for a pymatgen-indexed CP2K data dir
        "O": {"basis": "DZVP-MOLOPT-SR-GTH", "potential": "GTH-PBE-q6"},
        "H": {"basis": "DZVP-MOLOPT-SR-GTH", "potential": "GTH-PBE-q1"},
        "basis_filenames": ["BASIS_MOLOPT"],
        "potential_filename": "GTH_POTENTIALS",
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")      # "Unable to validate basis ..." is expected
        cs = StaticSet(structure=box, project_name="h2o", xc_functionals=["PBE"],
                       basis_and_potential=bp, cutoff=400, rel_cutoff=50)
    # Round-trip through text so we work on a plain Cp2kInput object
    ci = Cp2kInput.from_str(cs.get_str())

    ci.update({"FORCE_EVAL": {
        "STRESS_TENSOR": "NONE",
        "DFT": {"UKS": False, "MULTIPLICITY": 1,
                "POISSON": {"PERIODIC": "NONE", "PSOLVER": "MT"}},
        "SUBSYS": {"CELL": {"PERIODIC": "NONE"}},
    }})
    dft = ci["FORCE_EVAL"]["DFT"]
    if "PRINT" in dft.subsections:            # PDOS/cube printing: not wanted here
        dft.subsections.pop("PRINT")
    fe_print = ci["FORCE_EVAL"].subsections.get("PRINT")
    if fe_print is not None and "STRESS_TENSOR" in fe_print.subsections:
        fe_print.subsections.pop("STRESS_TENSOR")
    ci.write_file("h2o_static.inp")
    print(ci.get_str())
    print("written h2o_static.inp ; check it with: cp2k.psmp -c h2o_static.inp")


def parse(outfile):
    from pymatgen.io.cp2k.outputs import Cp2kOutput
    out = Cp2kOutput(outfile)
    print("CP2K version      :", out.cp2k_version)
    print("run type          :", out.run_type)
    print("completed         :", out.completed)
    print("final energy (eV) :", out.final_energy)   # pymatgen converts Ha -> eV
    print("SCF converged     :", "SCF run NOT converged" not in open(outfile).read())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parse", metavar="OUT")
    a = ap.parse_args()
    if a.parse:
        parse(a.parse)
    else:
        build_input()


if __name__ == "__main__":
    main()
