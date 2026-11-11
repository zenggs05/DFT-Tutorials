#!/usr/bin/env python3
"""Write VASP input sets for Si with pymatgen (relax -> static -> band/DOS -> HSE).

POTCARs are licensed: with --spec (default) only POTCAR.spec (symbol list)
is written, so the script runs on any laptop.  Use --potcar on a machine
where PMG_VASP_PSP_DIR is configured (pmg config ...) to write real POTCARs.

Usage:
  python pmg_write_inputset.py                 # Si from the built-in structure
  python pmg_write_inputset.py --poscar my.vasp --potcar
"""
import argparse

from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp.sets import (MP24RelaxSet, MPHSEBSSet, MPNonSCFSet,
                                   MPRelaxSet, MPStaticSet, MatPESStaticSet)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def silicon():
    conv = Structure.from_spacegroup("Fd-3m", Lattice.cubic(5.43), ["Si"], [[0, 0, 0]])
    return SpacegroupAnalyzer(conv).get_primitive_standard_structure()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--poscar", help="structure file (POSCAR/CIF/...); default: Si")
    p.add_argument("--potcar", action="store_true", help="write real POTCARs (needs PMG_VASP_PSP_DIR)")
    args = p.parse_args()
    s = Structure.from_file(args.poscar) if args.poscar else silicon()
    spec = not args.potcar

    # 1) Relaxation: MP PBE settings + our overrides (user_incar_settings wins)
    relax = MPRelaxSet(s, user_incar_settings={"ENCUT": 520, "EDIFFG": -0.01,
                                                "ISMEAR": 0, "SIGMA": 0.05,
                                                "NCORE": 4, "LWAVE": False},
                       user_kpoints_settings={"reciprocal_density": 200})
    relax.write_input("01_relax", potcar_spec=spec)

    # 2) Static SCF on the same structure (in practice: on 01_relax/CONTCAR)
    MPStaticSet(s, user_incar_settings={"ENCUT": 520, "NCORE": 4}).write_input(
        "02_static", potcar_spec=spec)

    # 3) Non-SCF band structure along the high-symmetry line and a dense-mesh DOS.
    #    In a real workflow use MPNonSCFSet.from_prev_calc("02_static", mode="line")
    #    so that ISPIN/NBANDS/MAGMOM are inherited and CHGCAR is copied.
    MPNonSCFSet(s, mode="line", kpoints_line_density=40,
                user_incar_settings={"ENCUT": 520}).write_input("03_band", potcar_spec=spec)
    MPNonSCFSet(s, mode="uniform", reciprocal_density=400, nedos=3001,
                user_incar_settings={"ENCUT": 520}).write_input("04_dos", potcar_spec=spec)

    # 4) HSE06 band structure with zero-weight k-points along the path
    MPHSEBSSet(s, mode="line", user_incar_settings={"ENCUT": 520, "NCORE": 4, "LASPH": True}).write_input(
        "05_hse_band", potcar_spec=spec)

    # 5) Newer standards: MP24 (r2SCAN-based) and MatPES (PBE, ML training data)
    MP24RelaxSet(s).write_input("06_mp24_relax", potcar_spec=spec)
    MatPESStaticSet(s).write_input("07_matpes_static", potcar_spec=spec)

    for d in ["01_relax", "03_band", "05_hse_band", "06_mp24_relax"]:
        with open(f"{d}/INCAR") as f:
            tags = dict(l.split(" = ", 1) for l in f.read().splitlines() if " = " in l)
        keys = ["ENCUT", "ISMEAR", "ICHARG", "LHFCALC", "METAGGA", "KSPACING", "NSW"]
        print(f"{d:16s}", {k: tags[k] for k in keys if k in tags})


if __name__ == "__main__":
    main()
