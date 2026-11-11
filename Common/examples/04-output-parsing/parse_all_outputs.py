#!/usr/bin/env python3
"""One entry point to read energy / forces / final structure from VASP,
Gaussian, CP2K and OpenMX outputs, returning ASE Atoms objects.

  VASP     OUTCAR or vasprun.xml        -> ase.io.read(format="vasp-out"/"vasp-xml")
  Gaussian *.log / *.out                -> ase.io.read(format="gaussian-out")
                                           + cclib (frequencies, HOMO/LUMO) if installed
  CP2K     PROJECT-pos-1.xyz + main out -> ase.io.read(format="xyz") + energy regex
  OpenMX   System.Name.out              -> ase.calculators.openmx.reader.read_openmx

Usage:
  python parse_all_outputs.py FILE [FILE ...] [--extxyz all.xyz]
Example:
  python parse_all_outputs.py relax/OUTCAR g16/h2o_opt.log cp2k/h2o.out omx/si.out
"""
import argparse
import os
import re

from ase.io import read, write

HA2EV = 27.211386245988


def detect(path):
    name = os.path.basename(path)
    if name.startswith("OUTCAR"):
        return "vasp-out"
    if name.startswith("vasprun") and ".xml" in name:
        return "vasp-xml"
    with open(path, errors="ignore") as f:
        head = f.read(20000)
    if "Entering Gaussian System" in head or "Gaussian, Inc." in head:
        return "gaussian"
    if "CP2K|" in head or "PROGRAM STARTED AT" in head:
        return "cp2k"
    if "OpenMX" in head:
        return "openmx"
    raise ValueError(f"cannot detect code of {path}")


def read_vasp(path, fmt):
    atoms = read(path, format=fmt, index=-1)          # last ionic step
    return atoms, {"energy_eV": atoms.get_potential_energy(),
                   "fmax_eV/A": float(abs(atoms.get_forces()).max())}


def read_gaussian(path):
    atoms = read(path, format="gaussian-out", index=-1)
    info = {"energy_eV": atoms.get_potential_energy()}
    try:
        import cclib  # stable 1.8.1 recommended; 2.0a1 has Gaussian parser bugs
        data = cclib.io.ccread(path)
        if hasattr(data, "vibfreqs"):
            info["nimag"] = int((data.vibfreqs < 0).sum())
            info["lowest_freq_cm-1"] = float(data.vibfreqs.min())
        if hasattr(data, "homos") and hasattr(data, "moenergies"):
            h = data.homos[0]
            info["HOMO_eV"] = float(data.moenergies[0][h])
            info["LUMO_eV"] = float(data.moenergies[0][h + 1])
    except ImportError:
        pass
    return atoms, info


def read_cp2k(path):
    with open(path, errors="ignore") as f:
        txt = f.read()
    # 2025.1+: "[hartree]"; 2024.x and older: "[a.u.]:"
    e = re.findall(r"ENERGY\| Total FORCE_EVAL \( \w+ \) energy \[(?:hartree|a\.u\.)\]:?\s+(-?\d+\.\d+)", txt)
    proj = re.search(r"GLOBAL\| Project name\s+(\S+)", txt)
    info = {"energy_eV": float(e[-1]) * HA2EV if e else None,
            "finished": "PROGRAM ENDED AT" in txt,
            "scf_failed": "SCF run NOT converged" in txt}
    atoms = None
    if proj:
        pos = os.path.join(os.path.dirname(path), f"{proj.group(1)}-pos-1.xyz")
        if os.path.exists(pos):
            atoms = read(pos, format="xyz", index=-1)   # geometry of last step
    return atoms, info


def read_omx(path):
    from ase.calculators.openmx.reader import read_openmx
    label = os.path.splitext(path)[0]          # reads label.out (+ .dat/.log if present)
    atoms = read_openmx(label)
    return atoms, {"energy_eV": atoms.get_potential_energy()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--extxyz", help="write all final structures (+energy) to an extxyz file")
    a = ap.parse_args()
    collected = []
    for p in a.files:
        code = detect(p)
        if code.startswith("vasp"):
            atoms, info = read_vasp(p, code)
        elif code == "gaussian":
            atoms, info = read_gaussian(p)
        elif code == "cp2k":
            atoms, info = read_cp2k(p)
        else:
            atoms, info = read_omx(p)
        formula = atoms.get_chemical_formula() if atoms is not None else "?"
        print(f"{p}: [{code}] {formula}  " + "  ".join(f"{k}={v}" for k, v in info.items()))
        if atoms is not None:
            atoms.info.update({"source": p, "code": code})
            collected.append(atoms)
    if a.extxyz and collected:
        write(a.extxyz, collected, format="extxyz")
        print(f"wrote {len(collected)} structures to {a.extxyz}")


if __name__ == "__main__":
    main()
