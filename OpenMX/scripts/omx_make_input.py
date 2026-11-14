#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_make_input.py -- generate an OpenMX 4.0 input file (.dat) from POSCAR / CIF / xyz.
                     从 POSCAR/CIF/xyz 等结构文件生成 OpenMX 输入：内置 DFT_DATA19 推荐基组表
                     (Quick/Standard/Precise) 与价电子数，自动写 Definition.of.Atomic.Species、
                     坐标行的上/下自旋初始占据（FM/AFM/MAGMOM）、按 k 间距设置 scf.Kgrid，
                     并可用 seekpath 生成 Band.kpath。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

What it writes
--------------
* System.CurrrentDirectory / System.Name / DATA.PATH / level.of.*
  (the keyword System.CurrrentDirectory really has three r's)
* Species.Number + <Definition.of.Atomic.Species: one species per element,
  PAO and VPS taken from Tables 1-2 of manual Sec. 12.3 (Database Ver. 2019)
  for the chosen accuracy level (quick | standard | precise).
* Atoms.SpeciesAndCoordinates (FRAC for periodic cells, Ang for molecules)
  with up/down initial electrons whose sum equals the VPS valence charge.
* Atoms.UnitVectors in Ang (molecules: box = extent + 2 x --vacuum).
* scf.* keywords: GGA-PBE (or LDA), spin, cutoff, solver (Band/Cluster),
  scf.Kgrid from a k-point spacing, scf.Generation.Kpoint GCenter (4.0 only),
  RMM-DIISK mixing, criterion.
* optional MD.Type (geometry / cell optimization), Band.kpath (seekpath),
  Dos.fileout.

Valence electrons
-----------------
The built-in table gives the valence charge of every *_PBE19 VPS listed in
manual Sec. 12.3.  With --data-path the value is read from the actual
VPS file (line "valence.electron") whenever that file exists, which is the
recommended way and is required for LDA (*_CA19) potentials.
Elements with two potentials (Fe, Co, Ni, Cu, Zn: *_PBE19H with semicore
3s states, *_PBE19S without) default to H; use --soft for S.
Elements missing from the table (e.g. Eu, Gd, Tb, Yb) need --basis.

Usage examples
--------------
  python omx_make_input.py POSCAR -o si.dat --data-path /opt/openmx4.0/DFT_DATA19
  python omx_make_input.py Fe.cif -o fe.dat --spin on --moment Fe=3.0 --etemp 1000
  python omx_make_input.py NiO_afm.vasp --spin on --moment Ni=1.8 --afm Ni
  python omx_make_input.py POSCAR --magmom "2*3.0 2*-3.0 4*0" --spin on
  python omx_make_input.py POSCAR --kpath --kpath-density 30 -o band.dat
  python omx_make_input.py h2o.xyz --vacuum 6 -o h2o.dat       # molecule -> Cluster
  python omx_make_input.py POSCAR --level precise --basis Si=Si7.0-s3p3d2
  python omx_make_input.py POSCAR --basis Fe=Fe6.0S-s3p2d2f1:Fe_PBE19S:14  # PAO:VPS:valence

Dependencies
------------
python >= 3.9, numpy, and ASE (preferred) or pymatgen for reading
structures; optional spglib (--primitive) and seekpath (--kpath).

Always review the generated file: cutoff, k-grid, basis level, spin and
mixing parameters must be checked by convergence tests (chapters 4-5).
"""
from __future__ import annotations

import argparse
import datetime
import math
import os
import re
import sys

# ---------------------------------------------------------------------------
# DFT_DATA19 recommended basis sets (OpenMX 4.0 manual, Sec. 12.3, Tables 1-2)
#   VPS name: (valence electrons, Quick, Standard, Precise)
# ---------------------------------------------------------------------------
DFT_DATA19 = {
    "H_PBE19": (1.0, "H5.0-s2", "H6.0-s2p1", "H7.0-s2p2d1"),
    "He_PBE19": (2.0, "He8.0-s1p1", "He8.0-s2p1", "He10.0-s2p2d1"),
    "Li_PBE19": (3.0, "Li8.0-s3p1", "Li8.0-s3p2", "Li8.0-s3p2d1"),
    "Be_PBE19": (2.0, "Be7.0-s2p1", "Be7.0-s2p2", "Be7.0-s3p2d1"),
    "B_PBE19": (3.0, "B7.0-s2p2", "B7.0-s2p2d1", "B7.0-s3p2d2"),
    "C_PBE19": (4.0, "C6.0-s2p2", "C6.0-s2p2d1", "C6.0-s3p2d2"),
    "N_PBE19": (5.0, "N6.0-s2p2", "N6.0-s2p2d1", "N6.0-s3p2d2"),
    "O_PBE19": (6.0, "O6.0-s2p2", "O6.0-s2p2d1", "O6.0-s3p2d2"),
    "F_PBE19": (7.0, "F6.0-s2p2", "F6.0-s2p2d1", "F6.0-s3p3d2f1"),
    "Ne_PBE19": (8.0, "Ne9.0-s2p2", "Ne9.0-s2p2d1", "Ne9.0-s3p2d2"),
    "Na_PBE19": (9.0, "Na9.0-s3p2", "Na9.0-s3p2d1", "Na9.0-s3p2d2"),
    "Mg_PBE19": (8.0, "Mg9.0-s2p2", "Mg9.0-s3p2d1", "Mg9.0-s3p2d2"),
    "Al_PBE19": (3.0, "Al7.0-s2p1d1", "Al7.0-s2p2d1", "Al7.0-s3p2d2"),
    "Si_PBE19": (4.0, "Si7.0-s2p1d1", "Si7.0-s2p2d1", "Si7.0-s3p3d2"),
    "P_PBE19": (5.0, "P7.0-s2p2d1", "P7.0-s2p2d1f1", "P7.0-s3p2d2f1"),
    "S_PBE19": (6.0, "S7.0-s2p2d1", "S7.0-s2p2d1f1", "S7.0-s3p2d2f1"),
    "Cl_PBE19": (7.0, "Cl7.0-s2p2d1", "Cl7.0-s2p2d1f1", "Cl7.0-s3p2d2f1"),
    "Ar_PBE19": (8.0, "Ar9.0-s2p2d1", "Ar9.0-s2p2d1f1", "Ar9.0-s3p2d2f1"),
    "K_PBE19": (9.0, "K10.0-s3p2", "K10.0-s3p2d1", "K10.0-s3p2d2"),
    "Ca_PBE19": (10.0, "Ca9.0-s3p2", "Ca9.0-s3p2d1", "Ca9.0-s3p2d2"),
    "Sc_PBE19": (11.0, "Sc9.0-s2p2d1", "Sc9.0-s3p2d1", "Sc9.0-s3p2d2"),
    "Ti_PBE19": (12.0, "Ti7.0-s2p2d1", "Ti7.0-s3p2d1", "Ti7.0-s3p2d2f1"),
    "V_PBE19": (13.0, "V6.0-s2p2d1", "V6.0-s3p2d1", "V6.0-s3p2d2f1"),
    "Cr_PBE19": (14.0, "Cr6.0-s2p2d1", "Cr6.0-s3p2d1", "Cr6.0-s3p2d2f1"),
    "Mn_PBE19": (15.0, "Mn6.0-s2p2d1", "Mn6.0-s3p2d1", "Mn6.0-s3p2d2f1"),
    "Fe_PBE19H": (16.0, "Fe5.5H-s2p2d1", "Fe5.5H-s3p2d1", "Fe5.5H-s3p2d2f1"),
    "Fe_PBE19S": (14.0, "Fe6.0S-s2p2d1", "Fe6.0S-s3p2d1", "Fe6.0S-s3p2d2f1"),
    "Co_PBE19H": (17.0, "Co6.0H-s2p2d1", "Co6.0H-s3p2d1", "Co6.0H-s3p2d2f1"),
    "Co_PBE19S": (15.0, "Co6.0S-s2p2d1", "Co6.0S-s3p2d1", "Co6.0S-s3p2d2f1"),
    "Ni_PBE19H": (18.0, "Ni6.0H-s2p2d1", "Ni6.0H-s3p2d1", "Ni6.0H-s3p2d2f1"),
    "Ni_PBE19S": (16.0, "Ni6.0S-s2p2d1", "Ni6.0S-s3p2d1", "Ni6.0S-s3p2d2f1"),
    "Cu_PBE19H": (19.0, "Cu6.0H-s2p2d1", "Cu6.0H-s3p2d1", "Cu6.0H-s3p2d2f1"),
    "Cu_PBE19S": (11.0, "Cu6.0S-s2p1d1", "Cu6.0S-s3p2d1", "Cu6.0S-s3p2d2f1"),
    "Zn_PBE19H": (20.0, "Zn6.0H-s2p2d1", "Zn6.0H-s3p2d1", "Zn6.0H-s3p2d2f1"),
    "Zn_PBE19S": (12.0, "Zn6.0S-s2p1d1", "Zn6.0S-s3p2d1", "Zn6.0S-s3p2d2f1"),
    "Ga_PBE19": (13.0, "Ga7.0-s2p2d1", "Ga7.0-s3p2d2", "Ga7.0-s3p2d2f1"),
    "Ge_PBE19": (4.0, "Ge7.0-s2p1d1", "Ge7.0-s3p2d2", "Ge7.0-s3p2d2f1"),
    "As_PBE19": (15.0, "As7.0-s3p2d1", "As7.0-s3p2d2", "As7.0-s3p2d2f1"),
    "Se_PBE19": (6.0, "Se7.0-s3p2d1", "Se7.0-s3p2d2", "Se7.0-s3p2d2f1"),
    "Br_PBE19": (7.0, "Br7.0-s3p2d1", "Br7.0-s3p2d2", "Br7.0-s3p2d2f1"),
    "Kr_PBE19": (8.0, "Kr10.0-s2p2d1", "Kr10.0-s3p2d2", "Kr10.0-s3p2d2f1"),
    "Rb_PBE19": (9.0, "Rb11.0-s2p2d1", "Rb11.0-s3p2d2", "Rb11.0-s3p2d2f1"),
    "Sr_PBE19": (10.0, "Sr10.0-s2p2d1", "Sr10.0-s3p2d2", "Sr10.0-s3p3d2f1"),
    "Y_PBE19": (11.0, "Y10.0-s3p2d1", "Y10.0-s3p2d2", "Y10.0-s3p3d2f1"),
    "Zr_PBE19": (12.0, "Zr7.0-s3p2d1", "Zr7.0-s3p2d2", "Zr7.0-s3p2d2f1"),
    "Nb_PBE19": (13.0, "Nb7.0-s3p2d1", "Nb7.0-s3p2d2", "Nb7.0-s3p2d2f1"),
    "Mo_PBE19": (14.0, "Mo7.0-s3p2d1", "Mo7.0-s3p2d2", "Mo7.0-s3p2d2f1"),
    "Tc_PBE19": (15.0, "Tc7.0-s3p2d1", "Tc7.0-s3p2d2", "Tc7.0-s3p2d2f1"),
    "Ru_PBE19": (14.0, "Ru7.0-s3p2d1", "Ru7.0-s3p2d2", "Ru7.0-s3p2d2f1"),
    "Rh_PBE19": (15.0, "Rh7.0-s3p2d1", "Rh7.0-s3p2d2", "Rh7.0-s3p2d2f1"),
    "Pd_PBE19": (16.0, "Pd7.0-s3p2d1", "Pd7.0-s3p2d2", "Pd7.0-s3p2d2f1"),
    "Ag_PBE19": (17.0, "Ag7.0-s3p2d1", "Ag7.0-s3p2d2", "Ag7.0-s3p2d2f1"),
    "Cd_PBE19": (12.0, "Cd7.0-s3p2d1", "Cd7.0-s3p2d2", "Cd7.0-s3p2d2f1"),
    "In_PBE19": (13.0, "In7.0-s3p2d1", "In7.0-s3p2d2", "In7.0-s3p2d2f1"),
    "Sn_PBE19": (14.0, "Sn7.0-s2p2d2", "Sn7.0-s3p2d3", "Sn7.0-s3p2d3f1"),
    "Sb_PBE19": (15.0, "Sb7.0-s3p2d1", "Sb7.0-s3p2d2", "Sb7.0-s3p2d2f1"),
    "Te_PBE19": (16.0, "Te7.0-s3p2d2", "Te7.0-s3p2d2f1", "Te7.0-s3p3d2f1"),
    "I_PBE19": (7.0, "I7.0-s3p2d2", "I7.0-s3p2d2f1", "I7.0-s3p3d2f1"),
    "Xe_PBE19": (8.0, "Xe11.0-s3p2d1", "Xe11.0-s3p2d2", "Xe11.0-s3p2d2f1"),
    "Cs_PBE19": (9.0, "Cs12.0-s3p2d1", "Cs12.0-s3p2d2", "Cs12.0-s3p2d2f1"),
    "Ba_PBE19": (10.0, "Ba10.0-s3p2d1", "Ba10.0-s3p2d2", "Ba10.0-s3p2d2f1"),
    "La_PBE19": (11.0, "La8.0-s3p2d1f1", "La8.0-s3p2d2f1", "La8.0-s3p3d2f1"),
    "Ce_PBE19": (12.0, "Ce8.0-s3p2d1f1", "Ce8.0-s3p2d2f1", "Ce8.0-s3p3d2f1"),
    "Pr_PBE19": (13.0, "Pr8.0-s3p2d1f1", "Pr8.0-s3p2d2f1", "Pr8.0-s3p3d2f1"),
    "Nd_PBE19": (14.0, "Nd8.0-s3p2d1f1", "Nd8.0-s3p2d2f1", "Nd8.0-s3p3d2f1"),
    "Pm_PBE19": (15.0, "Pm8.0-s3p2d1f1", "Pm8.0-s3p2d2f1", "Pm8.0-s3p3d2f1"),
    "Sm_PBE19": (16.0, "Sm8.0-s3p2d1f1", "Sm8.0-s3p2d2f1", "Sm8.0-s3p3d2f1"),
    "Dy_PBE19": (20.0, "Dy8.0-s3p2d1f1", "Dy8.0-s3p2d2f1", "Dy8.0-s3p3d2f1"),
    "Ho_PBE19": (21.0, "Ho8.0-s3p2d1f1", "Ho8.0-s3p2d2f1", "Ho8.0-s3p3d2f1"),
    "Lu_PBE19": (11.0, "Lu8.0-s3p2d2", "Lu8.0-s3p2d2f1", "Lu8.0-s3p3d2f1"),
    "Hf_PBE19": (12.0, "Hf9.0-s3p2d2", "Hf9.0-s3p2d2f1", "Hf9.0-s3p3d2f1"),
    "Ta_PBE19": (13.0, "Ta7.0-s3p2d2", "Ta7.0-s3p2d2f1", "Ta7.0-s3p3d2f1"),
    "W_PBE19": (12.0, "W7.0-s3p2d2", "W7.0-s3p2d2f1", "W7.0-s3p3d2f1"),
    "Re_PBE19": (15.0, "Re7.0-s3p2d2", "Re7.0-s3p2d2f1", "Re7.0-s3p3d2f1"),
    "Os_PBE19": (14.0, "Os7.0-s3p2d2", "Os7.0-s3p2d2f1", "Os7.0-s3p3d2f1"),
    "Ir_PBE19": (15.0, "Ir7.0-s3p2d2", "Ir7.0-s3p2d2f1", "Ir7.0-s3p3d2f1"),
    "Pt_PBE19": (16.0, "Pt7.0-s3p2d2", "Pt7.0-s3p2d2f1", "Pt7.0-s3p3d2f1"),
    "Au_PBE19": (17.0, "Au7.0-s3p2d2", "Au7.0-s3p2d2f1", "Au7.0-s3p3d2f1"),
    "Hg_PBE19": (18.0, "Hg8.0-s3p2d2", "Hg8.0-s3p2d2f1", "Hg8.0-s3p3d2f1"),
    "Tl_PBE19": (19.0, "Tl8.0-s3p2d2", "Tl8.0-s3p2d2f1", "Tl8.0-s3p3d2f1"),
    "Pb_PBE19": (14.0, "Pb8.0-s3p2d2", "Pb8.0-s3p2d2f1", "Pb8.0-s3p3d2f1"),
    "Bi_PBE19": (15.0, "Bi8.0-s3p2d2", "Bi8.0-s3p2d2f1", "Bi8.0-s3p3d2f1"),
}

LEVELS = {"quick": 1, "standard": 2, "precise": 3}

# elements with an H (semicore) and an S (soft) potential in the table
_TWO_VPS = sorted({k.split("_")[0] for k in DFT_DATA19 if k.endswith(("PBE19H", "PBE19S"))})


# ---------------------------------------------------------------------------
# helpers that other scripts of this folder import (omx_magconf.py, omx_parse_out.py,
# omx_conv_test.py); keep their signatures stable
# ---------------------------------------------------------------------------
def table_key(vps: str) -> str:
    """Map a VPS name to the key of DFT_DATA19 (CA19 -> PBE19, same valence assumed)."""
    return vps.replace("_CA19", "_PBE19")


def read_vps_valence(vps: str, data_path: str | None) -> float | None:
    """Return 'valence.electron' from DATA.PATH/VPS/<vps>.vps, or None if unavailable."""
    if not data_path:
        return None
    path = os.path.join(data_path, "VPS", vps + ".vps")
    if not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="ascii", errors="replace") as fh:
            for line in fh:
                tok = line.split()
                if len(tok) >= 2 and tok[0].lower() == "valence.electron":
                    return float(tok[1])
    except (OSError, ValueError):
        return None
    return None


def valence_of(vps: str, data_path: str | None = None) -> float | None:
    """Valence charge of a VPS: the VPS file if found, else the built-in table."""
    val = read_vps_valence(vps, data_path)
    if val is not None:
        return val
    entry = DFT_DATA19.get(table_key(vps))
    return float(entry[0]) if entry else None


def recommended_basis(vps: str, level: str = "standard") -> str | None:
    """Recommended PAO string (e.g. 'Si7.0-s2p2d1') for a VPS and accuracy level."""
    entry = DFT_DATA19.get(table_key(vps))
    if entry is None:
        return None
    return entry[LEVELS[level.lower()]]


def default_vps(elem: str, xc: str = "GGA-PBE", soft: bool = False) -> str | None:
    """Default VPS name of an element in DFT_DATA19 (None if not tabulated)."""
    if elem + "_PBE19" in DFT_DATA19:
        name = elem + "_PBE19"
    elif elem in _TWO_VPS:
        name = elem + ("_PBE19S" if soft else "_PBE19H")
    else:
        return None
    if xc.upper().startswith("LDA") or xc.upper().startswith("LSDA"):
        name = name.replace("_PBE19", "_CA19")
    return name


# ---------------------------------------------------------------------------
# small linear algebra without forcing numpy at import time
# ---------------------------------------------------------------------------
def _np():
    try:
        import numpy as np  # noqa: F401
    except ImportError:
        sys.exit("ERROR: numpy is required (pip install numpy)")
    return np


def reciprocal(cell):
    """Reciprocal vectors (rows) including the 2*pi factor, 1/Angstrom."""
    np = _np()
    a = np.asarray(cell, dtype=float)
    return 2.0 * math.pi * np.linalg.inv(a).T


def kgrid_from_spacing(cell, spacing: float, pbc=(True, True, True)):
    """n_i = max(1, ceil(|b_i| / spacing)), |b_i| including 2*pi (VASP KSPACING convention)."""
    np = _np()
    b = reciprocal(cell)
    grid = []
    for i in range(3):
        if not pbc[i]:
            grid.append(1)
            continue
        grid.append(max(1, int(math.ceil(np.linalg.norm(b[i]) / spacing - 1e-8))))
    return grid


# ---------------------------------------------------------------------------
# structure input
# ---------------------------------------------------------------------------
def read_structure(path: str, fmt: str | None = None) -> dict:
    """Return dict(symbols, positions[Ang, Cartesian], cell[Ang] or None, pbc)."""
    if not os.path.isfile(path):
        sys.exit("ERROR: structure file not found: %s" % path)
    np = _np()
    try:
        import ase.io  # type: ignore
    except ImportError:
        ase_io = None
    else:
        ase_io = ase.io
    if ase_io is not None:
        try:
            atoms = ase_io.read(path, format=fmt) if fmt else ase_io.read(path)
        except Exception as exc:  # ASE raises many different exception types
            sys.exit("ERROR: ASE could not read %s: %s" % (path, exc))
        pbc = tuple(bool(x) for x in atoms.get_pbc())
        cell = np.array(atoms.get_cell()) if any(pbc) else None
        if cell is not None and abs(np.linalg.det(cell)) < 1e-8:
            cell = None
            pbc = (False, False, False)
        return {"symbols": atoms.get_chemical_symbols(),
                "positions": np.array(atoms.get_positions()),
                "cell": cell, "pbc": pbc}
    try:
        from pymatgen.core import Molecule, Structure  # type: ignore
    except ImportError:
        sys.exit("ERROR: reading structures needs ASE (pip install ase) or pymatgen")
    try:
        st = Structure.from_file(path)
        return {"symbols": [s.specie.symbol for s in st],
                "positions": np.array(st.cart_coords),
                "cell": np.array(st.lattice.matrix), "pbc": (True, True, True)}
    except Exception:
        try:
            mol = Molecule.from_file(path)
        except Exception as exc:
            sys.exit("ERROR: pymatgen could not read %s: %s" % (path, exc))
        return {"symbols": [s.specie.symbol for s in mol],
                "positions": np.array(mol.cart_coords), "cell": None,
                "pbc": (False, False, False)}


def _atomic_numbers(symbols):
    try:
        from ase.data import atomic_numbers  # type: ignore
        return [atomic_numbers[s] for s in symbols]
    except ImportError:
        uniq = sorted(set(symbols))
        return [uniq.index(s) + 1 for s in symbols]   # only identity matters for spglib


def to_primitive(st: dict, symprec: float) -> dict:
    try:
        import spglib  # type: ignore
    except ImportError:
        sys.exit("ERROR: --primitive needs spglib (pip install spglib)")
    np = _np()
    cell = st["cell"]
    frac = np.linalg.solve(cell.T, st["positions"].T).T
    nums = _atomic_numbers(st["symbols"])
    res = spglib.find_primitive((cell, frac, nums), symprec=symprec)
    if res is None:
        print("WARNING: spglib found no primitive cell; the input cell is kept", file=sys.stderr)
        return st
    lat, pos, num = res
    back = {n: s for n, s in zip(nums, st["symbols"])}
    return {"symbols": [back[int(n)] for n in num], "positions": np.dot(pos, lat),
            "cell": np.array(lat), "pbc": st["pbc"]}


def seekpath_path(st: dict, symprec: float, density: float):
    """Standardized primitive cell + Band.kpath lines from seekpath (HPKOT convention)."""
    try:
        import seekpath  # type: ignore
    except ImportError:
        sys.exit("ERROR: --kpath needs seekpath (pip install seekpath)")
    np = _np()
    cell = st["cell"]
    frac = np.linalg.solve(cell.T, st["positions"].T).T
    nums = _atomic_numbers(st["symbols"])
    res = seekpath.get_path((cell, frac, nums), with_time_reversal=True, symprec=symprec)
    lat = np.array(res["primitive_lattice"])
    pos = np.array(res["primitive_positions"])
    back = {n: s for n, s in zip(nums, st["symbols"])}
    new = {"symbols": [back[int(n)] for n in res["primitive_types"]],
           "positions": np.dot(pos, lat), "cell": lat, "pbc": st["pbc"]}
    b = reciprocal(lat)
    lines = []
    for lab1, lab2 in res["path"]:
        k1 = np.array(res["point_coords"][lab1], dtype=float)
        k2 = np.array(res["point_coords"][lab2], dtype=float)
        length = np.linalg.norm(np.dot(k2 - k1, b))
        npts = max(2, int(round(length * density)))
        lines.append((npts, k1, k2, _label(lab1), _label(lab2)))
    return new, lines


def _label(lab: str) -> str:
    return "g" if lab.upper() in ("GAMMA", "G", "\u0393") else lab.replace(" ", "")


# ---------------------------------------------------------------------------
# magnetic moments
# ---------------------------------------------------------------------------
def expand_magmom(text: str, natoms: int):
    """VASP-like list: '2*3.0 -3 0' -> [3.0, 3.0, -3.0, 0.0]."""
    vals = []
    for tok in text.replace(",", " ").split():
        if "*" in tok:
            n, v = tok.split("*", 1)
            vals.extend([float(v)] * int(n))
        else:
            vals.append(float(tok))
    if len(vals) != natoms:
        sys.exit("ERROR: --magmom has %d values but the structure has %d atoms" % (len(vals), natoms))
    return vals


def parse_pairs(items, what: str):
    out = {}
    for it in items or []:
        if "=" not in it:
            sys.exit("ERROR: %s expects Element=value, got '%s'" % (what, it))
        k, v = it.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def occupation(valence: float, moment: float, label: str):
    """Split the valence charge into (up, down) with up - down = moment."""
    if abs(moment) > valence + 1e-9:
        sys.exit("ERROR: |moment| = %.3f exceeds the valence charge %.3f of %s"
                 % (abs(moment), valence, label))
    return 0.5 * (valence + moment), 0.5 * (valence - moment)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate an OpenMX 4.0 .dat file from a structure file (POSCAR, CIF, xyz, ...). "
                    "Basis sets and valence charges follow DFT_DATA19 (manual Sec. 12.3).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python omx_make_input.py POSCAR -o si.dat --data-path /opt/openmx4.0/DFT_DATA19\n"
               "  python omx_make_input.py Fe.cif --spin on --moment Fe=3.0 --etemp 1000\n"
               "  python omx_make_input.py NiO.vasp --spin on --moment Ni=1.8 --afm Ni\n"
               "  python omx_make_input.py POSCAR --kpath -o band.dat\n"
               "  python omx_make_input.py h2o.xyz --vacuum 6 -o h2o.dat")
    p.add_argument("structure", help="structure file readable by ASE/pymatgen (POSCAR, CIF, xyz, ...)")
    p.add_argument("--format", default=None, help="force the ASE file format (e.g. vasp, cif, xyz)")
    p.add_argument("-o", "--output", default=None, help="output .dat file (default: stdout)")
    p.add_argument("--name", default=None, help="System.Name (default: output file stem or 'omx')")
    p.add_argument("--data-path", default="../DFT_DATA19",
                   help="DATA.PATH: directory containing PAO/ and VPS/ (default ../DFT_DATA19). "
                        "If it exists locally, PAO/VPS files and valence charges are checked")
    p.add_argument("--level", choices=sorted(LEVELS), default="standard",
                   help="basis level from the DFT_DATA19 table (default standard)")
    p.add_argument("--soft", action="store_true",
                   help="use the *_PBE19S potentials for Fe, Co, Ni, Cu, Zn (default *_PBE19H)")
    p.add_argument("--basis", action="append", metavar="EL=PAO[:VPS[:VALENCE]]",
                   help="override the basis of one element, e.g. Si=Si7.0-s3p3d2 or "
                        "Fe=Fe6.0S-s3p2d2f1:Fe_PBE19S:14 (repeatable); elements missing from the "
                        "table need all three fields (take the valence from the VPS file)")
    p.add_argument("--xc", choices=["GGA-PBE", "LDA", "LSDA-CA", "LSDA-PW"], default="GGA-PBE",
                   help="scf.XcType (LDA/LSDA use *_CA19 potentials; default GGA-PBE)")
    p.add_argument("--spin", choices=["off", "on", "nc"], default=None,
                   help="scf.SpinPolarization (default: on if any moment is given, else off)")
    p.add_argument("--moment", action="append", metavar="EL=MU",
                   help="initial moment (muB, = up - down electrons) per element, e.g. Fe=3.0 (repeatable)")
    p.add_argument("--afm", action="append", metavar="EL",
                   help="alternate the sign of the moment of this element in file order (repeatable)")
    p.add_argument("--magmom", default=None,
                   help="per-atom initial moments, VASP style ('2*3.0 2*-3.0 4*0'); overrides --moment")
    p.add_argument("--ecut", type=float, default=220.0, help="scf.energycutoff in Ry (default 220)")
    p.add_argument("--kspacing", type=float, default=0.20,
                   help="k-point spacing in 1/Ang incl. 2*pi, for scf.Kgrid (default 0.20)")
    p.add_argument("--kgrid", type=int, nargs=3, default=None, metavar=("N1", "N2", "N3"),
                   help="explicit scf.Kgrid (overrides --kspacing)")
    p.add_argument("--etemp", type=float, default=300.0,
                   help="scf.ElectronicTemperature in K (default 300; metals often 500-1000)")
    p.add_argument("--criterion", type=float, default=1.0e-8, help="scf.criterion in Hartree (default 1e-8)")
    p.add_argument("--maxiter", type=int, default=200, help="scf.maxIter (default 200)")
    p.add_argument("--charge", type=float, default=0.0,
                   help="scf.system.charge (positive = holes); written only if nonzero")
    p.add_argument("--omx-version", choices=["4.0", "3.9"], default="4.0",
                   help="target version: 4.0 writes scf.Generation.Kpoint GCenter, 3.9 omits it")
    p.add_argument("--md-type", default="nomd",
                   help="MD.Type, e.g. nomd, RF, EF, BFGS, RFC5 (default nomd)")
    p.add_argument("--primitive", action="store_true", help="reduce to the primitive cell (spglib)")
    p.add_argument("--symprec", type=float, default=1e-3, help="symmetry tolerance for spglib/seekpath")
    p.add_argument("--kpath", action="store_true",
                   help="add Band.dispersion + Band.kpath from seekpath; the structure is replaced "
                        "by the seekpath standardized primitive cell so that labels are correct")
    p.add_argument("--kpath-density", type=float, default=25.0,
                   help="k points per 1/Ang along the band path (default 25)")
    p.add_argument("--dos", action="store_true",
                   help="add Dos.fileout on, Dos.Erange -20 20 and Dos.Kgrid = 2 x scf.Kgrid")
    p.add_argument("--vacuum", type=float, default=6.0,
                   help="molecules only: vacuum on each side of the box, Ang (default 6)")
    p.add_argument("--fileout", type=int, choices=[0, 1, 2], default=1,
                   help="level.of.fileout (default 1; 0 = no cube files)")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    np = _np()
    st = read_structure(args.structure, args.format)
    periodic = st["cell"] is not None

    if args.primitive:
        if not periodic:
            sys.exit("ERROR: --primitive needs a periodic structure")
        st = to_primitive(st, args.symprec)
    kpath_lines = None
    if args.kpath:
        if not periodic:
            sys.exit("ERROR: --kpath needs a periodic structure")
        st, kpath_lines = seekpath_path(st, args.symprec, args.kpath_density)
        print("NOTE: --kpath: the structure was converted to the seekpath standardized primitive "
              "cell (%d atoms); compare it with your input." % len(st["symbols"]), file=sys.stderr)

    symbols = list(st["symbols"])
    natoms = len(symbols)
    if natoms == 0:
        sys.exit("ERROR: no atoms found in %s" % args.structure)
    elements = []
    for s in symbols:
        if s not in elements:
            elements.append(s)

    data_local = os.path.isdir(args.data_path)
    overrides = parse_pairs(args.basis, "--basis")
    species = {}          # element -> (PAO, VPS, valence)
    for el in elements:
        if el in overrides:
            parts = overrides[el].split(":")
            pao = parts[0]
            vps = parts[1] if len(parts) > 1 and parts[1] else default_vps(el, args.xc, args.soft)
            if vps is None:
                sys.exit("ERROR: %s is not in the DFT_DATA19 table; give --basis %s=PAO:VPS:VALENCE"
                         % (el, el))
            val = float(parts[2]) if len(parts) > 2 and parts[2] else valence_of(vps, args.data_path)
        else:
            vps = default_vps(el, args.xc, args.soft)
            if vps is None:
                sys.exit("ERROR: element %s is not in the built-in DFT_DATA19 table; use "
                         "--basis %s=PAO:VPS:VALENCE (see manual Sec. 12.3 / the VPS file)" % (el, el))
            pao = recommended_basis(vps, args.level)
            val = valence_of(vps, args.data_path)
        if val is None:
            sys.exit("ERROR: unknown valence charge for %s (%s); add it as --basis %s=%s:%s:VALENCE"
                     % (el, vps, el, pao, vps))
        if vps.endswith("_CA19") and read_vps_valence(vps, args.data_path) is None:
            print("WARNING: valence of %s taken from the PBE19 table; verify with "
                  "'grep valence.electron %s/VPS/%s.vps'" % (vps, args.data_path, vps), file=sys.stderr)
        if data_local:
            pao_file = os.path.join(args.data_path, "PAO", pao.split("-")[0] + ".pao")
            vps_file = os.path.join(args.data_path, "VPS", vps + ".vps")
            for f in (pao_file, vps_file):
                if not os.path.isfile(f):
                    print("WARNING: %s not found" % f, file=sys.stderr)
        species[el] = (pao, vps, val)

    # ---------------- moments
    moments = [0.0] * natoms
    if args.magmom:
        moments = expand_magmom(args.magmom, natoms)
    else:
        per_el = {k: float(v) for k, v in parse_pairs(args.moment, "--moment").items()}
        for k in per_el:
            if k not in elements:
                sys.exit("ERROR: --moment for %s, but the structure has no %s atom" % (k, k))
        afm = set(args.afm or [])
        counter = {}
        for i, s in enumerate(symbols):
            m = per_el.get(s, 0.0)
            if s in afm:
                n = counter.get(s, 0)
                m = m if n % 2 == 0 else -m
                counter[s] = n + 1
            moments[i] = m
    spin = args.spin or ("on" if any(abs(m) > 1e-12 for m in moments) else "off")
    if spin == "off" and any(abs(m) > 1e-12 for m in moments):
        print("WARNING: moments given but --spin off; they are ignored (up = down)",
              file=sys.stderr)

    # ---------------- cell and coordinates
    if periodic:
        cell = np.asarray(st["cell"], dtype=float)
        frac = np.linalg.solve(cell.T, np.asarray(st["positions"]).T).T
        frac = frac - np.floor(frac + 1e-10)
        coords = frac
        unit = "FRAC"
        solver = "Band"
        kgrid = args.kgrid or kgrid_from_spacing(cell, args.kspacing, st["pbc"])
    else:
        pos = np.asarray(st["positions"], dtype=float)
        lo, hi = pos.min(axis=0), pos.max(axis=0)
        size = hi - lo + 2.0 * args.vacuum
        cell = np.diag(size)
        coords = pos - lo + args.vacuum
        unit = "Ang"
        solver = "Cluster"
        kgrid = None

    name = args.name
    if not name:
        name = os.path.splitext(os.path.basename(args.output))[0] if args.output else "omx"
    name = re.sub(r"\s+", "_", name)

    L = []
    w = L.append
    w("#")
    w("# Generated by omx_make_input.py (DFT-Tutorials) from %s on %s"
      % (os.path.basename(args.structure), datetime.date.today().isoformat()))
    w("# Target: OpenMX %s, database DFT_DATA19, basis level '%s'." % (args.omx_version, args.level))
    w("# Review every keyword and run convergence tests before production calculations.")
    w("#")
    w("")
    w("System.CurrrentDirectory    ./          # official spelling has three r's")
    w("System.Name                 %s" % name)
    w("DATA.PATH                   %s   # directory that contains PAO/ and VPS/" % args.data_path)
    w("level.of.stdout             1")
    w("level.of.fileout            %d" % args.fileout)
    w("")
    w("Species.Number              %d" % len(elements))
    w("<Definition.of.Atomic.Species")
    for el in elements:
        pao, vps, val = species[el]
        w("  %-4s %-20s %-12s # valence %g" % (el, pao, vps, val))
    w("Definition.of.Atomic.Species>")
    w("")
    w("Atoms.Number                %d" % natoms)
    w("Atoms.SpeciesAndCoordinates.Unit   %s" % unit)
    w("<Atoms.SpeciesAndCoordinates")
    tot_up = tot_dn = 0.0
    for i, s in enumerate(symbols):
        up, dn = occupation(species[s][2], moments[i] if spin != "off" else 0.0, "atom %d (%s)" % (i + 1, s))
        tot_up += up
        tot_dn += dn
        x, y, z = coords[i]
        w("  %4d  %-4s %14.9f %14.9f %14.9f   %8.4f %8.4f" % (i + 1, s, x, y, z, up, dn))
    w("Atoms.SpeciesAndCoordinates>")
    w("Atoms.UnitVectors.Unit      Ang")
    w("<Atoms.UnitVectors")
    for v in cell:
        w("  %16.10f %16.10f %16.10f" % tuple(v))
    w("Atoms.UnitVectors>")
    w("")
    w("scf.XcType                  %s" % args.xc)
    w("scf.SpinPolarization        %s" % spin)
    w("scf.ElectronicTemperature   %.1f      # K" % args.etemp)
    w("scf.energycutoff            %.1f      # Ry (real-space grid); test convergence" % args.ecut)
    w("scf.maxIter                 %d" % args.maxiter)
    w("scf.EigenvalueSolver        %s" % solver)
    if kgrid is not None:
        w("scf.Kgrid                   %d %d %d" % tuple(kgrid))
        if args.omx_version == "4.0":
            w("scf.Generation.Kpoint       GCenter     # 4.0 default; use Regular to reproduce 3.9")
    if abs(args.charge) > 0:
        w("scf.system.charge           %g     # positive = hole doping" % args.charge)
    w("scf.Mixing.Type             RMM-DIISK")
    w("scf.Init.Mixing.Weight      0.10")
    w("scf.Min.Mixing.Weight       0.001")
    w("scf.Max.Mixing.Weight       0.30")
    w("scf.Mixing.History          30")
    w("scf.Mixing.StartPulay       10")
    w("scf.Mixing.EveryPulay       1")
    w("scf.criterion               %.1e     # Hartree" % args.criterion)
    w("")
    md = args.md_type
    w("MD.Type                     %s" % md)
    if md.lower() != "nomd":
        w("MD.maxIter                  200")
        w("MD.Opt.criterion            3.0e-4      # Hartree/Bohr")
    else:
        w("MD.maxIter                  1")
    if kpath_lines:
        w("")
        w("Band.dispersion             on")
        w("Band.Nkpath                 %d" % len(kpath_lines))
        w("<Band.kpath")
        for npts, k1, k2, l1, l2 in kpath_lines:
            w("  %4d   %9.6f %9.6f %9.6f   %9.6f %9.6f %9.6f   %s %s"
              % (npts, k1[0], k1[1], k1[2], k2[0], k2[1], k2[2], l1, l2))
        w("Band.kpath>")
    if args.dos:
        w("")
        w("Dos.fileout                 on")
        w("Dos.Erange                  -20.0  20.0   # eV, relative to the chemical potential")
        dk = [2 * n for n in kgrid] if kgrid else [1, 1, 1]
        w("Dos.Kgrid                   %d %d %d" % tuple(dk))
    text = "\n".join(L) + "\n"

    if args.output:
        with open(args.output, "w", encoding="ascii", newline="\n") as fh:
            fh.write(text)
        print("wrote %s" % args.output, file=sys.stderr)
    else:
        sys.stdout.write(text)
    print("atoms %d, species %d, electrons %.2f (up %.2f / down %.2f), initial moment %.2f muB%s"
          % (natoms, len(elements), tot_up + tot_dn, tot_up, tot_dn, tot_up - tot_dn,
             (", scf.Kgrid %s" % " ".join(map(str, kgrid))) if kgrid else ", Cluster"),
          file=sys.stderr)
    if not data_local:
        print("NOTE: DATA.PATH '%s' does not exist here; PAO/VPS files were not checked." % args.data_path,
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
