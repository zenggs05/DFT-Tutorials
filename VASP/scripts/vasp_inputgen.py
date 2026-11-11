#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vasp_inputgen.py -- generate VASP input files from a POSCAR for common tasks.
                    根据 POSCAR 和任务类型生成带注释的 INCAR / KPOINTS / KPOINTS_OPT / POTCAR。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Purpose
-------
* INCAR      : task presets (relax/scf/band/dos/hse/dftu/soc/neb/phonon/md/mlff/
               optics/elastic) with a short comment on every tag (English by
               default = pure ASCII; ``--lang zh`` writes Chinese comments).
* KPOINTS    : automatic mesh from a k-point spacing (same N_i formula as the
               VASP KSPACING tag); Gamma-centred mesh is chosen automatically for
               hexagonal/trigonal cells, 2D slabs (vacuum => 1 k-point along that
               axis) and molecules (Gamma only).
* Line-mode  : high-symmetry path (ASE by default, seekpath optional) for band
               structures, written either as KPOINTS (ICHARG=11) or KPOINTS_OPT
               (VASP >= 6.3), or appended with zero weights to an IBZKPT for
               hybrid-functional band structures.
* POTCAR     : concatenated from ``$VASP_PP_PATH`` following the VASP
               recommended PAW list (October 2026 revision; ``--pp-flavor
               legacy`` for the earlier list, ``gw`` for _GW potentials).  The
               repository never ships POTCAR files.

Usage examples
--------------
    python vasp_inputgen.py relax                       # POSCAR in cwd
    python vasp_inputgen.py scf -p POSCAR -o scf --kspacing 0.2
    python vasp_inputgen.py band -o band                # ICHARG=11 + line-mode KPOINTS
    python vasp_inputgen.py band --kpoints-opt -o band  # SCF + KPOINTS_OPT (VASP >= 6.3)
    python vasp_inputgen.py dos --ispin 2 --magmom "2*3.0 4*0"
    python vasp_inputgen.py hse --band --ibzkpt ../pbe/IBZKPT   # zero-weight k-points
    python vasp_inputgen.py dftu --hubbard "Ni=6.2" --magmom "2 -2 2*0"
    python vasp_inputgen.py soc --saxis "0 0 1" --magmom "0 0 2.2"
    python vasp_inputgen.py neb --images 5 -p 00/POSCAR
    python vasp_inputgen.py md --temp 1000 --nsw 5000 --ensemble nvt
    python vasp_inputgen.py mlff --ml-mode train --temp 600
    python vasp_inputgen.py optics --optics-mode loptics
    python vasp_inputgen.py relax --stdout              # print, write nothing
    python vasp_inputgen.py -h

Dependencies
------------
python >= 3.9, numpy.  Optional: ase (k-path, Bravais lattice), seekpath
(``--path-backend seekpath``), spglib (crystal-system detection).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

__version__ = "1.0.0"

# ---------------------------------------------------------------------------
# Recommended PAW potentials (potpaw.64, PBE, "standard DFT" list valid since
# October 2026 according to https://www.vasp.at/wiki/index.php/Recommended_pseudopotentials).
# Values: element -> (POTCAR name, ZVAL, ENMAX [eV]).  ZVAL/ENMAX are only used
# as a fallback when no POTCAR can be read.
# ---------------------------------------------------------------------------
RECOMMENDED_PBE = {
    "H": ("H", 1, 250.0), "He": ("He", 2, 478.896), "Li": ("Li_sv", 3, 499.034),
    "Be": ("Be_sv", 4, 308.768), "B": ("B", 3, 318.614), "C": ("C", 4, 400.0),
    "N": ("N", 5, 400.0), "O": ("O", 6, 400.0), "F": ("F", 7, 400.0),
    "Ne": ("Ne", 8, 343.606), "Na": ("Na_pv", 7, 259.561), "Mg": ("Mg_pv", 8, 403.929),
    "Al": ("Al", 3, 240.3), "Si": ("Si", 4, 245.345), "P": ("P", 5, 255.04),
    "S": ("S", 6, 258.689), "Cl": ("Cl", 7, 262.472), "Ar": ("Ar", 8, 266.408),
    "K": ("K_sv", 9, 259.264), "Ca": ("Ca_sv", 10, 266.622), "Sc": ("Sc_sv", 11, 222.66),
    "Ti": ("Ti_sv", 12, 274.61), "V": ("V_sv", 13, 263.673), "Cr": ("Cr_pv", 12, 265.681),
    "Mn": ("Mn_pv", 13, 269.864), "Fe": ("Fe_pv", 14, 293.238), "Co": ("Co_pv", 15, 271.042),
    "Ni": ("Ni_pv", 16, 367.986), "Cu": ("Cu_pv", 17, 368.648), "Zn": ("Zn", 12, 276.723),
    "Ga": ("Ga_d", 13, 282.691), "Ge": ("Ge_d", 14, 310.294), "As": ("As", 5, 208.702),
    "Se": ("Se", 6, 211.555), "Br": ("Br", 7, 216.285), "Kr": ("Kr", 8, 185.331),
    "Rb": ("Rb_sv", 9, 220.112), "Sr": ("Sr_sv", 10, 229.353), "Y": ("Y_sv", 11, 202.626),
    "Zr": ("Zr_sv", 12, 229.898), "Nb": ("Nb_sv", 13, 293.235), "Mo": ("Mo_sv", 14, 242.676),
    "Tc": ("Tc_pv", 13, 263.523), "Ru": ("Ru_pv", 14, 240.049), "Rh": ("Rh_pv", 15, 247.408),
    "Pd": ("Pd", 10, 250.925), "Ag": ("Ag", 11, 249.844), "Cd": ("Cd", 12, 274.336),
    "In": ("In_d", 13, 239.211), "Sn": ("Sn_d", 14, 241.083), "Sb": ("Sb", 5, 172.069),
    "Te": ("Te", 6, 174.982), "I": ("I", 7, 175.647), "Xe": ("Xe", 8, 153.118),
    "Cs": ("Cs_sv", 9, 220.318), "Ba": ("Ba_sv", 10, 187.181), "La": ("La", 11, 219.292),
    "Ce": ("Ce", 12, 273.042), "Pr": ("Pr_3", 11, 181.719), "Nd": ("Nd_3", 11, 182.619),
    "Pm": ("Pm_3", 11, 176.959), "Sm": ("Sm_3", 11, 177.087), "Eu": ("Eu_3", 9, 129.057),
    "Gd": ("Gd_3", 9, 154.332), "Tb": ("Tb_3", 9, 155.613), "Dy": ("Dy_3", 9, 155.713),
    "Ho": ("Ho_3", 9, 154.137), "Er": ("Er_3", 9, 155.037), "Tm": ("Tm_3", 9, 149.221),
    "Yb": ("Yb_3", 9, 188.359), "Lu": ("Lu_3", 9, 154.992), "Hf": ("Hf_pv", 10, 220.334),
    "Ta": ("Ta_pv", 11, 223.667), "W": ("W_sv", 14, 223.057), "Re": ("Re_pv", 13, 226.216),
    "Os": ("Os_pv", 14, 228.022), "Ir": ("Ir", 9, 210.864), "Pt": ("Pt", 10, 230.283),
    "Au": ("Au", 11, 229.943), "Hg": ("Hg", 12, 233.204), "Tl": ("Tl_d", 13, 237.053),
    "Pb": ("Pb_d", 14, 237.835), "Bi": ("Bi_d", 15, 242.839), "Po": ("Po_d", 16, 264.565),
    "At": ("At", 7, 161.43), "Rn": ("Rn", 8, 151.497), "Fr": ("Fr_sv", 9, 214.54),
    "Ra": ("Ra_sv", 10, 237.367), "Ac": ("Ac", 11, 172.351), "Th": ("Th", 12, 247.306),
    "Pa": ("Pa", 13, 252.193), "U": ("U", 14, 252.502), "Np": ("Np", 15, 254.26),
    "Pu": ("Pu", 16, 254.353), "Am": ("Am", 17, 255.875), "Cm": ("Cm", 18, 257.953),
}

# Recommendation valid before October 2026 (only the entries that changed).
LEGACY_PBE = {
    "Be": ("Be", 2, 247.543), "Mg": ("Mg", 2, 200.0), "Fe": ("Fe", 8, 267.882),
    "Co": ("Co", 9, 267.968), "Ni": ("Ni", 10, 269.532), "Cu": ("Cu", 11, 295.446),
    "Eu": ("Eu_2", 8, 99.328), "Yb": ("Yb_2", 8, 112.578), "Re": ("Re", 7, 226.216),
    "Os": ("Os", 8, 228.022),
}

# Recommended _GW potentials (for optics / many unoccupied states); name, ENMAX.
GW_PBE = {
    "H": ("H_GW", 300.0), "He": ("He_GW", 405.78), "Li": ("Li_sv_GW", 433.699),
    "Be": ("Be_sv_GW", 537.454), "B": ("B_GW", 318.614), "C": ("C_GW", 413.992),
    "N": ("N_GW", 420.902), "O": ("O_GW", 414.635), "F": ("F_GW", 487.698),
    "Ne": ("Ne_GW", 432.275), "Na": ("Na_sv_GW", 372.853), "Mg": ("Mg_sv_GW", 429.893),
    "Al": ("Al_GW", 240.3), "Si": ("Si_GW", 245.345), "P": ("P_GW", 255.04),
    "S": ("S_GW", 258.689), "Cl": ("Cl_GW", 262.472), "Ar": ("Ar_GW", 290.599),
    "K": ("K_sv_GW", 248.998), "Ca": ("Ca_sv_GW", 281.43), "Sc": ("Sc_sv_GW", 378.961),
    "Ti": ("Ti_sv_GW", 383.774), "V": ("V_sv_GW", 382.321), "Cr": ("Cr_sv_GW", 384.932),
    "Mn": ("Mn_sv_GW", 384.627), "Fe": ("Fe_sv_GW", 387.837), "Co": ("Co_sv_GW", 387.491),
    "Ni": ("Ni_sv_GW", 389.645), "Cu": ("Cu_sv_GW", 467.331), "Zn": ("Zn_sv_GW", 401.665),
    "Ga": ("Ga_d_GW", 404.602), "Ge": ("Ge_d_GW", 375.434), "As": ("As_GW", 208.702),
    "Se": ("Se_GW", 211.555), "Br": ("Br_GW", 216.285), "Kr": ("Kr_GW", 252.232),
    "Rb": ("Rb_sv_GW", 221.197), "Sr": ("Sr_sv_GW", 224.817), "Y": ("Y_sv_GW", 339.758),
    "Zr": ("Zr_sv_GW", 346.364), "Nb": ("Nb_sv_GW", 353.872), "Mo": ("Mo_sv_GW", 344.914),
    "Tc": ("Tc_sv_GW", 351.044), "Ru": ("Ru_sv_GW", 348.106), "Rh": ("Rh_sv_GW", 351.206),
    "Pd": ("Pd_sv_GW", 356.093), "Ag": ("Ag_sv_GW", 354.43), "Cd": ("Cd_sv_GW", 361.806),
    "In": ("In_d_GW", 278.624), "Sn": ("Sn_d_GW", 260.066), "Sb": ("Sb_d_GW", 263.1),
    "Te": ("Te_GW", 174.982), "I": ("I_GW", 175.647), "Xe": ("Xe_GW", 179.547),
    "Cs": ("Cs_sv_GW", 198.101), "Ba": ("Ba_sv_GW", 267.02), "La": ("La_GW", 313.688),
    "Ce": ("Ce_GW", 304.625), "Hf": ("Hf_sv_GW", 309.037), "Ta": ("Ta_sv_GW", 286.008),
    "W": ("W_sv_GW", 317.132), "Re": ("Re_sv_GW", 317.012), "Os": ("Os_sv_GW", 319.773),
    "Ir": ("Ir_sv_GW", 319.843), "Pt": ("Pt_sv_GW", 323.669), "Au": ("Au_sv_GW", 306.658),
    "Hg": ("Hg_sv_GW", 312.028), "Tl": ("Tl_d_GW", 237.053), "Pb": ("Pb_d_GW", 237.809),
    "Bi": ("Bi_d_GW", 261.876), "Po": ("Po_d_GW", 267.847), "At": ("At_d_GW", 266.251),
    "Rn": ("Rn_d_GW", 267.347),
}

D_BLOCK = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
              "Hf Ta W Re Os Ir Pt Au Hg".split())
F_BLOCK = set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
              "Ac Th Pa U Np Pu Am Cm".split())

# Materials-Project-style U values (eV) for oxides/fluorides (Dudarev U_eff).
# Only used when --hubbard is not given for the dftu task; always check the
# literature for your own system.
MP_U = {"Co": 3.32, "Cr": 3.7, "Fe": 5.3, "Mn": 3.9, "Mo": 4.38,
        "Ni": 6.2, "V": 3.25, "W": 6.2}

# Starting magnetic moments (mu_B) used when ISPIN=2 and no --magmom is given:
# high-spin start for magnetic d/f elements, small value elsewhere.
DEFAULT_MAGMOM = {"Cr": 5.0, "Mn": 5.0, "Fe": 5.0, "Ni": 5.0, "V": 5.0, "Mo": 5.0,
                  "W": 5.0, "Ce": 5.0, "Eu": 10.0, "Gd": 7.0, "Co": 0.6}

FUNCTIONALS = {
    # name : (tag, value, VASP version note)
    "PBE": ("GGA", "PE", ""),
    "PBESOL": ("GGA", "PS", ""),
    "RPBE": ("GGA", "RP", ""),
    "REVPBE": ("GGA", "RE", ""),
    "PW91": ("GGA", "91", ""),
    "AM05": ("GGA", "AM", ""),
    "SCAN": ("METAGGA", "SCAN", "VASP >= 5.4.4"),
    "R2SCAN": ("METAGGA", "R2SCAN", "VASP >= 6.2"),
}

VDW = {  # name: (IVDW, note)
    "d2": (1, "DFT-D2 (Grimme 2006)"),
    "d3": (11, "DFT-D3 zero damping"),
    "d3bj": (12, "DFT-D3 Becke-Johnson damping"),
    "d4": (13, "DFT-D4, needs VASP built with the external dftd4 library"),
    "sd3": (15, "simple-DFT-D3 library, VASP >= 6.6.0 built with it; see wiki page simple-DFT-D3"),
    "ts": (20, "Tkatchenko-Scheffler"),
    "mbd": (202, "many-body dispersion MBD@rsSCS"),
}

TASKS = ("relax", "scf", "band", "dos", "hse", "dftu", "soc", "neb",
         "phonon", "md", "mlff", "optics", "elastic")

DEFAULT_KSPACING = {  # 1/Angstrom (VASP convention, 2*pi included)
    "relax": 0.25, "scf": 0.22, "band": 0.22, "dos": 0.15, "hse": 0.30,
    "dftu": 0.22, "soc": 0.20, "neb": 0.30, "phonon": 0.25, "md": 0.50,
    "mlff": 0.50, "optics": 0.15, "elastic": 0.18,
}


def die(msg: str, code: int = 1) -> None:
    """Print an error message and exit."""
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


# ---------------------------------------------------------------------------
# POSCAR handling
# ---------------------------------------------------------------------------
class Structure:
    """Minimal crystal structure container (lattice in Angstrom, fractional coords)."""

    def __init__(self, comment, lattice, species, counts, frac, selective=None):
        self.comment = comment
        self.lattice = np.asarray(lattice, dtype=float)
        self.species = list(species)
        self.counts = list(counts)
        self.frac = np.asarray(frac, dtype=float).reshape(-1, 3)
        self.selective = selective  # list of [bool,bool,bool] or None

    @property
    def natoms(self) -> int:
        return int(sum(self.counts))

    @property
    def symbols(self) -> list:
        out = []
        for s, n in zip(self.species, self.counts):
            out += [s] * n
        return out

    @property
    def cart(self) -> np.ndarray:
        return self.frac @ self.lattice


def _clean_symbol(tok: str) -> str:
    """'Fe_pv' -> 'Fe', 'Si/abc123' -> 'Si'."""
    m = re.match(r"([A-Z][a-z]?)", tok.strip())
    if not m:
        raise ValueError(f"cannot interpret element symbol '{tok}'")
    return m.group(1)


def read_poscar(path, elements=None) -> Structure:
    """Read a VASP 5/6 POSCAR/CONTCAR (VASP 4 files need ``elements``)."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"structure file '{path}' not found")
    lines = p.read_text(errors="replace").splitlines()
    if len(lines) < 8:
        raise ValueError(f"'{path}' is too short to be a POSCAR")
    comment = lines[0].strip()
    try:
        scale = [float(x) for x in lines[1].split()[:3]]
        lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"malformed scale factor / lattice in '{path}': {exc}") from None
    i = 5
    tok = lines[i].split()
    if not tok:
        raise ValueError(f"line 6 of '{path}' is empty")
    if all(re.fullmatch(r"\d+", t) for t in tok):
        species = None                    # VASP 4 format
    else:
        species = [_clean_symbol(t) for t in tok]
        i += 1
    try:
        counts = [int(x) for x in lines[i].split()]
    except ValueError:
        raise ValueError(f"cannot read atom counts on line {i + 1} of '{path}'") from None
    if species is None:
        if elements:
            species = list(elements)
        else:
            # try the comment line (ASE/VASP4 convention)
            cand = comment.split()
            if len(cand) >= len(counts) and all(re.fullmatch(r"[A-Z][a-z]?", c) for c in cand[:len(counts)]):
                species = cand[:len(counts)]
            else:
                raise ValueError("POSCAR has no element line (VASP 4 format); "
                                 "use --elements 'Si O' to specify them")
    if len(species) != len(counts):
        raise ValueError(f"{len(species)} element symbols but {len(counts)} counts in '{path}'")
    i += 1
    selective = False
    if lines[i].strip()[:1] in ("s", "S"):
        selective = True
        i += 1
    cartesian = lines[i].strip()[:1] in ("c", "C", "k", "K")
    i += 1
    nat = sum(counts)
    frac, flags = [], []
    for j in range(nat):
        try:
            t = lines[i + j].split()
            frac.append([float(x) for x in t[:3]])
        except (IndexError, ValueError):
            raise ValueError(f"cannot read coordinates of atom {j + 1} in '{path}'") from None
        if selective:
            f = [x.upper().startswith("T") for x in t[3:6]] if len(t) >= 6 else [True] * 3
            flags.append(f)
    frac = np.array(frac)
    # scaling
    if len(scale) == 1:
        if scale[0] < 0:
            vol = abs(np.linalg.det(lat))
            s = (abs(scale[0]) / vol) ** (1.0 / 3.0)
        else:
            s = scale[0]
        lat_s = lat * s
        if cartesian:
            frac = (frac * s) @ np.linalg.inv(lat_s)
    else:  # three scale factors act on Cartesian x, y, z
        sv = np.array(scale[:3])
        lat_s = lat * sv[None, :]
        if cartesian:
            frac = (frac * sv[None, :]) @ np.linalg.inv(lat_s)
    return Structure(comment, lat_s, species, counts, frac, flags if selective else None)


def write_poscar(path, s: Structure, comment=None) -> None:
    lines = [comment or s.comment or "generated by vasp_inputgen.py", "   1.0"]
    for v in s.lattice:
        lines.append("  " + "  ".join(f"{x:20.12f}" for x in v))
    lines.append("  " + "  ".join(f"{e:>4s}" for e in s.species))
    lines.append("  " + "  ".join(f"{n:>4d}" for n in s.counts))
    if s.selective is not None:
        lines.append("Selective dynamics")
    lines.append("Direct")
    for k, r in enumerate(s.frac):
        row = "  " + "  ".join(f"{x:18.12f}" for x in r)
        if s.selective is not None:
            row += "  " + " ".join("T" if f else "F" for f in s.selective[k])
        lines.append(row)
    Path(path).write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------
def reciprocal(lattice) -> np.ndarray:
    """Reciprocal vectors WITHOUT 2*pi (rows), b_i . a_j = delta_ij."""
    return np.linalg.inv(np.asarray(lattice)).T


def vacuum_axes(s: Structure, threshold: float = 8.0) -> list:
    """Return [bool]*3: True if the largest empty gap along axis i (measured
    perpendicular to the lattice plane) exceeds ``threshold`` Angstrom."""
    rec = reciprocal(s.lattice)
    out = []
    for i in range(3):
        d_perp = 1.0 / np.linalg.norm(rec[i])          # interplanar spacing
        x = np.sort(np.mod(s.frac[:, i], 1.0))
        if len(x) == 1:
            gap = 1.0
        else:
            gaps = np.diff(np.concatenate([x, [x[0] + 1.0]]))
            gap = gaps.max()
        out.append(gap * d_perp > threshold)
    return out


def kmesh(lattice, kspacing: float, vac=(False, False, False)) -> list:
    """N_i = max(1, ceil(|b_i| * 2*pi / KSPACING)), N_i = 1 along vacuum."""
    rec = reciprocal(lattice)
    mesh = []
    for i in range(3):
        if vac[i]:
            mesh.append(1)
        else:
            mesh.append(max(1, int(math.ceil(np.linalg.norm(rec[i]) * 2 * math.pi / kspacing - 1e-8))))
    return mesh


def lattice_family(s: Structure, symprec: float = 1e-3) -> str:
    """Crystal family: 'hexagonal' (incl. trigonal), 'cubic', ... or 'unknown'."""
    try:
        import spglib  # optional
        cell = (s.lattice, s.frac, [s.species.index(x) + 1 for x in s.symbols])
        ds = spglib.get_symmetry_dataset(cell, symprec=symprec)
        num = getattr(ds, "number", None) if ds is not None else None
        if num is None and isinstance(ds, dict):
            num = ds.get("number")
        if num:
            if 143 <= num <= 194:
                return "hexagonal"
            if num >= 195:
                return "cubic"
            if num >= 75:
                return "tetragonal"
            if num >= 16:
                return "orthorhombic"
            if num >= 3:
                return "monoclinic"
            return "triclinic"
    except ImportError:
        pass
    except Exception:  # spglib failure -> fall back
        pass
    try:
        from ase.cell import Cell
        vac = vacuum_axes(s)
        pbc = [True, True, False] if vac == [False, False, True] else True
        name = Cell(s.lattice).get_bravais_lattice(pbc=pbc, eps=max(symprec, 2e-4)).name
        if name in ("HEX", "RHL", "HEX2D"):
            return "hexagonal"
        if name in ("CUB", "FCC", "BCC"):
            return "cubic"
        return name.lower()
    except Exception:
        pass
    # geometric fallback: two equal vectors at 60/120 degrees
    a = s.lattice
    ln = np.linalg.norm(a, axis=1)
    for i, j in ((0, 1), (1, 2), (0, 2)):
        if abs(ln[i] - ln[j]) < 1e-3 * ln[i]:
            ang = math.degrees(math.acos(np.dot(a[i], a[j]) / (ln[i] * ln[j])))
            if abs(ang - 120) < 0.1 or abs(ang - 60) < 0.1:
                return "hexagonal"
    return "unknown"


# ---------------------------------------------------------------------------
# k-point files
# ---------------------------------------------------------------------------
def kpoints_auto_text(mesh, gamma: bool, comment: str) -> str:
    style = "Gamma" if gamma else "Monkhorst-Pack"
    return (f"{comment}\n0\n{style}\n  {mesh[0]} {mesh[1]} {mesh[2]}\n  0 0 0\n")


def _vasp_label(lab: str) -> str:
    return "GAMMA" if lab in ("G", "Gamma", "GAMMA", "\\Gamma") else lab


def kpath_segments(s: Structure, backend: str, pbc, symprec: float):
    """Return (segments, description, new_structure_or_None).

    segments: list of ((label1, k1), (label2, k2)) in fractional reciprocal
    coordinates of the structure that has to be used in the calculation."""
    if backend == "seekpath":
        try:
            import seekpath
        except ImportError:
            die("seekpath is not installed (pip install seekpath) - or use --path-backend ase")
        numbers = [s.species.index(x) + 1 for x in s.symbols]
        res = seekpath.get_path((s.lattice.tolist(), s.frac.tolist(), numbers),
                                with_time_reversal=True, symprec=symprec)
        pc = res["point_coords"]
        segs = [((a, np.array(pc[a])), (b, np.array(pc[b]))) for a, b in res["path"]]
        types = list(res["primitive_types"])
        uniq = []
        for t in types:
            if t not in uniq:
                uniq.append(t)
        # reorder atoms by species
        order = sorted(range(len(types)), key=lambda k: uniq.index(types[k]))
        frac = np.array(res["primitive_positions"])[order]
        species = [s.species[t - 1] for t in uniq]
        counts = [types.count(t) for t in uniq]
        prim = Structure(s.comment + " (seekpath standardized primitive)",
                         res["primitive_lattice"], species, counts, frac)
        desc = f"seekpath {res.get('bravais_lattice_extended', '')}"
        return segs, desc, prim
    # ASE backend (default)
    try:
        from ase.cell import Cell
        from ase.dft.kpoints import parse_path_string
    except ImportError:
        die("ASE is required for the k-path (pip install ase), or use --path-backend seekpath")
    cell = Cell(s.lattice)
    try:
        bp = cell.bandpath(npoints=0, pbc=pbc, eps=max(symprec, 2e-4))
    except Exception as exc:  # e.g. 2D with non-orthogonal third vector
        warn(f"ASE could not build a 2D path ({exc}); falling back to the 3D path")
        bp = cell.bandpath(npoints=0, eps=max(symprec, 2e-4))
    sp = bp.special_points
    segs = []
    for branch in parse_path_string(bp.path):
        for a, b in zip(branch[:-1], branch[1:]):
            segs.append(((a, np.array(sp[a])), (b, np.array(sp[b]))))
    try:
        lat_name = cell.get_bravais_lattice(pbc=pbc).name
    except Exception:
        lat_name = "?"
    return segs, f"ASE {lat_name} path {bp.path}", None


def kpoints_line_text(segs, npts: int, comment: str) -> str:
    out = [comment, str(npts), "Line-mode", "Reciprocal"]
    for (la, ka), (lb, kb) in segs:
        out.append("  " + " ".join(f"{x:12.8f}" for x in ka) + f"   ! {_vasp_label(la)}")
        out.append("  " + " ".join(f"{x:12.8f}" for x in kb) + f"   ! {_vasp_label(lb)}")
        out.append("")
    return "\n".join(out).rstrip() + "\n"


def read_ibzkpt(path):
    """Read an IBZKPT file -> (kpts[N,3], weights[N])."""
    p = Path(path)
    if not p.is_file():
        die(f"IBZKPT file '{path}' not found (run a PBE SCF with the same KPOINTS first)")
    lines = p.read_text().splitlines()
    try:
        n = int(lines[1].split()[0])
    except (IndexError, ValueError):
        die(f"'{path}' does not look like an IBZKPT file (line 2 must be the number of k-points)")
    if not lines[2].strip().lower().startswith("r"):
        die(f"'{path}': only IBZKPT files in reciprocal (fractional) coordinates are supported")
    k, w = [], []
    for ln in lines[3:3 + n]:
        t = ln.split()
        k.append([float(x) for x in t[:3]])
        w.append(float(t[3]))
    if len(k) != n:
        die(f"'{path}' is truncated: expected {n} k-points, found {len(k)}")
    return np.array(k), np.array(w)


def path_points(segs, npts: int):
    """Explicit points along the path (shared end points not duplicated)."""
    pts, labels = [], []
    last = None
    for (la, ka), (lb, kb) in segs:
        for t in np.linspace(0.0, 1.0, npts):
            k = ka + t * (kb - ka)
            if last is not None and np.allclose(k, last):
                continue
            pts.append(k)
            labels.append(la if t == 0 else (lb if t == 1 else ""))
            last = k
    return np.array(pts), labels


def kpoints_zero_weight_text(ibz_k, ibz_w, segs, npts: int) -> str:
    pk, labels = path_points(segs, npts)
    n = len(ibz_k) + len(pk)
    out = [f"IBZKPT ({len(ibz_k)} weighted) + band path ({len(pk)} zero-weight) by vasp_inputgen.py",
           str(n), "Reciprocal lattice"]
    for k, w in zip(ibz_k, ibz_w):
        out.append("  " + " ".join(f"{x:14.10f}" for x in k) + f"  {w:g}")
    for k, lab in zip(pk, labels):
        out.append("  " + " ".join(f"{x:14.10f}" for x in k) + "  0" + (f"   ! {_vasp_label(lab)}" if lab else ""))
    return "\n".join(out) + "\n"


# ---------------------------------------------------------------------------
# POTCAR
# ---------------------------------------------------------------------------
def parse_map(text: str) -> dict:
    """'Fe=Fe_pv,O=O_s' or 'Fe:5.3 Ni:6.2' -> dict."""
    out = {}
    if not text:
        return out
    for item in re.split(r"[,\s;]+", text.strip()):
        if not item:
            continue
        m = re.fullmatch(r"([A-Z][a-z]?)\s*[=:]\s*(\S+)", item)
        if not m:
            die(f"cannot parse '{item}' (expected Element=value, e.g. Fe=5.3)")
        out[m.group(1)] = m.group(2)
    return out


def potcar_names(species, flavor: str, overrides: dict) -> list:
    names = []
    for el in species:
        if el in overrides:
            names.append(overrides[el])
        elif flavor == "gw" and el in GW_PBE:
            names.append(GW_PBE[el][0])
        elif flavor == "legacy" and el in LEGACY_PBE:
            names.append(LEGACY_PBE[el][0])
        elif el in RECOMMENDED_PBE:
            names.append(RECOMMENDED_PBE[el][0])
        else:
            names.append(el)
    return names


def find_pp_dir(pp_path: str, pp_dir: str | None) -> Path:
    root = Path(pp_path).expanduser()
    if not root.is_dir():
        die(f"pseudopotential root '{root}' does not exist (check $VASP_PP_PATH / --pp-path)")
    cands = [pp_dir] if pp_dir else ["potpaw_PBE.64", "potpaw_PBE.54", "potpaw_PBE", "PBE.64", "PBE"]
    for c in cands:
        if c and (root / c).is_dir():
            return root / c
    if (root / "Si").is_dir() or any(root.glob("*/POTCAR")):
        return root  # VASP_PP_PATH points directly to the set
    die(f"none of {cands} found below '{root}'; use --pp-dir to name the PBE set directory")
    return root


def _read_potcar_file(f: Path) -> str:
    if f.name.endswith(".gz"):
        with gzip.open(f, "rt", errors="replace") as fh:
            return fh.read()
    return f.read_text(errors="replace")


def build_potcar(names, pp_dir: Path):
    """Concatenate POTCARs; return (text, info list of dicts)."""
    texts, info = [], []
    for nm in names:
        d = pp_dir / nm
        f = None
        for cand in ("POTCAR", "POTCAR.gz"):
            if (d / cand).is_file():
                f = d / cand
                break
        if f is None:
            if (d / "POTCAR.Z").is_file():
                die(f"{d}/POTCAR.Z is compressed with 'compress'; run 'uncompress {d}/POTCAR.Z' first")
            die(f"POTCAR for '{nm}' not found in {pp_dir} (override with --potcar-map {nm.split('_')[0]}=NAME)")
        t = _read_potcar_file(f)
        if not t.endswith("\n"):
            t += "\n"
        texts.append(t)
        info.append(potcar_header_info(t, nm))
    return "".join(texts), info


def potcar_header_info(t: str, name: str = "?") -> dict:
    titel = re.search(r"TITEL\s*=\s*(.+)", t)
    enmax = re.search(r"ENMAX\s*=\s*([-\d.]+)", t)
    zval = re.search(r"ZVAL\s*=\s*([-\d.]+)", t)
    return {"name": name,
            "titel": titel.group(1).strip() if titel else "?",
            "enmax": float(enmax.group(1)) if enmax else None,
            "zval": float(zval.group(1)) if zval else None}


def read_existing_potcar(path: Path) -> list:
    text = path.read_text(errors="replace")
    blocks = re.split(r"End of Dataset", text)
    info = []
    for b in blocks:
        if "TITEL" in b:
            info.append(potcar_header_info(b))
    return info


# ---------------------------------------------------------------------------
# MAGMOM helpers
# ---------------------------------------------------------------------------
def expand_magmom(text: str) -> list:
    """'2*0 2 -2' or '2*0_2_-2' -> [0,0,2,-2]."""
    vals = []
    for tok in re.split(r"[\s_,]+", text.strip()):
        if not tok:
            continue
        if "*" in tok:
            n, v = tok.split("*", 1)
            try:
                vals += [float(v)] * int(n)
            except ValueError:
                die(f"bad MAGMOM item '{tok}'")
        else:
            try:
                vals.append(float(tok))
            except ValueError:
                die(f"bad MAGMOM item '{tok}'")
    return vals


def compress_values(vals) -> str:
    """[0,0,2,-2] -> '2*0 2 -2' (VASP N*value syntax)."""
    out = []
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and abs(vals[j + 1] - vals[i]) < 1e-10:
            j += 1
        v = f"{vals[i]:g}"
        out.append(f"{j - i + 1}*{v}" if j > i else v)
        i = j + 1
    return " ".join(out)


# ---------------------------------------------------------------------------
# INCAR builder
# ---------------------------------------------------------------------------
SECTIONS = [
    ("start", "Start / system", "启动与体系"),
    ("elec", "Electronic minimisation", "电子步"),
    ("smear", "Smearing / BZ integration", "展宽与布里渊区积分"),
    ("xc", "Exchange-correlation", "交换关联泛函"),
    ("mag", "Spin / magnetism", "自旋与磁性"),
    ("ldau", "DFT+U (Dudarev)", "DFT+U（Dudarev）"),
    ("soc", "Spin-orbit coupling (run vasp_ncl)", "自旋轨道耦合（用 vasp_ncl）"),
    ("hf", "Hybrid functional", "杂化泛函"),
    ("vdw", "van der Waals correction", "范德华修正"),
    ("ion", "Ionic motion", "离子步"),
    ("neb", "Nudged elastic band", "NEB 过渡态"),
    ("md", "Molecular dynamics", "分子动力学"),
    ("ml", "Machine-learned force field (VASP >= 6.4 syntax)", "机器学习力场（VASP >= 6.4 写法）"),
    ("resp", "Linear response / optics", "线性响应与光学"),
    ("out", "Output", "输出控制"),
    ("kp", "k-points via KSPACING", "KSPACING 自动 k 点"),
    ("par", "Parallelisation", "并行设置"),
]


class Incar:
    def __init__(self):
        self.tags = {}

    def set(self, tag, value, en="", zh="", sec="elec"):
        self.tags[tag] = [str(value), en, zh or en, sec]

    def drop(self, tag):
        self.tags.pop(tag, None)

    def has(self, tag):
        return tag in self.tags

    def value(self, tag):
        return self.tags[tag][0] if tag in self.tags else None

    def render(self, lang: str, header) -> str:
        out = [f"# {h}" for h in header]
        width = max([len(t) for t in self.tags] + [6])
        for key, en_t, zh_t in SECTIONS:
            items = [(t, v) for t, v in self.tags.items() if v[3] == key]
            if not items:
                continue
            title = zh_t if lang == "zh" else en_t
            out.append("")
            out.append(f"# ---- {title} ----")
            for tag, (val, en, zh, _sec) in items:
                cmt = zh if lang == "zh" else en
                line = f"{tag:<{width}} = {val}"
                if cmt and lang != "none":
                    line = f"{line:<{width + 22}} # {cmt}"
                out.append(line)
        return "\n".join(out) + "\n"


def round_up(x: float, step: float = 10.0) -> int:
    return int(math.ceil(x / step - 1e-9) * step)


def build_incar(a, s: Structure, info: dict) -> Incar:
    """Assemble the INCAR for task ``a.task`` (``info`` holds derived data)."""
    t = a.task
    inc = Incar()
    nions = s.natoms
    mol = info["molecule"]
    is2d = info["is2d"]
    el = set(s.species)
    has_d = bool(el & D_BLOCK)
    has_f = bool(el & F_BLOCK)

    inc.set("SYSTEM", info["system"], "job title", "任务名称", "start")
    if t in ("band", "dos") and not (t == "band" and a.kpoints_opt):
        inc.set("ICHARG", 11, "non-SCF: read CHGCAR of the SCF run, keep density fixed",
                "非自洽：读取自洽步的 CHGCAR 并固定电荷密度", "start")
        inc.set("ISTART", 0, "do not read WAVECAR", "不读 WAVECAR", "start")
    if t == "hse":
        inc.set("ISTART", 1, "start from a PBE WAVECAR if present (faster, more stable)",
                "若有 PBE 的 WAVECAR 则从它开始（更快更稳）", "start")
    if t == "soc" and a.soc_from_chgcar:
        inc.set("ICHARG", 11, "read collinear CHGCAR, keep it fixed (force theorem / MAE)",
                "读取共线计算的 CHGCAR 并固定（力定理/MAE）", "start")

    # ---- electronic ----
    inc.set("PREC", "Accurate", "accurate FFT grids, avoids wrap-around errors",
            "精确 FFT 网格，避免混叠误差")
    enc_note = info["encut_note"]
    inc.set("ENCUT", info["encut"], enc_note, enc_note.replace("x max(ENMAX)", "x 最大 ENMAX"))
    ediff = {"relax": "1E-6", "scf": "1E-6", "band": "1E-6", "dos": "1E-6", "hse": "1E-6",
             "dftu": "1E-6", "soc": "1E-7", "neb": "1E-6", "phonon": "1E-8", "md": "1E-5",
             "mlff": "1E-5", "optics": "1E-7", "elastic": "1E-7"}[t]
    if a.ediff:
        ediff = a.ediff
    inc.set("EDIFF", ediff, "SCF convergence (eV); phonons/elastic need 1E-7..1E-8",
            "电子步收敛判据 (eV)；声子/弹性要 1E-7~1E-8")
    if t in ("hse", "dftu", "soc", "optics") or a.hubbard:
        inc.set("NELM", 120, "allow more SCF steps (default 60)", "允许更多电子步（默认 60）")
    if t in ("md", "mlff"):
        inc.set("NELMIN", 4, "at least 4 SCF steps per MD step", "每个 MD 步至少 4 个电子步")
    if t == "hse":
        if a.metal:
            inc.set("ALGO", "Normal", "blocked Davidson (robust for metals with HF)",
                    "分块 Davidson（金属杂化计算更稳）")
        else:
            inc.set("ALGO", "Damped", "damped velocity-friction algorithm for hybrids",
                    "杂化泛函常用的阻尼算法")
            inc.set("TIME", 0.4, "step for ALGO=Damped", "ALGO=Damped 的步长")
    elif t == "optics":
        inc.set("ALGO", "Normal", "Davidson: reliable for many empty bands", "Davidson：空带多时更可靠")
    lreal = ".FALSE."
    if nions > 30 and t in ("relax", "neb", "md", "mlff", "scf", "dftu"):
        lreal = "Auto"
    if a.lreal:
        lreal = a.lreal
    inc.set("LREAL", lreal,
            "real-space projection; Auto for >~30 atoms, .FALSE. for accurate energies",
            "实空间投影；约 30 原子以上用 Auto，精确能量用 .FALSE.")
    inc.set("LASPH", ".TRUE.", "aspherical PAW gradient corrections (d/f, meta-GGA, vdW-DF)",
            "PAW 球内非球形梯度修正（d/f 元素、meta-GGA、vdW-DF 必需）")
    if t in ("phonon", "elastic"):
        inc.set("ADDGRID", ".TRUE.", "extra support grid for augmentation charges (smoother forces)",
                "增强电荷的附加网格（受力更平滑）")

    # ---- smearing ----
    ismear, sigma = 0, 0.05
    if mol:
        ismear, sigma = 0, 0.01
    elif t == "dos":
        if info["nk_total"] >= 4:
            ismear, sigma = -5, None
        else:
            warn("fewer than 4 k-points: tetrahedron method (ISMEAR=-5) impossible, using ISMEAR=0")
    elif a.metal:
        ismear, sigma = (1, 0.2) if t in ("relax", "neb", "md", "mlff") else (1, 0.1)
    elif t == "phonon":
        sigma = 0.01
    if a.ismear is not None:
        ismear = a.ismear
    if a.sigma is not None:
        sigma = a.sigma
    note = {-5: "tetrahedron + Bloechl (DOS / accurate total energy; needs >= 4 k-points)",
            0: "Gaussian smearing (semiconductors, insulators, molecules)",
            1: "Methfessel-Paxton order 1 (metals); never for semiconductors",
            2: "Methfessel-Paxton order 2 (metals)",
            -1: "Fermi-Dirac smearing"}.get(ismear, "smearing method")
    note_zh = {-5: "四面体+Bloechl 修正（DOS/精确总能，需 >= 4 个 k 点）",
               0: "Gaussian 展宽（半导体、绝缘体、分子）",
               1: "Methfessel-Paxton 一阶（金属），半导体勿用",
               2: "Methfessel-Paxton 二阶（金属）",
               -1: "Fermi-Dirac 展宽"}.get(ismear, "展宽方法")
    inc.set("ISMEAR", ismear, note, note_zh, "smear")
    if sigma is not None and ismear != -5:
        inc.set("SIGMA", f"{sigma:g}", "smearing width (eV); check T*S/atom < 1 meV in OUTCAR",
                "展宽宽度 (eV)；检查 OUTCAR 中 T*S/原子 < 1 meV", "smear")

    # ---- exchange-correlation ----
    if t != "hse":
        f = a.functional.upper()
        tag, val, ver = FUNCTIONALS[f]
        if tag == "GGA":
            inc.set("GGA", val, f"{a.functional} (use PBE POTCARs)", f"{a.functional}（使用 PBE 赝势）", "xc")
        else:
            inc.set("METAGGA", val, f"{a.functional} meta-GGA, {ver}", f"{a.functional} meta-GGA，{ver}", "xc")
            inc.set("LMIXTAU", ".TRUE.", "also mix the kinetic-energy density (helps meta-GGA SCF)",
                    "同时混合动能密度（有助于 meta-GGA 收敛）", "xc")
    else:
        inc.set("GGA", "PE", "PBE base for HSE06", "HSE06 以 PBE 为基础", "xc")

    # ---- spin / magnetism ----
    soc = t == "soc"
    magmom_vals = expand_magmom(a.magmom) if a.magmom else None
    want_spin = (a.ispin == 2) or (magmom_vals is not None) or t == "dftu" or bool(a.hubbard and a.ispin != 1)
    if a.ispin == 1:
        want_spin = False
    if soc:
        if magmom_vals is None:
            magmom_vals = []
            for sym in s.symbols:
                magmom_vals += [0.0, 0.0, DEFAULT_MAGMOM.get(sym, 0.0)]
        elif len(magmom_vals) == nions:
            mv = []
            for m in magmom_vals:
                mv += [0.0, 0.0, m]
            magmom_vals = mv
        if len(magmom_vals) != 3 * nions:
            die(f"SOC needs 3 MAGMOM values per atom ({3 * nions}) or one per atom ({nions}); got {len(magmom_vals)}")
        inc.set("MAGMOM", compress_values(magmom_vals), "mx my mz per atom (in SAXIS frame)",
                "每个原子 mx my mz（SAXIS 坐标系）", "mag")
    elif want_spin:
        inc.set("ISPIN", 2, "spin-polarised", "自旋极化", "mag")
        if magmom_vals is None:
            magmom_vals = [DEFAULT_MAGMOM.get(sym, 0.6) for sym in s.symbols]
            warn("no --magmom given: ferromagnetic high-spin start " + compress_values(magmom_vals)
                 + " (give --magmom for AFM orderings)")
        if len(magmom_vals) != nions:
            die(f"MAGMOM must have one value per atom ({nions}); got {len(magmom_vals)}")
        inc.set("MAGMOM", compress_values(magmom_vals), "initial moments, one per atom (POSCAR order)",
                "初始磁矩，每个原子一个（与 POSCAR 顺序一致）", "mag")
    if a.nupdown is not None:
        inc.set("NUPDOWN", a.nupdown, "fixed total moment N_up - N_down", "固定总磁矩", "mag")
    if (has_d or has_f) and (t in ("scf", "band", "dos", "dftu", "soc", "hse") or a.hubbard):
        inc.set("LMAXMIX", 6 if has_f else 4, "l-cutoff of mixed/written density (4 for d, 6 for f)",
                "混合/写入 CHGCAR 的角动量上限（d 用 4，f 用 6）", "mag")

    # ---- DFT+U ----
    if t == "dftu" or a.hubbard:
        u = {k: float(v) for k, v in parse_map(a.hubbard).items()} if a.hubbard else {}
        if not u:
            u = {e: MP_U[e] for e in s.species if e in MP_U}
            if not u:
                die("no element of this structure has a default U; give --hubbard 'Cu=4.0'")
            warn("no --hubbard given: using Materials-Project-style U for oxides: "
                 + ", ".join(f"{k}={v}" for k, v in u.items()))
        bad = [k for k in u if k not in s.species]
        if bad:
            die(f"--hubbard elements {bad} are not in the POSCAR ({s.species})")
        ldaul = [(3 if e in F_BLOCK else 2) if u.get(e, 0) != 0 else -1 for e in s.species]
        inc.set("LDAU", ".TRUE.", "switch on DFT+U", "打开 DFT+U", "ldau")
        inc.set("LDAUTYPE", 2, "Dudarev: only U_eff = U - J matters", "Dudarev 方法，只有 U-J 起作用", "ldau")
        inc.set("LDAUL", " ".join(str(x) for x in ldaul), "l per species: 2=d, 3=f, -1=no U",
                "每种元素的 l：2=d，3=f，-1=不加 U", "ldau")
        inc.set("LDAUU", " ".join(f"{u.get(e, 0.0):g}" for e in s.species), "U (eV) per species",
                "每种元素的 U (eV)", "ldau")
        inc.set("LDAUJ", " ".join("0" for _ in s.species), "J (eV) per species", "每种元素的 J (eV)", "ldau")
        inc.set("LDAUPRINT", 1, "write occupation matrices to OUTCAR", "在 OUTCAR 中输出占据矩阵", "ldau")
        inc.set("LMAXMIX", 6 if any(x == 3 for x in ldaul) else 4,
                "required for DFT+U (4 for d, 6 for f)", "DFT+U 必须（d 用 4，f 用 6）", "mag")

    # ---- SOC ----
    if soc:
        inc.set("LSORBIT", ".TRUE.", "spin-orbit coupling (implies LNONCOLLINEAR)", "自旋轨道耦合（自动非共线）", "soc")
        inc.set("SAXIS", a.saxis, "spin quantisation axis", "自旋量子化轴", "soc")
        inc.set("GGA_COMPAT", ".FALSE.", "full gradient correction symmetry (important for MAE)",
                "保持梯度修正的对称性（MAE 计算重要）", "soc")
        inc.set("ISYM", -1, "no symmetry (magnetisation + SOC break it)", "关闭对称性（磁化+SOC 破坏对称）", "soc")
        inc.set("LORBMOM", ".TRUE.", "write orbital moments (needs LORBIT)", "输出轨道磁矩（需 LORBIT）", "soc")

    # ---- hybrid ----
    if t == "hse":
        inc.set("LHFCALC", ".TRUE.", "hybrid functional", "杂化泛函", "hf")
        inc.set("HFSCREEN", 0.2, "HSE06 screening parameter (1/Angstrom)", "HSE06 屏蔽参数 (1/埃)", "hf")
        inc.set("AEXX", 0.25, "fraction of exact exchange (default 0.25)", "精确交换比例（默认 0.25）", "hf")
        inc.set("PRECFOCK", "Normal", "FFT grid for exact exchange (Fast = cheaper)", "精确交换的 FFT 网格（Fast 更省）", "hf")
        if a.nkred:
            inc.set("NKRED", a.nkred, "down-sample k-mesh for the Fock operator", "Fock 算符 k 点降采样", "hf")

    # ---- vdW ----
    if a.vdw:
        iv, note = VDW[a.vdw]
        inc.set("IVDW", iv, note, note, "vdw")
        if a.vdw in ("d3", "d3bj", "d4", "sd3") and t == "phonon" and a.dfpt:
            warn("DFPT (IBRION=7/8) does not include IVDW corrections; use finite displacements")

    # ---- ionic part per task ----
    if t == "relax":
        isif = a.isif if a.isif is not None else (2 if (mol or is2d) else 3)
        inc.set("IBRION", 2, "conjugate gradient (robust); 1 = RMM-DIIS near minimum", "共轭梯度（稳健）；接近极小用 1", "ion")
        inc.set("ISIF", isif, {2: "relax ions only", 3: "relax ions + cell shape + volume",
                               4: "ions + shape, fixed volume", 7: "volume only"}.get(isif, "ISIF"),
                {2: "只优化原子", 3: "原子+晶胞形状+体积", 4: "原子+形状，体积固定", 7: "只变体积"}.get(isif, "ISIF"), "ion")
        inc.set("NSW", a.nsw or 100, "max ionic steps", "最大离子步数", "ion")
        inc.set("EDIFFG", a.ediffg or "-0.02", "stop when all forces < 0.02 eV/A", "所有原子受力 < 0.02 eV/埃时停止", "ion")
        if a.relax_2d_cell:
            inc.set("ISIF", 3, "relax ions + cell, constrained below", "优化原子+晶胞，下面约束 c 轴", "ion")
            inc.set("LATTICE_CONSTRAINTS", ".TRUE. .TRUE. .FALSE.",
                    "relax a1,a2 only (2D); VASP >= 6.4.3, >= 6.5.0 recommended",
                    "只优化面内 a1,a2；VASP >= 6.4.3（建议 >= 6.5.0）", "ion")
        if isif == 3 or a.relax_2d_cell:
            inc.set("ENCUT", info["encut"], enc_note + "; >= 1.3 x ENMAX against Pulay stress",
                    enc_note + "；变胞需 >= 1.3 x ENMAX 以抑制 Pulay 应力")
    elif t in ("scf", "band", "dos", "hse", "dftu", "soc", "optics"):
        inc.set("NSW", 0, "single point", "单点计算", "ion")
        inc.set("IBRION", -1, "no ionic update", "不移动原子", "ion")
    elif t == "phonon":
        if a.dfpt:
            inc.set("IBRION", 8, "DFPT with symmetry (7 = without); no IVDW, no NCORE",
                    "DFPT（利用对称性，7 为不用）；不支持 IVDW 与 NCORE", "ion")
        else:
            inc.set("IBRION", -1, "single point on a phonopy displaced supercell", "phonopy 位移超胞单点", "ion")
            inc.set("NSW", 0, "single point", "单点计算", "ion")
            inc.set("ISYM", 0, "symmetry off for displaced cells (optional)", "位移超胞关闭对称性（可选）", "ion")
    elif t == "elastic":
        inc.set("IBRION", 6, "finite differences: elastic tensor from strain-stress", "有限差分：应力-应变求弹性张量", "ion")
        inc.set("ISIF", 3, "needed to obtain the elastic tensor", "计算弹性张量必须", "ion")
        inc.set("NFREE", 4, "4 displacements per degree of freedom (2 is cheaper)", "每个自由度 4 个位移（2 更省）", "ion")
        inc.set("POTIM", 0.015, "displacement step (Angstrom)", "位移步长（埃）", "ion")
        inc.set("NSW", 1, "required > 0 for IBRION=6", "IBRION=6 需要 NSW > 0", "ion")
    elif t == "neb":
        if a.no_vtst:
            inc.set("IBRION", 1, "built-in VASP NEB with RMM-DIIS (no climbing image)", "VASP 自带 NEB（无爬坡）", "ion")
        else:
            inc.set("IBRION", 3, "VTST: let IOPT drive the optimiser", "VTST：由 IOPT 优化器控制", "ion")
            inc.set("POTIM", 0, "VTST: must be 0", "VTST：必须为 0", "ion")
            inc.set("IOPT", a.iopt, "1=LBFGS 2=CG 3=QuickMin 7=FIRE (VTST)", "1=LBFGS 2=CG 3=QuickMin 7=FIRE（VTST）", "ion")
        inc.set("ISIF", 2, "fixed cell", "晶胞固定", "ion")
        inc.set("NSW", a.nsw or 200, "max ionic steps", "最大离子步数", "ion")
        inc.set("EDIFFG", a.ediffg or "-0.03", "NEB force criterion (eV/A)", "NEB 受力收敛判据 (eV/埃)", "ion")
        inc.set("IMAGES", a.images, "number of intermediate images (dirs 01..NN)", "中间图像数（目录 01..NN）", "neb")
        inc.set("SPRING", -5, "spring constant (eV/A^2), negative = nudged", "弹簧常数（负值=nudged）", "neb")
        if not a.no_vtst:
            inc.set("ICHAIN", 0, "0 = NEB (VTST)", "0 = NEB（VTST）", "neb")
            inc.set("LCLIMB", ".TRUE." if a.climb else ".FALSE.",
                    "climbing image; first run a plain NEB, then switch on", "爬坡图像；建议先普通 NEB 再打开", "neb")
    elif t in ("md", "mlff"):
        has_h = "H" in el
        potim = a.potim or (0.5 if has_h else (1.0 if t == "md" else 1.5))
        inc.set("IBRION", 0, "molecular dynamics", "分子动力学", "ion")
        inc.set("NSW", a.nsw or (5000 if t == "md" else 20000), "number of MD steps", "MD 步数", "ion")
        inc.set("POTIM", potim, "time step (fs); 0.5-1 fs with H", "时间步长 (fs)；含 H 用 0.5-1 fs", "ion")
        inc.set("TEBEG", a.temp, "start temperature (K)", "起始温度 (K)", "md")
        inc.set("TEEND", a.temp_end if a.temp_end is not None else a.temp, "end temperature (K)", "终止温度 (K)", "md")
        ens = a.ensemble
        if ens == "nvt":
            inc.set("MDALGO", 2, "Nose-Hoover thermostat (NVT)", "Nose-Hoover 热浴（NVT）", "md")
            inc.set("SMASS", 0, "Nose mass giving a 40-step period", "Nose 质量（约 40 步周期）", "md")
            inc.set("ISIF", 2, "fixed cell", "晶胞固定", "md")
        elif ens == "nve":
            inc.set("MDALGO", 1, "Andersen thermostat ...", "Andersen 热浴……", "md")
            inc.set("ANDERSEN_PROB", 0.0, "... with zero collision probability = NVE", "……碰撞概率为 0 即 NVE", "md")
            inc.set("ISIF", 2, "fixed cell", "晶胞固定", "md")
        else:  # npt
            gam = " ".join("10.0" for _ in s.species)
            inc.set("MDALGO", 3, "Langevin thermostat + Parrinello-Rahman (NpT)", "Langevin 热浴 + Parrinello-Rahman（NpT）", "md")
            inc.set("ISIF", 3, "variable cell for NpT", "NpT 需要变胞", "md")
            inc.set("LANGEVIN_GAMMA", gam, "friction (1/ps) per species", "每种元素的摩擦系数 (1/ps)", "md")
            inc.set("LANGEVIN_GAMMA_L", 10.0, "lattice friction (1/ps)", "晶格摩擦系数 (1/ps)", "md")
            inc.set("PMASS", 1000, "fictitious lattice mass", "晶格虚拟质量", "md")
            inc.set("PSTRESS", a.pstress if a.pstress is not None else 0.0, "external pressure (kB)", "外压 (kB)", "md")
            inc.set("ENCUT", info["encut"], enc_note + "; raise ENCUT for NpT (Pulay stress)",
                    enc_note + "；NpT 需提高 ENCUT（Pulay 应力）")
        inc.set("ISYM", 0, "no symmetry in MD", "MD 关闭对称性", "md")
        if t == "md":
            inc.set("PREC", "Normal", "MD: Normal is usually enough", "MD 通常 Normal 足够")
        else:
            inc.set("PREC", "Accurate", "MLFF: accurate reference data", "MLFF：训练数据需要较高精度")
        if t == "mlff":
            inc.set("ML_LMLFF", ".TRUE.", "switch on machine-learned force field (VASP >= 6.3)",
                    "打开机器学习力场（VASP >= 6.3）", "ml")
            mode_note = {"train": "on-the-fly training (copy ML_ABN -> ML_AB to continue)",
                         "select": "re-select local configurations from ML_AB",
                         "refit": "refit ML_AB for fast production force field -> ML_FFN",
                         "run": "production run with ML_FF (copy ML_FFN -> ML_FF)"}[a.ml_mode]
            mode_zh = {"train": "在线训练（续训前 cp ML_ABN ML_AB）",
                       "select": "从 ML_AB 重新挑选局域构型",
                       "refit": "由 ML_AB 重拟合快速力场 -> ML_FFN",
                       "run": "用 ML_FF 做生产计算（cp ML_FFN ML_FF）"}[a.ml_mode]
            inc.set("ML_MODE", a.ml_mode, mode_note + "; VASP 6.3: use ML_ISTART instead",
                    mode_zh + "；VASP 6.3 请改用 ML_ISTART", "ml")
            if a.ml_mode in ("refit", "select"):
                inc.set("NSW", 1, "set automatically by ML_MODE=refit/select", "ML_MODE=refit/select 会自动设为 1", "ion")

    # ---- linear response / optics ----
    if t == "optics":
        if a.optics_mode == "loptics":
            inc.set("LOPTICS", ".TRUE.", "frequency-dependent dielectric function (IPA)", "独立粒子近似频率相关介电函数", "resp")
            inc.set("CSHIFT", 0.1, "complex shift (eV) in Kramers-Kronig", "Kramers-Kronig 复数位移 (eV)", "resp")
            inc.set("NEDOS", 2001, "energy grid points for epsilon(omega)", "介电函数能量网格点数", "resp")
            if info.get("nbands"):
                inc.set("NBANDS", info["nbands"], info["nbands_note"], info["nbands_note"], "resp")
        elif a.optics_mode == "lepsilon":
            inc.set("LEPSILON", ".TRUE.", "static dielectric tensor + Born charges by DFPT (no NCORE)",
                    "DFPT 静态介电张量与 Born 有效电荷（不要设 NCORE）", "resp")
            if a.ionic:
                inc.set("IBRION", 8, "add ionic contribution (DFPT phonons at Gamma)", "加入离子贡献（Gamma 点 DFPT）", "ion")
                inc.drop("NSW")
        else:
            inc.set("LCALCEPS", ".TRUE.", "dielectric tensor from finite electric field (works with hybrids)",
                    "有限电场法介电张量（可配合杂化泛函）", "resp")
            if a.ionic:
                inc.set("IBRION", 6, "add ionic contribution by finite differences", "有限差分加入离子贡献", "ion")
                inc.set("NFREE", 2, "2 displacements per degree of freedom", "每个自由度 2 个位移", "ion")
                inc.set("POTIM", 0.015, "displacement (Angstrom)", "位移 (埃)", "ion")
                inc.set("NSW", 1, "required > 0 for IBRION=6", "IBRION=6 需要 NSW > 0", "ion")
    if t == "phonon" and a.dfpt and a.born:
        inc.set("LEPSILON", ".TRUE.", "Born charges + epsilon_inf for NAC (phonopy-vasp-born)",
                "Born 电荷与 epsilon_inf，供非解析修正", "resp")

    # ---- output ----
    write_chg = t in ("scf", "dftu", "hse") or (t == "band" and a.kpoints_opt) or (t == "soc" and not a.soc_from_chgcar)
    write_wav = a.write_wavecar and t in ("scf", "hse", "dftu")
    inc.set("LCHARG", ".TRUE." if write_chg else ".FALSE.", "write CHGCAR" + (" (needed by band/dos)" if write_chg else ""),
            "写 CHGCAR" + ("（能带/DOS 需要）" if write_chg else ""), "out")
    inc.set("LWAVE", ".TRUE." if write_wav else ".FALSE.", "write WAVECAR", "写 WAVECAR", "out")
    if t in ("scf", "band", "dos", "hse", "dftu", "soc"):
        inc.set("LORBIT", 11, "site- and lm-projected DOS/PROCAR (needs no RWIGS)", "按原子和 lm 投影的 DOS/PROCAR", "out")
    if t == "dos":
        inc.set("NEDOS", 3001, "number of DOS grid points", "DOS 能量网格点数", "out")
    if a.bandgap:
        inc.set("BANDGAP", a.bandgap, "report band gap in OUTCAR (VASP >= 6.4.3)", "在 OUTCAR 中报告带隙（VASP >= 6.4.3）", "out")

    # ---- k-points via KSPACING ----
    if a.kspacing_in_incar:
        inc.set("KSPACING", f"{info['kspacing']:g}", "k-point spacing (1/A, 2*pi included); no KPOINTS file",
                "k 点间距（1/埃，含 2*pi）；不需要 KPOINTS 文件", "kp")
        inc.set("KGAMMA", ".TRUE." if info["gamma"] else ".FALSE.", ".TRUE. = Gamma-centred, .FALSE. = shifted (MP-like)",
                ".TRUE. 为 Gamma 中心，.FALSE. 为平移网格（类 MP）", "kp")

    # ---- parallel ----
    no_ncore = (t == "phonon" and a.dfpt) or (t == "optics" and a.optics_mode in ("lepsilon", "lcalceps"))
    if a.ncore and not no_ncore:
        inc.set("NCORE", a.ncore, "cores per band group; ~sqrt(cores per node), divisor of it",
                "每个能带组的核数；取每节点核数的因子、约其平方根", "par")
    if no_ncore:
        inc.set("NCORE", 1, "DFPT / finite-field routines need NCORE = 1 (default)",
                "DFPT/有限电场需要 NCORE = 1（默认）", "par")
    if a.kpar:
        inc.set("KPAR", a.kpar, "k-point groups; must divide the number of MPI ranks",
                "k 点并行组数；须整除 MPI 进程数", "par")
    return inc


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="vasp_inputgen.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Generate INCAR / KPOINTS / KPOINTS_OPT / POTCAR from a POSCAR.\n"
                     "根据 POSCAR 与任务类型生成 VASP 输入文件。\n\n"
                     "Tasks: " + ", ".join(TASKS)),
        epilog=("Examples / 示例:\n"
                "  python vasp_inputgen.py relax -p POSCAR -o relax\n"
                "  python vasp_inputgen.py band --kpoints-opt -o band\n"
                "  python vasp_inputgen.py hse --band --ibzkpt ../scf/IBZKPT\n"
                "  python vasp_inputgen.py dftu --hubbard 'Ni=6.2' --magmom '2 -2 2*0'\n"
                "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License"))
    p.add_argument("task", choices=TASKS, help="calculation type / 任务类型")
    g = p.add_argument_group("files / 文件")
    g.add_argument("-p", "--poscar", default="POSCAR", help="input structure (default POSCAR) / 输入结构")
    g.add_argument("-o", "--outdir", default=".", help="output directory (default .) / 输出目录")
    g.add_argument("--elements", help="element symbols for VASP4-style POSCARs, e.g. 'Ga As'")
    g.add_argument("-f", "--force", action="store_true", help="overwrite existing files / 覆盖已有文件")
    g.add_argument("--stdout", action="store_true", help="print files to screen, write nothing / 只打印不写文件")
    g.add_argument("--lang", choices=("en", "zh", "none"), default="en",
                   help="INCAR comment language: en (ASCII, default), zh (UTF-8 Chinese), none")
    g.add_argument("--system", help="SYSTEM tag (default: POSCAR comment line)")

    g = p.add_argument_group("electronic / 电子结构")
    g.add_argument("--encut", type=float, help="ENCUT in eV (default: factor x max ENMAX)")
    g.add_argument("--encut-factor", type=float, default=1.3, help="ENCUT = factor x max(ENMAX) (default 1.3)")
    g.add_argument("--ediff", help="override EDIFF, e.g. 1E-7")
    g.add_argument("--functional", default="PBE", type=str.upper, choices=list(FUNCTIONALS),
                   help="XC functional for non-HSE tasks (default PBE)")
    g.add_argument("--metal", action="store_true", help="metallic system: MP smearing / 金属体系")
    g.add_argument("--ismear", type=int, help="override ISMEAR")
    g.add_argument("--sigma", type=float, help="override SIGMA (eV)")
    g.add_argument("--lreal", choices=("Auto", ".FALSE.", ".TRUE.", "On"), help="override LREAL")
    g.add_argument("--vdw", choices=list(VDW), help="vdW correction: " + ", ".join(f"{k}={v[0]}" for k, v in VDW.items()))
    g.add_argument("--ispin", type=int, choices=(1, 2), help="force ISPIN (default: 2 if --magmom/dftu)")
    g.add_argument("--magmom", help="initial moments, '2*3.0 4*0' or '2*3_4*0' (SOC: 3 per atom or 1 per atom = z)")
    g.add_argument("--nupdown", type=float, help="NUPDOWN (fixed total moment)")
    g.add_argument("--hubbard", help="DFT+U values, e.g. 'Ni=6.2,Fe=5.3' (Dudarev U_eff)")
    g.add_argument("--saxis", default="0 0 1", help="SAXIS for soc (default '0 0 1')")
    g.add_argument("--soc-from-chgcar", action="store_true",
                   help="soc: ICHARG=11 non-SCF on a collinear CHGCAR (MAE by force theorem)")
    g.add_argument("--nkred", type=int, help="hse: NKRED down-sampling")
    g.add_argument("--bandgap", choices=("COMPACT", "WEIGHT", "KPOINT"),
                   help="BANDGAP tag (VASP >= 6.4.3)")
    g.add_argument("--nbands", type=int, help="NBANDS (optics: default ~3x VASP default)")
    g.add_argument("--write-wavecar", action="store_true", help="scf/hse: keep WAVECAR (LWAVE=.TRUE.)")

    g = p.add_argument_group("ionic / MD / NEB")
    g.add_argument("--isif", type=int, help="relax: ISIF (default 3, or 2 for slabs/molecules)")
    g.add_argument("--relax-2d-cell", action="store_true", help="relax: ISIF=3 + LATTICE_CONSTRAINTS = T T F (2D)")
    g.add_argument("--nsw", type=int, help="NSW")
    g.add_argument("--ediffg", help="EDIFFG, e.g. -0.01")
    g.add_argument("--images", type=int, default=5, help="neb: number of intermediate images (default 5)")
    g.add_argument("--iopt", type=int, default=3, choices=(1, 2, 3, 4, 7), help="neb: VTST IOPT (default 3 = QuickMin)")
    g.add_argument("--climb", action=argparse.BooleanOptionalAction, default=True, help="neb: LCLIMB (default on)")
    g.add_argument("--no-vtst", action="store_true", help="neb: VASP built-in NEB instead of VTST")
    g.add_argument("--temp", type=float, default=300.0, help="md/mlff: TEBEG (K)")
    g.add_argument("--temp-end", type=float, help="md/mlff: TEEND (K), default = --temp")
    g.add_argument("--ensemble", choices=("nvt", "npt", "nve"), default="nvt", help="md/mlff ensemble")
    g.add_argument("--potim", type=float, help="md/mlff: time step in fs")
    g.add_argument("--pstress", type=float, help="npt: PSTRESS in kB")
    g.add_argument("--ml-mode", choices=("train", "select", "refit", "run"), default="train", help="mlff: ML_MODE")
    g.add_argument("--dfpt", action="store_true", help="phonon: DFPT (IBRION=8) instead of finite displacement")
    g.add_argument("--born", action="store_true", help="phonon --dfpt: also LEPSILON for Born charges")
    g.add_argument("--optics-mode", choices=("loptics", "lepsilon", "lcalceps"), default="loptics", help="optics flavour")
    g.add_argument("--ionic", action="store_true", help="optics lepsilon/lcalceps: include ionic contribution")

    g = p.add_argument_group("k-points / k 点")
    g.add_argument("--kspacing", type=float, help="k spacing in 1/A (VASP convention, 2*pi included)")
    g.add_argument("--kmesh", help="explicit mesh 'n1 n2 n3' (overrides --kspacing)")
    g.add_argument("--kmode", choices=("auto", "gamma", "mp"), default="auto", help="mesh centring (default auto)")
    g.add_argument("--kspacing-in-incar", action="store_true", help="write KSPACING/KGAMMA in INCAR, no KPOINTS")
    g.add_argument("--vacuum", type=float, default=8.0, help="gap (A) treated as vacuum (default 8)")
    g.add_argument("--band", action="store_true", help="hse: also produce a band path (KPOINTS_OPT or zero-weight)")
    g.add_argument("--kpoints-opt", action="store_true", help="band: SCF + KPOINTS_OPT (VASP >= 6.3) instead of ICHARG=11")
    g.add_argument("--ibzkpt", help="hse --band: IBZKPT of a PBE run -> zero-weight KPOINTS")
    g.add_argument("--nkpath", type=int, default=40, help="points per path segment (default 40)")
    g.add_argument("--path-backend", choices=("ase", "seekpath"), default="ase", help="k-path generator")
    g.add_argument("--symprec", type=float, default=1e-3, help="symmetry tolerance for lattice detection/k-path")

    g = p.add_argument_group("POTCAR")
    g.add_argument("--potcar", choices=("auto", "skip"), default="auto",
                   help="auto: build from $VASP_PP_PATH if available; skip: never")
    g.add_argument("--pp-path", default=os.environ.get("VASP_PP_PATH"), help="PP root (default $VASP_PP_PATH)")
    g.add_argument("--pp-dir", help="PBE set directory below the root (default: first of potpaw_PBE.64, potpaw_PBE.54, potpaw_PBE)")
    g.add_argument("--pp-flavor", choices=("recommended", "legacy", "gw"), default="recommended",
                   help="recommended (VASP list Oct 2026), legacy (pre-Oct-2026 list), gw (_GW potentials)")
    g.add_argument("--potcar-map", help="override names, e.g. 'Fe=Fe,O=O_h'")

    g = p.add_argument_group("parallel / 并行")
    g.add_argument("--ncore", type=int, default=4, help="NCORE (default 4; 0 = do not write)")
    g.add_argument("--kpar", type=int, help="KPAR")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    if a.task == "hse" and a.ibzkpt:
        a.band = True
    elements = a.elements.split() if a.elements else None
    try:
        s = read_poscar(a.poscar, elements)
    except (FileNotFoundError, ValueError) as exc:
        die(str(exc))

    outdir = Path(a.outdir)
    files = {}          # name -> text
    notes = []

    # ---- POTCAR / ENMAX / NELECT ----
    overrides = parse_map(a.potcar_map)
    names = potcar_names(s.species, a.pp_flavor, overrides)
    pinfo = None
    existing = outdir / "POTCAR"
    if a.potcar == "auto" and a.pp_path:
        pp_dir = find_pp_dir(a.pp_path, a.pp_dir)
        text, pinfo = build_potcar(names, pp_dir)
        files["POTCAR"] = text
        notes.append(f"POTCAR from {pp_dir}: " + " ".join(names))
    elif existing.is_file():
        pinfo = read_existing_potcar(existing)
        notes.append(f"using ENMAX/ZVAL of the existing {existing}")
    if pinfo:
        tit_el = []
        for d in pinfo:
            parts = d["titel"].split()
            tit_el.append(_clean_symbol(parts[1]) if len(parts) > 1 else "?")
        if tit_el != s.species:
            msg = f"POTCAR elements {tit_el} do not match POSCAR elements {s.species}"
            if "POTCAR" in files:
                die(msg)
            warn(msg)
    # ENMAX
    enmax = None
    if pinfo and all(d["enmax"] for d in pinfo):
        enmax = max(d["enmax"] for d in pinfo)
        src = "POTCAR"
    else:
        vals = []
        for e, nm in zip(s.species, names):
            if a.pp_flavor == "gw" and e in GW_PBE and nm == GW_PBE[e][0]:
                vals.append(GW_PBE[e][1])
            elif a.pp_flavor == "legacy" and e in LEGACY_PBE and nm == LEGACY_PBE[e][0]:
                vals.append(LEGACY_PBE[e][2])
            elif e in RECOMMENDED_PBE and nm == RECOMMENDED_PBE[e][0]:
                vals.append(RECOMMENDED_PBE[e][2])
        if len(vals) == len(s.species):
            enmax = max(vals)
            src = "built-in table"
    if a.encut:
        encut = int(round(a.encut))
        enc_note = "user value (eV)"
    elif enmax:
        encut = round_up(a.encut_factor * enmax)
        enc_note = f"{a.encut_factor:g} x max(ENMAX) = {a.encut_factor:g} x {enmax:.1f} eV ({src})"
    else:
        encut = 520
        enc_note = "placeholder: ENMAX unknown, set ENCUT >= 1.3 x max ENMAX of your POTCAR"
        warn("ENMAX unknown (no POTCAR and custom names); ENCUT set to 520 eV - please check")
    # NELECT (for optics NBANDS)
    zvals = None
    if pinfo and all(d["zval"] for d in pinfo):
        zvals = [d["zval"] for d in pinfo]
    elif all(e in RECOMMENDED_PBE for e in s.species) and a.pp_flavor == "recommended" and not overrides:
        zvals = [RECOMMENDED_PBE[e][1] for e in s.species]
    nbands, nb_note = None, ""
    if a.nbands:
        nbands, nb_note = a.nbands, "user value"
    elif a.task == "optics" and zvals:
        nelect = sum(z * n for z, n in zip(zvals, s.counts))
        nb0 = max(0.5 * round(nelect + 2) + max(0.5 * s.natoms, 3), int(0.6 * nelect))
        nbands = int(math.ceil(3 * nb0 / 8.0) * 8)
        nb_note = f"~3x default ({int(nb0)}) for NELECT={nelect:g}; keep divisible by ranks/NCORE"

    # ---- k-points ----
    vac = vacuum_axes(s, a.vacuum)
    molecule = all(vac)
    is2d = (sum(vac) == 1)
    fam = lattice_family(s, a.symprec)
    kspacing = a.kspacing or DEFAULT_KSPACING[a.task]
    if a.kmesh:
        try:
            mesh = [int(x) for x in a.kmesh.replace("x", " ").split()]
            assert len(mesh) == 3 and min(mesh) >= 1
        except (ValueError, AssertionError):
            die("--kmesh needs three positive integers, e.g. '8 8 1'")
    elif molecule:
        mesh = [1, 1, 1]
    else:
        mesh = kmesh(s.lattice, kspacing, vac)
    if a.kmode == "gamma":
        gamma = True
    elif a.kmode == "mp":
        gamma = False
        if fam == "hexagonal" and any(m % 2 == 0 for m in mesh):
            warn("even Monkhorst-Pack mesh in a hexagonal cell breaks the symmetry (IBZKPT error) - prefer Gamma")
    else:
        gamma = (fam == "hexagonal" or is2d or molecule or
                 a.task in ("dos", "optics", "band", "hse", "md", "mlff", "phonon", "soc"))
    nk_total = int(np.prod(mesh))
    system = a.system or (s.comment.split()[0] if s.comment.split() else "".join(s.species))
    system = re.sub(r"[^\x20-\x7e]", "", system) or "".join(s.species)

    info = dict(encut=encut, encut_note=enc_note, molecule=molecule, is2d=is2d,
                nk_total=nk_total, system=system, nbands=nbands, nbands_note=nb_note,
                kspacing=kspacing, gamma=gamma)
    if a.ncore == 0:
        a.ncore = None
    inc = build_incar(a, s, info)

    header = [f"INCAR generated by vasp_inputgen.py {__version__} (task: {a.task}) on {_dt.date.today().isoformat()}",
              "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License",
              f"Structure: {a.poscar}  ({s.natoms} atoms: " + " ".join(f"{e}{n}" for e, n in zip(s.species, s.counts)) + ")",
              "Check every value before production runs."]
    files["INCAR"] = inc.render(a.lang, header)
    if a.lang == "zh":
        warn("--lang zh writes UTF-8 comments; VASP 6 normally ignores them but keep ASCII if your build complains")

    # ---- KPOINTS ----
    kcomment = (f"{'Gamma' if gamma else 'MP'} mesh {mesh[0]}x{mesh[1]}x{mesh[2]} from k-spacing "
                f"{kspacing:g} 1/A ({fam}{', 2D' if is2d else ''}{', molecule' if molecule else ''}) by vasp_inputgen.py")
    if a.kmesh:
        kcomment = f"{'Gamma' if gamma else 'MP'} mesh {mesh[0]}x{mesh[1]}x{mesh[2]} (user) by vasp_inputgen.py"
    need_path = (a.task == "band") or (a.task == "hse" and a.band)
    uniform = kpoints_auto_text(mesh, gamma, kcomment)
    if need_path:
        pbc = [not v for v in vac]
        if sum(pbc) == 2 and not pbc[2]:
            pbc_arg = pbc
        else:
            pbc_arg = True
            if is2d:
                warn("vacuum is not along the 3rd lattice vector; using the 3D k-path")
        segs, desc, prim = kpath_segments(s, a.path_backend, pbc_arg, a.symprec)
        if prim is not None:
            files["POSCAR_seekpath"] = None  # placeholder, written below
            notes.append("seekpath standardised the cell: run SCF AND band with POSCAR_seekpath "
                         "(the k-path refers to that cell)")
        line_text = kpoints_line_text(segs, a.nkpath, f"k-path {desc} by vasp_inputgen.py")
        if a.task == "band" and not a.kpoints_opt:
            files["KPOINTS"] = line_text
            notes.append("band (ICHARG=11): copy CHGCAR of the SCF run (same ENCUT/ISPIN/LMAXMIX/POTCAR) here")
        elif a.task == "hse" and a.ibzkpt:
            kk, ww = read_ibzkpt(a.ibzkpt)
            files["KPOINTS"] = kpoints_zero_weight_text(kk, ww, segs, max(8, a.nkpath // 3))
            notes.append(f"zero-weight KPOINTS built from {a.ibzkpt}; plot only the zero-weight points")
        else:
            files["KPOINTS"] = uniform
            files["KPOINTS_OPT"] = line_text
            notes.append("KPOINTS_OPT needs VASP >= 6.3 (ISPIN=2 PDOS bug fixed in 6.6.1); results in vaspout.h5 / PROCAR_OPT"
                         + (" (screened HSE needs no HFRCUT; unscreened hybrids: HFRCUT = -1)" if a.task == "hse" else ""))
    elif not a.kspacing_in_incar:
        files["KPOINTS"] = uniform
    if molecule and a.task not in ("band",):
        notes.append("isolated molecule detected: Gamma-only, consider vasp_gam")
    if a.task == "neb":
        notes.append(f"NEB: create 00..{a.images + 1:02d} with 'python neb_tools.py interp POSCAR_IS POSCAR_FS -n {a.images}'; "
                     f"MPI ranks must be a multiple of IMAGES={a.images}")
    if a.task == "soc":
        notes.append("SOC: run with vasp_ncl; NBANDS is doubled automatically")
    if a.task == "phonon" and not a.dfpt:
        notes.append("phonon: use this INCAR for every phonopy displacement (phonopy -d --dim ...)")

    # ---- output ----
    if a.stdout:
        for nm, txt in files.items():
            if nm == "POTCAR":
                print(f"===== POTCAR: {' '.join(names)} (content not shown) =====")
                continue
            if txt is None:
                continue
            print(f"===== {nm} =====")
            print(txt, end="")
    else:
        outdir.mkdir(parents=True, exist_ok=True)
        clash = [nm for nm in files if (outdir / nm).exists()]
        if clash and not a.force:
            die(f"{', '.join(clash)} already exist in {outdir} (use -f to overwrite)")
        for nm, txt in files.items():
            if nm == "POSCAR_seekpath":
                write_poscar(outdir / nm, prim)
                continue
            with open(outdir / nm, "w", newline="\n", encoding="utf-8") as fh:
                fh.write(txt)
        print(f"wrote {', '.join(files)} to {outdir}")
    if "POTCAR" not in files and not (outdir / "POTCAR").is_file():
        notes.append("POTCAR not written (set VASP_PP_PATH or --pp-path). Suggested: "
                     + " ".join(names) + "  e.g. cat $VASP_PP_PATH/potpaw_PBE/{"
                     + ",".join(names) + "}/POTCAR > POTCAR")
    print(f"structure: {s.natoms} atoms, {fam}{' (2D)' if is2d else ''}{' (molecule)' if molecule else ''}; "
          f"k-mesh {mesh} ({'Gamma' if gamma else 'MP'}); ENCUT = {encut} eV")
    for n in notes:
        print("note:", n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
