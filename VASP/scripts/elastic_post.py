#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
elastic_post.py -- post-process VASP elastic constants (IBRION = 6, ISIF >= 3).
                   解析 OUTCAR 弹性张量：kBar->GPa、Voigt 重排、Born 稳定性、VRH 模量、德拜温度。

Purpose
-------
* Parse the "ELASTIC MODULI ... (kBar)" blocks of OUTCAR (TOTAL, clamped-ion /
  symmetrized, ionic-relaxation contribution), convert kBar -> GPa (divide by 10)
  and reorder VASP's XX YY ZZ XY YZ ZX into Voigt order 11 22 33 23 13 12
  (i.e. swap rows/columns 4-6: VASP XY,YZ,ZX -> Voigt YZ,ZX,XY).
* Alternatively read a 6x6 matrix from a text file (--matrix).
* Born mechanical stability: eigenvalue test (always) plus the explicit
  crystal-system conditions of Mouhat & Coudert, PRB 90, 224104 (2014).
  The crystal system is taken from --system, from spglib (if installed) applied
  to POSCAR/CONTCAR, or guessed from the pattern of the matrix.
* Voigt / Reuss / Hill bulk and shear moduli, Young's modulus, Poisson ratio,
  Pugh ratio B/G, universal anisotropy A^U, Zener ratio and Cauchy pressure
  (cubic), Chen's empirical Vickers hardness.
* With POSCAR/CONTCAR (masses from ASE): density, longitudinal / transverse /
  mean sound velocities and the Debye temperature (Anderson's method).
* 2D materials (--2d): in-plane stiffness in N/m (C x c), 2D Young's moduli,
  Poisson ratios and 2D stability conditions.

Usage examples
--------------
  python elastic_post.py                                # ./OUTCAR (+ ./CONTCAR or ./POSCAR)
  python elastic_post.py -f OUTCAR --poscar POSCAR --system cubic
  python elastic_post.py --block symmetrized            # clamped-ion (rigid-ion) moduli
  python elastic_post.py --matrix Cij_GPa.txt --system hexagonal --density 1.74
  python elastic_post.py --2d --poscar CONTCAR          # monolayer, N/m
  python elastic_post.py --json elastic.json

Dependencies
------------
  Python >= 3.9, numpy; optional: ase (atomic masses), spglib (crystal system).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

VASP_ORDER = ["XX", "YY", "ZZ", "XY", "YZ", "ZX"]
TO_VOIGT = [0, 1, 2, 4, 5, 3]  # Voigt index (11,22,33,23,13,12) -> VASP row
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
H_PLANCK = 6.62607015e-34
K_B = 1.380649e-23
AMU = 1.66053906660e-27


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def check_file(path: str, what: str) -> None:
    if not os.path.isfile(path):
        die(f"{what} not found: {path}")
    if os.path.getsize(path) == 0:
        die(f"{what} is empty: {path}")


# --------------------------------------------------------------------------- input
def read_outcar_moduli(path: str) -> Dict[str, np.ndarray]:
    """Return {kind: 6x6 array in kBar, VASP order} for every ELASTIC MODULI block (last wins)."""
    check_file(path, "OUTCAR")
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    blocks: Dict[str, np.ndarray] = {}
    for i, line in enumerate(lines):
        if "ELASTIC MODULI" not in line.upper() or "KBAR" not in line.upper():
            continue
        up = line.upper()
        if "TOTAL" in up:
            kind = "total"
        elif "IONIC" in up:
            kind = "ionic"
        elif "SYMMETRIZED" in up:
            kind = "symmetrized"
        else:
            kind = "clamped"
        rows: Dict[str, List[float]] = {}
        for j in range(i + 1, min(i + 14, len(lines))):
            t = lines[j].split()
            if t and t[0] in VASP_ORDER:
                vals = [float(x) for x in FLOAT_RE.findall(lines[j][lines[j].index(t[0]) + 2:])]
                if len(vals) >= 6:
                    rows[t[0]] = vals[:6]
            if len(rows) == 6:
                break
        if len(rows) == 6:
            blocks[kind] = np.array([rows[k] for k in VASP_ORDER])
        else:
            warn(f"incomplete '{line.strip()}' block at line {i + 1} of {path}")
    return blocks


def vasp_to_voigt(c: np.ndarray) -> np.ndarray:
    return c[np.ix_(TO_VOIGT, TO_VOIGT)]


def read_matrix(path: str) -> np.ndarray:
    check_file(path, "matrix file")
    rows = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            body = re.split(r"[#!]", line)[0]
            vals = [float(x) for x in FLOAT_RE.findall(body)]
            if len(vals) >= 6:
                rows.append(vals[-6:])
    if len(rows) < 6:
        die(f"{path}: need 6 rows with 6 numbers each (found {len(rows)})")
    return np.array(rows[:6])


def read_poscar(path: str) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    check_file(path, "POSCAR/CONTCAR")
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    try:
        scale = [float(x) for x in lines[1].split()[:3]]
        lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
        if len(scale) == 3:
            fac = np.array(scale)
        elif scale[0] < 0:
            fac = np.full(3, (abs(scale[0]) / abs(np.linalg.det(lat))) ** (1 / 3))
        else:
            fac = np.full(3, scale[0])
        lat = lat * fac[None, :]
        toks = lines[5].split()
        if toks and not re.match(r"^\d+$", toks[0]):
            names = [re.match(r"[A-Z][a-z]?", re.match(r"[A-Za-z]+", t).group(0).capitalize()).group(0)
                     for t in toks]
            counts = [int(x) for x in lines[6].split()[:len(names)]]
            k = 7
        else:
            die(f"{path}: VASP 4 format without element names is not supported here")
        if lines[k].strip()[:1] in ("s", "S"):
            k += 1
        cart = lines[k].strip()[:1] in ("c", "C", "k", "K")
        nat = sum(counts)
        pos = np.array([[float(x) for x in lines[k + 1 + i].split()[:3]] for i in range(nat)])
        if cart:
            pos = (pos * fac[None, :]) @ np.linalg.inv(lat)
    except (IndexError, ValueError, AttributeError) as exc:
        die(f"cannot parse {path}: {exc}")
    symbols: List[str] = []
    for n, c in zip(names, counts):
        symbols += [n] * c
    return lat, pos, symbols


# --------------------------------------------------------------------------- symmetry
def system_from_spglib(lat, pos, symbols) -> Optional[Tuple[str, int]]:
    try:
        import spglib  # optional
    except ImportError:
        return None
    uniq = list(dict.fromkeys(symbols))
    nums = [uniq.index(s) + 1 for s in symbols]
    try:
        ds = spglib.get_symmetry_dataset((lat, pos, nums), symprec=1e-3)
    except Exception:  # noqa: BLE001 - spglib raises different errors across versions
        return None
    if ds is None:
        return None
    n = ds.number if hasattr(ds, "number") else ds["number"]
    if n <= 2:
        s = "triclinic"
    elif n <= 15:
        s = "monoclinic"
    elif n <= 74:
        s = "orthorhombic"
    elif n <= 88:
        s = "tetragonal2"
    elif n <= 142:
        s = "tetragonal"
    elif n <= 148:
        s = "trigonal2"
    elif n <= 167:
        s = "trigonal"
    elif n <= 194:
        s = "hexagonal"
    else:
        s = "cubic"
    return s, n


def guess_system(C: np.ndarray) -> str:
    tol = 0.02 * float(np.abs(C).max()) + 1e-6
    eq = lambda a, b: abs(a - b) <= tol  # noqa: E731
    z = lambda a: abs(a) <= tol  # noqa: E731
    c = C
    normal_shear = all(z(c[i, j]) for i in range(3) for j in range(3, 6))
    shear_off = all(z(c[i, j]) for i in range(3, 6) for j in range(3, 6) if i != j)
    if normal_shear and shear_off:
        if (eq(c[0, 0], c[1, 1]) and eq(c[1, 1], c[2, 2]) and eq(c[0, 1], c[0, 2]) and eq(c[0, 2], c[1, 2])
                and eq(c[3, 3], c[4, 4]) and eq(c[4, 4], c[5, 5])):
            return "cubic"
        if eq(c[0, 0], c[1, 1]) and eq(c[0, 2], c[1, 2]) and eq(c[3, 3], c[4, 4]):
            return "hexagonal" if eq(c[5, 5], (c[0, 0] - c[0, 1]) / 2) else "tetragonal"
        return "orthorhombic"
    if eq(c[0, 0], c[1, 1]) and eq(c[3, 3], c[4, 4]) and not z(c[0, 3]) and eq(c[0, 3], -c[1, 3]):
        return "trigonal"
    if eq(c[0, 0], c[1, 1]) and eq(c[3, 3], c[4, 4]) and not z(c[0, 5]) and eq(c[0, 5], -c[1, 5]):
        return "tetragonal2"
    nz = sum(1 for i in range(6) for j in range(i, 6) if not z(c[i, j]))
    return "monoclinic" if nz <= 13 else "triclinic"


# --------------------------------------------------------------------------- stability
def born_conditions(C: np.ndarray, system: str) -> List[Tuple[str, bool]]:
    c = lambda i, j: C[i - 1, j - 1]  # noqa: E731  (1-based Voigt indices)
    conds: List[Tuple[str, bool]] = []
    if system == "cubic":
        conds = [("C11 - C12 > 0", c(1, 1) - c(1, 2) > 0),
                 ("C11 + 2 C12 > 0", c(1, 1) + 2 * c(1, 2) > 0),
                 ("C44 > 0", c(4, 4) > 0)]
    elif system == "hexagonal":
        conds = [("C11 > |C12|", c(1, 1) > abs(c(1, 2))),
                 ("2 C13^2 < C33 (C11 + C12)", 2 * c(1, 3) ** 2 < c(3, 3) * (c(1, 1) + c(1, 2))),
                 ("C44 > 0", c(4, 4) > 0),
                 ("C66 = (C11 - C12)/2 > 0", (c(1, 1) - c(1, 2)) / 2 > 0)]
    elif system in ("tetragonal", "tetragonal2"):
        conds = [("C11 > |C12|", c(1, 1) > abs(c(1, 2))),
                 ("2 C13^2 < C33 (C11 + C12)", 2 * c(1, 3) ** 2 < c(3, 3) * (c(1, 1) + c(1, 2))),
                 ("C44 > 0", c(4, 4) > 0),
                 ("C66 > 0", c(6, 6) > 0)]
        if system == "tetragonal2":
            conds.append(("2 C16^2 < C66 (C11 - C12)", 2 * c(1, 6) ** 2 < c(6, 6) * (c(1, 1) - c(1, 2))))
    elif system in ("trigonal", "trigonal2"):
        conds = [("C11 > |C12|", c(1, 1) > abs(c(1, 2))),
                 ("C44 > 0", c(4, 4) > 0),
                 ("C13^2 < C33 (C11 + C12)/2", c(1, 3) ** 2 < 0.5 * c(3, 3) * (c(1, 1) + c(1, 2)))]
        if system == "trigonal":
            conds.append(("C14^2 < C44 (C11 - C12)/2", c(1, 4) ** 2 < 0.5 * c(4, 4) * (c(1, 1) - c(1, 2))))
        else:
            conds.append(("C14^2 + C15^2 < C44 (C11 - C12)/2",
                          c(1, 4) ** 2 + c(1, 5) ** 2 < 0.5 * c(4, 4) * (c(1, 1) - c(1, 2))))
    elif system == "orthorhombic":
        det3 = (c(1, 1) * c(2, 2) * c(3, 3) + 2 * c(1, 2) * c(1, 3) * c(2, 3)
                - c(1, 1) * c(2, 3) ** 2 - c(2, 2) * c(1, 3) ** 2 - c(3, 3) * c(1, 2) ** 2)
        conds = [("C11 > 0", c(1, 1) > 0),
                 ("C11 C22 > C12^2", c(1, 1) * c(2, 2) > c(1, 2) ** 2),
                 ("det(C_ij, i,j<=3) > 0", det3 > 0),
                 ("C44 > 0", c(4, 4) > 0), ("C55 > 0", c(5, 5) > 0), ("C66 > 0", c(6, 6) > 0)]
    return conds


# --------------------------------------------------------------------------- moduli
def vrh(C: np.ndarray) -> Dict[str, float]:
    S = np.linalg.inv(C)
    kv = (C[0, 0] + C[1, 1] + C[2, 2] + 2 * (C[0, 1] + C[0, 2] + C[1, 2])) / 9.0
    gv = (C[0, 0] + C[1, 1] + C[2, 2] - (C[0, 1] + C[0, 2] + C[1, 2]) + 3 * (C[3, 3] + C[4, 4] + C[5, 5])) / 15.0
    kr = 1.0 / (S[0, 0] + S[1, 1] + S[2, 2] + 2 * (S[0, 1] + S[0, 2] + S[1, 2]))
    gr = 15.0 / (4 * (S[0, 0] + S[1, 1] + S[2, 2]) - 4 * (S[0, 1] + S[0, 2] + S[1, 2])
                 + 3 * (S[3, 3] + S[4, 4] + S[5, 5]))
    out = {"K_V": kv, "K_R": kr, "K_H": (kv + kr) / 2, "G_V": gv, "G_R": gr, "G_H": (gv + gr) / 2}
    for s in ("V", "R", "H"):
        k, g = out[f"K_{s}"], out[f"G_{s}"]
        out[f"E_{s}"] = 9 * k * g / (3 * k + g)
        out[f"nu_{s}"] = (3 * k - 2 * g) / (2 * (3 * k + g))
    out["A_U"] = 5 * gv / gr + kv / kr - 6
    out["Pugh_B/G"] = out["K_H"] / out["G_H"]
    if out["K_H"] > 0 and out["G_H"] > 0:
        kk = out["G_H"] / out["K_H"]
        out["Hv_Chen"] = 2 * (kk * kk * out["G_H"]) ** 0.585 - 3
    else:
        out["Hv_Chen"] = float("nan")
    out["S"] = S
    return out


def masses_for(symbols: List[str]) -> Optional[np.ndarray]:
    try:
        from ase.data import atomic_masses, atomic_numbers  # optional
    except ImportError:
        return None
    try:
        return np.array([atomic_masses[atomic_numbers[s]] for s in symbols])
    except KeyError as exc:
        warn(f"unknown element {exc} for the mass table")
        return None


def debye(K: float, G: float, rho_kg_m3: float, n_atoms: int, vol_A3: float) -> Dict[str, float]:
    vl = math.sqrt((K + 4 * G / 3) * 1e9 / rho_kg_m3)
    vt = math.sqrt(G * 1e9 / rho_kg_m3)
    vm = (((2 / vt ** 3) + (1 / vl ** 3)) / 3) ** (-1 / 3)
    theta = H_PLANCK / K_B * (3 * n_atoms / (4 * math.pi * vol_A3 * 1e-30)) ** (1 / 3) * vm
    return {"v_l": vl, "v_t": vt, "v_m": vm, "Theta_D": theta}


# --------------------------------------------------------------------------- output helpers
def print_matrix(C: np.ndarray, title: str, unit: str) -> None:
    print(title)
    print("        " + "".join(f"{j + 1:>10d}" for j in range(6)))
    for i in range(6):
        print(f"  {i + 1:>4d}  " + "".join(f"{C[i, j]:10.2f}" for j in range(6)))
    print(f"  (unit: {unit}; Voigt order 1=xx 2=yy 3=zz 4=yz 5=zx 6=xy)")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Elastic constants from VASP OUTCAR (IBRION = 6, ISIF >= 3): unit conversion, Voigt "
                    "reordering, Born stability, VRH averages, sound velocities and Debye temperature.\n"
                    "解析 OUTCAR 弹性张量（kBar -> GPa，重排为 Voigt 顺序），判断 Born 稳定性，计算 VRH 模量等。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples / 示例:\n"
               "  python elastic_post.py -f OUTCAR --poscar CONTCAR\n"
               "  python elastic_post.py --matrix Cij.txt --system hexagonal --density 1.74\n"
               "  python elastic_post.py --2d --poscar CONTCAR\n"
               "Note: OUTCAR values are kBar in the order XX YY ZZ XY YZ ZX; 1 GPa = 10 kBar.")
    p.add_argument("-f", "--file", default="OUTCAR", help="OUTCAR (default ./OUTCAR)")
    p.add_argument("--block", choices=["total", "symmetrized", "clamped", "ionic"], default="total",
                   help="which OUTCAR block to use (default TOTAL ELASTIC MODULI, i.e. incl. ionic relaxation)")
    p.add_argument("--matrix", default=None,
                   help="read a 6x6 matrix (Voigt order, GPa) from a text file instead of OUTCAR")
    p.add_argument("--matrix-units", choices=["GPa", "kBar"], default="GPa", help="units of --matrix")
    p.add_argument("--matrix-vasp-order", action="store_true",
                   help="--matrix file is in VASP order XX YY ZZ XY YZ ZX (will be reordered)")
    p.add_argument("--poscar", default=None,
                   help="structure for density/Debye/crystal system (default CONTCAR/POSCAR next to OUTCAR)")
    p.add_argument("--system", default="auto",
                   choices=["auto", "cubic", "hexagonal", "tetragonal", "tetragonal2", "trigonal", "trigonal2",
                            "orthorhombic", "monoclinic", "triclinic"],
                   help="crystal system for the Born criteria (tetragonal2 = classes 4, -4, 4/m; "
                        "trigonal2 = classes 3, -3); default auto")
    p.add_argument("--density", type=float, default=None, help="density in g/cm^3 (overrides POSCAR + ASE masses)")
    p.add_argument("--2d", dest="two_d", action="store_true",
                   help="2D material in the xy plane: report in-plane constants in N/m")
    p.add_argument("--c-length", type=float, default=None,
                   help="cell length along z in Angstrom for --2d (default from POSCAR)")
    p.add_argument("--json", default=None, help="write all results to a JSON file")
    p.add_argument("--dat", default=None, help="write the Voigt matrix (GPa) to a text file")
    return p


def main(argv=None) -> None:
    try:
        sys.stdout.reconfigure(errors="backslashreplace")  # type: ignore[attr-defined]
    except (AttributeError, ValueError):
        pass
    args = build_parser().parse_args(argv)
    results: Dict = {}
    folder = os.path.dirname(os.path.abspath(args.file))

    if args.matrix:
        C = read_matrix(args.matrix)
        if args.matrix_vasp_order:
            C = vasp_to_voigt(C)
        if args.matrix_units == "kBar":
            C = C / 10.0
        source = args.matrix
        folder = os.path.dirname(os.path.abspath(args.matrix))
    else:
        blocks = read_outcar_moduli(args.file)
        if not blocks:
            die(f"no 'ELASTIC MODULI (kBar)' block in {args.file}: was it run with IBRION = 6 and ISIF >= 3 "
                "and did it finish?")
        print(f"OUTCAR blocks found: {', '.join(sorted(blocks))}")
        kind = args.block
        if kind not in blocks:
            alt = next((k for k in ("total", "symmetrized", "clamped") if k in blocks), None)
            if alt is None:
                die(f"block '{kind}' not found in {args.file}")
            warn(f"block '{kind}' not in OUTCAR, using '{alt}'")
            kind = alt
        if kind != "total":
            print(f"NOTE: using the '{kind}' block - ionic-relaxation contributions are "
                  f"{'the only part' if kind == 'ionic' else 'NOT included'}.")
        C_vasp = blocks[kind]
        print_matrix(C_vasp, f"Raw '{kind}' block from OUTCAR (VASP order XX YY ZZ XY YZ ZX):", "kBar")
        C = vasp_to_voigt(C_vasp) / 10.0
        source = f"{args.file} [{kind}]"
        if "ionic" in blocks and kind == "total":
            ion = vasp_to_voigt(blocks["ionic"]) / 10.0
            print(f"Largest ionic-relaxation contribution: {np.abs(ion).max():.2f} GPa")
    asym = float(np.abs(C - C.T).max())
    if asym > 1e-3 * max(1.0, float(np.abs(C).max())):
        warn(f"matrix is not symmetric (max |Cij - Cji| = {asym:.3f} GPa); it is symmetrised")
    C = 0.5 * (C + C.T)
    print()
    print_matrix(C, f"Elastic constants C_ij from {source}:", "GPa")
    results["C_GPa"] = C.tolist()
    results["source"] = source

    # structure
    pos_path = args.poscar
    if pos_path is None:
        for cand in ("CONTCAR", "POSCAR"):
            c = os.path.join(folder, cand)
            if os.path.isfile(c) and os.path.getsize(c) > 0:
                pos_path = c
                break
    struct = read_poscar(pos_path) if pos_path else None

    # crystal system
    system = args.system
    sys_src = "command line"
    if system == "auto":
        sg = system_from_spglib(*struct) if struct is not None else None
        if sg:
            system, num = sg
            sys_src = f"spglib on {pos_path}, space group #{num}"
        else:
            system = guess_system(C)
            sys_src = "pattern of the C_ij matrix (install spglib and give POSCAR for a reliable answer)"
    print(f"\nCrystal system: {system} ({sys_src})")
    results["system"] = system

    # stability
    eig = np.linalg.eigvalsh(C)
    print("\nBorn stability")
    print(f"  eigenvalues of C (GPa): {' '.join(f'{e:.2f}' for e in eig)}")
    stable = bool(eig.min() > 0)
    print(f"  [{'PASS' if stable else 'FAIL'}] all eigenvalues positive (necessary and sufficient)")
    for desc, ok in born_conditions(C, system):
        print(f"  [{'PASS' if ok else 'FAIL'}] {desc}")
        stable &= ok
    print(f"  => mechanically {'STABLE' if stable else 'UNSTABLE'}")
    results["eigenvalues_GPa"] = eig.tolist()
    results["stable"] = stable

    if args.two_d:
        clen = args.c_length or (float(np.linalg.norm(struct[0][2])) if struct is not None else None)
        if clen is None:
            die("--2d needs the cell length along z (--c-length or POSCAR)")
        f = clen / 10.0  # GPa * Angstrom -> N/m
        c11, c22, c12, c66 = C[0, 0] * f, C[1, 1] * f, C[0, 1] * f, C[5, 5] * f
        c16, c26 = C[0, 5] * f, C[1, 5] * f
        print(f"\n2D in-plane stiffness (c = {clen:.3f} A; C_2D = C_3D x c):")
        print(f"  C11 = {c11:.2f}  C22 = {c22:.2f}  C12 = {c12:.2f}  C66 = {c66:.2f}  C16 = {c16:.2f}  C26 = {c26:.2f} N/m")
        m2 = np.array([[c11, c12, c16], [c12, c22, c26], [c16, c26, c66]])
        ok2 = bool(np.linalg.eigvalsh(m2).min() > 0)
        print(f"  [{'PASS' if c11 > 0 else 'FAIL'}] C11 > 0   [{'PASS' if c66 > 0 else 'FAIL'}] C66 > 0   "
              f"[{'PASS' if c11 * c22 > c12 ** 2 else 'FAIL'}] C11 C22 > C12^2   "
              f"[{'PASS' if ok2 else 'FAIL'}] 3x3 eigenvalues > 0")
        det = c11 * c22 - c12 ** 2
        if c11 > 0 and c22 > 0:
            print(f"  Y_x = {det / c22:.2f} N/m, Y_y = {det / c11:.2f} N/m, nu_xy = {c12 / c22:.3f}, "
                  f"nu_yx = {c12 / c11:.3f}")
        results["2D_N_per_m"] = {"C11": c11, "C22": c22, "C12": c12, "C66": c66,
                                 "Y_x": det / c22 if c22 else None, "Y_y": det / c11 if c11 else None}
        print("  (3D averages below are meaningless for a slab with vacuum.)")

    # VRH
    try:
        m = vrh(C)
    except np.linalg.LinAlgError:
        die("C is singular; cannot compute compliances")
    print("\nPolycrystalline averages (GPa)")
    print(f"  {'':10s}{'Voigt':>10s}{'Reuss':>10s}{'Hill':>10s}")
    for key, lab in (("K", "Bulk B"), ("G", "Shear G"), ("E", "Young E")):
        print(f"  {lab:10s}{m[key + '_V']:10.2f}{m[key + '_R']:10.2f}{m[key + '_H']:10.2f}")
    print(f"  {'Poisson nu':10s}{m['nu_V']:10.3f}{m['nu_R']:10.3f}{m['nu_H']:10.3f}")
    print(f"  Pugh ratio B/G = {m['Pugh_B/G']:.3f} ({'ductile' if m['Pugh_B/G'] > 1.75 else 'brittle'}"
          " by the 1.75 rule of thumb)")
    print(f"  Universal anisotropy A^U = {m['A_U']:.4f} (0 = isotropic)")
    if system == "cubic":
        zener = 2 * C[3, 3] / (C[0, 0] - C[0, 1]) if C[0, 0] != C[0, 1] else float("inf")
        print(f"  Zener ratio A = 2 C44/(C11 - C12) = {zener:.4f};  Cauchy pressure C12 - C44 = "
              f"{C[0, 1] - C[3, 3]:.2f} GPa")
        results["Zener"] = zener
    if np.isfinite(m["Hv_Chen"]) and m["Hv_Chen"] > 0:
        print(f"  Vickers hardness (Chen et al. 2011 empirical model) ~ {m['Hv_Chen']:.1f} GPa")
    results["VRH"] = {k: v for k, v in m.items() if k != "S"}
    results["S_1_per_GPa"] = m["S"].tolist()

    # density / Debye
    rho = None
    if args.density:
        rho = args.density
        print(f"\nDensity (command line): {rho:.4f} g/cm^3")
    elif struct is not None:
        mass = masses_for(struct[2])
        vol = abs(np.linalg.det(struct[0]))
        if mass is None:
            warn("ase not available: no atomic masses -> give --density to get sound velocities")
        else:
            rho = mass.sum() * AMU / (vol * 1e-30) / 1000.0
            print(f"\nDensity from {pos_path}: {rho:.4f} g/cm^3 ({len(struct[2])} atoms, V = {vol:.3f} A^3)")
    if rho and not args.two_d:
        if struct is not None:
            nat, vol = len(struct[2]), abs(np.linalg.det(struct[0]))
        else:
            nat = vol = None
        if m["K_H"] > 0 and m["G_H"] > 0:
            v_l = math.sqrt((m["K_H"] + 4 * m["G_H"] / 3) * 1e9 / (rho * 1000))
            v_t = math.sqrt(m["G_H"] * 1e9 / (rho * 1000))
            v_m = (((2 / v_t ** 3) + (1 / v_l ** 3)) / 3) ** (-1 / 3)
            print(f"  Sound velocities (Hill): v_l = {v_l:.0f} m/s, v_t = {v_t:.0f} m/s, v_m = {v_m:.0f} m/s")
            results["sound_m_per_s"] = {"v_l": v_l, "v_t": v_t, "v_m": v_m}
            if nat:
                d = debye(m["K_H"], m["G_H"], rho * 1000, nat, vol)
                print(f"  Debye temperature (Anderson) Theta_D = {d['Theta_D']:.1f} K")
                results["Debye_K"] = d["Theta_D"]
            else:
                print("  (Debye temperature needs the structure: give --poscar)")
        else:
            warn("non-positive Hill moduli: sound velocities skipped")
        results["density_g_cm3"] = rho

    if args.dat:
        np.savetxt(args.dat, C, fmt="%12.4f", header="Elastic constants C_ij (GPa), Voigt order 11 22 33 23 13 12")
        print(f"\nMatrix written: {args.dat}")
    if args.json:
        def clean(o):
            if isinstance(o, dict):
                return {k: clean(v) for k, v in o.items()}
            if isinstance(o, (list, tuple)):
                return [clean(v) for v in o]
            if isinstance(o, (float, np.floating)):
                return float(o) if np.isfinite(o) else None
            if isinstance(o, np.bool_):
                return bool(o)
            return o
        with open(args.json, "w") as fh:
            json.dump(clean(results), fh, indent=2)
        print(f"Results written: {args.json}")


if __name__ == "__main__":
    main()
