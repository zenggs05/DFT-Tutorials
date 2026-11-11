#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_dos.py -- plot total / projected DOS from a VASP DOSCAR and compute d-band centres.
              读取 DOSCAR 画总态密度与投影态密度，并计算 d 带中心与带宽。

Purpose
-------
* Read DOSCAR written with LORBIT = 10 (s, p, d[, f]) or LORBIT = 11/12 (lm-resolved),
  for ISPIN = 1, ISPIN = 2 and non-collinear runs (only the total part of the
  (tot, mx, my, mz) PDOS is used for non-collinear runs).
* Select projections by element, atom index and orbital:
  ``-p "Ni:d"``, ``-p "O:p"``, ``-p "1-4,9:eg"``, ``-p "Pt:dz2"``; or ``--by-element``.
  Orbital groups: s, p, d, f, t2g (dxy+dyz+dxz), eg (dz2+x2-y2), all, or single lm names.
* Spin-polarised DOS is drawn mirrored (spin down negative).
* d-band centre (first moment), RMS width (second moment) and filling for any
  selection (``--dband "Pt"`` or ``--dband "1-9:d"``), with an optional energy window.
* Optional Gaussian broadening (``--smear``) and column-data export (``--dat``).

Usage examples
--------------
  python plot_dos.py                                   # ./DOSCAR, PDOS per element
  python plot_dos.py -f DOSCAR -p "Ni:d" -p "O:p" --xlim -8 4 -o NiO_pdos.png
  python plot_dos.py -p "Ni:t2g" -p "Ni:eg" --smear 0.05
  python plot_dos.py --dband "Pt" --dband-window -10 2 --no-plot
  python plot_dos.py --by-element --dat dos.dat

Dependencies
------------
  Python >= 3.9, numpy, matplotlib.

Notes
-----
  E-fermi is read from line 6 of DOSCAR (override with --efermi).  The element names
  come from POSCAR/CONTCAR next to DOSCAR (or --poscar).  DOSCAR does not tell
  whether a run was non-collinear; this is detected from OUTCAR/INCAR next to it,
  or forced with --ncl.

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK = "#1f1f1f"
MUTED = "#8a8986"
TDOS_FILL = "#dcdbd6"
TDOS_LINE = "#6f6e6a"
LM = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2",
      "fy3x2", "fxyz", "fyz2", "fz3", "fxz2", "fzx2", "fx3"]
ORB_SETS = {1: ["s"], 3: ["s", "p", "d"], 4: ["s", "p", "d", "f"], 9: LM[:9], 16: LM}
FIX_EXP = re.compile(r"(\d)([+-]\d{2,3})(?=\s|$)")


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


def trapz(y, x):
    f = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return f(y, x)


def parse_block(lines: List[str], what: str) -> np.ndarray:
    """Convert a list of numeric text lines to a 2D float array."""
    text = " ".join(lines)
    try:
        flat = np.array(text.split(), dtype=float)
    except ValueError:
        # Fortran drops the 'E' for 3-digit exponents: 0.1234-100
        flat = np.array(FIX_EXP.sub(r"\1E\2", text).split(), dtype=float)
    ncol = len(lines[0].split())
    if flat.size != ncol * len(lines):
        die(f"DOSCAR {what}: rows have different numbers of columns (file corrupted?)")
    return flat.reshape(len(lines), ncol)


# --------------------------------------------------------------------------- readers
def read_poscar_symbols(path: str) -> List[str]:
    check_file(path, "POSCAR/CONTCAR")
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    try:
        toks = lines[5].split()
        if toks and not re.match(r"^\d+$", toks[0]):
            names = [re.match(r"[A-Z][a-z]?", re.match(r"[A-Za-z]+", t).group(0).capitalize()).group(0)
                     for t in toks]
            counts = [int(t) for t in lines[6].split()[:len(names)]]
        else:
            counts = [int(t) for t in toks]
            names = [f"X{i + 1}" for i in range(len(counts))]
            warn(f"{path}: VASP 4 format without element names; species are called X1, X2, ...")
    except (IndexError, ValueError, AttributeError) as exc:
        die(f"cannot read element names from {path}: {exc}")
    out: List[str] = []
    for n, c in zip(names, counts):
        out += [n] * c
    return out


def detect_ncl(folder: str) -> bool:
    """True if OUTCAR/INCAR next to DOSCAR indicate a non-collinear run."""
    out = os.path.join(folder, "OUTCAR")
    if os.path.isfile(out):
        with open(out, "r", errors="replace") as fh:
            for i, line in enumerate(fh):
                if "LNONCOLLINEAR" in line and "=" in line:
                    return line.split("=")[1].split()[0].upper().startswith("T")
                if i > 5000:
                    break
    inc = os.path.join(folder, "INCAR")
    if os.path.isfile(inc):
        with open(inc, "r", errors="replace") as fh:
            for line in fh:
                body = re.split(r"[#!]", line)[0]
                m = re.search(r"\b(LSORBIT|LNONCOLLINEAR)\s*=\s*\.?([TtFf])", body)
                if m and m.group(2).upper() == "T":
                    return True
    return False


def read_doscar(path: str, ncl_flag: Optional[bool]) -> Dict:
    check_file(path, "DOSCAR")
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    try:
        natoms = int(lines[0].split()[0])
        h = [float(x) for x in FIX_EXP.sub(r"\1E\2", lines[5]).split()[:4]]
        emax, emin, nedos, efermi = h[0], h[1], int(round(h[2])), h[3]
    except (IndexError, ValueError):
        die(f"{path}: header is not in DOSCAR format")
    if len(lines) < 6 + nedos:
        die(f"{path}: truncated (expected {nedos} TDOS lines)")
    tdos = parse_block(lines[6:6 + nedos], "total DOS")
    energies = tdos[:, 0]
    if tdos.shape[1] == 5:
        ispin = 2
        tot = tdos[:, [1, 2]].T
        itot = tdos[:, [3, 4]].T
    elif tdos.shape[1] == 3:
        ispin = 1
        tot = tdos[:, [1]].T
        itot = tdos[:, [2]].T
    else:
        die(f"{path}: unexpected number of TDOS columns ({tdos.shape[1]})")

    folder = os.path.dirname(os.path.abspath(path))
    ncl = detect_ncl(folder) if ncl_flag is None else ncl_flag
    pdos = []
    pos = 6 + nedos
    for ia in range(natoms):
        if pos >= len(lines) or not lines[pos].strip():
            break
        if pos + 1 + nedos > len(lines):
            warn(f"{path}: PDOS of atom {ia + 1} is truncated; ignoring the remaining atoms")
            break
        blk = parse_block(lines[pos + 1: pos + 1 + nedos], f"PDOS of atom {ia + 1}")
        pos += 1 + nedos
        ncol = blk.shape[1] - 1
        if ispin == 2:
            ncomp = 2
        elif ncl and ncol % 4 == 0 and (ncol // 4) in ORB_SETS:
            ncomp = 4
        elif ncol in ORB_SETS:
            ncomp = 1
        elif ncol % 4 == 0 and (ncol // 4) in ORB_SETS:
            ncomp = 4
            if not ncl:
                warn("PDOS column count only fits a non-collinear run; treating it as such")
                ncl = True
        else:
            die(f"{path}: cannot interpret {ncol} PDOS columns for atom {ia + 1}")
        norb = ncol // ncomp
        names = ORB_SETS.get(norb) or [f"orb{i + 1}" for i in range(norb)]
        data = blk[:, 1:].reshape(nedos, norb, ncomp)
        # store (nspin_out, nedos, norb); non-collinear -> total component only
        arr = np.moveaxis(data[:, :, :2] if ncomp == 2 else data[:, :, :1], 2, 0)
        pdos.append({"names": names, "data": arr})
    if natoms and not pdos:
        warn("DOSCAR contains no projected DOS (set LORBIT = 10 or 11)")
    elif pdos and len(pdos) != natoms:
        warn(f"PDOS found for {len(pdos)} of {natoms} atoms")
    return {"energies": energies, "tdos": tot, "itdos": itot, "efermi": efermi, "ispin": ispin,
            "ncl": ncl, "pdos": pdos, "natoms": natoms, "nedos": nedos}


# --------------------------------------------------------------------------- selections
def orbital_indices(token: str, names: List[str]) -> List[int]:
    t = token.strip().lower()
    low = [n.lower() for n in names]
    if t in ("all", "tot", "total", "*", ""):
        return list(range(len(names)))
    groups = {
        "s": [i for i, n in enumerate(low) if n == "s"],
        "p": [i for i, n in enumerate(low) if n.startswith("p")],
        "d": [i for i, n in enumerate(low) if n.startswith("d") or n == "x2-y2"],
        "f": [i for i, n in enumerate(low) if n.startswith("f")],
        "t2g": [i for i, n in enumerate(low) if n in ("dxy", "dyz", "dxz")],
        "eg": [i for i, n in enumerate(low) if n in ("dz2", "x2-y2")],
    }
    if t in groups:
        if t in ("t2g", "eg") and not groups[t]:
            die(f"'{t}' needs lm-resolved PDOS (LORBIT = 11); this DOSCAR only has {names}")
        return groups[t]
    t = {"dx2-y2": "x2-y2", "dx2": "x2-y2"}.get(t, t)
    if t in low:
        return [low.index(t)]
    die(f"unknown orbital '{token}'. Available: {', '.join(names)} (+ s,p,d,f,t2g,eg,all)")
    return []


def atom_indices(spec: str, symbols: Optional[List[str]], nat: int) -> List[int]:
    out: List[int] = []
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok.lower() in ("all", "*"):
            out += list(range(nat))
        elif re.match(r"^\d+(-\d+)?$", tok):
            a, _, b = tok.partition("-")
            lo, hi = int(a), int(b or a)
            if lo < 1 or hi > nat or lo > hi:
                die(f"atom range '{tok}' outside 1..{nat}")
            out += list(range(lo - 1, hi))
        else:
            if symbols is None:
                die("element names need POSCAR/CONTCAR next to DOSCAR (or --poscar)")
            idx = [i for i, s in enumerate(symbols) if s.lower() == tok.lower()]
            if not idx:
                die(f"element '{tok}' not found; elements present: {sorted(set(symbols))}")
            out += idx
    return sorted(set(out))


def select(dos: Dict, spec: str, symbols: Optional[List[str]], default_orb: str = "all") -> Tuple[np.ndarray, str]:
    """Sum PDOS over the selected atoms/orbitals. Returns (nspin, nedos) array."""
    if not dos["pdos"]:
        die("no PDOS in DOSCAR (LORBIT = 10/11 required for projections)")
    atoms_s, _, orbs_s = spec.partition(":")
    orbs_s = orbs_s or default_orb
    ia = atom_indices(atoms_s, symbols, len(dos["pdos"]))
    total = None
    for i in ia:
        names = dos["pdos"][i]["names"]
        io: List[int] = []
        for t in orbs_s.split(","):
            io += orbital_indices(t, names)
        part = dos["pdos"][i]["data"][:, :, sorted(set(io))].sum(axis=2)
        total = part if total is None else total + part
    return total, spec


def gaussian_smear(y: np.ndarray, e: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0 or len(e) < 3:
        return y
    de = float(np.mean(np.diff(e)))
    half = int(np.ceil(4 * sigma / de))
    xk = np.arange(-half, half + 1) * de
    k = np.exp(-0.5 * (xk / sigma) ** 2)
    k /= k.sum()
    return np.array([np.convolve(row, k, mode="same") for row in np.atleast_2d(y)])


# --------------------------------------------------------------------------- band moments
def band_moments(e_rel: np.ndarray, rho: np.ndarray, window: Optional[Tuple[float, float]]) -> Dict:
    lo, hi = (window if window else (e_rel[0], e_rel[-1]))
    m = (e_rel >= lo) & (e_rel <= hi)
    if m.sum() < 3:
        return {"ok": False}
    e, r = e_rel[m], np.clip(rho[m], 0, None)
    n0 = trapz(r, e)
    if n0 <= 1e-12:
        return {"ok": False}
    c = trapz(e * r, e) / n0
    w = np.sqrt(trapz((e - c) ** 2 * r, e) / n0)
    occ = e <= 0.0
    filling = trapz(r[occ], e[occ]) / n0 if occ.sum() > 1 else 0.0
    res = {"ok": True, "center": c, "width": w, "filling": filling, "states": n0}
    if occ.sum() > 2:
        n1 = trapz(r[occ], e[occ])
        if n1 > 1e-12:
            res["center_occ"] = trapz(e[occ] * r[occ], e[occ]) / n1
    return res


def report_dband(dos: Dict, spec: str, symbols, window, smear: float) -> List[Tuple[str, float]]:
    rho, _ = select(dos, spec if ":" in spec else spec + ":d", symbols)
    e_rel = dos["energies"] - dos["efermi_used"]
    rho = gaussian_smear(rho, e_rel, smear)
    wtxt = f"[{window[0]:g}, {window[1]:g}] eV" if window else "full energy range of DOSCAR"
    print(f"d-band analysis for '{spec if ':' in spec else spec + ':d'}' (E relative to E-fermi, window {wtxt})")
    rows = []
    chans = [("total", rho.sum(axis=0))]
    if rho.shape[0] == 2:
        chans += [("spin up", rho[0]), ("spin down", rho[1])]
    centers = []
    for name, r in chans:
        res = band_moments(e_rel, r, window)
        if not res["ok"]:
            print(f"  {name:9s}: no states in the window")
            continue
        extra = f", occupied-part centre {res['center_occ']:.3f} eV" if "center_occ" in res else ""
        print(f"  {name:9s}: centre = {res['center']:7.3f} eV, RMS width = {res['width']:6.3f} eV, "
              f"filling = {res['filling']:.3f}, states = {res['states']:.3f}{extra}")
        rows.append((name, res["center"]))
        if name == "total":
            centers.append((spec, res["center"]))
    return centers


# --------------------------------------------------------------------------- plotting
def setup_matplotlib(show: bool):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot(args, e_rel, tdos, curves, centers, plt) -> None:
    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    nspin = tdos.shape[0]
    sign = [1.0, -1.0]
    xlo, xhi = args.xlim
    vis = (e_rel >= xlo) & (e_rel <= xhi)
    ymax = 0.0
    if not args.no_tdos:
        for s in range(nspin):
            y = sign[s] * tdos[s] * args.tdos_scale
            ax.fill_between(e_rel, 0, y, color=TDOS_FILL, lw=0, zorder=1)
            lab = "Total" + (f" (x{args.tdos_scale:g})" if args.tdos_scale != 1 else "")
            ax.plot(e_rel, y, color=TDOS_LINE, lw=0.8, zorder=2, label=lab if s == 0 else None)
            ymax = max(ymax, float(np.abs(y[vis]).max()) if vis.any() else 0)
    for j, (lab, y2) in enumerate(curves):
        c = SERIES[j % len(SERIES)]
        for s in range(nspin):
            y = sign[s] * y2[s]
            ax.plot(e_rel, y, color=c, lw=1.4, zorder=3, label=lab if s == 0 else None)
            ymax = max(ymax, float(np.abs(y[vis]).max()) if vis.any() else 0)
    for j, (lab, c0) in enumerate(centers):
        ax.axvline(c0, color=INK, lw=0.9, ls=":", zorder=4)
        ax.annotate(f"$\\varepsilon_d$ = {c0:.2f} eV", xy=(c0, 0.92 if nspin == 1 else 0.96),
                    xycoords=("data", "axes fraction"), ha="right", va="top", fontsize=8,
                    rotation=90, color=INK,
                    bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.0))
    ax.axvline(0.0, color=MUTED, lw=0.8, ls="--", zorder=2)
    if nspin == 2:
        ax.axhline(0.0, color=INK, lw=0.6, zorder=2)
    ax.set_xlim(xlo, xhi)
    if args.ylim:
        ax.set_ylim(*args.ylim)
    elif ymax > 0:
        ax.set_ylim(-1.08 * ymax if nspin == 2 else 0, 1.08 * ymax)
    ax.set_xlabel("$E - E_\\mathrm{F}$ (eV)")
    ax.set_ylabel("DOS (states/eV)")
    ax.tick_params(direction="in", length=3)
    if nspin == 2:
        ax.text(0.01, 0.98, "spin up", transform=ax.transAxes, ha="left", va="top", fontsize=8, color=MUTED)
        ax.text(0.01, 0.02, "spin down", transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color=MUTED)
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) >= 2:
        ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.9)
    if args.title:
        ax.set_title(args.title)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)
    print(f"Figure written: {args.output}")
    if args.show:
        plt.show()


def write_dat(path: str, e_rel, tdos, curves, efermi: float) -> None:
    nspin = tdos.shape[0]
    cols = [e_rel]
    names = ["E-Ef(eV)"]
    sfx = ["_up", "_down"] if nspin == 2 else [""]
    for s in range(nspin):
        cols.append(tdos[s])
        names.append("TDOS" + sfx[s])
    for lab, y in curves:
        safe = re.sub(r"[^A-Za-z0-9_+-]", "_", lab)
        for s in range(nspin):
            cols.append(y[s])
            names.append(safe + sfx[s])
    with open(path, "w") as fh:
        fh.write(f"# DOS exported by plot_dos.py, E-fermi = {efermi:.6f} eV (energies relative to it)\n")
        fh.write(" ".join(names) + "\n")
        np.savetxt(fh, np.column_stack(cols), fmt="%12.6f")
    print(f"Data written: {path}")


# --------------------------------------------------------------------------- main
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot TDOS/PDOS from a VASP DOSCAR (LORBIT = 10/11) and compute d-band centres.\n"
                    "读取 DOSCAR 画总/投影态密度（按元素、原子、轨道选择），计算 d 带中心与带宽。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Selection syntax / 选择语法:  ATOMS:ORBITALS\n"
               "  ATOMS    = element (Ni), indices (1-4,9) or all\n"
               "  ORBITALS = s,p,d,f,t2g,eg,all or lm names (px, dz2, x2-y2 ...); default all\n"
               "Examples / 示例:\n"
               "  python plot_dos.py -p 'Ni:d' -p 'O:p' --xlim -8 4\n"
               "  python plot_dos.py -p 'Ni:t2g' -p 'Ni:eg' --smear 0.05\n"
               "  python plot_dos.py --dband Pt --dband-window -10 2 --no-plot\n")
    p.add_argument("-f", "--file", default="DOSCAR", help="DOSCAR file (default ./DOSCAR)")
    p.add_argument("--poscar", default=None,
                   help="POSCAR/CONTCAR for element names (default: CONTCAR/POSCAR next to DOSCAR)")
    p.add_argument("-p", "--pdos", action="append", default=[],
                   help="projection to plot, repeatable, e.g. -p 'Ni:d' -p '1-4:eg' / 投影选择")
    p.add_argument("--by-element", action="store_true",
                   help="plot the PDOS of every element (default when no -p is given and PDOS exists)")
    p.add_argument("--by-orbital", default=None, metavar="ORBS",
                   help="with --by-element: split each element into these orbitals, e.g. 'spd' or 's,p,d'")
    p.add_argument("--efermi", type=float, default=None, help="override E-fermi (eV) / 覆盖费米能")
    p.add_argument("--ncl", action="store_true", default=None,
                   help="force non-collinear interpretation of PDOS columns (auto-detected from OUTCAR/INCAR)")
    p.add_argument("--smear", type=float, default=0.0, help="extra Gaussian broadening sigma in eV (default 0)")
    p.add_argument("--dband", action="append", default=[],
                   help="compute d-band centre/width for a selection, e.g. --dband Pt or --dband '1-9:d'")
    p.add_argument("--dband-window", type=float, nargs=2, default=None, metavar=("EMIN", "EMAX"),
                   help="integration window relative to E-fermi (default: whole DOSCAR range)")
    p.add_argument("--xlim", type=float, nargs=2, default=[-8.0, 4.0], metavar=("EMIN", "EMAX"))
    p.add_argument("--ylim", type=float, nargs=2, default=None, metavar=("YMIN", "YMAX"))
    p.add_argument("--no-tdos", action="store_true", help="do not draw the total DOS")
    p.add_argument("--tdos-scale", type=float, default=1.0, help="scale factor for the total DOS curve")
    p.add_argument("--figsize", type=float, nargs=2, default=[5.5, 3.8], metavar=("W", "H"))
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--title", default=None)
    p.add_argument("-o", "--output", default="dos.png", help="figure file (default dos.png)")
    p.add_argument("--dat", default=None, help="export columns (E-Ef, TDOS, PDOS...) to this file")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--show", action="store_true", help="open an interactive window")
    return p


def main(argv=None) -> None:
    try:
        sys.stdout.reconfigure(errors="backslashreplace")  # type: ignore[attr-defined]
    except (AttributeError, ValueError):
        pass
    args = build_parser().parse_args(argv)
    dos = read_doscar(args.file, True if args.ncl else None)
    folder = os.path.dirname(os.path.abspath(args.file))
    ef = dos["efermi"] if args.efermi is None else args.efermi
    dos["efermi_used"] = ef
    e_rel = dos["energies"] - ef

    symbols = None
    pos = args.poscar
    if pos is None:
        for cand in ("CONTCAR", "POSCAR"):
            c = os.path.join(folder, cand)
            if os.path.isfile(c) and os.path.getsize(c) > 0:
                pos = c
                break
    if pos:
        symbols = read_poscar_symbols(pos)
        if dos["pdos"] and len(symbols) != len(dos["pdos"]):
            warn(f"{pos} has {len(symbols)} atoms but DOSCAR has PDOS for {len(dos['pdos'])}")

    kind = "non-collinear" if dos["ncl"] else ("spin-polarised" if dos["ispin"] == 2 else "non-spin-polarised")
    orbs = dos["pdos"][0]["names"] if dos["pdos"] else []
    print(f"DOSCAR: {args.file}  NEDOS = {dos['nedos']}  atoms = {dos['natoms']}  {kind}")
    print(f"E-fermi = {ef:.4f} eV ({'command line' if args.efermi is not None else 'DOSCAR line 6'}); "
          f"PDOS orbitals: {' '.join(orbs) if orbs else 'none'}")
    de = float(np.mean(np.diff(dos["energies"])))
    print(f"Energy grid: {dos['energies'][0]:.3f} .. {dos['energies'][-1]:.3f} eV, step {de:.4f} eV")
    i0 = int(np.argmin(np.abs(e_rel)))
    if dos["ispin"] == 2:
        print(f"Integrated DOS at E-fermi: up {dos['itdos'][0][i0]:.3f}, down {dos['itdos'][1][i0]:.3f} "
              f"-> moment ~ {dos['itdos'][0][i0] - dos['itdos'][1][i0]:.3f} muB")
    else:
        print(f"Integrated DOS at E-fermi: {dos['itdos'][0][i0]:.3f} electrons")

    tdos = gaussian_smear(dos["tdos"], e_rel, args.smear)
    curves: List[Tuple[str, np.ndarray]] = []
    specs = list(args.pdos)
    if (args.by_element or not specs) and dos["pdos"]:
        if symbols is None:
            if args.by_element:
                die("--by-element needs POSCAR/CONTCAR (use --poscar)")
        else:
            uniq = list(dict.fromkeys(symbols))
            if args.by_orbital:
                orb_list = [o for o in re.split(r"[,\s]+", args.by_orbital) if o]
                if len(orb_list) == 1 and len(orb_list[0]) > 1 and orb_list[0].lower() not in ("t2g", "eg", "all"):
                    orb_list = list(orb_list[0])
                specs += [f"{el}:{o}" for el in uniq for o in orb_list]
            else:
                specs += [f"{el}:all" for el in uniq]
    for spec in specs:
        y, lab = select(dos, spec, symbols)
        curves.append((lab.replace(":all", ""), gaussian_smear(y, e_rel, args.smear)))
    if len(curves) > len(SERIES):
        warn(f"{len(curves)} curves requested; colours repeat after {len(SERIES)} - consider fewer projections")

    centers = []
    for spec in args.dband:
        centers += report_dband(dos, spec, symbols,
                                tuple(args.dband_window) if args.dband_window else None, args.smear)

    if args.dat:
        write_dat(args.dat, e_rel, tdos, curves, ef)
    if not args.no_plot:
        plt = setup_matplotlib(args.show)
        plot(args, e_rel, tdos, curves, centers, plt)


if __name__ == "__main__":
    main()
