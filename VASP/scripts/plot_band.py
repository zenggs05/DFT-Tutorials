#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_band.py -- plot VASP band structures and report the band gap.
               绘制 VASP 能带并报告带隙 / VBM / CBM。

Purpose
-------
* Read eigenvalues from ``vasprun.xml`` (including the ``<eigenvalues_kpoints_opt>``
  block written by VASP >= 6.3 when a KPOINTS_OPT file is present) or from
  ``EIGENVAL`` (+ POSCAR/CONTCAR for the reciprocal lattice).
* Spin-polarised (ISPIN = 2) and non-collinear runs are supported.
* Hybrid-functional (HSE) band runs with zero-weight k-points: the weighted
  SCF mesh points are removed automatically.
* High-symmetry labels are taken from a line-mode KPOINTS / KPOINTS_OPT file,
  from labels in an explicit k-point list (e.g. "0.5 0.0 0.5 0 X"), from a
  vaspkit KPATH.in file, or from ``--labels``.
* Fat bands from PROCAR / PROCAR_OPT (LORBIT = 11/12), e.g. ``--proj "Ga:s"``.
* Band gap, VBM, CBM (per spin channel and overall), direct vs. indirect.
* Export of Origin/gnuplot friendly column data (``--dat``).

Usage examples
--------------
  python plot_band.py                                  # ./vasprun.xml (+ ./KPOINTS labels)
  python plot_band.py -f band/vasprun.xml --efermi scf/OUTCAR -o Si_band.png
  python plot_band.py -f EIGENVAL --kpoints KPOINTS --poscar POSCAR --efermi 5.76
  python plot_band.py -f vasprun.xml --opt              # KPOINTS_OPT band (VASP >= 6.3)
  python plot_band.py -f vasprun.xml --kpath KPATH.in   # HSE zero-weight run (labels from vaspkit)
  python plot_band.py --procar PROCAR --proj "Mo:d" --proj "S:p" --ylim -3 3
  python plot_band.py --labels "L G X U|K G" --zero vbm --dat band.dat

Dependencies
------------
  Python >= 3.9, numpy, matplotlib (xml.etree from the standard library).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple

import numpy as np

FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?")
# Categorical colours (fixed order) used for projections.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK = "#1f1f1f"
MUTED = "#8a8986"
GRID = "#c9c8c4"
LM_NAMES = ["s", "py", "pz", "px", "dxy", "dyz", "dz2", "dxz", "x2-y2",
            "fy3x2", "fxyz", "fyz2", "fz3", "fxz2", "fzx2", "fx3"]


# --------------------------------------------------------------------------- utilities
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


def floats_in(line: str) -> List[float]:
    """Extract all floats of a line (robust to glued numbers like -0.5-0.5)."""
    return [float(x.replace("D", "E").replace("d", "e")) for x in FLOAT_RE.findall(line)]


def to_float(tok: str) -> float:
    try:
        return float(tok)
    except ValueError:
        return float("nan")  # e.g. '********' overflow in Fortran output


def same_k(a, b, tol: float = 1e-4) -> bool:
    return bool(np.all(np.abs(np.asarray(a) - np.asarray(b)) < tol))


def nice_label(lab: Optional[str]) -> Optional[str]:
    """Normalise a k-point label for plotting (G/GAMMA/\\Gamma -> Greek Gamma)."""
    if lab is None:
        return None
    s = lab.strip().strip("$").strip()
    if not s:
        return None
    s = s.replace("\\mathrm", "").replace("{", "").replace("}", "")
    parts = []
    for p in s.split("|"):
        q = p.strip()
        if q.upper() in ("G", "GAMMA", "\\GAMMA", "GM", "Γ"):
            q = "Γ"
        q = q.replace("\\Gamma", "Γ").replace("\\Sigma", "Σ").replace("\\Delta", "Δ")
        parts.append(q)
    s = "|".join(parts)
    if "_" in s:
        s = "$\\mathregular{" + s.replace("Γ", "\\Gamma") + "}$"
    return s


def ascii_label(lab: Optional[str]) -> str:
    """Label usable in console output (avoid encoding problems)."""
    if not lab:
        return ""
    return (lab.replace("Γ", "Gamma").replace("$\\mathregular{", "")
            .replace("}$", "").replace("\\Gamma", "Gamma"))


def reciprocal(lattice: np.ndarray) -> np.ndarray:
    """Rows are b1, b2, b3 including the 2*pi factor (1/Angstrom)."""
    return 2.0 * np.pi * np.linalg.inv(lattice).T


# --------------------------------------------------------------------------- POSCAR
def read_poscar(path: str) -> Tuple[np.ndarray, List[str]]:
    """Return (lattice[3x3] in Angstrom, list of element symbols per atom)."""
    check_file(path, "POSCAR/CONTCAR")
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    try:
        scale = floats_in(lines[1])
        lat = np.array([floats_in(lines[i])[:3] for i in (2, 3, 4)], dtype=float)
        if len(scale) == 3:
            lat = lat * np.array(scale)[None, :]
        elif scale[0] < 0:
            lat = lat * (abs(scale[0]) / abs(np.linalg.det(lat))) ** (1.0 / 3.0)
        else:
            lat = lat * scale[0]
        toks = lines[5].split()
        if toks and not re.match(r"^\d+$", toks[0]):
            names = [re.match(r"[A-Za-z]+", t).group(0) for t in toks]
            names = [re.match(r"[A-Z][a-z]?", n.capitalize()).group(0) for n in names]
            counts = [int(t) for t in lines[6].split()[:len(names)]]
        else:
            counts = [int(t) for t in toks]
            names = [f"X{i + 1}" for i in range(len(counts))]
            warn(f"{path}: no element line (VASP 4 format); species named X1, X2, ...")
    except (IndexError, ValueError, AttributeError) as exc:
        die(f"cannot parse {path} as POSCAR: {exc}")
    symbols: List[str] = []
    for n, c in zip(names, counts):
        symbols += [n] * c
    return lat, symbols


# --------------------------------------------------------------------------- vasprun.xml
def _parse_kpoints_elem(elem) -> Dict:
    out: Dict = {"kpoints": None, "weights": None, "divisions": None,
                 "generation": None, "vertices": []}
    gen = elem.find("generation")
    if gen is not None:
        out["generation"] = gen.get("param")
        for child in gen:
            if child.get("name") == "divisions" and child.text:
                vals = child.text.split()
                if len(vals) == 1:
                    out["divisions"] = int(float(vals[0]))
            elif child.tag == "v" and child.get("name") is None and child.text:
                out["vertices"].append([float(x) for x in child.text.split()[:3]])
    for va in elem.findall("varray"):
        rows = [[to_float(x) for x in v.text.split()] for v in va.findall("v")]
        if va.get("name") == "kpointlist":
            out["kpoints"] = np.array(rows, dtype=float)
        elif va.get("name") == "weights":
            out["weights"] = np.array(rows, dtype=float).reshape(-1)
    return out


def _parse_eigen_elem(elem) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    arr = elem.find("array")
    if arr is None:
        return None, None
    top = arr.find("set")
    if top is None:
        return None, None
    eig_s, occ_s = [], []
    for sp in top.findall("set"):
        e_k, o_k = [], []
        for kp in sp.findall("set"):
            rows = [[to_float(x) for x in r.text.split()] for r in kp.findall("r")]
            a = np.array(rows, dtype=float)
            e_k.append(a[:, 0])
            o_k.append(a[:, 1] if a.shape[1] > 1 else np.full(a.shape[0], np.nan))
        eig_s.append(e_k)
        occ_s.append(o_k)
    try:
        return np.array(eig_s, dtype=float), np.array(occ_s, dtype=float)
    except ValueError:
        die("inconsistent number of bands/k-points in vasprun.xml eigenvalue block")
    return None, None


def parse_vasprun(path: str) -> Dict:
    """Stream-parse vasprun.xml and keep only what band plotting needs."""
    check_file(path, "vasprun.xml")
    d: Dict = {"main": None, "opt": None, "lattice": None, "efermi": None,
               "efermi_opt": None, "symbols": None, "eig": None, "occ": None,
               "nelect": None, "lsorbit": False}
    stack: List[str] = []
    try:
        for event, elem in ET.iterparse(path, events=("start", "end")):
            if event == "start":
                stack.append(elem.tag)
                continue
            tag = elem.tag
            stack.pop()
            parent = stack[-1] if stack else None
            if tag == "kpoints" and parent == "modeling":
                d["main"] = _parse_kpoints_elem(elem)
            elif tag == "i" and elem.get("name") in ("NELECT",) and elem.text:
                d["nelect"] = to_float(elem.text)
            elif tag == "i" and elem.get("name") in ("LSORBIT", "LNONCOLLINEAR") and elem.text:
                if elem.text.strip().upper().startswith("T"):
                    d["lsorbit"] = True
            elif tag == "atominfo":
                arr = elem.find("array[@name='atoms']")
                if arr is not None:
                    d["symbols"] = [rc.find("c").text.strip() for rc in arr.iter("rc")]
            elif tag == "structure":
                basis = elem.find("crystal/varray[@name='basis']")
                if basis is not None:
                    d["lattice"] = np.array([[float(x) for x in v.text.split()]
                                             for v in basis.findall("v")])
            elif tag == "eigenvalues" and parent == "calculation":
                d["eig"], d["occ"] = _parse_eigen_elem(elem)
                elem.clear()
            elif tag == "eigenvalues_kpoints_opt":
                kp = elem.find("kpoints")
                ev = elem.find("eigenvalues")
                opt = _parse_kpoints_elem(kp) if kp is not None else {}
                if ev is not None:
                    opt["eig"], opt["occ"] = _parse_eigen_elem(ev)
                d["opt"] = opt
                elem.clear()
            elif tag == "dos":
                ef = elem.find("i[@name='efermi']")
                if ef is not None and ef.text:
                    if elem.get("comment") == "kpoints_opt":
                        d["efermi_opt"] = to_float(ef.text)
                    else:
                        d["efermi"] = to_float(ef.text)
                elem.clear()
            elif tag in ("projected", "projected_kpoints_opt", "dielectricfunction",
                         "dynmat", "calculation"):
                elem.clear()
    except ET.ParseError as exc:
        warn(f"{path} is incomplete or corrupted ({exc}); using the data read so far.")
    if d["eig"] is None and (d["opt"] is None or d["opt"].get("eig") is None):
        die(f"no eigenvalues found in {path} (job not finished?)")
    return d


# --------------------------------------------------------------------------- EIGENVAL
def parse_eigenval(path: str) -> Dict:
    check_file(path, "EIGENVAL")
    with open(path, "r", errors="replace") as fh:
        lines = fh.read().splitlines()
    try:
        ispin = int(lines[0].split()[3])
        h = lines[5].split()
        nelect, nk, nb = float(h[0]), int(h[1]), int(h[2])
    except (IndexError, ValueError):
        die(f"{path}: header is not in EIGENVAL format")
    kpts = np.zeros((nk, 3))
    wts = np.zeros(nk)
    eig = np.zeros((ispin, nk, nb))
    occ = np.full((ispin, nk, nb), np.nan)
    i = 6
    for ik in range(nk):
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i + nb >= len(lines) + 1:
            die(f"{path}: truncated at k-point {ik + 1}/{nk}")
        v = floats_in(lines[i])
        kpts[ik] = v[:3]
        wts[ik] = v[3] if len(v) > 3 else np.nan
        i += 1
        for ib in range(nb):
            t = [to_float(x) for x in lines[i].split()[1:]]
            i += 1
            if ispin == 1:
                eig[0, ik, ib] = t[0]
                if len(t) > 1:
                    occ[0, ik, ib] = t[1]
            else:
                eig[0, ik, ib], eig[1, ik, ib] = t[0], t[1]
                if len(t) >= 4:
                    occ[0, ik, ib], occ[1, ik, ib] = t[2], t[3]
    return {"kpoints": kpts, "weights": wts, "eig": eig, "occ": occ, "nelect": nelect}


# --------------------------------------------------------------------------- KPOINTS
def read_kpoints_file(path: str) -> Optional[Dict]:
    """Parse a KPOINTS-like file. Returns dict with mode 'line'/'explicit'/'auto'."""
    if not path or not os.path.isfile(path):
        return None
    with open(path, "r", errors="replace") as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    if len(lines) < 3:
        return None

    def point(line: str):
        label = None
        body = line
        if "!" in line:
            body, lab = line.split("!", 1)
            label = lab.strip() or None
        nums, words = [], []
        for t in body.split():
            try:
                nums.append(float(t))
            except ValueError:
                words.append(t)
        if label is None and words:
            label = " ".join(words)
        return nums, label

    mode = lines[2].strip().lower()
    if mode.startswith("l"):
        try:
            nper = int(lines[1].split()[0])
        except (IndexError, ValueError):
            return None
        cart = lines[3].strip().lower()[:1] in ("c", "k")
        pts = []
        for ln in lines[4:]:
            if not ln.strip():
                continue
            nums, lab = point(ln)
            if len(nums) >= 3:
                pts.append((nums[:3], lab))
        segs = [(pts[2 * i], pts[2 * i + 1]) for i in range(len(pts) // 2)]
        return {"mode": "line", "N": nper, "cart": cart, "segments": segs}
    try:
        n = int(lines[1].split()[0])
    except (IndexError, ValueError):
        return None
    if n <= 0:
        return {"mode": "auto"}
    cart = lines[2].strip().lower()[:1] in ("c", "k")
    pts = []
    for ln in lines[3:3 + n]:
        nums, lab = point(ln)
        if len(nums) >= 3:
            pts.append((nums[:3], lab))
    return {"mode": "explicit", "cart": cart, "points": pts}


def vertices_from_kfile(info: Optional[Dict]) -> List[Tuple[List[float], str]]:
    if not info or info.get("cart"):
        return []
    out = []
    if info["mode"] == "line":
        for a, b in info["segments"]:
            for c, lab in (a, b):
                if lab:
                    out.append((c, lab))
    elif info["mode"] == "explicit":
        out = [(c, lab) for c, lab in info["points"] if lab]
    return out


# --------------------------------------------------------------------------- PROCAR
def read_procar(path: str) -> Dict:
    """Read PROCAR / PROCAR_OPT (LORBIT = 10/11/12; collinear, spin, non-collinear).

    Returns {'proj': array (nspin, nk, nb, nions, norb), 'orbitals': names}.
    For non-collinear runs only the total (first) block is kept.
    """
    check_file(path, "PROCAR")
    spins: List[np.ndarray] = []
    names: Optional[List[str]] = None
    nk = nb = nions = 0
    cur = None
    ik = ib = -1
    need = False
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("# of k-points"):
                v = [int(x) for x in re.findall(r"\d+", s)]
                nk, nb, nions = v[0], v[1], v[2]
                cur = None
                spins.append(None)  # placeholder until the orbital header is known
                need = False
            elif s.startswith("k-point"):
                m = re.match(r"k-point\s+(\d+)", s)
                ik = int(m.group(1)) - 1 if m else ik + 1
            elif s.startswith("band"):
                m = re.match(r"band\s+(\d+)", s)
                ib = int(m.group(1)) - 1 if m else ib + 1
                need = True
            elif s.startswith("ion") and need:
                cols = s.split()[1:]
                if cols and cols[-1] == "tot":
                    cols = cols[:-1]
                if names is None:
                    names = cols
                if cur is None:
                    if not spins:
                        die(f"{path}: missing '# of k-points' header")
                    cur = np.zeros((nk, nb, nions, len(names)), dtype=np.float32)
                    spins[-1] = cur
                norb = len(names)
                for _ in range(nions):
                    row = next(fh).split()
                    try:
                        ia = int(row[0]) - 1
                        cur[ik, ib, ia, :] = [float(x) for x in row[1:1 + norb]]
                    except (ValueError, IndexError):
                        die(f"{path}: malformed projection line at k-point {ik + 1}, band {ib + 1}")
                need = False
    spins = [s for s in spins if s is not None]
    if not spins or names is None:
        die(f"{path}: no projections found (was LORBIT >= 10 set?)")
    return {"proj": np.array(spins), "orbitals": names}


def orbital_indices(token: str, names: List[str]) -> List[int]:
    t = token.strip().lower()
    low = [n.lower() for n in names]
    if t in ("all", "tot", "total", "*"):
        return list(range(len(names)))
    groups = {
        "s": [i for i, n in enumerate(low) if n == "s"],
        "p": [i for i, n in enumerate(low) if n.startswith("p")],
        "d": [i for i, n in enumerate(low) if n.startswith("d") or n in ("x2-y2",)],
        "f": [i for i, n in enumerate(low) if n.startswith("f")],
        "t2g": [i for i, n in enumerate(low) if n in ("dxy", "dyz", "dxz")],
        "eg": [i for i, n in enumerate(low) if n in ("dz2", "x2-y2", "dx2", "dx2-y2")],
    }
    if t in groups:
        if t in ("t2g", "eg") and len(groups[t]) == 0:
            die(f"orbital group '{t}' needs lm-resolved projections (LORBIT = 11/12)")
        return groups[t]
    aliases = {"dx2-y2": "x2-y2", "dx2": "x2-y2"}
    t = aliases.get(t, t)
    if t in low:
        return [low.index(t)]
    if t == "x2-y2" and "dx2" in low:
        return [low.index("dx2")]
    die(f"unknown orbital '{token}'. Available: {', '.join(names)} (+ s,p,d,f,t2g,eg,all)")
    return []


def atom_indices(token_str: str, symbols: Optional[List[str]], nions: int) -> List[int]:
    out: List[int] = []
    for tok in token_str.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok.lower() in ("all", "*"):
            out += list(range(nions))
        elif re.match(r"^\d+(-\d+)?$", tok):
            a, _, b = tok.partition("-")
            lo, hi = int(a), int(b or a)
            if lo < 1 or hi > nions or lo > hi:
                die(f"atom range '{tok}' outside 1..{nions}")
            out += list(range(lo - 1, hi))
        else:
            if symbols is None:
                die("element names in --proj need POSCAR/CONTCAR or vasprun.xml (use --poscar)")
            idx = [i for i, s in enumerate(symbols) if s.lower() == tok.lower()]
            if not idx:
                die(f"element '{tok}' not found; elements present: {sorted(set(symbols))}")
            out += idx
    return sorted(set(out))


# --------------------------------------------------------------------------- k-path geometry
def build_kpath(kfrac: np.ndarray, lattice: Optional[np.ndarray], kinfo: Optional[Dict],
                nper: Optional[int], extra_vertices, manual_labels: Optional[str],
                jump_factor: float = 5.0):
    """Return x (cumulative distance), ticks [(x, label)], boundaries."""
    nk = len(kfrac)
    if lattice is None:
        warn("no lattice available: k-distances are in fractional units (use --poscar)")
        kc = kfrac.copy()
    else:
        kc = kfrac @ reciprocal(lattice)
    steps = np.linalg.norm(np.diff(kc, axis=0), axis=1) if nk > 1 else np.zeros(0)
    jump = np.zeros(max(nk - 1, 0), dtype=bool)
    dup = steps < 1e-8
    labels_at: Dict[int, str] = {}
    boundaries = {0, nk - 1}

    use_line = (nper is not None and nper > 1 and nk % nper == 0)
    if use_line:
        nseg = nk // nper
        segs = kinfo["segments"] if (kinfo and kinfo.get("mode") == "line") else None
        for s in range(nseg):
            a, b = s * nper, s * nper + nper - 1
            boundaries.update((a, b))
            if s > 0 and not same_k(kfrac[a - 1], kfrac[a]):
                jump[a - 1] = True
            if segs is not None and s < len(segs) and not kinfo.get("cart"):
                if segs[s][0][1]:
                    labels_at[a] = segs[s][0][1]
                if segs[s][1][1]:
                    labels_at[b] = segs[s][1][1]
    else:
        verts = list(extra_vertices or [])
        for i, k in enumerate(kfrac):
            for c, lab in verts:
                if same_k(k, c):
                    labels_at[i] = lab
                    break
        for i in range(nk - 1):
            if not dup[i] and i in labels_at and (i + 1) in labels_at:
                jump[i] = True
        good = steps[~dup] if nk > 1 else np.zeros(0)
        if good.size:
            med = float(np.median(good))
            jump |= steps > jump_factor * med
        for i in range(nk - 1):
            if dup[i] or jump[i]:
                boundaries.update((i, i + 1))
        boundaries.update(labels_at.keys())

    steps = steps.copy()
    steps[jump] = 0.0
    x = np.concatenate([[0.0], np.cumsum(steps)])

    # Group labels by x position.
    ticks: List[Tuple[float, str]] = []
    if manual_labels:
        labs = manual_labels.split()
        bx = sorted({round(float(x[i]), 6) for i in boundaries})
        if len(labs) != len(bx):
            warn(f"--labels has {len(labs)} entries but the path has {len(bx)} "
                 f"segment boundaries at x = {', '.join(f'{v:.3f}' for v in bx)}; labels ignored")
        else:
            ticks = [(v, nice_label(lab) or "") for v, lab in zip(bx, labs)]
    if not ticks:
        byx: Dict[float, List[str]] = {}
        for i in sorted(boundaries):
            key = round(float(x[i]), 6)
            byx.setdefault(key, [])
            lab = labels_at.get(i)
            if lab:
                lab = nice_label(lab)
                if lab not in byx[key]:
                    byx[key].append(lab)
        ticks = [(k, "|".join(v)) for k, v in sorted(byx.items())]
    return x, ticks, labels_at


# --------------------------------------------------------------------------- gap analysis
def gap_report(eig: np.ndarray, occ: Optional[np.ndarray], kfrac: np.ndarray, efermi: Optional[float],
               method: str, nocc: Optional[int], tol: float, labels_at: Dict[int, str]) -> Dict:
    """Determine VBM/CBM/gap. eig: (nspin, nk, nb)."""
    nspin, nk, nb = eig.shape
    if nocc is not None:
        if not 0 < nocc < nb:
            die(f"--nocc {nocc} must be between 1 and NBANDS-1 = {nb - 1}")
        occmask = np.zeros_like(eig, dtype=bool)
        occmask[:, :, :nocc] = True
        used = f"first {nocc} bands occupied (--nocc)"
    elif method == "occ":
        occmask = occ > 0.5
        used = "occupations (occ > 0.5)"
    else:
        if efermi is None:
            return {"ok": False, "why": "E-fermi unknown (give --efermi or --nocc)"}
        occmask = eig <= efermi + tol
        used = f"E-fermi = {efermi:.4f} eV (tolerance {tol} eV)"

    def one(mask_e, m_occ, spins):
        res = {"metal": False}
        # band crossing: a band partly occupied, partly empty
        for s in spins:
            for b in range(nb):
                col = m_occ[s, :, b]
                if col.any() and (~col).any():
                    res["metal"] = True
                    res["cross"] = (s, b)
                    return res
        ev = np.where(m_occ[spins], mask_e[spins], -np.inf)
        ec = np.where(~m_occ[spins], mask_e[spins], np.inf)
        if not np.isfinite(ev.max()) or not np.isfinite(ec.min()):
            res["metal"] = True
            return res
        # among (near-)degenerate VBM states report the highest band index
        cand = np.argwhere(ev >= ev.max() - 1e-4)
        iv = tuple(cand[np.argmax(cand[:, 2])])
        ic = np.unravel_index(np.argmin(ec), ec.shape)
        res["vbm"] = (float(ev[iv]), spins[iv[0]], int(iv[1]), int(iv[2]))
        res["cbm"] = (float(ec[ic]), spins[ic[0]], int(ic[1]), int(ic[2]))
        res["gap"] = res["cbm"][0] - res["vbm"][0]
        dk = ec.min(axis=(0, 2)) - ev.max(axis=(0, 2))
        kd = int(np.argmin(dk))
        res["direct"] = (float(dk[kd]), kd)
        res["is_direct"] = same_k(kfrac[iv[1]], kfrac[ic[1]])
        return res

    out = {"ok": True, "used": used, "channels": []}
    if nspin == 2:
        for s in (0, 1):
            out["channels"].append((("spin up", "spin down")[s], one(eig, occmask, [s])))
    out["overall"] = one(eig, occmask, list(range(nspin)))
    return out


def print_gap(rep: Dict, kfrac: np.ndarray, labels_at: Dict[int, str], efermi: Optional[float]) -> None:
    print("-" * 64)
    if not rep.get("ok"):
        print(f"Band gap: not analysed ({rep.get('why')})")
        return
    print(f"Occupied states determined from: {rep['used']}")

    def kdesc(ik):
        lab = ascii_label(nice_label(labels_at.get(ik))) if ik in labels_at else ""
        k = kfrac[ik]
        return f"k#{ik + 1} ({k[0]:.4f} {k[1]:.4f} {k[2]:.4f}){' [' + lab + ']' if lab else ''}"

    def show(name, r):
        if r["metal"]:
            extra = ""
            if "cross" in r:
                extra = f" (band {r['cross'][1] + 1} crosses the Fermi level)"
            print(f"[{name}] metallic / no gap{extra}")
            return
        ev, sv, kv, bv = r["vbm"]
        ec, sc, kc, bc = r["cbm"]
        sp = ("up", "down")
        print(f"[{name}] VBM = {ev:9.4f} eV  band {bv + 1:4d} spin {sp[sv]:4s} {kdesc(kv)}")
        print(f"[{name}] CBM = {ec:9.4f} eV  band {bc + 1:4d} spin {sp[sc]:4s} {kdesc(kc)}")
        kind = "direct" if r["is_direct"] else "indirect"
        print(f"[{name}] Band gap = {r['gap']:.4f} eV ({kind})")
        print(f"[{name}] Smallest direct gap = {r['direct'][0]:.4f} eV at {kdesc(r['direct'][1])}")

    for name, r in rep["channels"]:
        show(name, r)
    show("overall", rep["overall"])
    print("-" * 64)


# --------------------------------------------------------------------------- Fermi level
def efermi_from_file(path: str) -> Optional[float]:
    check_file(path, "Fermi-level source")
    base = os.path.basename(path).lower()
    if base.endswith(".xml"):
        ef = None
        try:
            for _, elem in ET.iterparse(path, events=("end",)):
                if elem.tag == "dos" and elem.get("comment") != "kpoints_opt":
                    i = elem.find("i[@name='efermi']")
                    if i is not None:
                        ef = float(i.text)
                    elem.clear()
                elif elem.tag in ("calculation", "projected"):
                    elem.clear()
        except ET.ParseError:
            pass
        return ef
    if base.startswith("doscar"):
        with open(path, "r", errors="replace") as fh:
            for _ in range(5):
                fh.readline()
            v = floats_in(fh.readline())
            return v[3] if len(v) >= 4 else None
    ef = None
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if "E-fermi" in line:
                v = floats_in(line.split(":", 1)[1]) if ":" in line else floats_in(line)
                if v:
                    ef = v[0]
            elif "Fermi energy:" in line and ef is None:
                v = floats_in(line.split(":", 1)[1])
                if v:
                    ef = v[0]
    return ef


def resolve_efermi(arg: Optional[str]) -> Tuple[Optional[float], str]:
    if arg is None:
        return None, ""
    try:
        return float(arg), "command line"
    except ValueError:
        ef = efermi_from_file(arg)
        if ef is None:
            die(f"could not find the Fermi energy in {arg}")
        return ef, arg


# --------------------------------------------------------------------------- output
def write_dat(path: str, x: np.ndarray, eig: np.ndarray, ticks, zero: float, zero_name: str,
              proj_list) -> None:
    nspin = eig.shape[0]
    root, ext = os.path.splitext(path)
    ext = ext or ".dat"
    written = []
    for s in range(nspin):
        fname = path if nspin == 1 else f"{root}_{('up', 'down')[s]}{ext}"
        nb = eig.shape[2]
        with open(fname, "w") as fh:
            fh.write(f"# Band structure exported by plot_band.py; energies in eV relative to "
                     f"{zero_name} = {zero:.6f} eV\n")
            fh.write("# High-symmetry ticks: " +
                     "; ".join(f"{t:.6f} {ascii_label(lab) or '-'}" for t, lab in ticks) + "\n")
            fh.write("k_dist(1/A) " + " ".join(f"band_{b + 1}" for b in range(nb)) + "\n")
            data = np.column_stack([x, eig[s] - zero])
            np.savetxt(fh, data, fmt="%12.6f")
        written.append(fname)
        for j, (lab, w) in enumerate(proj_list):
            pname = f"{root}_proj{j + 1}{'' if nspin == 1 else '_' + ('up', 'down')[s]}{ext}"
            with open(pname, "w") as fh:
                fh.write(f"# Projection weights for '{lab}' (rows = k-points, columns = bands)\n")
                fh.write("k_dist(1/A) " + " ".join(f"w_band_{b + 1}" for b in range(nb)) + "\n")
                np.savetxt(fh, np.column_stack([x, w[s]]), fmt="%12.6f")
            written.append(pname)
    tname = f"{root}_ticks{ext}"
    with open(tname, "w") as fh:
        fh.write("# x_position label\n")
        for t, lab in ticks:
            fh.write(f"{t:12.6f} {ascii_label(lab) or '-'}\n")
    written.append(tname)
    print("Data written: " + ", ".join(written))


def setup_matplotlib(show: bool):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def make_plot(args, x, eig_rel, ticks, proj_list, zero_name, plt) -> None:
    nspin = eig_rel.shape[0]
    split = bool(proj_list) and nspin == 2
    npan = 2 if split else 1
    w, h = args.figsize
    fig, axes = plt.subplots(1, npan, figsize=(w * npan if split else w, h), sharey=True,
                             squeeze=False)
    axes = axes[0]
    ymin, ymax = args.ylim
    for p in range(npan):
        ax = axes[p]
        spins = [p] if split else list(range(nspin))
        for s in spins:
            if proj_list:
                color, lw, ls = "#9c9b97", 0.8, "-"
            elif nspin == 2:
                color, lw, ls = (SERIES[0], SERIES[1])[s], 1.2, ("-", "--")[s]
            else:
                color, lw, ls = INK, 1.2, "-"
            for b in range(eig_rel.shape[2]):
                lab = None
                if nspin == 2 and not proj_list and b == 0:
                    lab = ("spin up", "spin down")[s]
                ax.plot(x, eig_rel[s, :, b], color=color, lw=lw, ls=ls, label=lab, zorder=2)
            for j, (plab, wgt) in enumerate(proj_list):
                c = SERIES[j % len(SERIES)]
                sel = (eig_rel[s] > ymin - 1) & (eig_rel[s] < ymax + 1)
                xx = np.broadcast_to(x[:, None], eig_rel[s].shape)
                ax.scatter(xx[sel], eig_rel[s][sel], s=np.clip(wgt[s][sel], 0, None) * args.fat_scale,
                           color=c, alpha=0.55, edgecolors="none", zorder=3,
                           label=plab if (s == spins[0]) else None)
        for t, _ in ticks:
            ax.axvline(t, color=GRID, lw=0.8, zorder=1)
        if zero_name != "none":
            ax.axhline(0.0, color=MUTED, lw=0.8, ls="--", zorder=1)
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(ymin, ymax)
        ax.set_xticks([t for t, _ in ticks])
        ax.set_xticklabels([lab for _, lab in ticks])
        ax.tick_params(direction="in", length=3)
        if split:
            ax.set_title(("Spin up", "Spin down")[p], fontsize=10)
        if p == 0:
            ylab = {"efermi": "$E - E_\\mathrm{F}$ (eV)", "vbm": "$E - E_\\mathrm{VBM}$ (eV)",
                    "none": "Energy (eV)"}[zero_name]
            ax.set_ylabel(ylab)
        handles, labels = ax.get_legend_handles_labels()
        if labels and (p == 0):
            leg = ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.9,
                            markerscale=0.8 if proj_list else 1.0)
            if proj_list:
                for hnd in leg.legend_handles if hasattr(leg, "legend_handles") else leg.legendHandles:
                    try:
                        hnd.set_sizes([30])
                    except AttributeError:
                        pass
    if args.title:
        fig.suptitle(args.title)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)
    print(f"Figure written: {args.output}")
    if args.show:
        plt.show()


# --------------------------------------------------------------------------- main
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot a VASP band structure from vasprun.xml (incl. KPOINTS_OPT) or EIGENVAL, "
                    "report the band gap and optionally draw PROCAR fat bands.\n"
                    "从 vasprun.xml（含 KPOINTS_OPT）或 EIGENVAL 画能带，报告带隙/VBM/CBM，可画 PROCAR 投影能带。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples / 示例:\n"
               "  python plot_band.py -f vasprun.xml --efermi ../scf/OUTCAR\n"
               "  python plot_band.py -f EIGENVAL --kpoints KPOINTS --poscar POSCAR\n"
               "  python plot_band.py --opt                       # KPOINTS_OPT (VASP >= 6.3)\n"
               "  python plot_band.py --procar PROCAR --proj 'Mo:d' --proj 'S:p'\n"
               "  python plot_band.py --labels 'L G X U|K G' --zero vbm --dat band.dat\n")
    p.add_argument("-f", "--file", default=None,
                   help="vasprun.xml or EIGENVAL (default: ./vasprun.xml, else ./EIGENVAL) / 输入文件")
    p.add_argument("--kpoints", default=None,
                   help="KPOINTS file used for labels (default: KPOINTS or KPOINTS_OPT next to the input; 'none' disables) / k 点文件")
    p.add_argument("--kpath", default=None,
                   help="extra line-mode file with labelled vertices, e.g. vaspkit KPATH.in (for HSE zero-weight runs)")
    p.add_argument("--labels", default=None,
                   help="manual labels for the segment boundaries, e.g. 'L G X U|K G' / 手动高对称点标签")
    p.add_argument("--poscar", default=None,
                   help="POSCAR/CONTCAR for lattice and element names (EIGENVAL mode; default: CONTCAR/POSCAR next to input)")
    p.add_argument("--efermi", default=None,
                   help="Fermi energy in eV, or a file to read it from (OUTCAR, vasprun.xml, DOSCAR). "
                        "For ICHARG=11 runs pass the SCF OUTCAR / 费米能（数值或文件）")
    p.add_argument("--zero", choices=["efermi", "vbm", "none"], default="efermi",
                   help="energy zero: Fermi level (default), VBM, or raw eigenvalues / 能量零点")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--opt", action="store_true", help="force use of the KPOINTS_OPT block of vasprun.xml")
    g.add_argument("--no-opt", action="store_true", help="ignore the KPOINTS_OPT block even if present")
    p.add_argument("--keep-weighted", action="store_true",
                   help="keep weighted k-points even when zero-weight points exist (disable HSE filtering)")
    p.add_argument("--gap-method", choices=["auto", "occ", "efermi"], default="auto",
                   help="how occupied states are identified (auto: occupations if available, "
                        "else E-fermi) / 判断占据态的方法")
    p.add_argument("--gap-tol", type=float, default=0.05,
                   help="tolerance (eV) above E-fermi still counted as occupied (default 0.05)")
    p.add_argument("--nocc", type=int, default=None,
                   help="number of occupied bands per spin channel (insulators; overrides --gap-method)")
    p.add_argument("--procar", default=None, help="PROCAR or PROCAR_OPT for fat bands / 投影能带文件")
    p.add_argument("--proj", action="append", default=[],
                   help="projection 'ATOMS:ORBITALS', ATOMS = element|1-4,7|all, ORBITALS = s,p,d,f,t2g,eg,"
                        "px,dz2,...|all. Repeatable, e.g. --proj 'Ga:s' --proj 'As:p'")
    p.add_argument("--fat-scale", type=float, default=60.0, help="marker size for weight = 1 (default 60)")
    p.add_argument("--ylim", type=float, nargs=2, default=[-6.0, 6.0], metavar=("EMIN", "EMAX"),
                   help="energy window in eV relative to the zero (default -6 6)")
    p.add_argument("--figsize", type=float, nargs=2, default=[4.5, 5.0], metavar=("W", "H"))
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--title", default=None)
    p.add_argument("-o", "--output", default="band.png", help="figure file (png/pdf/svg; default band.png)")
    p.add_argument("--dat", default=None,
                   help="also export column data (Origin/gnuplot), e.g. band.dat / 导出数据")
    p.add_argument("--no-plot", action="store_true", help="only print the gap report")
    p.add_argument("--show", action="store_true", help="open an interactive window")
    return p


def main(argv=None) -> None:
    try:
        sys.stdout.reconfigure(errors="backslashreplace")  # type: ignore[attr-defined]
    except (AttributeError, ValueError):
        pass
    args = build_parser().parse_args(argv)

    path = args.file
    if path is None:
        path = "vasprun.xml" if os.path.isfile("vasprun.xml") else "EIGENVAL"
    check_file(path, "input file")
    folder = os.path.dirname(os.path.abspath(path))
    with open(path, "rb") as fh:
        is_xml = path.lower().endswith(".xml") or fh.read(200).lstrip().startswith(b"<?xml")

    lattice = None
    symbols = None
    efermi_src = None
    efermi = None
    nper = None
    gen_vertices: List = []
    use_opt = False
    occ = None

    if is_xml:
        d = parse_vasprun(path)
        lattice, symbols = d["lattice"], d["symbols"]
        has_opt = d["opt"] is not None and d["opt"].get("eig") is not None
        use_opt = has_opt and not args.no_opt
        if args.opt and not has_opt:
            die("--opt given but vasprun.xml has no <eigenvalues_kpoints_opt> block (VASP >= 6.3 + KPOINTS_OPT)")
        if use_opt:
            blk = d["opt"]
            print("Using the KPOINTS_OPT eigenvalues (<eigenvalues_kpoints_opt>) of vasprun.xml")
        else:
            blk = dict(d["main"] or {})
            blk["eig"], blk["occ"] = d["eig"], d["occ"]
            if d["eig"] is None:
                die("vasprun.xml only contains KPOINTS_OPT eigenvalues; rerun without --no-opt")
        kfrac, weights, eig, occ = blk.get("kpoints"), blk.get("weights"), blk["eig"], blk.get("occ")
        if kfrac is None:
            die("k-point list not found in vasprun.xml")
        if blk.get("generation") and blk["generation"].lower().startswith("list"):
            nper = blk.get("divisions")
            gen_vertices = blk.get("vertices", [])
        if d["efermi"] is not None:
            efermi, efermi_src = d["efermi"], "vasprun.xml <dos>"
        if occ is not None and np.all(np.isnan(occ)):
            occ = None
    else:
        e = parse_eigenval(path)
        kfrac, weights, eig, occ = e["kpoints"], e["weights"], e["eig"], e["occ"]
        if np.all(np.isnan(occ)):
            occ = None
        pos = args.poscar
        if pos is None:
            for cand in ("CONTCAR", "POSCAR"):
                c = os.path.join(folder, cand)
                if os.path.isfile(c) and os.path.getsize(c) > 0:
                    pos = c
                    break
        if pos:
            lattice, symbols = read_poscar(pos)
        out = os.path.join(folder, "OUTCAR")
        if os.path.isfile(out):
            efermi = efermi_from_file(out)
            efermi_src = "OUTCAR next to EIGENVAL"
    if args.poscar and is_xml:
        lattice, symbols = read_poscar(args.poscar)

    if eig.shape[1] != len(kfrac):
        die(f"number of k-points in eigenvalues ({eig.shape[1]}) != k-point list ({len(kfrac)})")

    # Fermi level from the command line overrides everything.
    if args.efermi is not None:
        efermi, efermi_src = resolve_efermi(args.efermi)
    elif efermi is not None:
        print(f"E-fermi = {efermi:.4f} eV taken from {efermi_src}. For a non-self-consistent band run "
              "(ICHARG = 11) pass the SCF value with --efermi.")
    else:
        warn("Fermi energy unknown; use --efermi (energies are plotted unshifted unless --zero vbm works)")

    # Zero-weight (HSE) filtering.
    full_nk = len(kfrac)
    kmask = np.ones(full_nk, dtype=bool)
    if weights is not None and not args.keep_weighted and not use_opt:
        w = np.asarray(weights)
        if np.any(np.abs(w) < 1e-12) and np.any(np.abs(w) > 1e-12):
            kmask = np.abs(w) < 1e-12
            print(f"Zero-weight k-points detected: keeping {kmask.sum()} zero-weight points, "
                  f"dropping {full_nk - kmask.sum()} weighted SCF points (hybrid-functional band run).")
            nper = None
    kfrac = kfrac[kmask]
    eig = eig[:, kmask, :]
    if occ is not None:
        occ = occ[:, kmask, :]
    if len(kfrac) < 2:
        die("fewer than 2 k-points to plot")

    # Labels.
    kfile = args.kpoints
    if kfile is None:
        names = ["KPOINTS_OPT", "KPOINTS"] if use_opt else ["KPOINTS"]
        for n in names:
            c = os.path.join(folder, n)
            if os.path.isfile(c):
                kfile = c
                break
    if kfile is not None and kfile.lower() == "none":
        kfile = None
    kinfo = read_kpoints_file(kfile) if kfile else None
    if kfile and kinfo is None:
        warn(f"could not interpret {kfile} as a KPOINTS file; no labels from it")
    if kinfo and kinfo["mode"] == "line" and nper is None and kmask.all():
        if len(kfrac) % kinfo["N"] == 0:
            nper = kinfo["N"]
    if kinfo and kinfo["mode"] == "line" and nper and kinfo["N"] != nper:
        warn(f"{kfile} has {kinfo['N']} points per segment but the calculation used {nper}; labels may be wrong")
    vertices = vertices_from_kfile(kinfo)
    if args.kpath:
        vertices += vertices_from_kfile(read_kpoints_file(args.kpath))
    if nper and kinfo is None and gen_vertices:
        pseudo = [((v, None)) for v in gen_vertices]
        kinfo = {"mode": "line", "N": nper, "cart": False,
                 "segments": [(pseudo[2 * i], pseudo[2 * i + 1]) for i in range(len(pseudo) // 2)]}
    x, ticks, labels_at = build_kpath(kfrac, lattice, kinfo, nper, vertices, args.labels)

    # Gap.
    method = args.gap_method
    if method == "auto":
        method = "occ" if (occ is not None and args.efermi is None) else "efermi"
    if method == "occ" and (occ is None or not np.any(occ > 0.5)):
        if args.gap_method == "occ":
            warn("no usable occupations available; falling back to the E-fermi criterion")
        method = "efermi"
    rep = gap_report(eig, occ, kfrac, efermi, method, args.nocc, args.gap_tol, labels_at)
    nspin = eig.shape[0]
    print(f"Input: {path}   spin channels: {nspin}   k-points plotted: {len(kfrac)}   bands: {eig.shape[2]}")
    if efermi is not None:
        print(f"E-fermi = {efermi:.4f} eV ({efermi_src})")
    print_gap(rep, kfrac, labels_at, efermi)

    # Energy zero.
    zero_name = args.zero
    zero = 0.0
    if zero_name == "efermi":
        if efermi is None:
            warn("no Fermi energy: plotting raw eigenvalues")
            zero_name = "none"
        else:
            zero = efermi
    elif zero_name == "vbm":
        ov = rep.get("overall", {}) if rep.get("ok") else {}
        if ov.get("vbm"):
            zero = ov["vbm"][0]
        else:
            warn("VBM not defined (metal or unknown E-fermi); using E-fermi as zero")
            zero_name = "efermi" if efermi is not None else "none"
            zero = efermi or 0.0
    eig_rel = eig - zero

    # Projections.
    proj_list = []
    if args.proj and not args.procar:
        cand = os.path.join(folder, "PROCAR_OPT" if use_opt else "PROCAR")
        if os.path.isfile(cand):
            args.procar = cand
        else:
            die("--proj needs --procar PROCAR (LORBIT = 11)")
    if args.procar:
        pr = read_procar(args.procar)
        P = pr["proj"]  # (nspin, nk, nb, nions, norb)
        if P.shape[0] != nspin:
            die(f"PROCAR has {P.shape[0]} spin block(s) but eigenvalues have {nspin}")
        if P.shape[1] == full_nk:
            P = P[:, kmask]
        elif P.shape[1] != len(kfrac):
            die(f"PROCAR has {P.shape[1]} k-points, eigenvalues have {full_nk}; use matching files "
                "(PROCAR_OPT for KPOINTS_OPT runs)")
        nb_use = min(P.shape[2], eig.shape[2])
        if P.shape[2] != eig.shape[2]:
            warn("PROCAR and eigenvalue files have different NBANDS; using the common bands")
            eig_rel = eig_rel[:, :, :nb_use]
        if not args.proj:
            args.proj = ["all:all"]
        for spec in args.proj:
            atoms_s, _, orbs_s = spec.partition(":")
            orbs_s = orbs_s or "all"
            ia = atom_indices(atoms_s, symbols, P.shape[3])
            io: List[int] = []
            for t in orbs_s.split(","):
                io += orbital_indices(t, pr["orbitals"])
            io = sorted(set(io))
            wgt = P[:, :, :nb_use][:, :, :, ia][..., io].sum(axis=(3, 4))
            proj_list.append((spec, wgt))
            print(f"Projection '{spec}': {len(ia)} atom(s), orbitals {[pr['orbitals'][i] for i in io]}, "
                  f"max weight {float(wgt.max()):.3f}")

    if args.dat:
        write_dat(args.dat, x, eig_rel + zero, ticks, zero, zero_name, proj_list)
    if not args.no_plot:
        plt = setup_matplotlib(args.show)
        make_plot(args, x, eig_rel, ticks, proj_list, zero_name, plt)


if __name__ == "__main__":
    main()
