#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
md_analysis.py -- analyse VASP molecular-dynamics runs (AIMD / MLFF).
                  分析 VASP 分子动力学：温度/能量/体积曲线、RDF、MSD 与扩散系数、ML_LOGFILE 误差。

Sub-commands
------------
  thermo : temperature, energies (OSZICAR) and cell volume (XDATCAR) vs. time,
           averages after equilibration and drift of the conserved energy E.
  rdf    : radial distribution functions g(r) for element pairs + running
           coordination number (XDATCAR, constant or variable cell).
  msd    : mean-square displacement per element (unwrapped, centre-of-mass drift
           removed, all time origins via FFT) and diffusion coefficient
           D = slope / (2 d) (Einstein relation).
  mlff   : parse ML_LOGFILE (VASP >= 6.3 machine-learned force fields): plot the
           ERR (training RMSE) and BEEF (Bayesian error / threshold) lines and
           count the STATUS states (accurate / threshold / learning / critical ...).
  all    : run every analysis whose input files exist in the directory.

Time step: POTIM is read from --potim, else OUTCAR, else INCAR (default 1 fs).
XDATCAR frames are written every NBLOCK steps; the frame spacing is derived from
the "Direct configuration=" numbers and NBLOCK (override with --dt).

Usage examples
--------------
  python md_analysis.py thermo                          # ./OSZICAR (+ ./XDATCAR)
  python md_analysis.py rdf --pair Si-Si --skip-ps 2 --rmax 6
  python md_analysis.py rdf --pair O-H --pair O-O -d water_npt/
  python md_analysis.py msd --skip-ps 2 --fit-range 1 5 --species Li
  python md_analysis.py mlff -f ML_LOGFILE
  python md_analysis.py all -d run1/ --skip-ps 1

Dependencies
------------
  Python >= 3.9, numpy, matplotlib; optional: ase (atomic masses for the COM correction).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

from __future__ import annotations

import argparse
import itertools
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK = "#1f1f1f"
MUTED = "#8a8986"
KV_RE = re.compile(r"([A-Za-z0-9]+)=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?)")


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


def element(tok: str) -> str:
    m = re.match(r"[A-Za-z]+", tok)
    s = m.group(0).capitalize() if m else tok
    m2 = re.match(r"[A-Z][a-z]?", s)
    return m2.group(0) if m2 else s


# --------------------------------------------------------------------------- run parameters
def read_run_params(folder: str, potim: Optional[float]) -> Dict:
    out = {"potim": None, "nblock": None, "src": None}
    oc = os.path.join(folder, "OUTCAR")
    if os.path.isfile(oc):
        with open(oc, "r", errors="replace") as fh:
            for i, line in enumerate(fh):
                if out["potim"] is None and re.match(r"^\s*POTIM\s*=", line):
                    try:
                        out["potim"] = float(line.split("=")[1].split()[0])
                        out["src"] = "OUTCAR"
                    except (ValueError, IndexError):
                        pass
                if out["nblock"] is None and "NBLOCK" in line and "=" in line:
                    m = re.search(r"NBLOCK\s*=\s*(\d+)", line)
                    if m:
                        out["nblock"] = int(m.group(1))
                if (out["potim"] is not None and out["nblock"] is not None) or i > 20000:
                    break
    inc = os.path.join(folder, "INCAR")
    if os.path.isfile(inc):
        with open(inc, "r", errors="replace") as fh:
            for line in fh:
                body = re.split(r"[#!]", line)[0]
                for part in body.split(";"):
                    m = re.match(r"\s*(\w+)\s*=\s*(\S+)", part)
                    if not m:
                        continue
                    k, v = m.group(1).upper(), m.group(2)
                    try:
                        if k == "POTIM" and out["potim"] is None:
                            out["potim"], out["src"] = float(v), "INCAR"
                        elif k == "NBLOCK" and out["nblock"] is None:
                            out["nblock"] = int(float(v))
                    except ValueError:
                        pass
    if potim is not None:
        out["potim"], out["src"] = potim, "command line"
    if out["potim"] is None:
        out["potim"], out["src"] = 1.0, "default (no OUTCAR/INCAR found)"
        warn("POTIM unknown: assuming 1 fs (use --potim)")
    if out["nblock"] is None:
        out["nblock"] = 1
    return out


# --------------------------------------------------------------------------- XDATCAR
def read_xdatcar(path: str, skip_frames: int = 0, stride: int = 1) -> Dict:
    """Read XDATCAR (constant or variable cell). Positions are unwrapped on the fly
    using every frame; only frames idx >= skip_frames with (idx-skip) % stride == 0 are stored."""
    check_file(path, "XDATCAR")
    frac_l, cart_l, lat_l, step_l = [], [], [], []
    with open(path, "r", errors="replace") as fh:

        def read_header(first_line: str):
            scale = [float(x) for x in fh.readline().split()[:3]]
            lat = np.array([[float(x) for x in fh.readline().split()[:3]] for _ in range(3)])
            if len(scale) == 3:
                lat = lat * np.array(scale)[None, :]
            elif scale[0] < 0:
                lat = lat * (abs(scale[0]) / abs(np.linalg.det(lat))) ** (1 / 3)
            else:
                lat = lat * scale[0]
            line = fh.readline()
            toks = line.split()
            if toks and not re.match(r"^\d+$", toks[0]):
                names = [element(t) for t in toks]
                counts = [int(x) for x in fh.readline().split()[:len(names)]]
            else:
                counts = [int(x) for x in toks]
                names = [element(t) for t in first_line.split()[:len(counts)]]
                if len(names) != len(counts):
                    names = [f"X{i + 1}" for i in range(len(counts))]
            return lat, names, counts

        try:
            first = fh.readline()
            lat, names, counts = read_header(first)
        except (ValueError, IndexError) as exc:
            die(f"{path}: cannot parse the header ({exc})")
        nat = sum(counts)
        symbols = sum(([n] * c for n, c in zip(names, counts)), [])
        prev = None
        unwrapped = None
        idx = 0
        variable = False
        while True:
            line = fh.readline()
            if not line:
                break
            s = line.strip()
            if not s:
                continue
            low = s.lower()
            if "configuration" in low or low.startswith(("direct", "cartesian")):
                nums = re.findall(r"\d+", s)
                step = int(nums[-1]) if nums else idx + 1
                rows = []
                for _ in range(nat):
                    ln = fh.readline()
                    if not ln:
                        break
                    rows.append(ln.split()[:3])
                if len(rows) < nat:
                    warn(f"{path}: last frame (configuration {step}) is incomplete and was ignored")
                    break
                try:
                    fr = np.array(rows, dtype=float)
                except ValueError:
                    warn(f"{path}: malformed coordinates in configuration {step}; stopping here")
                    break
                if low.startswith("cartesian"):
                    fr = fr @ np.linalg.inv(lat)
                if prev is None:
                    unwrapped = fr @ lat
                else:
                    d = fr - prev
                    d -= np.round(d)
                    unwrapped = unwrapped + d @ lat
                prev = fr
                if idx >= skip_frames and (idx - skip_frames) % stride == 0:
                    frac_l.append(fr % 1.0)
                    cart_l.append(unwrapped.copy())
                    lat_l.append(lat.copy())
                    step_l.append(step)
                idx += 1
            else:
                # new header block -> variable-cell XDATCAR (ISIF = 3 MD)
                try:
                    lat, nm, ct = read_header(line)
                except (ValueError, IndexError):
                    warn(f"{path}: unexpected line '{s[:40]}' - stopping here")
                    break
                variable = True
                if sum(ct) != nat:
                    die(f"{path}: number of atoms changes between frames")
    if not frac_l:
        die(f"{path}: no complete frames after skipping {skip_frames}")
    return {"symbols": symbols, "frac": np.array(frac_l), "cart": np.array(cart_l), "lat": np.array(lat_l),
            "steps": np.array(step_l), "variable": variable, "nframes_total": idx}


def frame_dt(xd: Dict, prm: Dict, stride: int, dt_override: Optional[float]) -> float:
    """Time between two stored frames in fs."""
    if dt_override:
        return dt_override * stride
    steps = xd["steps"]
    inc = int(np.median(np.diff(steps))) // stride if len(steps) > 1 else 1
    inc = max(inc, 1)
    nb = prm["nblock"]
    if inc == 1 and nb > 1:
        print(f"Note: configurations are numbered consecutively; using NBLOCK = {nb} steps per frame")
        inc = nb
    return prm["potim"] * inc * stride


def perp_widths(lat: np.ndarray) -> np.ndarray:
    vol = abs(np.linalg.det(lat))
    return np.array([vol / np.linalg.norm(np.cross(lat[(i + 1) % 3], lat[(i + 2) % 3])) for i in range(3)])


def skip_to_frames(args, dt_raw: float) -> int:
    if args.skip_ps is not None:
        return int(round(args.skip_ps * 1000.0 / dt_raw))
    return args.skip


# --------------------------------------------------------------------------- OSZICAR
def read_oszicar(path: str) -> Dict[str, np.ndarray]:
    check_file(path, "OSZICAR")
    rows: List[Dict[str, float]] = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if "T=" not in line or "E=" not in line:
                continue
            t = line.split()
            if not t or not t[0].isdigit():
                continue
            d = {"step": float(t[0])}
            for k, v in KV_RE.findall(line):
                try:
                    d[k] = float(v)
                except ValueError:
                    pass
            rows.append(d)
    if not rows:
        die(f"{path}: no MD lines ('N T= ... E= ...') found - is this an IBRION = 0 run?")
    keys = sorted({k for r in rows for k in r})
    return {k: np.array([r.get(k, np.nan) for r in rows]) for k in keys}


# --------------------------------------------------------------------------- ML_LOGFILE
ML_DEFAULT_COLS = {
    "ERR": ["nstep", "rmse_energy", "rmse_force", "rmse_stress"],
    "BEEF": ["nstep", "bee_energy", "bee_max_force", "bee_ave_force", "threshold", "bee_max_stress",
             "bee_ave_stress"],
    "STATUS": ["nstep", "state", "is", "doabin", "iff", "nsample", "ngenff"],
}


def read_ml_logfile(path: str) -> Dict:
    check_file(path, "ML_LOGFILE")
    cols: Dict[str, List[str]] = {}
    data: Dict[str, List[List[str]]] = {}
    word = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                t = s[1:].split()
                if len(t) > 2 and t[1] == "nstep" and all(word.match(x) for x in t[1:]):
                    cols[t[0]] = t[1:]
                continue
            t = s.split()
            if t[0] in ("ERR", "BEEF", "STATUS", "THRUPD", "THRHIST", "LCONF", "SPRSC", "REGR", "REGRF"):
                data.setdefault(t[0], []).append(t[1:])
    for k, v in ML_DEFAULT_COLS.items():
        cols.setdefault(k, v)
    out: Dict = {"cols": cols, "raw": data}
    for tag in ("ERR", "BEEF"):
        rows = data.get(tag, [])
        if rows:
            n = min(len(r) for r in rows)
            arr = np.array([[float(x) for x in r[:n]] for r in rows])
            out[tag] = {name: arr[:, i] for i, name in enumerate(cols[tag][:n])}
    return out


# --------------------------------------------------------------------------- analysis helpers
def species_masses(symbols: List[str]) -> Tuple[np.ndarray, bool]:
    try:
        from ase.data import atomic_masses, atomic_numbers  # optional
        return np.array([atomic_masses[atomic_numbers[s]] for s in symbols]), True
    except (ImportError, KeyError):
        return np.ones(len(symbols)), False


def msd_fft(x: np.ndarray) -> np.ndarray:
    """MSD for all lags, all time origins. x: (nframes, natoms, 3). Returns (nframes, natoms)."""
    n = x.shape[0]
    d = np.square(x).sum(axis=2)  # (n, N)
    d = np.vstack([d, np.zeros((1, d.shape[1]))])
    q = 2 * d[:n].sum(axis=0)
    s1 = np.zeros((n, x.shape[1]))
    for m in range(n):
        q = q - d[m - 1] - d[n - m]
        s1[m] = q / (n - m)
    f = np.fft.rfft(x, n=2 * n, axis=0)
    ac = np.fft.irfft(f * np.conj(f), axis=0)[:n].sum(axis=2)  # (n, N)
    s2 = ac / (n - np.arange(n))[:, None]
    return s1 - 2 * s2


def rdf_pair(frac: np.ndarray, lats: np.ndarray, ia: np.ndarray, ib: np.ndarray, rmax: float, dr: float,
             same: bool) -> Tuple[np.ndarray, np.ndarray, float]:
    edges = np.arange(0.0, rmax + dr / 2, dr)
    hist = np.zeros(len(edges) - 1)
    rho_sum = 0.0
    shifts_f = np.array(list(itertools.product((-1, 0, 1), repeat=3)), dtype=float)
    nf = frac.shape[0]
    for f in range(nf):
        lat = lats[f]
        vol = abs(np.linalg.det(lat))
        w = perp_widths(lat)
        ortho = np.allclose(lat - np.diag(np.diag(lat)), 0.0, atol=1e-6)
        sh = np.zeros((1, 3)) if (rmax <= 0.5 * w.min() and ortho) else shifts_f @ lat
        A = frac[f, ia]
        B = frac[f, ib]
        nb_eff = len(ib) - 1 if same else len(ib)
        rho_b = nb_eff / vol
        rho_sum += rho_b
        h = np.zeros_like(hist)
        chunk = max(1, int(2_000_000 // max(1, len(ib) * len(sh))))
        for c0 in range(0, len(A), chunk):
            d = B[None, :, :] - A[c0:c0 + chunk, None, :]
            d -= np.round(d)
            dc = d @ lat
            r = np.linalg.norm(dc[:, :, None, :] + sh[None, None, :, :], axis=-1).ravel()
            r = r[(r > 1e-6) & (r < rmax)]
            h += np.histogram(r, bins=edges)[0]
        shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
        hist += h / (len(ia) * rho_b * shell)
    g = hist / nf
    r = 0.5 * (edges[1:] + edges[:-1])
    return r, g, rho_sum / nf


def first_shell(r: np.ndarray, g: np.ndarray, cn: np.ndarray) -> Optional[Tuple[float, float, float, float]]:
    """First peak and first minimum of a (lightly smoothed) g(r); CN at the minimum."""
    if len(g) < 9 or g.max() <= 0:
        return None
    w = 5
    gs = np.convolve(g, np.ones(w) / w, mode="same")
    above = np.where(gs > max(1.0, 0.5 * gs.max()))[0]
    if above.size == 0:
        return None
    k = int(above[0])
    while k + 1 < len(gs) and gs[k + 1] >= gs[k]:
        k += 1
    j = k
    while j + 1 < len(gs) and gs[j + 1] <= gs[j]:
        j += 1
    if j == k or j >= len(gs) - 2:
        return r[k], gs[k], float("nan"), float("nan")
    return r[k], gs[k], r[j], cn[j]


def setup_matplotlib(show: bool):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def save(fig, path: str, args, plt) -> None:
    fig.tight_layout()
    fig.savefig(path, dpi=args.dpi)
    print(f"Figure written: {path}")
    if args.show:
        plt.show()


# --------------------------------------------------------------------------- commands
def cmd_thermo(args, folder: str, plt) -> None:
    prm = read_run_params(folder, args.potim)
    osz = read_oszicar(os.path.join(folder, args.oszicar))
    t = osz["step"] * prm["potim"] / 1000.0  # ps
    print(f"OSZICAR: {len(t)} MD steps, POTIM = {prm['potim']} fs ({prm['src']}), total {t[-1]:.3f} ps")
    t0 = args.skip_ps if args.skip_ps is not None else args.skip * prm["potim"] / 1000.0
    m = t >= t0
    if m.sum() < 2:
        die(f"nothing left after skipping {t0} ps")
    print(f"Averages over t >= {t0:.3f} ps ({m.sum()} steps):")
    for key, lab, unit in (("T", "Temperature", "K"), ("E", "Total (conserved) E", "eV"),
                           ("F", "Free energy F", "eV"), ("E0", "E0 (sigma->0)", "eV"),
                           ("EK", "Ionic kinetic EK", "eV"), ("mag", "Magnetisation", "muB")):
        if key in osz and np.isfinite(osz[key][m]).any():
            v = osz[key][m]
            print(f"  {lab:22s}: {np.nanmean(v):14.4f} +- {np.nanstd(v):.4f} {unit}")
    if "E" in osz and m.sum() > 10:
        slope = np.polyfit(t[m], osz["E"][m], 1)[0]
        print(f"  Drift of E (conserved quantity): {slope * 1000:.2f} meV/ps "
              "(large drift -> smaller POTIM or tighter EDIFF)")
    vol = None
    xpath = os.path.join(folder, args.xdatcar)
    if os.path.isfile(xpath):
        xd = read_xdatcar(xpath)
        dt = frame_dt(xd, prm, 1, args.dt)
        tv = np.arange(1, len(xd["steps"]) + 1) * dt / 1000.0  # ps
        vols = np.abs(np.linalg.det(xd["lat"]))
        if xd["variable"]:
            mv = tv >= t0
            abc = np.linalg.norm(xd["lat"][mv], axis=2).mean(axis=0)
            print(f"  Volume (XDATCAR, variable cell): {vols[mv].mean():.3f} +- {vols[mv].std():.3f} A^3; "
                  f"<a,b,c> = {abc[0]:.4f} {abc[1]:.4f} {abc[2]:.4f} A")
            vol = (tv, vols)
        else:
            print(f"  Volume (constant cell): {vols[0]:.3f} A^3")
    if args.no_plot:
        return
    panels = [("T", "T (K)")] + [(k, f"{k} (eV)") for k in ("F", "E") if k in osz]
    n = len(panels) + (1 if vol else 0)
    fig, axs = plt.subplots(n, 1, figsize=(args.figsize[0], 1.9 * n + 0.4), sharex=True, squeeze=False)
    axs = axs[:, 0]
    for ax, (k, lab) in zip(axs, panels):
        ax.plot(t, osz[k], color=SERIES[0], lw=0.9)
        ax.axhline(np.nanmean(osz[k][m]), color=MUTED, lw=0.8, ls="--")
        ax.set_ylabel(lab)
    if vol:
        axs[-1].plot(vol[0], vol[1], color=SERIES[0], lw=0.9)
        axs[-1].set_ylabel("V (Å$^3$)")
    for ax in axs:
        if t0 > 0:
            ax.axvspan(t[0], t0, color="#eeeeec", lw=0)
        ax.tick_params(direction="in", length=3)
    axs[-1].set_xlabel("Time (ps)")
    save(fig, os.path.join(folder, args.prefix + "thermo.png"), args, plt)


def load_traj(args, folder: str):
    prm = read_run_params(folder, args.potim)
    xpath = os.path.join(folder, args.xdatcar)
    # first pass with stride 1 just to know the frame spacing cheaply is avoided: estimate from header numbers
    xd_full = read_xdatcar(xpath, 0, 1) if (args.skip_ps is not None) else None
    if xd_full is not None:
        dt_raw = frame_dt(xd_full, prm, 1, args.dt)
        skip = skip_to_frames(args, dt_raw)
        sel = slice(skip, None, args.stride)
        xd = {k: (v[sel] if isinstance(v, np.ndarray) and v.ndim >= 1 and len(v) == len(xd_full["steps"]) else v)
              for k, v in xd_full.items()}
        if len(xd["steps"]) == 0:
            die(f"no frames left after skipping {args.skip_ps} ps")
    else:
        xd = read_xdatcar(xpath, args.skip, args.stride)
        dt_raw = frame_dt(xd, prm, args.stride, args.dt) / args.stride
        skip = args.skip
    dt = dt_raw * args.stride
    print(f"XDATCAR: {xd['nframes_total']} frames, {len(xd['symbols'])} atoms "
          f"({', '.join(f'{s}{xd['symbols'].count(s)}' for s in dict.fromkeys(xd['symbols']))}), "
          f"{'variable' if xd['variable'] else 'constant'} cell; POTIM = {prm['potim']} fs ({prm['src']}); "
          f"frame spacing {dt_raw:.3f} fs; using {len(xd['steps'])} frames after skipping {skip} (stride {args.stride})")
    return xd, dt


def parse_pairs(pairs: List[str], symbols: List[str]) -> List[Tuple[str, str]]:
    uniq = list(dict.fromkeys(symbols))
    if not pairs:
        if len(uniq) > 3:
            return [("all", "all")]
        return [(a, b) for i, a in enumerate(uniq) for b in uniq[i:]]
    out = []
    for p in pairs:
        if p.lower() == "all":
            out.append(("all", "all"))
            continue
        a, sep, b = p.partition("-")
        if not sep:
            die(f"pair '{p}' must look like A-B (e.g. Si-Si, O-H) or 'all'")
        a, b = element(a), element(b)
        for x in (a, b):
            if x not in uniq:
                die(f"element '{x}' not in XDATCAR ({uniq})")
        out.append((a, b))
    return out


def cmd_rdf(args, folder: str, plt) -> None:
    xd, dt = load_traj(args, folder)
    sym = np.array(xd["symbols"])
    wmin = min(perp_widths(l).min() for l in xd["lat"])
    rmax = args.rmax if args.rmax else min(0.5 * wmin, 10.0)
    if rmax > wmin:
        warn(f"--rmax {rmax} exceeds the smallest cell width {wmin:.2f} A; truncated")
        rmax = wmin
    elif rmax > 0.5 * wmin:
        print(f"Note: rmax {rmax:.2f} A > half cell width {0.5 * wmin:.2f} A: periodic images are included")
    results = []
    for a, b in parse_pairs(args.pair, xd["symbols"]):
        ia = np.arange(len(sym)) if a == "all" else np.where(sym == a)[0]
        ib = np.arange(len(sym)) if b == "all" else np.where(sym == b)[0]
        same = (a == b)
        r, g, rho = rdf_pair(xd["frac"], xd["lat"], ia, ib, rmax, args.dr, same)
        cn = np.cumsum(4 * np.pi * r ** 2 * rho * g * args.dr)
        fs = first_shell(r, g, cn)
        lab = f"{a}-{b}"
        if fs:
            print(f"g({lab}): first peak r = {fs[0]:.3f} A (g = {fs[1]:.2f}); first minimum r = {fs[2]:.3f} A; "
                  f"coordination number up to the minimum = {fs[3]:.2f}")
        results.append((lab, r, g, cn))
    if args.dat:
        cols = [results[0][1]] + [x for _, _, g, cn in results for x in (g, cn)]
        names = ["r(A)"] + [n for lab, *_ in results for n in (f"g_{lab}", f"CN_{lab}")]
        with open(os.path.join(folder, args.dat) if not os.path.isabs(args.dat) else args.dat, "w") as fh:
            fh.write(f"# RDF from {args.xdatcar}, {len(xd['steps'])} frames\n" + " ".join(names) + "\n")
            np.savetxt(fh, np.column_stack(cols), fmt="%12.5f")
        print(f"Data written: {args.dat}")
    if args.no_plot:
        return
    fig, axs = plt.subplots(2, 1, figsize=(args.figsize[0], args.figsize[1] * 1.3), sharex=True,
                            gridspec_kw={"height_ratios": [2, 1]})
    for k, (lab, r, g, cn) in enumerate(results):
        c = SERIES[k % len(SERIES)]
        axs[0].plot(r, g, color=c, lw=1.4, label=lab)
        axs[1].plot(r, cn, color=c, lw=1.4, label=lab)
    axs[0].axhline(1.0, color=MUTED, lw=0.8, ls="--")
    axs[0].set_ylabel("g(r)")
    axs[1].set_ylabel("n(r)")
    axs[1].set_xlabel("r (Å)")
    axs[1].set_xlim(0, rmax)
    for ax in axs:
        ax.tick_params(direction="in", length=3)
    if len(results) > 1:
        axs[0].legend(fontsize=8)
    save(fig, os.path.join(folder, args.prefix + "rdf.png"), args, plt)


def cmd_msd(args, folder: str, plt) -> None:
    xd, dt = load_traj(args, folder)
    sym = np.array(xd["symbols"])
    pos = xd["cart"].copy()  # (nf, N, 3) unwrapped
    nf = pos.shape[0]
    if nf < 10:
        die(f"only {nf} frames: MSD needs a longer trajectory")
    if not args.keep_com:
        mass, real = species_masses(xd["symbols"])
        if not real:
            warn("ase not available: centre-of-mass drift removed with equal masses")
        com = (pos * mass[None, :, None]).sum(axis=1) / mass.sum()
        pos = pos - (com - com[0])[:, None, :]
    t = np.arange(nf) * dt / 1000.0  # ps
    total = t[-1]
    if total < 5.0:
        warn(f"trajectory after skipping is only {total:.2f} ps: diffusion coefficients will be rough")
    lo, hi = args.fit_range if args.fit_range else (0.1 * total, 0.5 * total)
    fit = (t >= lo) & (t <= hi)
    if fit.sum() < 3:
        die(f"fit range {lo}-{hi} ps contains fewer than 3 points")
    species = args.species or list(dict.fromkeys(xd["symbols"]))
    curves = []
    print(f"MSD fit window: {lo:.3f}-{hi:.3f} ps; D = slope / (2 x {args.dim})")
    for s in species:
        idx = np.where(sym == element(s))[0]
        if idx.size == 0:
            die(f"species '{s}' not in XDATCAR")
        msd = msd_fft(pos[:, idx, :]).mean(axis=1)
        a, b = np.polyfit(t[fit], msd[fit], 1)
        dval = a / (2 * args.dim)
        print(f"  {s:3s}: slope {a:.4f} A^2/ps -> D = {dval:.4e} A^2/ps = {dval * 1e-4:.3e} cm^2/s "
              f"= {dval * 1e-8:.3e} m^2/s")
        curves.append((s, msd, a, b))
    if args.dat:
        with open(os.path.join(folder, args.dat) if not os.path.isabs(args.dat) else args.dat, "w") as fh:
            fh.write(f"# MSD (A^2) vs lag time, all time origins, {nf} frames\n")
            fh.write("t(ps) " + " ".join(f"MSD_{s}" for s, *_ in curves) + "\n")
            np.savetxt(fh, np.column_stack([t] + [c[1] for c in curves]), fmt="%14.6f")
        print(f"Data written: {args.dat}")
    if args.no_plot:
        return
    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    nshow = max(3, int(0.8 * nf))  # the last lags have too few time origins
    for k, (s, msd, a, b) in enumerate(curves):
        c = SERIES[k % len(SERIES)]
        ax.plot(t[:nshow], msd[:nshow], color=c, lw=1.4, label=f"{s}  D = {a / (2 * args.dim) * 1e-4:.2e} cm$^2$/s")
        ax.plot(t[fit], a * t[fit] + b, color=INK, lw=0.9, ls="--")
    ax.axvspan(lo, hi, color="#eeeeec", lw=0, zorder=0)
    ax.set_xlabel("Lag time (ps)")
    ax.set_ylabel("MSD (Å$^2$)")
    ax.tick_params(direction="in", length=3)
    ax.legend(fontsize=8)
    save(fig, os.path.join(folder, args.prefix + "msd.png"), args, plt)


def cmd_mlff(args, folder: str, plt) -> None:
    path = args.file if args.file else os.path.join(folder, "ML_LOGFILE")
    ml = read_ml_logfile(path)
    raw = ml["raw"]
    print(f"ML_LOGFILE: {path}")
    if "STATUS" in raw:
        states: Dict[str, int] = {}
        for r in raw["STATUS"]:
            if len(r) > 1:
                states[r[1]] = states.get(r[1], 0) + 1
        print("  STATUS counts: " + ", ".join(f"{k} {v}" for k, v in sorted(states.items(), key=lambda x: -x[1])))
        n_abinitio = states.get("threshold", 0) + states.get("learning", 0) + states.get("critical", 0)
        if sum(states.values()):
            print(f"  steps with ab-initio sampling/learning: {n_abinitio} of {sum(states.values())} "
                  f"({100.0 * n_abinitio / sum(states.values()):.1f} %)")
    for tag in ("ERR", "BEEF"):
        if tag in ml:
            d = ml[tag]
            last = {k: v[-1] for k, v in d.items()}
            print(f"  last {tag}: " + ", ".join(f"{k} = {v:.4g}" for k, v in last.items()))
        else:
            print(f"  no {tag} lines")
    if args.no_plot or ("ERR" not in ml and "BEEF" not in ml):
        return
    panels = []
    if "BEEF" in ml:
        panels.append("beef")
    if "ERR" in ml:
        panels += ["err_e", "err_f", "err_s"]
    fig, axs = plt.subplots(len(panels), 1, figsize=(args.figsize[0], 1.8 * len(panels) + 0.4), sharex=True,
                            squeeze=False)
    axs = axs[:, 0]
    for ax, p in zip(axs, panels):
        if p == "beef":
            d = ml["BEEF"]
            ax.plot(d["nstep"], d.get("bee_max_force", d[list(d)[2]]), color=SERIES[0], lw=0.9,
                    label="Bayesian max force error")
            if "threshold" in d:
                ax.plot(d["nstep"], d["threshold"], color=SERIES[1], lw=1.2, ls="--", label="threshold (ML_CTIFOR)")
            ax.set_ylabel("eV/Å")
            ax.set_yscale("log")
            ax.legend(fontsize=8)
        else:
            d = ml["ERR"]
            key, lab, scale = {"err_e": ("rmse_energy", "RMSE E (meV/atom)", 1000.0),
                               "err_f": ("rmse_force", "RMSE F (eV/Å)", 1.0),
                               "err_s": ("rmse_stress", "RMSE stress (kB)", 1.0)}[p]
            if key in d:
                ax.plot(d["nstep"], d[key] * scale, color=SERIES[0], lw=1.2, marker="o", ms=2.5)
            ax.set_ylabel(lab)
        ax.tick_params(direction="in", length=3)
    axs[-1].set_xlabel("MD step (nstep)")
    save(fig, os.path.join(folder, args.prefix + "mlff.png"), args, plt)


# --------------------------------------------------------------------------- CLI
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Analyse VASP MD runs: thermodynamics (OSZICAR), RDF and MSD/diffusion (XDATCAR), "
                    "MLFF errors (ML_LOGFILE).\n"
                    "VASP 分子动力学分析：温度/能量/体积、RDF、MSD/扩散系数、ML_LOGFILE 误差曲线。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples / 示例:\n"
               "  python md_analysis.py thermo --skip-ps 1\n"
               "  python md_analysis.py rdf --pair Si-Si --skip-ps 2\n"
               "  python md_analysis.py msd --species Li --fit-range 1 5\n"
               "  python md_analysis.py mlff\n"
               "  python md_analysis.py all -d run1 --skip-ps 1\n")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("-d", "--dir", default=".", help="calculation directory (default .)")
        sp.add_argument("--xdatcar", default="XDATCAR", help="XDATCAR name inside --dir")
        sp.add_argument("--oszicar", default="OSZICAR", help="OSZICAR name inside --dir")
        sp.add_argument("--potim", type=float, default=None, help="time step in fs (default from OUTCAR/INCAR)")
        sp.add_argument("--dt", type=float, default=None,
                        help="time between XDATCAR frames in fs (overrides POTIM x NBLOCK)")
        g = sp.add_mutually_exclusive_group()
        g.add_argument("--skip", type=int, default=0, help="frames/steps to discard as equilibration")
        g.add_argument("--skip-ps", type=float, default=None, help="equilibration time to discard (ps)")
        sp.add_argument("--stride", type=int, default=1, help="use every n-th XDATCAR frame (rdf/msd)")
        sp.add_argument("--prefix", default="", help="prefix for output figure names")
        sp.add_argument("--figsize", type=float, nargs=2, default=[5.5, 3.8], metavar=("W", "H"))
        sp.add_argument("--dpi", type=int, default=300)
        sp.add_argument("--no-plot", action="store_true")
        sp.add_argument("--show", action="store_true", help="open interactive windows")

    th = sub.add_parser("thermo", help="T, energies, volume vs time (OSZICAR + XDATCAR)")
    common(th)
    rd = sub.add_parser("rdf", help="radial distribution function g(r) and coordination number")
    common(rd)
    rd.add_argument("--pair", action="append", default=[], help="element pair, e.g. Si-Si, O-H, or 'all' (repeatable)")
    rd.add_argument("--rmax", type=float, default=None, help="max r in Angstrom (default half the cell width, <= 10)")
    rd.add_argument("--dr", type=float, default=0.02, help="bin width in Angstrom (default 0.02)")
    rd.add_argument("--dat", default=None, help="export r, g(r), n(r)")
    ms = sub.add_parser("msd", help="mean-square displacement and diffusion coefficient")
    common(ms)
    ms.add_argument("--species", nargs="+", default=None, help="elements to analyse (default all)")
    ms.add_argument("--fit-range", type=float, nargs=2, default=None, metavar=("T0", "T1"),
                    help="linear-fit window in ps (default 10%%-50%% of the trajectory)")
    ms.add_argument("--dim", type=int, choices=[1, 2, 3], default=3, help="dimensionality d in D = slope/(2d)")
    ms.add_argument("--keep-com", action="store_true", help="do not remove the centre-of-mass drift")
    ms.add_argument("--dat", default=None, help="export t, MSD per species")
    ml = sub.add_parser("mlff", help="ML_LOGFILE: ERR / BEEF curves and STATUS statistics")
    common(ml)
    ml.add_argument("-f", "--file", default=None, help="ML_LOGFILE path (default <dir>/ML_LOGFILE)")
    al = sub.add_parser("all", help="run all analyses for which files exist")
    common(al)
    al.add_argument("--pair", action="append", default=[])
    al.add_argument("--rmax", type=float, default=None)
    al.add_argument("--dr", type=float, default=0.02)
    al.add_argument("--species", nargs="+", default=None)
    al.add_argument("--fit-range", type=float, nargs=2, default=None)
    al.add_argument("--dim", type=int, choices=[1, 2, 3], default=3)
    al.add_argument("--keep-com", action="store_true")
    al.add_argument("-f", "--file", default=None)
    al.set_defaults(dat=None)
    return p


def main(argv=None) -> None:
    try:
        sys.stdout.reconfigure(errors="backslashreplace")  # type: ignore[attr-defined]
    except (AttributeError, ValueError):
        pass
    args = build_parser().parse_args(argv)
    folder = args.dir
    if not os.path.isdir(folder):
        die(f"directory not found: {folder}")
    if args.stride < 1:
        die("--stride must be >= 1")
    plt = None if args.no_plot else setup_matplotlib(args.show)
    if args.cmd == "thermo":
        cmd_thermo(args, folder, plt)
    elif args.cmd == "rdf":
        cmd_rdf(args, folder, plt)
    elif args.cmd == "msd":
        cmd_msd(args, folder, plt)
    elif args.cmd == "mlff":
        cmd_mlff(args, folder, plt)
    else:
        ran = False
        if os.path.isfile(os.path.join(folder, args.oszicar)):
            print("=== thermo ===")
            cmd_thermo(args, folder, plt)
            ran = True
        if os.path.isfile(os.path.join(folder, args.xdatcar)):
            print("=== rdf ===")
            cmd_rdf(args, folder, plt)
            print("=== msd ===")
            cmd_msd(args, folder, plt)
            ran = True
        if os.path.isfile(args.file or os.path.join(folder, "ML_LOGFILE")):
            print("=== mlff ===")
            cmd_mlff(args, folder, plt)
            ran = True
        if not ran:
            die(f"no OSZICAR, XDATCAR or ML_LOGFILE in {folder}")


if __name__ == "__main__":
    main()
