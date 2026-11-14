#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
md_analyze.py -- post-processing of CP2K AIMD runs (&MOTION/&MD).
                 CP2K 从头算分子动力学后处理：.ener 温度/势能/守恒量漂移，
                 RDF 与配位数，MSD 与扩散系数，VACF 与振动态密度，偶极轨迹 IR 光谱。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Sub-commands
------------
  ener  <PROJECT>-1.ener
        Columns (CP2K): Step Nr., Time[fs], Kin.[a.u.], Temp[K], Pot.[a.u.],
        Cons Qty[a.u.], UsedTime[s].  Prints mean/std of T and E_pot, the drift
        of the conserved quantity (linear fit and end points; per atom with
        --natoms) and plots T, E_pot and Cons Qty versus time.
  rdf   <PROJECT>-pos-1.xyz --pairs O-O O-H  (+ cell)
        Partial radial distribution functions g(r) with the minimum-image
        convention for any (also triclinic, also NPT) cell, plus the running
        coordination number n(r); first peak / first minimum are reported.
  msd   <PROJECT>-pos-1.xyz  (+ cell)
        Mean-square displacement per element (multiple time origins), after
        unwrapping jumps across the cell and removing the centre-of-mass drift;
        D from a linear fit MSD = 6 D t + c in a chosen time window.
  vacf  <PROJECT>-vel-1.xyz    (or a position trajectory with --from-positions)
        Velocity autocorrelation function (optionally mass-weighted) and the
        vibrational density of states (Fourier transform with a Hann window).
  ir    <dipole file of &DFT/&PRINT/&MOMENTS (PERIODIC T)> --dt-fs DT
        IR line shape from the autocorrelation of the dipole time derivative;
        Berry-phase dipoles are unwrapped with the printed cell quanta.

Cell for rdf/msd (positions in CP2K xyz trajectories carry no cell):
  --cell <PROJECT>-1.cell   per-frame cell (NPT; matched by step number)
  --abc A B C [--angles ALPHA BETA GAMMA]   constant cell (Angstrom, degrees)
  --lattice 9 numbers       constant cell vectors a b c (Angstrom)
Frame times are taken from the xyz titles (" i = ..., time = ..."); use
--dt-fs if the titles do not carry a time.

Usage examples
--------------
    python md_analyze.py ener water32_nve-1.ener --natoms 96 -o nve_ener.png
    python md_analyze.py rdf water32_nvt-pos-1.xyz --abc 9.8654 9.8654 9.8654 --pairs O-O O-H --skip-fs 1000
    python md_analyze.py rdf water32_npt-pos-1.xyz --cell water32_npt-1.cell --pairs O-O
    python md_analyze.py msd water32_nvt-pos-1.xyz --abc 9.8654 9.8654 9.8654 --elements O --fit 1 4
    python md_analyze.py vacf water32_nve-vel-1.xyz --elements H O --max-lag-fs 1000
    python md_analyze.py ir water32_nve_dipole.dat --dt-fs 1.0 --temperature 300 --qcf harmonic

Notes: the analysis assumes equally spaced frames; statistics from a few ps of
a 32-molecule box are noisy and D needs a finite-size correction
(Yeh & Hummer, J. Phys. Chem. B 108, 15873 (2004)).

Dependencies: python >= 3.9, numpy, matplotlib (ase / MDAnalysis are not needed).
"""
import argparse
import math
import os
import re
import sys

NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
HA2EV = 27.211386245988
FS_INV_TO_CM = 1.0e15 / 2.99792458e10      # 1/fs -> cm^-1
HBAR_EVS = 6.582119569e-16
KB_EV = 8.617333262e-5
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]

MASSES = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867,
    "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95, "Tc": 98.0,
    "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33, "La": 138.91,
    "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 145.0, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25,
    "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08,
    "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209.0, "At": 210.0,
    "Rn": 222.0, "D": 2.014,
}


def die(msg):
    sys.exit("md_analyze.py: error: " + msg)


def need_numpy():
    try:
        import numpy as np
    except ImportError:
        die("numpy is required (pip install numpy)")
    return np


def get_plt(show):
    try:
        import matplotlib
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        die("matplotlib is required for plotting (pip install matplotlib); use -o none")
    return plt


def element(sym):
    m = re.match(r"([A-Z][a-z]?)", sym)
    return m.group(1) if m else sym


def mass_of(sym):
    el = element(sym)
    if el not in MASSES:
        die(f"no atomic mass for element '{sym}'")
    return MASSES[el]


def finish(fig, plt, args):
    fig.tight_layout()
    if args.output.lower() != "none":
        fig.savefig(args.output, dpi=args.dpi)
        print(f"wrote {args.output}")
    if args.show:
        plt.show()


def style(ax):
    ax.grid(alpha=0.25, lw=0.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


# --------------------------------------------------------------------------- readers
def read_ener(path):
    np = need_numpy()
    rows = []
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                t = line.split()
                if len(t) < 6:
                    continue
                try:
                    rows.append([float(x) for x in t[:7]] + [0.0] * (7 - len(t[:7])))
                except ValueError:
                    print(f"warning: skipping malformed line: {line.strip()[:60]}")
    except OSError as exc:
        die(f"cannot read {path}: {exc}")
    if not rows:
        die(f"{path}: no data - expected the CP2K <PROJECT>-1.ener format")
    d = np.array(rows)
    # restarted runs append to the same file: keep the last entry of every step
    steps = d[:, 0].astype(int)
    _, last = np.unique(steps[::-1], return_index=True)
    keep = np.sort(len(steps) - 1 - last)
    if keep.size < len(steps):
        print(f"note: {len(steps) - keep.size} repeated step(s) (restarts) removed")
    d = d[keep]
    return d[np.argsort(d[:, 0])]


def read_xyz(path, stride=1, max_frames=None):
    """Return symbols, coords (nf, na, 3), steps (nf), times (nf, NaN if absent)."""
    np = need_numpy()
    syms, frames, steps, times = None, [], [], []
    try:
        fh = open(path, "r", errors="replace")
    except OSError as exc:
        die(f"cannot read {path}: {exc}")
    with fh:
        iframe = 0
        while True:
            line = fh.readline()
            if not line:
                break
            if not line.strip():
                continue
            try:
                na = int(line.split()[0])
            except ValueError:
                die(f"{path}: expected the atom count, got: {line.strip()[:60]}")
            title = fh.readline()
            block = [fh.readline() for _ in range(na)]
            if len(block) < na or not block[-1]:
                print(f"warning: {path}: incomplete last frame ignored")
                break
            if iframe % stride == 0:
                s, xyz = [], []
                for ln in block:
                    t = ln.split()
                    if len(t) < 4:
                        die(f"{path}: malformed atom line: {ln.strip()[:60]}")
                    s.append(t[0])
                    xyz.append([float(t[1]), float(t[2]), float(t[3])])
                if syms is None:
                    syms = s
                elif len(s) != len(syms):
                    die(f"{path}: atom count changes between frames")
                frames.append(xyz)
                m = re.search(r"\bi\s*=\s*(\d+)", title)
                steps.append(int(m.group(1)) if m else iframe)
                m = re.search(r"time\s*=\s*(" + NUM + ")", title)
                times.append(float(m.group(1)) if m else float("nan"))
                if max_frames and len(frames) >= max_frames:
                    break
            iframe += 1
    if not frames:
        die(f"{path}: no frames read")
    return syms, np.array(frames), np.array(steps), np.array(times)


def read_cell_file(path):
    np = need_numpy()
    steps, cells = [], []
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                t = line.split()
                if len(t) < 11:
                    continue
                steps.append(int(float(t[0])))
                cells.append([float(x) for x in t[2:11]])
    except OSError as exc:
        die(f"cannot read {path}: {exc}")
    if not steps:
        die(f"{path}: no data - expected the CP2K <PROJECT>-1.cell format")
    return np.array(steps), np.array(cells).reshape(-1, 3, 3)   # rows: A, B, C


def cell_from_abc(np, abc, angles):
    a, b, c = abc
    al, be, ga = [math.radians(x) for x in angles]
    va = [a, 0.0, 0.0]
    vb = [b * math.cos(ga), b * math.sin(ga), 0.0]
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz2 = c * c - cx * cx - cy * cy
    if cz2 <= 0:
        die("--angles give an impossible cell")
    return np.array([va, vb, [cx, cy, math.sqrt(cz2)]])


def frame_cells(np, args, steps):
    """Array (nf, 3, 3) of cell rows per frame, or None."""
    if args.cell:
        cs, cc = read_cell_file(args.cell)
        idx = np.searchsorted(cs, steps, side="right") - 1
        if (idx < 0).any():
            print("warning: some frames precede the first .cell entry; using the first cell for them")
        idx = np.clip(idx, 0, len(cs) - 1)
        return cc[idx]
    if args.lattice:
        h = np.array(args.lattice, dtype=float).reshape(3, 3)
    elif args.abc:
        h = cell_from_abc(np, args.abc, args.angles)
    else:
        return None
    if abs(np.linalg.det(h)) < 1e-6:
        die("cell vectors are linearly dependent")
    return np.repeat(h[None], len(steps), axis=0)


def frame_dt(np, args, times):
    """Time between frames in fs."""
    if args.dt_fs:
        return args.dt_fs
    if len(times) < 2 or np.isnan(times).any():
        die("no 'time =' in the xyz titles; give --dt-fs (time between stored frames)")
    dts = np.diff(times)
    dt = float(np.median(dts))
    if dt <= 0:
        die("non-increasing times in the trajectory (restarted run appended?); clean the file or use --dt-fs")
    if np.abs(dts - dt).max() > 1e-3 * max(1.0, dt):
        print("warning: frames are not equally spaced in time (appended restarts?); results assume dt = "
              f"{dt:.4f} fs")
    return dt


def skip_frames(np, args, times, dt):
    if not args.skip_fs:
        return 0
    if np.isnan(times).any():
        if not dt:
            die("--skip-fs: no 'time =' in the xyz titles; give --dt-fs as well")
        return int(round(args.skip_fs / dt))
    return int(np.searchsorted(times, times[0] + args.skip_fs - 1e-9))


def select_atoms(np, syms, elements):
    els = [element(s) for s in syms]
    if not elements:
        return list(dict.fromkeys(els)), els
    missing = [e for e in elements if e not in els]
    if missing:
        die(f"elements {missing} not in the trajectory (found {sorted(set(els))})")
    return elements, els


def unwrap(np, X, cells):
    """Remove jumps across periodic boundaries using minimum-image steps."""
    if cells is None:
        return X
    U = X.copy()
    for t in range(1, len(X)):
        h = cells[t]
        d = X[t] - X[t - 1]
        f = d @ np.linalg.inv(h)
        f -= np.round(f)
        U[t] = U[t - 1] + f @ h
    return U


# --------------------------------------------------------------------------- ener
def cmd_ener(args):
    np = need_numpy()
    d = read_ener(args.file)
    t = d[:, 1]
    sel = t >= t[0] + args.skip_fs
    if sel.sum() < 2:
        die("not enough data after --skip-fs")
    d = d[sel]
    step, t, kin, temp, pot, cons = (d[:, i] for i in range(6))
    tps = (t - t[0]) / 1000.0
    print(f"file            : {args.file}")
    print(f"steps           : {int(step[0])} .. {int(step[-1])} ({len(step)} lines), time {t[0]:.1f} .. {t[-1]:.1f} fs")
    print(f"temperature     : {temp.mean():.2f} +- {temp.std():.2f} K   (min {temp.min():.1f}, max {temp.max():.1f})")
    print(f"potential energy: {pot.mean():.8f} +- {pot.std():.8f} Ha")
    span = tps[-1] - tps[0]
    if span > 0:
        slope, icpt = np.polyfit(tps, cons, 1)
        endp = (cons[-1] - cons[0]) / span
        print(f"Cons Qty drift  : fit {slope:+.3e} Ha/ps, end points {endp:+.3e} Ha/ps, "
              f"max |dev| {np.abs(cons - cons[0]).max():.3e} Ha")
        if args.natoms:
            mev = slope * HA2EV * 1000.0 / args.natoms
            print(f"                  = {mev:+.4f} meV/atom/ps ({slope / args.natoms:+.3e} Ha/atom/ps)")
            print("                  (rule of thumb for NVE: |drift| < ~1 meV/atom/ps; NVT/NPT Cons Qty "
                  "includes the thermostat/barostat energy)")
    if args.output.lower() == "none" and not args.show:
        return
    plt = get_plt(args.show)
    fig, axs = plt.subplots(3, 1, figsize=(6.0, 6.4), sharex=True)
    axs[0].plot(tps, temp, color=COLORS[0], lw=1.0)
    axs[0].axhline(temp.mean(), color="#8a8a86", lw=0.8, ls="--")
    axs[0].set_ylabel("T (K)")
    scale, unit = HA2EV * 1000.0, "meV"
    if args.natoms:
        scale, unit = scale / args.natoms, "meV/atom"
    axs[1].plot(tps, (pot - pot[0]) * scale, color=COLORS[1], lw=1.0)
    axs[1].set_ylabel(f"E_pot - E_pot(0) ({unit})")
    axs[2].plot(tps, (cons - cons[0]) * scale, color=COLORS[2], lw=1.0)
    axs[2].set_ylabel(f"Cons Qty drift ({unit})")
    axs[2].set_xlabel("time since first analysed step (ps)")
    for ax in axs:
        style(ax)
    finish(fig, plt, args)


# --------------------------------------------------------------------------- rdf
def cmd_rdf(args):
    np = need_numpy()
    syms, X, steps, times = read_xyz(args.file, args.stride)
    cells = frame_cells(np, args, steps)
    if cells is None:
        die("rdf needs a cell: --cell FILE, --abc A B C or --lattice ...")
    i0 = skip_frames(np, args, times, args.dt_fs)
    X, cells = X[i0:], cells[i0:]
    if len(X) == 0:
        die("no frames left after --skip-fs")
    els = [element(s) for s in syms]

    vols = np.abs(np.linalg.det(cells))
    widths = []
    for h in cells:
        v = abs(np.linalg.det(h))
        widths.append(min(v / np.linalg.norm(np.cross(h[1], h[2])), v / np.linalg.norm(np.cross(h[2], h[0])),
                          v / np.linalg.norm(np.cross(h[0], h[1]))))
    rlim = 0.5 * min(widths)
    rmax = args.rmax or rlim
    if rmax > rlim + 1e-9:
        print(f"warning: --rmax {rmax:.3f} exceeds half the smallest cell width ({rlim:.3f} A); clipped")
        rmax = rlim
    edges = np.arange(0.0, rmax + 0.5 * args.dr, args.dr)
    if edges[-1] > rmax:
        edges = edges[:-1]
    r = 0.5 * (edges[1:] + edges[:-1])
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)

    results = {}
    for pair in args.pairs:
        try:
            a, b = pair.split("-")
        except ValueError:
            die(f"--pairs: '{pair}' is not of the form A-B")
        ia = [i for i, e in enumerate(els) if e == a]
        ib = [i for i, e in enumerate(els) if e == b]
        if not ia or not ib:
            die(f"--pairs {pair}: element not found (found {sorted(set(els))})")
        same = a == b
        nb_eff = len(ib) - 1 if same else len(ib)
        if nb_eff <= 0:
            die(f"--pairs {pair}: need at least two {a} atoms")
        g = np.zeros(len(r))
        for f in range(len(X)):
            h = cells[f]
            hinv = np.linalg.inv(h)
            sa = X[f, ia] @ hinv
            sb = X[f, ib] @ hinv
            hist = np.zeros(len(r))
            for c0 in range(0, len(ia), 256):
                ds = sb[None, :, :] - sa[c0:c0 + 256, None, :]
                ds -= np.round(ds)
                dist = np.linalg.norm(ds @ h, axis=2)
                if same:
                    dist[dist < 1e-8] = np.inf
                hist += np.histogram(dist.ravel(), bins=edges)[0]
            rho_b = nb_eff / vols[f]
            g += hist / (len(ia) * rho_b * shell)
        g /= len(X)
        rho_mean = nb_eff / vols.mean()
        cn = np.cumsum(4.0 * np.pi * rho_mean * g * r * r * args.dr)
        results[pair] = (g, cn)
        ip = int(np.argmax(g))
        txt = f"{pair:7s}: first peak r = {r[ip]:.3f} A, g = {g[ip]:.3f}"
        after = g[ip:]
        if after.size > 3:
            imin = ip + int(np.argmin(after[: max(3, int(1.5 / args.dr))]))
            txt += f"; first minimum r = {r[imin]:.3f} A, g = {g[imin]:.3f}, n(r_min) = {cn[imin]:.2f}"
        print(txt)
    print(f"{len(X)} frames, r_max = {rmax:.3f} A, dr = {args.dr} A, <V> = {vols.mean():.3f} A^3")

    if args.csv:
        cols, heads = [r], ["r_A"]
        for p, (g, cn) in results.items():
            cols += [g, cn]
            heads += [f"g_{p}", f"n_{p}"]
        np.savetxt(args.csv, np.column_stack(cols), delimiter=",", header=",".join(heads), comments="", fmt="%.6f")
        print(f"wrote {args.csv}")
    if args.output.lower() == "none" and not args.show:
        return
    plt = get_plt(args.show)
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    for i, (p, (g, cn)) in enumerate(results.items()):
        ax.plot(r, g, color=COLORS[i % 8], lw=1.5, label=p)
    ax.axhline(1.0, color="#8a8a86", lw=0.8, ls=":")
    ax.set_xlabel(r"r ($\mathrm{\AA}$)")
    ax.set_ylabel("g(r)")
    ax.set_xlim(0, rmax)
    if len(results) > 1:
        ax.legend(frameon=False)
    style(ax)
    finish(fig, plt, args)


# --------------------------------------------------------------------------- msd
def cmd_msd(args):
    np = need_numpy()
    syms, X, steps, times = read_xyz(args.file, args.stride)
    dt = frame_dt(np, args, times)
    i0 = skip_frames(np, args, times, dt)
    cells = frame_cells(np, args, steps)
    X = X[i0:]
    cells = cells[i0:] if cells is not None else None
    if len(X) < 4:
        die("too few frames for an MSD")
    if cells is None:
        print("note: no cell given - coordinates are used as written (CP2K does not wrap MD "
              "trajectories by default); give --abc/--cell if your file is wrapped")
    U = unwrap(np, X, cells)
    if not args.keep_com:
        m = np.array([mass_of(s) for s in syms])
        com = (U * m[None, :, None]).sum(axis=1) / m.sum()
        U = U - (com - com[0])[:, None, :]
    wanted, els = select_atoms(np, syms, args.elements)
    nf = len(U)
    maxlag = nf // 2 if args.max_lag_fs is None else min(nf - 1, int(round(args.max_lag_fs / dt)))
    lags = np.arange(0, maxlag + 1)
    tlag = lags * dt / 1000.0       # ps
    res = {}
    for el in wanted:
        idx = [i for i, e in enumerate(els) if e == el]
        Y = U[:, idx, :]
        msd = np.zeros(len(lags))
        for k, L in enumerate(lags):
            if L == 0:
                continue
            o = np.arange(0, nf - L, args.origin_stride)
            d = Y[o + L] - Y[o]
            msd[k] = (d * d).sum(axis=2).mean()
        res[el] = msd
    tmax = tlag[-1]
    t1, t2 = args.fit if args.fit else (0.2 * tmax, 0.8 * tmax)
    fsel = (tlag >= t1) & (tlag <= t2)
    if fsel.sum() < 3:
        die(f"fit window {t1:.3f}-{t2:.3f} ps contains < 3 points (max lag {tmax:.3f} ps)")
    print(f"{nf} frames, dt = {dt:.3f} fs, max lag = {tmax:.3f} ps, fit window {t1:.3f}-{t2:.3f} ps")
    fits = {}
    for el, msd in res.items():
        slope, icpt = np.polyfit(tlag[fsel], msd[fsel], 1)
        D = slope / 6.0                       # A^2/ps
        fits[el] = (slope, icpt)
        print(f"{el:3s}: D = {D:.4f} A^2/ps = {D * 1e-4:.3e} cm^2/s = {D * 10:.3f} x 1e-5 cm^2/s"
              f"   (MSD at max lag {msd[-1]:.3f} A^2)")
    print("note: finite-size (Yeh-Hummer) correction not applied; check that the fit window lies in the "
          "linear (diffusive) regime on a log-log plot")
    if args.csv:
        cols, heads = [tlag], ["t_ps"]
        for el, msd in res.items():
            cols.append(msd)
            heads.append(f"msd_{el}_A2")
        np.savetxt(args.csv, np.column_stack(cols), delimiter=",", header=",".join(heads), comments="", fmt="%.6f")
        print(f"wrote {args.csv}")
    if args.output.lower() == "none" and not args.show:
        return
    plt = get_plt(args.show)
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    for i, (el, msd) in enumerate(res.items()):
        ax.plot(tlag, msd, color=COLORS[i % 8], lw=1.5, label=el)
        s, c = fits[el]
        ax.plot(tlag[fsel], s * tlag[fsel] + c, color="#4a4a48", lw=0.8, ls="--")
    ax.axvspan(t1, t2, color="#d9d8d2", alpha=0.3, lw=0)
    ax.set_xlabel("lag time (ps)")
    ax.set_ylabel(r"MSD ($\mathrm{\AA}^2$)")
    if len(res) > 1:
        ax.legend(frameon=False)
    style(ax)
    finish(fig, plt, args)


# --------------------------------------------------------------------------- correlation helpers
def autocorr(np, x):
    """Unbiased autocorrelation of x (n, ...) along axis 0, summed over the other axes."""
    n = x.shape[0]
    y = x.reshape(n, -1)
    nfft = 1 << int(math.ceil(math.log2(2 * n)))
    f = np.fft.rfft(y, n=nfft, axis=0)
    ac = np.fft.irfft(f * np.conj(f), n=nfft, axis=0)[:n].sum(axis=1)
    return ac / (n - np.arange(n))


def spectrum(np, c, dt_fs, pad):
    """One-sided cosine transform of a windowed correlation function; returns (cm^-1, intensity)."""
    n = len(c)
    w = np.hanning(2 * n)[n:]          # decaying half of a Hann window, 1 at t = 0
    y = c * w
    full = np.concatenate([y, y[-2:0:-1]])   # even extension
    nfft = len(full) * pad
    spec = np.fft.rfft(full, n=nfft).real * dt_fs
    freq = np.fft.rfftfreq(nfft, d=dt_fs) * FS_INV_TO_CM
    return freq, spec


# --------------------------------------------------------------------------- vacf
def cmd_vacf(args):
    np = need_numpy()
    syms, X, steps, times = read_xyz(args.file, args.stride)
    dt = frame_dt(np, args, times)
    i0 = skip_frames(np, args, times, dt)
    X = X[i0:]
    if args.from_positions:
        cells = frame_cells(np, args, steps)
        cells = cells[i0:] if cells is not None else None
        U = unwrap(np, X, cells)
        V = (U[2:] - U[:-2]) / (2.0 * dt)        # A/fs, central differences
        print("velocities from central finite differences of the positions")
    else:
        V = X                                      # CP2K -vel-1.xyz: bohr/au_time (units cancel after normalisation)
    if len(V) < 8:
        die("too few frames for a VACF")
    wanted, els = select_atoms(np, syms, args.elements)
    nf = len(V)
    maxlag = min(nf - 1, nf // 2 if args.max_lag_fs is None else int(round(args.max_lag_fs / dt)))
    res = {}
    for el in wanted:
        idx = [i for i, e in enumerate(els) if e == el]
        w = np.array([mass_of(syms[i]) for i in idx]) if args.mass_weighted else np.ones(len(idx))
        Vw = V[:, idx, :] * np.sqrt(w)[None, :, None]
        c = autocorr(np, Vw)[: maxlag + 1]
        res[el] = c
    tot = sum(res.values())
    tlag = np.arange(maxlag + 1) * dt
    freq, vdos = spectrum(np, tot / tot[0], dt, args.pad)
    parts = {el: spectrum(np, c / tot[0], dt, args.pad)[1] for el, c in res.items()}
    fsel = freq <= args.max_freq
    print(f"{nf} velocity frames, dt = {dt:.3f} fs, max lag = {tlag[-1]:.1f} fs, "
          f"resolution ~ {FS_INV_TO_CM / (2 * tlag[-1]):.1f} cm^-1 (Hann window)"
          + (", mass-weighted" if args.mass_weighted else ""))
    ip = np.argmax(np.where(fsel & (freq > 50), vdos, -np.inf))
    print(f"strongest VDOS peak above 50 cm^-1: {freq[ip]:.0f} cm^-1")
    if args.csv:
        cols = [freq[fsel], vdos[fsel]] + [parts[el][fsel] for el in res]
        heads = ["freq_cm-1", "vdos_total"] + [f"vdos_{el}" for el in res]
        np.savetxt(args.csv, np.column_stack(cols), delimiter=",", header=",".join(heads), comments="", fmt="%.6e")
        root, ext = os.path.splitext(args.csv)
        np.savetxt(root + "_vacf" + (ext or ".csv"), np.column_stack([tlag, tot / tot[0]]), delimiter=",",
                   header="t_fs,vacf_normalised", comments="", fmt="%.6e")
        print(f"wrote {args.csv} and {root}_vacf{ext or '.csv'}")
    if args.output.lower() == "none" and not args.show:
        return
    plt = get_plt(args.show)
    fig, axs = plt.subplots(2, 1, figsize=(6.0, 5.6))
    axs[0].plot(tlag, tot / tot[0], color=COLORS[0], lw=1.2)
    axs[0].axhline(0, color="#8a8a86", lw=0.8)
    axs[0].set_xlabel("t (fs)")
    axs[0].set_ylabel("VACF (normalised)")
    if len(res) > 1:
        for i, el in enumerate(res):
            axs[1].plot(freq[fsel], parts[el][fsel], color=COLORS[(i + 1) % 8], lw=1.0, label=el)
    axs[1].plot(freq[fsel], vdos[fsel], color="#4a4a48", lw=1.4, label="total")
    axs[1].set_xlabel(r"wavenumber (cm$^{-1}$)")
    axs[1].set_ylabel("VDOS (arb. units)")
    axs[1].set_xlim(0, args.max_freq)
    if len(res) > 1:
        axs[1].legend(frameon=False)
    for ax in axs:
        style(ax)
    finish(fig, plt, args)


# --------------------------------------------------------------------------- ir
def read_dipoles(path):
    """Dipole vectors of all 'Dipole moment [...]' blocks and the Berry-phase cell quanta (first block)."""
    np = need_numpy()
    dips, qrows = [], []
    want = False
    vec = r"\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")"
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                s = line.strip()
                if len(qrows) < 3:
                    m = re.match(r"\[([XYZ])\]\s*=?\s*\[" + vec, s)
                    if m and m.group(1) == "XYZ"[len(qrows)]:
                        qrows.append([float(m.group(i).replace("D", "E")) for i in (2, 3, 4)])
                        continue
                    if m:
                        qrows = []
                if s.startswith("Dipole moment ["):
                    want = True
                    continue
                if want:
                    m = re.search(r"X=\s*(" + NUM + r")\s+Y=\s*(" + NUM + r")\s+Z=\s*(" + NUM + r")", s)
                    if m:
                        dips.append([float(m.group(i).replace("D", "E")) for i in (1, 2, 3)])
                    want = False
    except OSError as exc:
        die(f"cannot read {path}: {exc}")
    if len(dips) < 8:
        die(f"{path}: found {len(dips)} 'Dipole moment' blocks - need a dipole trajectory "
            "(&DFT/&PRINT/&MOMENTS with &EACH MD n)")
    quanta = np.array(qrows) if len(qrows) == 3 else None
    return np.array(dips), quanta


def cmd_ir(args):
    np = need_numpy()
    M, quanta = read_dipoles(args.file)
    dt = args.dt_fs
    i0 = int(round(args.skip_fs / dt)) if args.skip_fs else 0
    M = M[i0:]
    if len(M) < 8:
        die("too few dipole frames after --skip-fs")
    if quanta is not None and abs(np.linalg.det(quanta)) > 1e-8:
        qinv = np.linalg.inv(quanta)
        jumps = 0
        for t in range(1, len(M)):
            d = M[t] - M[t - 1]
            n = np.round(qinv @ d)
            if np.any(n != 0):
                jumps += 1
                M[t:] -= quanta @ n
        print(f"Berry-phase dipole unwrapped with the printed cell quanta ({jumps} jump(s) removed)")
    else:
        print("note: no cell quanta found (non-periodic MOMENTS or old format); dipoles used as written")
    Md = (M[2:] - M[:-2]) / (2.0 * dt)          # Debye/fs
    nf = len(Md)
    maxlag = min(nf - 1, nf // 2 if args.max_lag_fs is None else int(round(args.max_lag_fs / dt)))
    c = autocorr(np, Md - Md.mean(axis=0))[: maxlag + 1]
    freq, I = spectrum(np, c / c[0], dt, args.pad)
    I = np.clip(I, 0.0, None)
    if args.qcf == "harmonic":
        x = (freq / FS_INV_TO_CM) * 1e15 * 2 * np.pi * HBAR_EVS / (KB_EV * args.temperature)   # hbar*omega/kT
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(x > 1e-8, x / (1.0 - np.exp(-x)), 1.0)
        I = I * f
    fsel = (freq <= args.max_freq)
    I = I / I[fsel & (freq > 50)].max() if (fsel & (freq > 50)).any() else I
    print(f"{len(M)} dipole frames, dt = {dt} fs, max lag = {maxlag * dt:.1f} fs, "
          f"resolution ~ {FS_INV_TO_CM / (2 * maxlag * dt):.1f} cm^-1")
    peaks = [i for i in range(1, len(I) - 1) if fsel[i] and freq[i] > 200 and I[i] > I[i - 1]
             and I[i] >= I[i + 1] and I[i] > 0.2]
    print("main bands (> 20 % of max): " + ", ".join(f"{freq[i]:.0f}" for i in peaks[:12]) + " cm^-1")
    print("note: classical line shape (|FT dmu/dt|^2); intensities are relative, not km/mol")
    if args.csv:
        np.savetxt(args.csv, np.column_stack([freq[fsel], I[fsel]]), delimiter=",",
                   header="freq_cm-1,intensity_rel", comments="", fmt="%.6e")
        print(f"wrote {args.csv}")
    if args.output.lower() == "none" and not args.show:
        return
    plt = get_plt(args.show)
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    ax.plot(freq[fsel], I[fsel], color=COLORS[0], lw=1.3)
    ax.set_xlim(0, args.max_freq)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"wavenumber (cm$^{-1}$)")
    ax.set_ylabel("IR intensity (normalised)")
    style(ax)
    finish(fig, plt, args)


# --------------------------------------------------------------------------- CLI
def build_parser():
    ap = argparse.ArgumentParser(
        prog="md_analyze.py",
        description="Post-processing of CP2K AIMD output: ener / rdf / msd / vacf / ir. "
                    "CP2K AIMD 后处理（守恒量漂移、RDF、MSD、VACF/VDOS、IR）。",
        epilog="Run 'python md_analyze.py <command> -h' for the options of each command.\n"
               "Examples:\n"
               "  python md_analyze.py ener water32_nve-1.ener --natoms 96\n"
               "  python md_analyze.py rdf water32_nvt-pos-1.xyz --abc 9.8654 9.8654 9.8654 --pairs O-O O-H\n"
               "  python md_analyze.py msd water32_nvt-pos-1.xyz --abc 9.8654 9.8654 9.8654 --elements O\n"
               "  python md_analyze.py vacf water32_nve-vel-1.xyz\n"
               "  python md_analyze.py ir water32_nve_dipole.dat --dt-fs 1.0\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common_out(p, default):
        p.add_argument("-o", "--output", default=default, help=f"figure file (default {default}; 'none' = no figure)")
        p.add_argument("--csv", metavar="FILE", help="write the data to a CSV file")
        p.add_argument("--dpi", type=int, default=200, help="figure resolution (default 200)")
        p.add_argument("--show", action="store_true", help="open an interactive window")

    def common_traj(p, cell=True):
        p.add_argument("--stride", type=int, default=1, help="use every N-th frame of the file (default 1)")
        p.add_argument("--skip-fs", type=float, default=0.0, help="discard the first T fs (equilibration)")
        p.add_argument("--dt-fs", type=float, help="time between stored frames in fs (default: from the xyz titles)")
        if cell:
            g = p.add_mutually_exclusive_group()
            g.add_argument("--cell", metavar="FILE", help="CP2K <PROJECT>-1.cell file (per-frame cell, NPT)")
            g.add_argument("--abc", nargs=3, type=float, metavar=("A", "B", "C"), help="constant cell lengths (A)")
            g.add_argument("--lattice", nargs=9, type=float, metavar="X", help="constant cell vectors a b c (A)")
            p.add_argument("--angles", nargs=3, type=float, default=[90.0, 90.0, 90.0],
                           metavar=("ALPHA", "BETA", "GAMMA"), help="cell angles for --abc (default 90 90 90)")

    p = sub.add_parser("ener", help="temperature, potential energy and conserved-quantity drift from .ener",
                       description="Analyse a CP2K <PROJECT>-1.ener file.")
    p.add_argument("file", help="<PROJECT>-1.ener")
    p.add_argument("--natoms", type=int, help="number of atoms (drift per atom)")
    p.add_argument("--skip-fs", type=float, default=0.0, help="discard the first T fs")
    common_out(p, "ener.png")
    p.set_defaults(func=cmd_ener)

    p = sub.add_parser("rdf", help="partial radial distribution functions g(r) and coordination numbers",
                       description="Partial g(r) from a CP2K position trajectory (Angstrom).")
    p.add_argument("file", help="<PROJECT>-pos-1.xyz")
    p.add_argument("--pairs", nargs="+", default=["O-O"], metavar="A-B", help="element pairs (default O-O)")
    p.add_argument("--rmax", type=float, help="maximum r in A (default: half the smallest cell width)")
    p.add_argument("--dr", type=float, default=0.02, help="bin width in A (default 0.02)")
    common_traj(p)
    common_out(p, "rdf.png")
    p.set_defaults(func=cmd_rdf)

    p = sub.add_parser("msd", help="mean-square displacement and diffusion coefficient",
                       description="MSD per element with multiple time origins; D from MSD = 6 D t + c.")
    p.add_argument("file", help="<PROJECT>-pos-1.xyz")
    p.add_argument("--elements", nargs="+", metavar="EL", help="elements to analyse (default: all)")
    p.add_argument("--max-lag-fs", type=float, help="longest lag time in fs (default: half the trajectory)")
    p.add_argument("--origin-stride", type=int, default=1, help="use every N-th frame as a time origin (default 1)")
    p.add_argument("--fit", nargs=2, type=float, metavar=("T1", "T2"),
                   help="fit window in ps (default: 20%%-80%% of the max lag)")
    p.add_argument("--keep-com", action="store_true", help="do not remove the centre-of-mass drift")
    common_traj(p)
    common_out(p, "msd.png")
    p.set_defaults(func=cmd_msd)

    p = sub.add_parser("vacf", help="velocity autocorrelation and vibrational density of states",
                       description="VACF and VDOS from a CP2K velocity trajectory (or positions).")
    p.add_argument("file", help="<PROJECT>-vel-1.xyz (or -pos-1.xyz with --from-positions)")
    p.add_argument("--from-positions", action="store_true", help="the file holds positions; differentiate them")
    p.add_argument("--elements", nargs="+", metavar="EL", help="elements to include (default: all)")
    p.add_argument("--no-mass", dest="mass_weighted", action="store_false",
                   help="do not mass-weight the VACF (default: mass-weighted)")
    p.add_argument("--max-lag-fs", type=float, help="correlation length in fs (default: half the trajectory)")
    p.add_argument("--max-freq", type=float, default=4500.0, help="upper wavenumber in cm^-1 (default 4500)")
    p.add_argument("--pad", type=int, default=4, help="zero-padding factor of the FFT (default 4)")
    common_traj(p)
    common_out(p, "vdos.png")
    p.set_defaults(func=cmd_vacf)

    p = sub.add_parser("ir", help="IR spectrum from a dipole trajectory (&MOMENTS)",
                       description="IR line shape from the dipole-derivative autocorrelation function.")
    p.add_argument("file", help="file written by &DFT/&PRINT/&MOMENTS (or a main output containing it)")
    p.add_argument("--dt-fs", type=float, required=True, help="time between dipole entries in fs "
                   "(MD TIMESTEP x the &EACH MD value)")
    p.add_argument("--skip-fs", type=float, default=0.0, help="discard the first T fs")
    p.add_argument("--max-lag-fs", type=float, help="correlation length in fs (default: half the trajectory)")
    p.add_argument("--max-freq", type=float, default=4500.0, help="upper wavenumber in cm^-1 (default 4500)")
    p.add_argument("--pad", type=int, default=4, help="zero-padding factor of the FFT (default 4)")
    p.add_argument("--qcf", choices=["none", "harmonic"], default="none",
                   help="quantum correction factor applied to the classical line shape (default none)")
    p.add_argument("--temperature", type=float, default=300.0, help="temperature in K for --qcf (default 300)")
    common_out(p, "ir.png")
    p.set_defaults(func=cmd_ir)
    return ap


def main():
    args = build_parser().parse_args()
    for name in ("stride", "origin_stride", "pad"):
        if getattr(args, name, 1) < 1:
            die(f"--{name.replace('_', '-')} must be >= 1")
    if getattr(args, "dr", 1.0) <= 0:
        die("--dr must be positive")
    if getattr(args, "dt_fs", None) is not None and args.dt_fs <= 0:
        die("--dt-fs must be positive")
    if not os.path.isfile(args.file):
        die(f"file not found: {args.file}")
    args.func(args)


if __name__ == "__main__":
    main()
