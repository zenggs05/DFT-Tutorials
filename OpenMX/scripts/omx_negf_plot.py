#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_negf_plot.py -- post-process OpenMX NEGF Step-3 output: transmission, conductance, I-V.
                    解析 System.Name.tran<i2>_<i3>（k 分辨透射）、.current、.conductance：
                    k 平均透射谱、Landauer 电导与 Seebeck 系数、I-V 曲线、k 分辨透射/电导图。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

File formats (checked against OpenMX 4.0 TRAN_Main_Analysis.c / TRAN_Main_Analysis_NC.c)
---------------------------------------------------------------------------------------
System.Name.tran<i2>_<i3>   one file per (k2, k3) point of NEGF.tran.Kgrid, header lines
    "# k2= ... k3= ...", "# SpinP_switch= n",
    "# Chemical potential (eV)      Left, Right Leads= muL muR".
    collinear / non-spin, 9 columns:
        iw  ReE(Ha) ImE(Ha)  ReE(eV) ImE(eV)  ReT_up ImT_up  ReT_down ImT_down
    non-collinear, 7 columns:   iw ReE ImE ReE(eV) ImE(eV) ReT ImT
    E is measured from the chemical potential of the LEFT lead.  For non-spin
    runs the "down" columns are a copy of "up" (T per spin channel), so the
    total transmission is T_up + T_down in all collinear cases and T for
    non-collinear runs.  Files of the full k grid are written (k and -k), so a
    plain average over all files is the correct k average.
    Note: <prefix>.traneval*/.tranevec*/.tranec* (eigenchannel output) are
    NOT transmission files and are skipped automatically.
System.Name.current      "# Corresponding bias voltage (Volt)= V",
    "# average current for up and down spins (ampere)= I_up I_down" (collinear) or
    "# average current (ampere)= I" (non-collinear); then i2 i3 k2 k3 I_up [I_down].
System.Name.conductance  "# average conductance for up and down spins (G0)= G_up G_down",
    G0 = e^2/h per spin channel (non-collinear: one value); then per-k lines.

Physics used by "trans"
-----------------------
    L_n = int T_tot(E) (E - mu)^n (-df/dE) dE            (E in eV)
    G   = (e^2/h) L_0                                      (T_tot already summed over spin)
    S   = -L_1 / (e T L_0)                                 (Seebeck, V/K)
    kappa_el = (L_2 - L_1^2 / L_0) / (h T)                 (electronic thermal conductance)
    I   = (e/h) int T_tot(E) [f_L(E) - f_R(E)] dE          (only meaningful with the T(E)
                                                            of the SAME bias calculation)
With --temp 0 the conductance is simply (e^2/h) T_tot(mu).

Sub-commands and usage examples
-------------------------------
    python omx_negf_plot.py trans negf-chain                  # all negf-chain.tran*_* files
    python omx_negf_plot.py trans negf-graphene.tran*_* --temp 300 --mu -0.2 0 0.2 -o T.png --csv T.csv
    python omx_negf_plot.py kmap negf-graphene --energy 0.0 -o Tk.png
    python omx_negf_plot.py cond negf-graphene.conductance -o Gk.png
    python omx_negf_plot.py iv negf-chain-V*.current -o iv.png --csv iv.csv

Dependencies: Python >= 3.9, numpy; matplotlib for plots.
"""

import argparse
import glob
import os
import re
import sys

KB_EV = 8.617333262e-5          # eV/K
E_CHARGE = 1.602176634e-19      # C
H_PLANCK = 6.62607015e-34       # J s
G_E2H = E_CHARGE ** 2 / H_PLANCK  # S  (= 3.874e-5 S)
TRAN_RE = re.compile(r"\.tran(\d+)_(\d+)$")


def die(msg):
    sys.exit("omx_negf_plot.py: error: " + msg)


def expand(paths):
    out = []
    for p in paths:
        hits = glob.glob(p)
        out.extend(hits if hits else [p])
    return sorted(set(out))


def tran_files(args_files):
    files = [f for f in expand(args_files) if TRAN_RE.search(f) and os.path.isfile(f)]
    if not files and len(args_files) == 1:
        files = [f for f in glob.glob(args_files[0] + ".tran*_*") if TRAN_RE.search(f)]
    if not files:
        die("no System.Name.tran<i2>_<i3> files found for %s" % " ".join(args_files))
    return sorted(files, key=lambda f: tuple(int(v) for v in TRAN_RE.search(f).groups()))


def header_values(fname):
    info = {"k2": None, "k3": None, "spin": None, "muL": None, "muR": None, "bias": None}
    with open(fname, errors="replace") as f:
        for line in f:
            if not line.startswith("#"):
                if line.strip():
                    break
                continue
            m = re.search(r"k2=\s*([-+0-9.eE]+)\s+k3=\s*([-+0-9.eE]+)", line)
            if m:
                info["k2"], info["k3"] = float(m.group(1)), float(m.group(2))
            m = re.search(r"SpinP_switch=\s*(\d+)", line)
            if m:
                info["spin"] = int(m.group(1))
            m = re.search(r"Chemical potential \(eV\)\s+Left, Right Leads=\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)", line)
            if m:
                info["muL"], info["muR"] = float(m.group(1)), float(m.group(2))
            m = re.search(r"Corresponding bias voltage \(Volt\)=\s*([-+0-9.eE]+)", line)
            if m:
                info["bias"] = float(m.group(1))
    return info


def load_tran(fname):
    """-> E (eV), T_up, T_down, T_tot, noncollinear flag."""
    import numpy as np
    try:
        d = np.loadtxt(fname, comments="#", ndmin=2)
    except (OSError, ValueError) as exc:
        die("cannot read %s (%s)" % (fname, exc))
    if d.shape[1] == 9:
        return d[:, 3], d[:, 5], d[:, 7], d[:, 5] + d[:, 7], False
    if d.shape[1] == 7:
        nan = np.full(len(d), np.nan)
        return d[:, 3], nan, nan, d[:, 5], True
    die("%s has %d columns; expected 9 (collinear) or 7 (non-collinear)" % (fname, d.shape[1]))


def k_average(files):
    import numpy as np
    E, tu, td, tt, nc = load_tran(files[0])
    tu, td, tt = tu.copy(), td.copy(), tt.copy()
    for fn in files[1:]:
        e2, a, b, c, nc2 = load_tran(fn)
        if nc2 != nc or len(e2) != len(E) or not np.allclose(e2, E, atol=1e-6):
            die("energy grid / format of %s differs from %s" % (fn, files[0]))
        tu += a
        td += b
        tt += c
    n = len(files)
    return E, tu / n, td / n, tt / n, nc


def check_grid(files):
    idx = [tuple(int(v) for v in TRAN_RE.search(f).groups()) for f in files]
    n2 = max(i for i, _ in idx) + 1
    n3 = max(j for _, j in idx) + 1
    if len(set(idx)) != n2 * n3:
        print("WARNING: %d files but the index range suggests a %d x %d grid - files missing? "
              "The k average will be biased." % (len(idx), n2, n3))
    return n2, n3, idx


def fermi_window(E, mu, kt):
    import numpy as np
    x = np.clip((E - mu) / kt, -60.0, 60.0)
    return 1.0 / (4.0 * kt * np.cosh(x / 2.0) ** 2)


def fermi(E, mu, kt):
    import numpy as np
    if kt <= 0:
        return (E < mu).astype(float)
    return 1.0 / (np.exp(np.clip((E - mu) / kt, -60.0, 60.0)) + 1.0)


def integrate(y, x):
    import numpy as np
    f = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return float(f(y, x))


# ----------------------------------------------------------------------------
def cmd_trans(a):
    import numpy as np
    files = tran_files(a.files)
    n2, n3, _ = check_grid(files)
    E, tu, td, tt, nc = k_average(files)
    info = header_values(files[0])
    print("# %d transmission file(s) (k grid %d x %d) averaged; %s; E = %.3f .. %.3f eV, %d points"
          % (len(files), n2, n3, "non-collinear" if nc else "collinear/non-spin (T_tot = T_up + T_down)",
             E[0], E[-1], len(E)))
    if info["muL"] is not None:
        dmu = info["muR"] - info["muL"]
        print("# mu_L = %.4f eV, mu_R = %.4f eV (E = 0 at mu_L); bias window [%.4f, %.4f] eV"
              % (info["muL"], info["muR"], min(0.0, dmu), max(0.0, dmu)))
    kt = KB_EV * a.temp
    print("# T = %.1f K" % a.temp)
    print("#  mu-mu_L(eV)  T_tot(mu)   G(e^2/h)   G(uS)      S(uV/K)    kappa_el(pW/K)  L/L0")
    L0_SOMMERFELD = (np.pi ** 2 / 3.0) * KB_EV ** 2       # V^2/K^2
    for mu in a.mu:
        t_at = float(np.interp(mu, E, tt))
        if a.temp <= 0:
            print("  %11.4f  %9.4f  %9.4f  %9.3f" % (mu, t_at, t_at, t_at * G_E2H * 1e6))
            continue
        if mu - 10 * kt < E[0] or mu + 10 * kt > E[-1]:
            print("  WARNING: mu = %.3f eV +- 10 kT exceeds the energy range of the T(E) data" % mu)
        w = fermi_window(E, mu, kt)
        l0 = integrate(tt * w, E)
        l1 = integrate(tt * (E - mu) * w, E)
        l2 = integrate(tt * (E - mu) ** 2 * w, E)
        if l0 < 1e-12:
            print("  %11.4f  %9.4f  %9.3e  (no transmission in the Fermi window)" % (mu, t_at, l0))
            continue
        seebeck = -(l1 / l0) / a.temp                          # V/K (energies in eV)
        kappa = (l2 - l1 ** 2 / l0) * E_CHARGE ** 2 / (H_PLANCK * a.temp)   # W/K
        lorenz = (l2 - l1 ** 2 / l0) / (l0 * a.temp ** 2) / L0_SOMMERFELD
        print("  %11.4f  %9.4f  %9.4f  %9.3f  %10.3f  %12.4f  %6.3f"
              % (mu, t_at, l0, l0 * G_E2H * 1e6, seebeck * 1e6, kappa * 1e12, lorenz))
    if info["muL"] is not None and abs(info["muR"] - info["muL"]) > 1e-6:
        dmu = info["muR"] - info["muL"]
        fl, fr = fermi(E, 0.0, kt), fermi(E, dmu, kt)
        lo, hi = min(0.0, dmu) - 10 * kt, max(0.0, dmu) + 10 * kt
        if lo < E[0] or hi > E[-1]:
            print("# WARNING: the T(E) energy range does not cover the bias window +- 10 kT")
        cur = G_E2H * integrate(tt * (fl - fr), E)          # A (integral in eV == V)
        print("# Landauer current from this T(E): %.4e A (compare with System.Name.current; OpenMX uses "
              "its own finer energy mesh NEGF.current.energy.step)" % cur)
    if a.csv:
        np.savetxt(a.csv, np.column_stack([E, tu, td, tt]), fmt="%.6e", delimiter=",",
                   header="E_minus_muL_eV,T_up,T_down,T_total", comments="")
        print("k-averaged transmission written to %s" % a.csv)
    if a.no_plot:
        return
    plt = get_plt(a.show)
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    if nc:
        ax.plot(E, tt, lw=1.2, color="#1f4e79", label="T (non-collinear)")
    elif np.max(np.abs(tu - td)) > 1e-8:
        ax.plot(E, tu, lw=1.2, color="#1f4e79", label="spin up")
        ax.plot(E, td, lw=1.2, ls="--", color="#c0504d", label="spin down")
        if a.total:
            ax.plot(E, tt, lw=1.0, color="0.3", label="total")
    else:
        if a.total:
            ax.plot(E, tt, lw=1.2, color="#1f4e79", label="total (both spins)")
        else:
            ax.plot(E, tu, lw=1.2, color="#1f4e79", label="per spin")
    if info["muL"] is not None and abs(info["muR"] - info["muL"]) > 1e-6:
        dmu = info["muR"] - info["muL"]
        ax.axvspan(min(0, dmu), max(0, dmu), color="#f2c14e", alpha=0.3, lw=0, label="bias window")
    ax.axvline(0.0, color="0.5", lw=0.8, ls=":")
    ax.set_xlabel(r"$E-\mu_L$ (eV)")
    ax.set_ylabel("transmission")
    ax.set_xlim(a.emin if a.emin is not None else E[0], a.emax if a.emax is not None else E[-1])
    if a.log:
        ax.set_yscale("log")
    else:
        ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=8)
    save(fig, plt, a, "transmission.png")


def cmd_kmap(a):
    import numpy as np
    files = tran_files(a.files)
    n2, n3, idx = check_grid(files)
    first = load_tran(files[0])
    E = first[0]
    ie = int(np.argmin(np.abs(E - a.energy)))
    tmap = np.full((n2, n3), np.nan)
    tek = {}
    k2v, k3v = np.full(n2, np.nan), np.full(n3, np.nan)
    for fn, (i2, i3) in zip(files, idx):
        e, tu, td, tt, nc = load_tran(fn)
        if len(e) != len(E):
            die("energy grid of %s differs" % fn)
        tmap[i2, i3] = tt[ie]
        tek[(i2, i3)] = tt
        h = header_values(fn)
        if h["k2"] is not None:
            k2v[i2], k3v[i3] = h["k2"], h["k3"]
    print("k-resolved T_tot at E - mu_L = %.4f eV (nearest grid point to %.4f): mean %.4f, min %.4f, max %.4f"
          % (E[ie], a.energy, np.nanmean(tmap), np.nanmin(tmap), np.nanmax(tmap)))
    if a.csv:
        with open(a.csv, "w") as f:
            f.write("i2,i3,k2,k3,T_total\n")
            for i2 in range(n2):
                for i3 in range(n3):
                    f.write("%d,%d,%.4f,%.4f,%.6e\n" % (i2, i3, k2v[i2], k3v[i3], tmap[i2, i3]))
        print("map written to %s" % a.csv)
    if a.no_plot:
        return
    plt = get_plt(a.show)
    if n2 > 1 and n3 > 1:
        fig, ax = plt.subplots(figsize=(4.8, 4.0))
        im = ax.imshow(tmap.T, origin="lower", aspect="auto", cmap="viridis",
                       extent=[np.nanmin(k2v), np.nanmax(k2v), np.nanmin(k3v), np.nanmax(k3v)])
        fig.colorbar(im, ax=ax, label="T_tot")
        ax.set_xlabel(r"$k_2$ (2$\pi$/b)")
        ax.set_ylabel(r"$k_3$ (2$\pi$/c)")
        ax.set_title(r"$E-\mu_L$ = %.3f eV" % E[ie])
    else:
        # one-dimensional k grid: show the full T(k, E) map
        line_k = k2v if n2 > 1 else k3v
        order = np.argsort(line_k)
        rows = [tek[(i, 0)] if n2 > 1 else tek[(0, i)] for i in range(max(n2, n3))]
        z = np.array(rows)[order]
        fig, ax = plt.subplots(figsize=(4.8, 4.0))
        im = ax.pcolormesh(line_k[order], E, z.T, shading="auto", cmap="viridis")
        fig.colorbar(im, ax=ax, label="T_tot")
        ax.axhline(E[ie], color="w", lw=0.8, ls="--")
        ax.set_xlabel(r"$k_%d$ (reciprocal-lattice units)" % (2 if n2 > 1 else 3))
        ax.set_ylabel(r"$E-\mu_L$ (eV)")
    save(fig, plt, a, "transmission_k.png")


def cmd_cond(a):
    import numpy as np
    if not os.path.isfile(a.file):
        die("file not found: %s" % a.file)
    avg = None
    with open(a.file, errors="replace") as f:
        for line in f:
            m = re.search(r"average conductance.*\(G0\)=\s*(.*)$", line)
            if m:
                avg = [float(v) for v in m.group(1).split()]
    try:
        d = np.loadtxt(a.file, comments="#", ndmin=2)
    except ValueError as exc:
        die("cannot read %s (%s)" % (a.file, exc))
    if d.shape[1] < 5:
        die("%s: expected 'i2 i3 k2 k3 G_up [G_down]' lines" % a.file)
    info = header_values(a.file)
    collinear = d.shape[1] >= 6
    gk = d[:, 4] + d[:, 5] if collinear else d[:, 4]
    if avg:
        tot = sum(avg) if collinear else avg[0]
        print("average conductance from header: %s G0 -> total %.5f e^2/h = %.4f uS"
              % (" ".join("%.5f" % v for v in avg), tot, tot * G_E2H * 1e6))
    print("k points: %d; k-resolved total G: min %.4f, max %.4f, mean %.5f e^2/h"
          % (len(d), gk.min(), gk.max(), gk.mean()))
    if info["bias"] is not None:
        print("bias voltage: %.4f V" % info["bias"])
    if a.no_plot:
        return
    n2, n3 = int(d[:, 0].max()) + 1, int(d[:, 1].max()) + 1
    plt = get_plt(a.show)
    fig, ax = plt.subplots(figsize=(4.8, 4.0))
    if n2 > 1 and n3 > 1:
        g = np.full((n2, n3), np.nan)
        g[d[:, 0].astype(int), d[:, 1].astype(int)] = gk
        im = ax.imshow(g.T, origin="lower", aspect="auto", cmap="magma",
                       extent=[d[:, 2].min(), d[:, 2].max(), d[:, 3].min(), d[:, 3].max()])
        fig.colorbar(im, ax=ax, label=r"G ($e^2/h$)")
        ax.set_xlabel(r"$k_2$")
        ax.set_ylabel(r"$k_3$")
    else:
        kk = d[:, 2] if n2 > 1 else d[:, 3]
        o = np.argsort(kk)
        ax.plot(kk[o], gk[o], marker="o", ms=3, lw=1.0)
        ax.set_xlabel(r"$k_%d$" % (2 if n2 > 1 else 3))
        ax.set_ylabel(r"G ($e^2/h$)")
    save(fig, plt, a, "conductance_k.png")


def cmd_iv(a):
    import numpy as np
    files = [f for f in expand(a.files) if os.path.isfile(f)]
    if not files:
        die("no .current files found")
    rows = []
    for fn in files:
        v, cur = None, None
        with open(fn, errors="replace") as f:
            for line in f:
                m = re.search(r"Corresponding bias voltage \(Volt\)=\s*([-+0-9.eE]+)", line)
                if m:
                    v = float(m.group(1))
                m = re.search(r"average current for up and down spins \(ampere\)=\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)", line)
                if m:
                    cur = (float(m.group(1)), float(m.group(2)))
                m = re.search(r"average current \(ampere\)=\s*([-+0-9.eE]+)", line)
                if m:
                    cur = (float(m.group(1)), 0.0)
        if v is None or cur is None:
            print("warning: %s has no bias/current header lines; skipped" % fn)
            continue
        rows.append((v, cur[0], cur[1], cur[0] + cur[1], fn))
    if not rows:
        die("no usable .current file")
    rows.sort()
    V = np.array([r[0] for r in rows])
    I = np.array([r[3] for r in rows])
    print("#  V(Volt)      I_up(A)       I_down(A)     I_tot(A)      file")
    for r in rows:
        print("  %8.4f  % .5e  % .5e  % .5e  %s" % r)
    print("# I_tot = I_up + I_down (for non-spin runs OpenMX writes the same value in both columns);")
    print("# sign: bias = mu_L - mu_R as printed in the file header")
    if len(V) > 2 and len(np.unique(V)) == len(V):
        dIdV = np.gradient(I, V)
        print("# differential conductance dI/dV (e^2/h): " +
              " ".join("%.3f" % (g / G_E2H) for g in dIdV))
    else:
        dIdV = None
    if a.csv:
        with open(a.csv, "w") as f:
            f.write("V,I_up_A,I_down_A,I_total_A,dIdV_e2h,file\n")
            for i, r in enumerate(rows):
                f.write("%.6f,%.6e,%.6e,%.6e,%s,%s\n" % (r[0], r[1], r[2], r[3],
                        "%.6f" % (dIdV[i] / G_E2H) if dIdV is not None else "", r[4]))
        print("I-V table written to %s" % a.csv)
    if a.no_plot:
        return
    plt = get_plt(a.show)
    n = 2 if dIdV is not None else 1
    fig, axes = plt.subplots(n, 1, figsize=(5.0, 2.8 * n), sharex=True)
    axes = np.atleast_1d(axes)
    axes[0].plot(V, I * 1e6, marker="o", ms=4, lw=1.2, color="#1f4e79")
    axes[0].axhline(0, color="0.6", lw=0.6)
    axes[0].set_ylabel(r"I ($\mu$A)")
    if dIdV is not None:
        axes[1].plot(V, dIdV / G_E2H, marker="s", ms=3, lw=1.0, color="#c0504d")
        axes[1].set_ylabel(r"dI/dV ($e^2/h$)")
    axes[-1].set_xlabel("bias (V)")
    save(fig, plt, a, "iv.png")


# ----------------------------------------------------------------------------
def get_plt(show):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def save(fig, plt, a, default):
    fig.tight_layout()
    if a.show:
        plt.show()
    else:
        out = a.out or default
        fig.savefig(out, dpi=300)
        print("figure written to %s" % out)


def main():
    ap = argparse.ArgumentParser(
        prog="omx_negf_plot.py",
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", metavar="COMMAND")
    sub.required = True

    def common(p, default_out):
        p.add_argument("-o", "--out", default="", help="figure file (default %s)" % default_out)
        p.add_argument("--no-plot", action="store_true", help="print/CSV only")
        p.add_argument("--show", action="store_true", help="interactive window instead of Agg")

    p = sub.add_parser("trans", help="k-averaged transmission, Landauer G, Seebeck, current")
    p.add_argument("files", nargs="+", help="System.Name (prefix) or System.Name.tran*_* files")
    p.add_argument("--temp", type=float, default=300.0, help="temperature in K (0 = zero-temperature G)")
    p.add_argument("--mu", type=float, nargs="+", default=[0.0],
                   help="chemical potential(s) in eV relative to mu_L (default 0)")
    p.add_argument("--total", action="store_true", help="plot T summed over spins instead of per spin")
    p.add_argument("--emin", type=float, default=None, help="plot range (eV)")
    p.add_argument("--emax", type=float, default=None, help="plot range (eV)")
    p.add_argument("--log", action="store_true", help="logarithmic transmission axis")
    p.add_argument("--csv", default="", help="write E, T_up, T_down, T_total")
    common(p, "transmission.png")
    p.set_defaults(func=cmd_trans)

    p = sub.add_parser("kmap", help="k-resolved transmission at one energy (or T(k,E) for 1D k grids)")
    p.add_argument("files", nargs="+", help="System.Name (prefix) or System.Name.tran*_* files")
    p.add_argument("--energy", type=float, default=0.0, help="E - mu_L in eV (default 0)")
    p.add_argument("--csv", default="", help="write i2, i3, k2, k3, T_total")
    common(p, "transmission_k.png")
    p.set_defaults(func=cmd_kmap)

    p = sub.add_parser("cond", help="read System.Name.conductance and plot the k-resolved conductance")
    p.add_argument("file", help="System.Name.conductance")
    common(p, "conductance_k.png")
    p.set_defaults(func=cmd_cond)

    p = sub.add_parser("iv", help="I-V curve from several System.Name.current files")
    p.add_argument("files", nargs="+", help=".current files of the bias scan (globs allowed)")
    p.add_argument("--csv", default="", help="write V, I_up, I_down, I_total, dI/dV")
    common(p, "iv.png")
    p.set_defaults(func=cmd_iv)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
