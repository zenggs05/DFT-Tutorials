#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spectrum_uvvis_ecd.py -- UV-Vis absorption and ECD spectra from Gaussian
                         TD-DFT / TDA / CIS logs.
                         解析 Gaussian TD-DFT 激发态（eV、nm、振子强度 f、
                         旋转强度 R），在能量 (eV) 空间做高斯展宽，输出 CSV
                         并绘制 UV-Vis 与 ECD 谱；支持能量平移、对映体翻转、
                         多构象 Boltzmann 加和、叠加实验谱。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Parsing
-------
* The LAST "Excitation energies and oscillator strengths:" block is used, so
  Opt TD(Root=n), --Link1-- jobs and TD=(50-50) outputs work.  Lines look like
      Excited State   1:      Singlet-BU     5.3351 eV  232.39 nm  f=0.1707  <S**2>=0.000
  (open-shell states carry labels such as "2.031-A").
* Rotatory strengths are read from the tables printed just before that block:
      Rotatory Strengths (R) in cgs (10**-40 erg-esu-cm/Gauss)
            state   XX   YY   ZZ   R(velocity)   E-M Angle
  R(velocity) (origin independent, default) or R(length) (--rot length).

Band shapes (Gaussian white paper "Creating UV/Visible Plots", GaussView
default sigma = 0.4 eV; sigma is the half-width at 1/e of the maximum)
-------------------------------------------------------------------------
  eps(E)  = sum_i 1.3062974e8 * f_i / sigma[cm-1] * exp(-((E - E_i)/sigma)^2)
  deps(E) = sum_i E_i * R_i / (2.296e-39 * sqrt(pi) * sigma) * exp(-((E - E_i)/sigma)^2)
  (eps, deps in L mol-1 cm-1; R_i in cgs; E_i/sigma is unit-free).
  FWHM = 2*sqrt(ln 2)*sigma = 1.665*sigma  (use --fwhm to give a FWHM instead).

Usage examples
--------------
  python spectrum_uvvis_ecd.py coumarin_td.log
  python spectrum_uvvis_ecd.py td.log --sigma 0.3 --shift -0.2 --xmin 200 --xmax 500
  python spectrum_uvvis_ecd.py conf1_td.log conf2_td.log --energies c1_freq.log c2_freq.log \
         --flip --exp-ecd exp_ecd.txt -o methyloxirane
  python spectrum_uvvis_ecd.py a.log b.log --weights 0.7 0.3 --xunit eV --no-ecd

Dependencies: Python >= 3.9, numpy, matplotlib
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import sys

FLOAT = r"[-+]?\d+\.\d*(?:[EeDd][-+]?\d+)?"
EV_TO_CM = 8065.543937          # 1 eV in cm-1
EV_NM = 1239.841984             # E(eV) * lambda(nm)
HARTREE_KCAL = 627.509474
HARTREE_KJ = 2625.499639
HARTREE_EV = 27.211386246
KB_HARTREE = 3.166811563e-6     # Boltzmann constant in Hartree/K
UV_CONST = 1.3062974e8          # Gaussian white paper constant
ECD_CONST = 22.96               # 2.296e-39 expressed in units of 1e-40 cgs

STATE_RE = re.compile(
    rf"^\s*Excited State\s+(\d+):\s+(\S+)\s+({FLOAT})\s+eV\s+({FLOAT})\s+nm\s+f=\s*({FLOAT})"
    rf"(?:\s+<S\*\*2>=\s*({FLOAT}))?")
ROT_ROW_RE = re.compile(rf"^\s*(\d+)\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})")
STOP_MARKERS = ("Population analysis", "Leave Link", "SavETr", "Normal termination",
                "Error termination", "Excitation energies and oscillator strengths")


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def safe_stdout() -> None:
    """Never crash on consoles that cannot encode Chinese help text."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass


def to_float(text: str) -> float:
    """Gaussian numbers may use a Fortran 'D' exponent."""
    return float(text.replace("D", "E").replace("d", "e"))


def read_lines(path: str) -> list[str]:
    """Read a text file; clear error if it does not exist."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"file not found: {path}")
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read().splitlines()


# ----------------------------------------------------------------------------
# parsing
# ----------------------------------------------------------------------------
def _rot_table(lines: list[str], start: int, stop: int, column: str) -> dict:
    """Last rotatory-strength table with header containing `column` in [start, stop)."""
    found = None
    for i in range(start, stop - 1):
        if "Rotatory Strengths (R)" in lines[i] and column in lines[i + 1]:
            found = i + 2
    if found is None:
        return {}
    table = {}
    for line in lines[found:stop]:
        m = ROT_ROW_RE.match(line)
        if not m:
            break
        table[int(m.group(1))] = to_float(m.group(5))
    return table


def parse_td_log(path: str) -> dict:
    """Excited states and rotatory strengths of the last TD block of a log."""
    lines = read_lines(path)
    headers = [i for i, ln in enumerate(lines) if "Excitation energies and oscillator strengths" in ln]
    if not headers:
        raise ValueError(f"{path}: no 'Excitation energies and oscillator strengths' block "
                         "(not a TD/CIS/TDA job, or it crashed before printing)")
    last = headers[-1]
    prev = headers[-2] if len(headers) > 1 else 0
    states = []
    for line in lines[last + 1:]:
        if any(marker in line for marker in STOP_MARKERS):
            break
        m = STATE_RE.match(line)
        if not m:
            continue
        num = int(m.group(1))
        if states and num <= states[-1]["state"]:
            break
        states.append({"state": num, "label": m.group(2), "E_eV": to_float(m.group(3)),
                       "nm": to_float(m.group(4)), "f": to_float(m.group(5)),
                       "S2": to_float(m.group(6)) if m.group(6) else None})
    if not states:
        raise ValueError(f"{path}: TD header found but no 'Excited State' lines could be parsed")
    rvel = _rot_table(lines, prev, last, "R(velocity)")
    rlen = _rot_table(lines, prev, last, "R(length)")
    for st in states:
        st["R_vel"] = rvel.get(st["state"])
        st["R_len"] = rlen.get(st["state"])
    tail = [ln for ln in lines[-15:] if ln.strip()]
    return {"file": path, "states": states, "has_rot": bool(rvel or rlen),
            "normal_term": any("Normal termination" in ln for ln in tail[-3:])}


def energy_from_token(token: str, unit: str) -> tuple[float, str]:
    """Boltzmann energy: a number (in `unit`) or a log (G if Freq was run, else E)."""
    if os.path.isfile(token):
        lines = read_lines(token)
        gvals = [ln for ln in lines if "Sum of electronic and thermal Free Energies=" in ln]
        if gvals:
            return to_float(gvals[-1].split("=")[-1]), "G(freq log)"
        scf = [re.search(rf"=\s*({FLOAT})", ln) for ln in lines if "SCF Done:" in ln]
        scf = [m for m in scf if m]
        if scf:
            return to_float(scf[-1].group(1)), "SCF energy (no Freq data!)"
        raise ValueError(f"{token}: no free energy or SCF energy found")
    try:
        value = float(token)
    except ValueError:
        raise ValueError(f"'{token}' is neither an existing file nor a number") from None
    factor = {"hartree": 1.0, "kcal": 1.0 / HARTREE_KCAL, "kj": 1.0 / HARTREE_KJ,
              "ev": 1.0 / HARTREE_EV}[unit]
    return value * factor, f"user value ({unit})"


def boltzmann_weights(energies_hartree: list[float], temp: float) -> list[float]:
    """Normalised Boltzmann weights."""
    emin = min(energies_hartree)
    terms = [math.exp(-(e - emin) / (KB_HARTREE * temp)) for e in energies_hartree]
    total = sum(terms)
    return [t / total for t in terms]


# ----------------------------------------------------------------------------
# spectra
# ----------------------------------------------------------------------------
def broaden(grid_ev, states: list[dict], sigma_uv: float, sigma_ecd: float, rot_key: str):
    """Gaussian-broadened eps(E) and delta-eps(E) on an eV grid."""
    import numpy as np
    eps = np.zeros_like(grid_ev)
    deps = np.zeros_like(grid_ev)
    sigma_cm = sigma_uv * EV_TO_CM
    for st in states:
        e_i = st["E_eV"]
        if e_i <= 0:
            continue
        eps += UV_CONST * st["f"] / sigma_cm * np.exp(-((grid_ev - e_i) / sigma_uv) ** 2)
        r_i = st.get(rot_key)
        if r_i:
            deps += (e_i / sigma_ecd) * r_i / (ECD_CONST * math.sqrt(math.pi)) * \
                np.exp(-((grid_ev - e_i) / sigma_ecd) ** 2)
    return eps, deps


def to_unit(e_ev, unit: str):
    """Convert energies in eV to the requested x unit."""
    import numpy as np
    e_ev = np.asarray(e_ev, dtype=float)
    if unit == "nm":
        return EV_NM / e_ev
    if unit == "cm-1":
        return e_ev * EV_TO_CM
    return e_ev


def from_unit(value: float, unit: str) -> float:
    """Convert a value in the x unit back to eV."""
    if unit == "nm":
        return EV_NM / value
    if unit == "cm-1":
        return value / EV_TO_CM
    return value


def read_xy(path: str) -> tuple[list[float], list[float]]:
    """Two-column text/CSV file (comments with #, header lines skipped)."""
    xs, ys = [], []
    for line in read_lines(path):
        parts = re.split(r"[,\s;]+", line.strip())
        if len(parts) < 2 or line.lstrip().startswith("#"):
            continue
        try:
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
        except ValueError:
            continue
    if not xs:
        raise ValueError(f"{path}: no numeric two-column data")
    return xs, ys


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=("UV-Vis and ECD spectra from Gaussian TD-DFT logs. "
                     "解析 TD-DFT 激发态并展宽为 UV-Vis/ECD 谱（eV 空间高斯线型）。"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=("examples:\n"
                "  spectrum_uvvis_ecd.py td.log\n"
                "  spectrum_uvvis_ecd.py td.log --sigma 0.3 --shift -0.2 --xmin 200 --xmax 500\n"
                "  spectrum_uvvis_ecd.py c1.log c2.log --energies c1_freq.log c2_freq.log --flip\n"
                "  spectrum_uvvis_ecd.py c1.log c2.log --weights 0.7 0.3 --xunit eV\n"
                "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License"))
    p.add_argument("logs", nargs="+", help="TD-DFT log files (conformers) / TD 计算输出")
    width = p.add_mutually_exclusive_group()
    width.add_argument("--sigma", type=float, default=0.4,
                       help="Gaussian half-width at 1/e in eV (default 0.4, GaussView default)")
    width.add_argument("--fwhm", type=float, default=None, help="full width at half maximum in eV")
    p.add_argument("--sigma-ecd", type=float, default=None,
                   help="separate width (1/e half-width, eV) for ECD (default = UV width)")
    p.add_argument("--shift", type=float, default=0.0,
                   help="shift all excitation energies by this many eV (e.g. -0.2) / 能量平移")
    p.add_argument("--flip", action="store_true",
                   help="flip the sign of R (mirror image / other enantiomer) / 对映体翻转")
    p.add_argument("--rot", choices=["velocity", "length"], default="velocity",
                   help="rotatory strength gauge (default velocity)")
    boltz = p.add_mutually_exclusive_group()
    boltz.add_argument("--weights", nargs="+", type=float, default=None,
                       help="explicit conformer weights (normalised automatically)")
    boltz.add_argument("--energies", nargs="+", default=None,
                       help="conformer energies: numbers or logs (G if Freq present) / 构象能量")
    p.add_argument("--energy-unit", choices=["hartree", "kcal", "kj", "ev"], default="hartree",
                   help="unit of numeric --energies (default hartree)")
    p.add_argument("-t", "--temp", type=float, default=298.15, help="temperature for Boltzmann (K)")
    p.add_argument("--xunit", choices=["nm", "eV", "cm-1"], default="nm", help="x axis unit")
    p.add_argument("--xmin", type=float, default=None, help="x range minimum (in --xunit)")
    p.add_argument("--xmax", type=float, default=None, help="x range maximum (in --xunit)")
    p.add_argument("--npts", type=int, default=2000, help="grid points (default 2000)")
    p.add_argument("--no-ecd", action="store_true", help="do not plot/write ECD")
    p.add_argument("--exp-uv", default=None, help="experimental UV file (x y), scaled to max")
    p.add_argument("--exp-ecd", default=None, help="experimental ECD file (x y), scaled to max")
    p.add_argument("--exp-unit", choices=["nm", "eV", "cm-1"], default="nm",
                   help="x unit of experimental files (default nm)")
    p.add_argument("-o", "--out", default=None, help="output prefix (default: first log name)")
    p.add_argument("--dpi", type=int, default=300, help="figure resolution")
    p.add_argument("--show", action="store_true", help="show the figure interactively")
    return p


def main(argv: list[str] | None = None) -> int:
    safe_stdout()
    args = build_parser().parse_args(argv)
    try:
        import numpy as np
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        sys.exit(f"ERROR: {exc}. Install numpy and matplotlib.")

    nconf = len(args.logs)
    sigma = args.fwhm / (2.0 * math.sqrt(math.log(2.0))) if args.fwhm else args.sigma
    sigma_ecd = args.sigma_ecd if args.sigma_ecd else sigma
    if sigma <= 0 or sigma_ecd <= 0:
        sys.exit("ERROR: widths must be positive")
    rot_key = "R_vel" if args.rot == "velocity" else "R_len"
    if args.xmin is not None and args.xmax is not None and args.xmin >= args.xmax:
        sys.exit("ERROR: --xmin must be smaller than --xmax")
    for value in (args.xmin, args.xmax):
        if value is not None and value <= 0:
            sys.exit("ERROR: x limits must be positive")
    if args.npts < 10:
        sys.exit("ERROR: --npts must be >= 10")

    try:
        exp_uv = read_xy(args.exp_uv) if args.exp_uv else None
        exp_ecd = read_xy(args.exp_ecd) if args.exp_ecd else None
        confs = [parse_td_log(path) for path in args.logs]
        if args.weights is not None:
            if len(args.weights) != nconf or min(args.weights) < 0 or sum(args.weights) <= 0:
                sys.exit(f"ERROR: --weights needs {nconf} non-negative numbers")
            weights = [w / sum(args.weights) for w in args.weights]
            wsrc = "user weights"
        elif args.energies is not None:
            if len(args.energies) != nconf:
                sys.exit(f"ERROR: --energies needs {nconf} values (one per log)")
            pairs = [energy_from_token(tok, args.energy_unit) for tok in args.energies]
            weights = boltzmann_weights([e for e, _ in pairs], args.temp)
            wsrc = f"Boltzmann {args.temp:.2f} K from " + ", ".join(sorted({s for _, s in pairs}))
        else:
            weights = [1.0 / nconf] * nconf
            wsrc = "equal weights" if nconf > 1 else "single structure"
    except (OSError, ValueError) as exc:
        sys.exit(f"ERROR: {exc}")

    # apply shift / flip
    for conf in confs:
        for st in conf["states"]:
            st["E_eV"] += args.shift
            st["nm"] = EV_NM / st["E_eV"] if st["E_eV"] > 0 else float("nan")
            if args.flip:
                for key in ("R_vel", "R_len"):
                    if st[key] is not None:
                        st[key] = -st[key]
    do_ecd = not args.no_ecd and any(c["has_rot"] for c in confs)
    if not args.no_ecd and not do_ecd:
        print("NOTE: no rotatory strengths in the log(s) -> ECD skipped.")
    if do_ecd and not all(c["has_rot"] for c in confs):
        print("WARNING: some logs lack rotatory strengths; they contribute 0 to ECD.")

    # report
    print(f"# width: sigma = {sigma:.4f} eV (FWHM {sigma * 2 * math.sqrt(math.log(2)):.4f} eV)"
          + (f", ECD sigma = {sigma_ecd:.4f} eV" if do_ecd else "")
          + f"; shift = {args.shift:+.3f} eV; flip = {args.flip}; weights: {wsrc}")
    for conf, w in zip(confs, weights):
        sts = conf["states"]
        bright = max(sts, key=lambda s: s["f"])
        print(f"\n== {os.path.basename(conf['file'])}: {len(sts)} states, weight {w:.4f}"
              + ("" if conf["normal_term"] else "  [WARNING: no Normal termination]"))
        print(f"   brightest: state {bright['state']} {bright['E_eV']:.4f} eV "
              f"{bright['nm']:.1f} nm f={bright['f']:.4f}")
        print("   state  label          E/eV     nm        f        R_vel     R_len")
        for st in sts:
            rv = "" if st["R_vel"] is None else f"{st['R_vel']:9.4f}"
            rl = "" if st["R_len"] is None else f"{st['R_len']:9.4f}"
            print(f"   {st['state']:5d}  {st['label']:12s} {st['E_eV']:8.4f} {st['nm']:8.2f} "
                  f"{st['f']:8.4f} {rv:>9s} {rl:>9s}")
        neg = [st for st in sts if st["E_eV"] <= 0]
        if neg:
            print(f"   WARNING: {len(neg)} state(s) with E <= 0 eV ignored (unstable reference?)")
        if do_ecd and conf["has_rot"] and all(abs(st.get(rot_key) or 0) < 1e-4 for st in sts):
            print("   NOTE: all rotatory strengths are ~0 (achiral or symmetric structure).")

    # grid in eV
    all_e = [st["E_eV"] for c in confs for st in c["states"] if st["E_eV"] > 0]
    if not all_e:
        sys.exit("ERROR: no state with positive excitation energy")
    smax = max(sigma, sigma_ecd)
    lo_ev = max(0.05, min(all_e) - 5 * smax)
    hi_ev = max(all_e) + 5 * smax
    if args.xunit == "nm":            # short wavelength = high energy
        if args.xmin is not None:
            hi_ev = EV_NM / args.xmin
        if args.xmax is not None:
            lo_ev = EV_NM / args.xmax
    else:
        if args.xmin is not None:
            lo_ev = from_unit(args.xmin, args.xunit)
        if args.xmax is not None:
            hi_ev = from_unit(args.xmax, args.xunit)
    if hi_ev <= lo_ev:
        sys.exit("ERROR: empty x range - check --xmin/--xmax against the state energies")
    grid = np.linspace(lo_ev, hi_ev, args.npts)
    eps_tot = np.zeros_like(grid)
    deps_tot = np.zeros_like(grid)
    per_conf = []
    for conf, w in zip(confs, weights):
        eps, deps = broaden(grid, conf["states"], sigma, sigma_ecd, rot_key)
        per_conf.append((eps, deps))
        eps_tot += w * eps
        deps_tot += w * deps

    prefix = args.out or os.path.splitext(os.path.basename(args.logs[0]))[0]
    states_csv = prefix + "_states.csv"
    with open(states_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["file", "weight", "state", "label", "E_eV", "lambda_nm", "f",
                         "R_velocity_1e-40cgs", "R_length_1e-40cgs", "S2"])
        for conf, w in zip(confs, weights):
            for st in conf["states"]:
                writer.writerow([os.path.basename(conf["file"]), f"{w:.6f}", st["state"], st["label"],
                                 f"{st['E_eV']:.5f}", f"{st['nm']:.3f}", st["f"],
                                 "" if st["R_vel"] is None else st["R_vel"],
                                 "" if st["R_len"] is None else st["R_len"],
                                 "" if st["S2"] is None else st["S2"]])
    spec_csv = prefix + "_spectrum.csv"
    with open(spec_csv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        head = ["E_eV", "lambda_nm", "wavenumber_cm-1", "epsilon_L/(mol*cm)"]
        if do_ecd:
            head.append("delta_epsilon_L/(mol*cm)")
        if nconf > 1:
            for conf in confs:
                base = os.path.splitext(os.path.basename(conf["file"]))[0]
                head.append(f"epsilon_{base}")
                if do_ecd:
                    head.append(f"delta_epsilon_{base}")
        writer.writerow(head)
        for k, e in enumerate(grid):
            row = [f"{e:.5f}", f"{EV_NM / e:.3f}", f"{e * EV_TO_CM:.2f}", f"{eps_tot[k]:.6g}"]
            if do_ecd:
                row.append(f"{deps_tot[k]:.6g}")
            if nconf > 1:
                for eps, deps in per_conf:
                    row.append(f"{eps[k]:.6g}")
                    if do_ecd:
                        row.append(f"{deps[k]:.6g}")
            writer.writerow(row)

    # plot
    x = to_unit(grid, args.xunit)
    nrows = 2 if do_ecd else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(6.4, 3.0 * nrows + 0.6), sharex=True, squeeze=False)
    ax_uv = axes[0, 0]
    ax_uv.plot(x, eps_tot, color="#1f4e9c", lw=1.6, label="calc.")
    stick_ax = ax_uv.twinx()
    for conf, w in zip(confs, weights):
        es = [st["E_eV"] for st in conf["states"] if st["E_eV"] > 0]
        fs = [st["f"] * w for st in conf["states"] if st["E_eV"] > 0]
        stick_ax.vlines(to_unit(es, args.xunit), 0, fs, color="#7f7f7f", lw=0.9)
    stick_ax.set_ylabel("oscillator strength f" + (" (weighted)" if nconf > 1 else ""))
    stick_ax.set_ylim(bottom=0)
    ax_uv.set_ylabel(r"$\varepsilon$ (L mol$^{-1}$ cm$^{-1}$)")
    ax_uv.set_ylim(bottom=0)
    if exp_uv:
        ex, ey = exp_uv
        ex_ev = [from_unit(v, args.exp_unit) for v in ex]
        scale = eps_tot.max() / max(abs(v) for v in ey) if max(abs(v) for v in ey) > 0 else 1.0
        ax_uv.plot(to_unit(ex_ev, args.xunit), [v * scale for v in ey], "k--", lw=1.0,
                   label="exp. (scaled)")
    ax_uv.legend(loc="upper right", frameon=False)
    ax_uv.set_title(f"UV-Vis (sigma = {sigma:.2f} eV, shift {args.shift:+.2f} eV)", fontsize=10)
    if do_ecd:
        ax_ecd = axes[1, 0]
        ax_ecd.axhline(0, color="k", lw=0.6)
        ax_ecd.plot(x, deps_tot, color="#b3261e", lw=1.6, label="calc." + (" (flipped)" if args.flip else ""))
        stick2 = ax_ecd.twinx()
        for conf, w in zip(confs, weights):
            es = [st["E_eV"] for st in conf["states"] if st["E_eV"] > 0 and st.get(rot_key) is not None]
            rs = [st[rot_key] * w for st in conf["states"] if st["E_eV"] > 0 and st.get(rot_key) is not None]
            if es:
                stick2.vlines(to_unit(es, args.xunit), 0, rs, color="#7f7f7f", lw=0.9)
        stick2.set_ylabel(r"R (10$^{-40}$ cgs)")
        lim = max(abs(stick2.get_ylim()[0]), abs(stick2.get_ylim()[1]))
        stick2.set_ylim(-lim, lim)
        ax_ecd.set_ylabel(r"$\Delta\varepsilon$ (L mol$^{-1}$ cm$^{-1}$)")
        if exp_ecd:
            ex, ey = exp_ecd
            ex_ev = [from_unit(v, args.exp_unit) for v in ex]
            ymax = max(abs(v) for v in ey)
            scale = np.abs(deps_tot).max() / ymax if ymax > 0 else 1.0
            ax_ecd.plot(to_unit(ex_ev, args.xunit), [v * scale for v in ey], "k--", lw=1.0,
                        label="exp. (scaled)")
        ylim = np.abs(deps_tot).max() * 1.1 or 1.0
        ax_ecd.set_ylim(-ylim, ylim)
        ax_ecd.legend(loc="upper right", frameon=False)
        ax_ecd.set_title(f"ECD ({args.rot} gauge, sigma = {sigma_ecd:.2f} eV)", fontsize=10)
    xlabel = {"nm": "wavelength (nm)", "eV": "excitation energy (eV)", "cm-1": r"wavenumber (cm$^{-1}$)"}
    axes[-1, 0].set_xlabel(xlabel[args.xunit])
    xs = sorted([float(x[0]), float(x[-1])])
    if args.xmin is not None:
        xs[0] = args.xmin
    if args.xmax is not None:
        xs[1] = args.xmax
    axes[-1, 0].set_xlim(xs[0], xs[1])
    fig.tight_layout()
    png = prefix + ".png"
    fig.savefig(png, dpi=args.dpi)
    print(f"\nwritten: {states_csv}, {spec_csv}, {png}")
    if args.show:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
