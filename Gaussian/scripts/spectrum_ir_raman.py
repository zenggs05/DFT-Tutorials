#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
spectrum_ir_raman.py -- IR and Raman spectra from Gaussian harmonic frequency
                        calculations (Freq, Freq=Raman, CPHF=RdFreq).
                        解析 Gaussian 频率计算的谐振频率、IR 强度和 Raman 活性，
                        做频率缩放与洛伦兹（或高斯）展宽并作图；可按激光波长和
                        温度把 Raman 活性换算为 Raman 强度。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Parsing
-------
The LAST "Harmonic frequencies (cm**-1)" block in normal (not HPModes) format
is read:  Frequencies --, IR Inten --, Raman Activ --, Depolar (P)/(U) --, and
for frequency-dependent (pre-resonance) Raman jobs the lines
"RamAct Fr= n--" / "Dep-P Fr= n--" (select with --raman-fr n).  Values that
overflowed to '*****' become NaN.  Imaginary modes are reported and dropped.

Formulas
--------
* scaled wavenumber   nu_i' = s * nu_i                  (--scale s)
* line shapes are area-normalised:  Lorentzian FWHM g (default 10 cm-1) or
  Gaussian FWHM g (--shape gaussian)
* IR: Gaussian prints Napierian integrated intensities A_i (km/mol);
  1 km/mol = 100 L mol-1 cm-2, hence the decadic molar absorption coefficient
      eps(nu) = (100 / ln 10) * sum_i A_i * L(nu - nu_i')   [L mol-1 cm-1]
* Raman activity S_i (A^4/amu) -> intensity for laser wavenumber nu0 at T:
      I_i = C (nu0 - nu_i)^4 S_i / ( nu_i [1 - exp(-h c nu_i / k T)] ),  C = 1e-12
  (relative units; ref. e.g. Krishnakumar et al., Spectrochim. Acta A 2004)

Usage examples
--------------
  python spectrum_ir_raman.py h2o_freq.log
  python spectrum_ir_raman.py mol.log --scale 0.967 --fwhm 8 --xmin 400 --xmax 4000
  python spectrum_ir_raman.py mol.log --raman-intensity --laser 785 --temp 298.15
  python spectrum_ir_raman.py h2o.log d2o.log --labels H2O D2O --mode ir -o isotopes
  python spectrum_ir_raman.py prer.log --raman-fr 2 --raman-intensity --laser 532

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
TOKEN_RE = re.compile(r"\*{2,}|" + FLOAT)
PROP_RE = re.compile(
    r"^\s*(Frequencies|Red\. masses|Frc consts|IR Inten|Raman Activ|Depolar \(P\)|Depolar \(U\)"
    r"|RamAct Fr=\s*\d+|Dep-P Fr=\s*\d+|Dep-U Fr=\s*\d+)\s*--(?!-)(.*)$")
H_PLANCK = 6.62607015e-34
K_BOLTZ = 1.380649e-23
C_CM = 2.99792458e10
IR_FACTOR = 100.0 / math.log(10.0)
RAMAN_C = 1.0e-12


def safe_stdout() -> None:
    """Never crash on consoles that cannot encode Chinese help text."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass


def read_lines(path: str) -> list[str]:
    """Read a text file; clear error if it does not exist."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"file not found: {path}")
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read().splitlines()


def values(text: str) -> list[float]:
    """Numbers of a property line; '****' overflow -> NaN."""
    out = []
    for tok in TOKEN_RE.findall(text):
        out.append(float("nan") if tok.startswith("*") else float(tok.replace("D", "E")))
    return out


def _parse_block(lines: list[str], start: int) -> list[dict]:
    """Parse one harmonic block starting after header line `start`."""
    modes: list[dict] = []
    group: list[dict] = []
    for i in range(start + 1, len(lines)):
        line = lines[i]
        if "- Thermochemistry -" in line or "Harmonic frequencies (cm**-1)" in line:
            break
        m = PROP_RE.match(line)
        if not m:
            continue
        key = re.sub(r"\s+", "", m.group(1))
        vals = values(m.group(2))
        if key == "Frequencies":
            prev = lines[i - 1].split()
            syms = prev if len(prev) == len(vals) and "--" not in lines[i - 1] else [""] * len(vals)
            group = [{"mode": len(modes) + k + 1, "sym": syms[k], "freq": v, "fr": {}}
                     for k, v in enumerate(vals)]
            modes.extend(group)
            continue
        if len(vals) != len(group):
            continue          # defensive: malformed or wrapped line
        for mode, v in zip(group, vals):
            if key.startswith("RamActFr="):
                mode["fr"].setdefault(int(key.split("=")[1]), {})["raman"] = v
            elif key.startswith("Dep-PFr=") or key.startswith("Dep-UFr="):
                continue
            else:
                mode[key] = v
    return modes


def parse_freq_log(path: str) -> dict:
    """Return the modes of the last normal-format harmonic block of a log."""
    lines = read_lines(path)
    headers = [i for i, ln in enumerate(lines) if "Harmonic frequencies (cm**-1)" in ln]
    if not headers:
        raise ValueError(f"{path}: no 'Harmonic frequencies' block (was Freq run?)")
    for hdr in reversed(headers):
        modes = _parse_block(lines, hdr)
        if modes:
            break
    else:
        raise ValueError(f"{path}: harmonic block without 'Frequencies --' lines")
    tail = [ln for ln in lines[-15:] if ln.strip()]
    return {"file": path, "modes": modes,
            "has_ir": any("IRInten" in m for m in modes),
            "has_raman": any("RamanActiv" in m for m in modes),
            "fr_indices": sorted({k for m in modes for k in m["fr"]}),
            "normal_term": any("Normal termination" in ln for ln in tail[-3:])}


def lineshape(x, x0: float, fwhm: float, shape: str):
    """Area-normalised Lorentzian or Gaussian (1/cm-1)."""
    import numpy as np
    if shape == "lorentzian":
        hw = 0.5 * fwhm
        return hw / math.pi / ((x - x0) ** 2 + hw ** 2)
    sig = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    return np.exp(-0.5 * ((x - x0) / sig) ** 2) / (sig * math.sqrt(2.0 * math.pi))


def raman_intensity(nu: float, activity: float, laser_nm: float, temp: float) -> float:
    """Convert a Raman activity to a (relative) Raman intensity."""
    nu0 = 1.0e7 / laser_nm
    if nu <= 0 or nu0 <= nu or math.isnan(activity):
        return float("nan")
    boltz = 1.0 - math.exp(-H_PLANCK * C_CM * nu / (K_BOLTZ * temp))
    return RAMAN_C * (nu0 - nu) ** 4 * activity / (nu * boltz)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=("IR/Raman spectra from Gaussian Freq logs. "
                     "解析频率、IR 强度、Raman 活性，缩放并展宽作图。"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=("examples:\n"
                "  spectrum_ir_raman.py freq.log\n"
                "  spectrum_ir_raman.py freq.log --scale 0.967 --fwhm 8 --xmin 400 --xmax 4000\n"
                "  spectrum_ir_raman.py freq.log --raman-intensity --laser 785\n"
                "  spectrum_ir_raman.py h2o.log d2o.log --labels H2O D2O --mode ir\n"
                "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License"))
    p.add_argument("logs", nargs="+", help="Gaussian Freq log files / 频率计算输出")
    p.add_argument("--labels", nargs="+", default=None, help="legend labels (one per log)")
    p.add_argument("-s", "--scale", type=float, default=1.0,
                   help="frequency scaling factor (default 1.0) / 频率缩放因子")
    p.add_argument("-w", "--fwhm", type=float, default=10.0,
                   help="full width at half maximum in cm-1 (default 10)")
    p.add_argument("--shape", choices=["lorentzian", "gaussian"], default="lorentzian")
    p.add_argument("--mode", choices=["auto", "ir", "raman", "both"], default="auto",
                   help="which spectra to plot (auto: Raman too if present)")
    p.add_argument("--raman-intensity", action="store_true",
                   help="convert Raman activity to intensity / 活性换算为强度")
    p.add_argument("--laser", type=float, default=532.0,
                   help="laser wavelength in nm for --raman-intensity (default 532)")
    p.add_argument("-t", "--temp", type=float, default=298.15, help="temperature in K")
    p.add_argument("--raman-fr", type=int, default=None,
                   help="use frequency-dependent Raman block 'RamAct Fr= N' (CPHF=RdFreq)")
    p.add_argument("--xmin", type=float, default=0.0, help="lower wavenumber (default 0)")
    p.add_argument("--xmax", type=float, default=None, help="upper wavenumber (default auto)")
    p.add_argument("--step", type=float, default=1.0, help="grid spacing in cm-1 (default 1)")
    p.add_argument("--normalize", action="store_true", help="normalise each curve to max = 1")
    p.add_argument("--sticks", action="store_true", help="draw stick spectra even for overlays")
    p.add_argument("--ir-ascending", action="store_true",
                   help="plot IR with increasing wavenumber (default: decreasing)")
    p.add_argument("-o", "--out", default=None, help="output prefix (default: first log name)")
    p.add_argument("--dpi", type=int, default=300)
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
    if args.fwhm <= 0 or args.step <= 0 or args.scale <= 0:
        sys.exit("ERROR: --fwhm, --step and --scale must be positive")
    if args.labels and len(args.labels) != len(args.logs):
        sys.exit(f"ERROR: --labels needs {len(args.logs)} entries")
    if args.laser <= 0 or args.temp <= 0:
        sys.exit("ERROR: --laser and --temp must be positive")

    try:
        datasets = [parse_freq_log(path) for path in args.logs]
    except (OSError, ValueError) as exc:
        sys.exit(f"ERROR: {exc}")
    labels = args.labels or [os.path.splitext(os.path.basename(p))[0] for p in args.logs]

    for ds, lab in zip(datasets, labels):
        modes = ds["modes"]
        imag = [m for m in modes if m["freq"] < 0]
        ds["real"] = [m for m in modes if m["freq"] > 0]
        print(f"== {lab}: {len(modes)} modes, {len(imag)} imaginary"
              + ("" if ds["normal_term"] else "  [WARNING: no Normal termination]"))
        if imag:
            print("   imaginary modes dropped: " + ", ".join(f"{m['freq']:.1f}" for m in imag))
        raman_key = None
        if args.raman_fr is not None:
            if args.raman_fr not in ds["fr_indices"]:
                sys.exit(f"ERROR: {ds['file']}: no 'RamAct Fr= {args.raman_fr}' block "
                         f"(available: {ds['fr_indices'] or 'none'})")
            raman_key = "fr"
        elif ds["has_raman"]:
            raman_key = "RamanActiv"
        ds["raman_key"] = raman_key
        if args.mode in ("raman", "both") and raman_key is None:
            sys.exit(f"ERROR: {ds['file']}: no Raman activities - run Freq=Raman "
                     "(DFT does not compute them by default)")
        for m in ds["real"]:
            m["scaled"] = m["freq"] * args.scale
            if raman_key == "fr":
                m["raman"] = m["fr"].get(args.raman_fr, {}).get("raman", float("nan"))
            elif raman_key:
                m["raman"] = m.get("RamanActiv", float("nan"))
            else:
                m["raman"] = float("nan")
            m["raman_I"] = raman_intensity(m["scaled"], m["raman"], args.laser, args.temp) \
                if raman_key else float("nan")
        if raman_key and args.raman_intensity:
            nu0 = 1.0e7 / args.laser
            above = [m for m in ds["real"] if m["scaled"] >= nu0]
            if above:
                print(f"   WARNING: {len(above)} mode(s) above the laser wavenumber "
                      f"{nu0:.0f} cm-1 get no Raman intensity")
        print("   mode  sym     freq   scaled   IR(km/mol)  Raman act.(A^4/amu)"
              + ("  Raman I (rel.)" if raman_key and args.raman_intensity else ""))
        for m in ds["real"]:
            ir = m.get("IRInten", float("nan"))
            extra = f"  {m['raman_I']:14.4g}" if raman_key and args.raman_intensity else ""
            print(f"   {m['mode']:4d}  {m['sym']:5s} {m['freq']:8.1f} {m['scaled']:8.1f} "
                  f"{ir:10.4f}  {m['raman']:12.4f}{extra}")

    do_raman = args.mode != "ir" and any(ds["raman_key"] for ds in datasets)
    do_ir = args.mode != "raman"

    all_nu = [m["scaled"] for ds in datasets for m in ds["real"]]
    if not all_nu:
        sys.exit("ERROR: no real vibrational modes found")
    xmax = args.xmax if args.xmax is not None else math.ceil((max(all_nu) + 300) / 100.0) * 100
    if xmax <= args.xmin:
        sys.exit("ERROR: --xmax must be larger than --xmin")
    grid = np.arange(args.xmin, xmax + 0.5 * args.step, args.step)

    curves = []
    for ds in datasets:
        ir = np.zeros_like(grid)
        ram = np.zeros_like(grid)
        for m in ds["real"]:
            shape = lineshape(grid, m["scaled"], args.fwhm, args.shape)
            a_ir = m.get("IRInten", 0.0)
            if not math.isnan(a_ir):
                ir += IR_FACTOR * a_ir * shape
            y = m["raman_I"] if args.raman_intensity else m["raman"]
            if ds["raman_key"] and not math.isnan(y):
                ram += y * shape
        if args.normalize:
            ir = ir / ir.max() if ir.max() > 0 else ir
            ram = ram / ram.max() if ram.max() > 0 else ram
        curves.append((ir, ram))

    prefix = args.out or os.path.splitext(os.path.basename(args.logs[0]))[0]
    modes_csv = prefix + "_modes.csv"
    with open(modes_csv, "w", newline="", encoding="utf-8") as handle:
        wr = csv.writer(handle)
        wr.writerow(["label", "mode", "sym", "freq_cm-1", "scaled_cm-1", "IR_km/mol",
                     "Raman_activity_A4/amu", "Raman_intensity_rel", "depolar_P", "depolar_U"])
        for ds, lab in zip(datasets, labels):
            for m in ds["real"]:
                wr.writerow([lab, m["mode"], m["sym"], m["freq"], f"{m['scaled']:.2f}",
                             m.get("IRInten", ""), "" if math.isnan(m["raman"]) else m["raman"],
                             "" if math.isnan(m["raman_I"]) else f"{m['raman_I']:.6g}",
                             m.get("Depolar(P)", ""), m.get("Depolar(U)", "")])
    spec_csv = prefix + "_spectrum.csv"
    ram_name = "Raman_intensity" if args.raman_intensity else "Raman_activity"
    with open(spec_csv, "w", newline="", encoding="utf-8") as handle:
        wr = csv.writer(handle)
        head = ["wavenumber_cm-1"]
        for lab in labels:
            if do_ir:
                head.append(f"IR_eps_{lab}" if not args.normalize else f"IR_norm_{lab}")
            if do_raman:
                head.append(f"{ram_name}_{lab}")
        wr.writerow(head)
        for k, x in enumerate(grid):
            row = [f"{x:.2f}"]
            for ir, ram in curves:
                if do_ir:
                    row.append(f"{ir[k]:.6g}")
                if do_raman:
                    row.append(f"{ram[k]:.6g}")
            wr.writerow(row)

    panels = [p for p, flag in (("ir", do_ir), ("raman", do_raman)) if flag]
    fig, axes = plt.subplots(len(panels), 1, figsize=(6.4, 2.9 * len(panels) + 0.6), squeeze=False)
    colors = ["#1f4e9c", "#b3261e", "#2e7d32", "#6a1b9a", "#ef6c00", "#00838f", "#5d4037"]
    show_sticks = args.sticks or len(datasets) == 1
    for ax, panel in zip(axes[:, 0], panels):
        stick_ax = ax.twinx() if show_sticks else None
        for k, ((ir, ram), ds, lab) in enumerate(zip(curves, datasets, labels)):
            col = colors[k % len(colors)]
            ycurve = ir if panel == "ir" else ram
            if panel == "raman" and not ds["raman_key"]:
                continue
            ax.plot(grid, ycurve, color=col, lw=1.3, label=lab)
            if stick_ax is not None:
                xs = [m["scaled"] for m in ds["real"]]
                if panel == "ir":
                    ys = [m.get("IRInten", 0.0) for m in ds["real"]]
                else:
                    ys = [m["raman_I"] if args.raman_intensity else m["raman"] for m in ds["real"]]
                ys = [0.0 if math.isnan(y) else y for y in ys]
                stick_ax.vlines(xs, 0, ys, color=col, alpha=0.45, lw=0.9)
        if stick_ax is not None:
            stick_ax.set_ylim(bottom=0)
            stick_ax.set_ylabel("IR intensity (km/mol)" if panel == "ir" else
                                ("Raman intensity (rel.)" if args.raman_intensity
                                 else r"Raman activity ($\AA^4$/amu)"))
        ax.set_ylim(bottom=0)
        if panel == "ir":
            ax.set_ylabel("normalised" if args.normalize else r"$\varepsilon$ (L mol$^{-1}$ cm$^{-1}$)")
            ax.set_title(f"IR (scale {args.scale:g}, {args.shape} FWHM {args.fwhm:g} cm$^{{-1}}$)",
                         fontsize=10)
            ax.set_xlim((args.xmin, xmax) if args.ir_ascending else (xmax, args.xmin))
            ax.set_xlabel(r"wavenumber (cm$^{-1}$)")
        else:
            ylab = "Raman intensity (rel.)" if args.raman_intensity else "Raman activity (broadened)"
            ax.set_ylabel("normalised" if args.normalize else ylab)
            extra = f", laser {args.laser:g} nm, {args.temp:g} K" if args.raman_intensity else ""
            ax.set_title(f"Raman (scale {args.scale:g}{extra})", fontsize=10)
            ax.set_xlim(args.xmin, xmax)
            ax.set_xlabel(r"Raman shift (cm$^{-1}$)")
        if len(datasets) > 1:
            ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    png = prefix + ".png"
    fig.savefig(png, dpi=args.dpi)
    print(f"\nwritten: {modes_csv}, {spec_csv}, {png}")
    if args.show:
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
