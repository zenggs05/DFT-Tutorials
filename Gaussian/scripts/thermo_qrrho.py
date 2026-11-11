#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
thermo_qrrho.py -- teaching-sized quasi-RRHO thermochemistry from Gaussian 16
                   frequency logs (Grimme qRRHO entropy, optional Head-Gordon
                   enthalpy, Truhlar quasi-harmonic, 1 atm -> 1 M correction,
                   frequency scaling, high-level single points, Boltzmann
                   weighting of conformers).
                   从 Gaussian 频率计算 log 出发，教学版重新实现热化学：
                   Grimme qRRHO 熵、可选 Head-Gordon 焓校正、Truhlar 准谐振、
                   1 atm -> 1 mol/L 标准态校正、频率缩放、高精度单点组合、
                   多构象 Boltzmann 加权。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

What it does
------------
1. Reads the LAST harmonic-frequency block and the thermochemistry section of
   a Gaussian 09/16 log (temperature, pressure, molecular mass, rotational
   constants, rotational symmetry number, multiplicity, frequencies).
2. Recomputes the ideal-gas RRHO thermochemistry from scratch (translation,
   rotation, vibration, electronic) and, when the settings equal those of the
   log, compares with Gaussian's own numbers (a self-check you can show to
   students).
3. Applies, on request:
     --qs grimme   Grimme qRRHO entropy: S = sum w*S_harm + (1-w)*S_free-rotor,
                   w = 1/(1+(nu0/nu)^alpha), nu0 = 100 cm-1, alpha = 4,
                   B_av = 1e-44 kg m^2 (Chem. Eur. J. 2012, 18, 9955).
     --qs truhlar  raise all frequencies below the cutoff to the cutoff for the
                   vibrational entropy (J. Phys. Chem. B 2011, 115, 14556).
     --qh          Head-Gordon-type interpolation of the vibrational energy
                   between harmonic oscillator and free rotor (1/2 RT)
                   (J. Phys. Chem. C 2015, 119, 1840).  Here the full
                   vibrational energy of each mode (ZPE + thermal) is
                   interpolated; other programs may differ in detail.
     --freq-scale / --zpe-scale   frequency scaling factors.
     --conc C      standard-state correction to C mol/L:
                   dG = RT ln(C*R*T/P)  (1.89 kcal/mol for 1 M at 298.15 K).
     --sp          replace the electronic energy by a high-level single point
                   (float in Hartree or a log file; the final archive energy is
                   used, e.g. CCSD(T) > MP2/double hybrid > HF/DFT).
     --boltz       Boltzmann weights and ensemble free energy of all logs
                   (treated as conformers of ONE species).
4. Prints a table (Hartree) and optionally writes CSV/JSON.

This is NOT a replacement for GoodVibes or Shermo -- use those for research
and cite them.  Different qRRHO implementations differ slightly (cut-off
function, B_av, treatment of the enthalpy); never mix them within one study.

Usage examples
--------------
  python thermo_qrrho.py h2o_opt_freq.log
  python thermo_qrrho.py ts.log --qs grimme --qh --conc 1.0
  python thermo_qrrho.py react.log --sp react_dsdpbep86.log --freq-scale 0.986
  python thermo_qrrho.py conf1.log conf2.log conf3.log --qs grimme --boltz \
         --csv confs.csv
  python thermo_qrrho.py a.log b.log --sp -382.31234 -382.30998 --temp 333.15
  python thermo_qrrho.py mol.log --qs grimme -v      # per-mode entropy table

Dependencies: Python >= 3.9, numpy
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys

# ----------------------------------------------------------------------------
# Physical constants (CODATA 2018, exact SI values where defined).  Gaussian 16
# uses CODATA 2010 internally; the differences are far below 0.01 kcal/mol.
# ----------------------------------------------------------------------------
H_PLANCK = 6.62607015e-34        # J s
K_BOLTZ = 1.380649e-23           # J/K
N_AVOG = 6.02214076e23           # 1/mol
C_CM = 2.99792458e10             # speed of light, cm/s
AMU_KG = 1.66053906660e-27       # kg
R_GAS = K_BOLTZ * N_AVOG         # J/(mol K)
HARTREE_J = 4.3597447222071e-18  # J
HARTREE_JMOL = HARTREE_J * N_AVOG
CAL = 4.184                      # J
HARTREE_KCAL = HARTREE_JMOL / CAL / 1000.0   # 627.5095 kcal/mol
ATM_PA = 101325.0
R_LATM = 0.0820573661            # L atm /(mol K)
B_AV = 1.0e-44                   # kg m^2, Grimme's average moment of inertia

FLOAT = r"[-+]?\d+\.\d*(?:[EeDd][-+]?\d+)?"


# ----------------------------------------------------------------------------
# Generic helpers
# ----------------------------------------------------------------------------
def _np():
    """Import numpy lazily so that '-h' works even without numpy."""
    try:
        import numpy as np  # noqa: WPS433
    except ImportError:
        sys.exit("ERROR: numpy is required (pip install numpy).")
    return np


def to_float(text: str) -> float:
    """Convert a Gaussian number (may use Fortran 'D' exponent) to float."""
    return float(text.replace("D", "E").replace("d", "e"))


def read_lines(path: str) -> list[str]:
    """Read a text file into lines; give a clear error for missing files."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"file not found: {path}")
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return handle.read().splitlines()


def numbers_in(text: str) -> list[float]:
    """Return all floats (incl. D exponents) found in a string."""
    return [to_float(x) for x in re.findall(FLOAT, text)]


def safe_stdout() -> None:
    """Never crash on consoles that cannot encode Chinese help text."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass


def normal_termination(lines: list[str]) -> bool:
    """True if the last non-empty lines contain 'Normal termination'."""
    tail = [ln for ln in lines[-15:] if ln.strip()]
    return any("Normal termination" in ln for ln in tail[-3:])


# ----------------------------------------------------------------------------
# Parsing of the frequency/thermochemistry part of a Gaussian log
# ----------------------------------------------------------------------------
FREQ_RE = re.compile(r"^\s*Frequencies --(?!-)(.*)$")


def _freqs_after(lines: list[str], start: int) -> tuple[list[float], int]:
    """Collect 'Frequencies --' values after a header until thermochemistry."""
    freqs: list[float] = []
    i = start + 1
    while i < len(lines):
        line = lines[i]
        if "- Thermochemistry -" in line or "Harmonic frequencies (cm**-1)" in line:
            break
        match = FREQ_RE.match(line)
        if match:
            freqs.extend(numbers_in(match.group(1)))
        i += 1
    return freqs, i


def parse_freq_log(path: str) -> dict:
    """Parse frequencies and thermochemistry input data from a Gaussian log.

    The last harmonic block (normal, not HPModes, format) and the
    thermochemistry section that follows it are used, so Opt+Freq,
    --Link1-- jobs and Freq=ReadIsotopes (last isotopologue) all work.
    """
    lines = read_lines(path)
    thermo_all = [i for i, ln in enumerate(lines) if "- Thermochemistry -" in ln]
    if not thermo_all:
        raise ValueError(f"{path}: no '- Thermochemistry -' section - was Freq run "
                         "and did it finish?")
    thermo_start = thermo_all[-1]
    headers = [i for i, ln in enumerate(lines[:thermo_start])
               if "Harmonic frequencies (cm**-1)" in ln]
    freqs: list[float] = []
    for hdr in reversed(headers):
        freqs, _ = _freqs_after(lines, hdr)
        if freqs:
            break
    if headers and not freqs:
        raise ValueError(f"{path}: harmonic block found but no 'Frequencies --' lines")
    # no harmonic block at all -> single atom (Gaussian prints thermo only)

    data: dict = {"file": path, "freqs": freqs, "normal_term": normal_termination(lines)}
    # Multiplicity: last 'Charge = .. Multiplicity = ..' before the thermo section
    mult = None
    for ln in lines[:thermo_start]:
        m = re.search(r"Charge\s*=\s*-?\d+\s+Multiplicity\s*=\s*(\d+)", ln)
        if m:
            mult = int(m.group(1))
    data["mult"] = mult

    patterns = {
        "temp_press": re.compile(rf"Temperature\s+({FLOAT})\s+Kelvin\.\s+Pressure\s+({FLOAT})\s+Atm"),
        "mass": re.compile(rf"Molecular mass:\s+({FLOAT})\s+amu"),
        "sigma": re.compile(r"Rotational symmetry number\s+(\d+)"),
        "rot_ghz": re.compile(r"Rotational constants? \(GHZ\):(.*)$"),
        "rot_k": re.compile(r"Rotational temperatures? \(Kelvin\)(.*)$"),
        "zpe": re.compile(rf"Zero-point correction=\s+({FLOAT})"),
        "e_corr": re.compile(rf"Thermal correction to Energy=\s+({FLOAT})"),
        "h_corr": re.compile(rf"Thermal correction to Enthalpy=\s+({FLOAT})"),
        "g_corr": re.compile(rf"Thermal correction to Gibbs Free Energy=\s+({FLOAT})"),
        "sum_zpe": re.compile(rf"Sum of electronic and zero-point Energies=\s+({FLOAT})"),
        "sum_g": re.compile(rf"Sum of electronic and thermal Free Energies=\s+({FLOAT})"),
        "s_total": re.compile(rf"^\s*Total\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$"),
    }
    end = len(lines)
    for i in range(thermo_start + 1, len(lines)):
        if "Harmonic frequencies (cm**-1)" in lines[i] or "- Thermochemistry -" in lines[i]:
            end = i
            break
    for ln in lines[thermo_start:end]:
        for key, pat in patterns.items():
            if key in data and key not in ("rot_ghz", "rot_k"):
                continue
            m = pat.search(ln)
            if not m:
                continue
            if key == "temp_press":
                data["T_log"] = to_float(m.group(1))
                data["P_log"] = to_float(m.group(2))
                data[key] = True
            elif key == "sigma":
                data[key] = int(m.group(1))
            elif key in ("rot_ghz", "rot_k"):
                if key not in data:   # first occurrence inside the thermo section
                    data[key] = [abs(v) for v in numbers_in(m.group(1))]
            elif key == "s_total":
                data[key] = to_float(m.group(3))
            else:
                data[key] = to_float(m.group(1))
    for key in ("mass", "temp_press"):
        if key not in data:
            raise ValueError(f"{path}: could not read '{key}' from the thermochemistry section")
    data.setdefault("sigma", 1)

    # electronic energy used by Gaussian = (E + ZPE) - ZPE (works for DFT,
    # MP2, double hybrids and TD excited states alike)
    if "sum_zpe" in data and "zpe" in data:
        data["E_elec"] = data["sum_zpe"] - data["zpe"]
    else:
        scf = [ln for ln in lines if "SCF Done:" in ln]
        if not scf:
            raise ValueError(f"{path}: no electronic energy found")
        data["E_elec"] = to_float(re.search(rf"=\s*({FLOAT})", scf[-1]).group(1))
    return data


# ----------------------------------------------------------------------------
# Final (single-point) energies
# ----------------------------------------------------------------------------
ARCHIVE_PRIORITY = ["CCSD(T)", "QCISD(T)", "CCSD", "QCISD", "MP4SDTQ", "MP4SDQ",
                    "MP4DQ", "MP3", "MP2", "HF"]


def archive_fields(lines: list[str]) -> dict:
    """Return the key=value fields of the LAST archive entry (\\ or | separated)."""
    start = None
    for i, ln in enumerate(lines):
        s = ln.strip()
        if s.startswith("1\\1\\") or s.startswith("1|1|"):
            start = i
    if start is None:
        return {}
    chunk = []
    for ln in lines[start:]:
        chunk.append(ln.strip())
        if "@" in ln:
            break
    text = "".join(chunk)
    sep = "\\" if text.startswith("1\\") else "|"
    fields = {}
    for item in text.split(sep):
        if "=" in item:
            key, _, val = item.partition("=")
            fields[key.strip()] = val.strip().rstrip("@")
    return fields


def final_energy(path: str, key: str = "auto") -> tuple[float, str]:
    """Final electronic energy of a single-point log (Hartree) and its label."""
    lines = read_lines(path)
    fields = archive_fields(lines)
    wanted = ARCHIVE_PRIORITY if key == "auto" else [key]
    for name in wanted:
        if name in fields:
            try:
                val = fields[name].split(",")[-1]   # scans: comma-separated list
                return to_float(val), f"archive {name}"
            except ValueError:
                continue
    if key != "auto":
        raise ValueError(f"{path}: archive field '{key}' not found")
    # fall-back: text lines of the log (no archive, e.g. #t or crashed later)
    pats = [
        ("CCSD(T)", re.compile(rf"^\s*CCSD\(T\)=\s*({FLOAT})")),
        ("double hybrid", re.compile(rf"E2\(\S+\)\s*=\s*{FLOAT}\s+E\(\S+\)\s*=\s*({FLOAT})")),
        ("MP2", re.compile(rf"EUMP2\s*=\s*({FLOAT})")),
        ("SCF", re.compile(rf"SCF Done:\s+E\(\S+\)\s*=\s*({FLOAT})")),
    ]
    for label, pat in pats:
        vals = [pat.search(ln) for ln in lines]
        vals = [m for m in vals if m]
        if vals:
            return to_float(vals[-1].group(1)), label
    raise ValueError(f"{path}: no final energy found")


def energy_token(token: str, key: str) -> tuple[float, str]:
    """A token is either a number (Hartree) or a log file name."""
    if os.path.isfile(token):
        return final_energy(token, key)
    try:
        return float(token), "user value"
    except ValueError:
        raise ValueError(f"'{token}' is neither an existing file nor a number") from None


# ----------------------------------------------------------------------------
# Thermochemistry
# ----------------------------------------------------------------------------
def rotational_temperatures(data: dict) -> list[float]:
    """Rotational temperatures (K): from GHz constants if possible (more digits)."""
    ghz = [v for v in data.get("rot_ghz", []) if v > 0]
    if ghz:
        return [H_PLANCK * v * 1.0e9 / K_BOLTZ for v in ghz]
    return [v for v in data.get("rot_k", []) if v > 0]


def harmonic_terms(nu, temp: float):
    """Per-mode harmonic thermal energy (no ZPE) and entropy, SI per mol."""
    np = _np()
    theta = H_PLANCK * C_CM * nu / K_BOLTZ
    x = theta / temp
    emx = np.exp(-x)
    denom = -np.expm1(-x)                         # 1 - exp(-x), stable
    e_th = R_GAS * theta * emx / denom            # J/mol
    s_v = R_GAS * (x * emx / denom - np.log(denom))  # J/(mol K)
    return theta, e_th, s_v


def free_rotor_entropy(nu, temp: float):
    """Grimme free-rotor entropy for each mode (J/(mol K))."""
    np = _np()
    mu = H_PLANCK / (8.0 * math.pi ** 2 * C_CM * nu)         # kg m^2
    mu_eff = mu * B_AV / (mu + B_AV)
    arg = 8.0 * math.pi ** 3 * mu_eff * K_BOLTZ * temp / H_PLANCK ** 2
    return R_GAS * (0.5 + 0.5 * np.log(arg))


def compute_thermo(data: dict, temp: float, p_atm: float, opts) -> dict:
    """Recompute RRHO and quasi-RRHO thermochemistry (results in Hartree)."""
    np = _np()
    mult = opts.mult or data.get("mult") or 1
    rt = R_GAS * temp
    # translation
    mass = data["mass"] * AMU_KG
    q_t = (2.0 * math.pi * mass * K_BOLTZ * temp / H_PLANCK ** 2) ** 1.5 * K_BOLTZ * temp / (p_atm * ATM_PA)
    s_trans = R_GAS * (math.log(q_t) + 2.5)
    e_trans = 1.5 * rt
    # rotation
    thetas = rotational_temperatures(data)
    sigma = data["sigma"]
    if len(thetas) == 0:                 # atom
        s_rot, e_rot, shape = 0.0, 0.0, "atom"
    elif len(thetas) == 1:               # linear
        q_r = temp / (sigma * thetas[0])
        s_rot, e_rot, shape = R_GAS * (math.log(q_r) + 1.0), rt, "linear"
    else:
        q_r = math.sqrt(math.pi) / sigma * temp ** 1.5 / math.sqrt(np.prod(thetas[:3]))
        s_rot, e_rot, shape = R_GAS * (math.log(q_r) + 1.5), 1.5 * rt, "nonlinear"
    # electronic
    s_el = R_GAS * math.log(mult)
    # vibration
    freqs = np.array(data["freqs"], dtype=float)
    imag = freqs[freqs < 0]
    real = freqs[freqs > 0]
    nu = real * opts.freq_scale
    zpe_scale = opts.zpe_scale if opts.zpe_scale is not None else opts.freq_scale
    zpe_modes = 0.5 * R_GAS * H_PLANCK * C_CM * real * zpe_scale / K_BOLTZ
    _, e_th, s_v = harmonic_terms(nu, temp)
    s_r = free_rotor_entropy(nu, temp)
    weight = 1.0 / (1.0 + (opts.cutoff / nu) ** opts.alpha)
    if opts.qs == "grimme":
        s_vib_q = weight * s_v + (1.0 - weight) * s_r
    elif opts.qs == "truhlar":
        _, _, s_v_raised = harmonic_terms(np.maximum(nu, opts.cutoff), temp)
        s_vib_q = s_v_raised
    else:
        s_vib_q = s_v
    u_harm = zpe_modes + e_th
    u_q = weight * u_harm + (1.0 - weight) * 0.5 * rt if opts.qh else u_harm

    zpe = zpe_modes.sum()
    e_corr = e_trans + e_rot + zpe + e_th.sum()
    h_corr = e_corr + rt
    h_corr_q = e_trans + e_rot + u_q.sum() + rt
    s_tot = s_trans + s_rot + s_el + s_v.sum()
    s_tot_q = s_trans + s_rot + s_el + s_vib_q.sum()
    d_conc = rt * math.log(opts.conc * R_LATM * temp / p_atm) if opts.conc else 0.0

    conv = 1.0 / HARTREE_JMOL
    out = {
        "T": temp, "P_atm": p_atm, "mult": mult, "sigma": sigma, "shape": shape,
        "n_imag": int(imag.size), "imag": imag.tolist(), "n_vib": int(real.size),
        "lowest": float(real.min()) if real.size else float("nan"),
        "ZPE": zpe * conv, "E_corr": e_corr * conv, "H_corr": h_corr * conv,
        "qh_H_corr": h_corr_q * conv,
        "TS": temp * s_tot * conv, "qh_TS": temp * s_tot_q * conv,
        "S_cal": {"trans": s_trans / CAL, "rot": s_rot / CAL, "elec": s_el / CAL,
                  "vib": s_v.sum() / CAL, "vib_q": s_vib_q.sum() / CAL,
                  "total": s_tot / CAL, "total_q": s_tot_q / CAL},
        "dG_conc": d_conc * conv,
        "modes": {"nu": nu.tolist(), "w": weight.tolist(), "S_V": (s_v / CAL).tolist(),
                  "S_R": (s_r / CAL).tolist(), "S_q": (s_vib_q / CAL).tolist()},
    }
    out["G_corr"] = out["H_corr"] - out["TS"] + out["dG_conc"]
    out["qh_G_corr"] = out["qh_H_corr"] - out["qh_TS"] + out["dG_conc"]
    return out


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------
def print_mode_table(res: dict) -> None:
    """Per-mode entropy contributions (cal/mol/K)."""
    m = res["modes"]
    print("   mode   nu(scaled)/cm-1   weight w    S_harm    S_rotor   S_used   [cal/(mol K)]")
    for k, (nu, w, sv, sr, sq) in enumerate(zip(m["nu"], m["w"], m["S_V"], m["S_R"], m["S_q"]), 1):
        print(f"   {k:4d}   {nu:14.2f}   {w:8.4f}  {sv:8.3f}  {sr:8.3f}  {sq:8.3f}")


def compare_with_gaussian(data: dict, res: dict) -> None:
    """Show the RRHO self-check against the numbers printed by Gaussian."""
    rows = [("ZPE", "zpe", "ZPE"), ("H_corr", "h_corr", "H_corr"), ("G_corr", "g_corr", "G_corr")]
    if not any(key in data for _, key, _ in rows) and "s_total" not in data:
        return
    print("   self-check vs Gaussian (RRHO, same T/P, no scaling):")
    for label, gkey, rkey in rows:
        if gkey in data:
            mine = res[rkey] - (res["dG_conc"] if rkey == "G_corr" else 0.0)
            diff = (mine - data[gkey]) * HARTREE_KCAL
            print(f"     {label:7s} Gaussian {data[gkey]:12.6f}  this script {mine:12.6f}"
                  f"  diff {diff:+.4f} kcal/mol")
    if "s_total" in data:
        print(f"     S_total Gaussian {data['s_total']:10.3f}  this script "
              f"{res['S_cal']['total']:10.3f} cal/(mol K)")


def boltzmann(values_hartree: list[float], temp: float) -> tuple[list[float], float]:
    """Boltzmann weights and ensemble free energy -RT ln sum exp(-G/RT)."""
    kt = R_GAS * temp / HARTREE_JMOL
    gmin = min(values_hartree)
    terms = [math.exp(-(g - gmin) / kt) for g in values_hartree]
    total = sum(terms)
    return [t / total for t in terms], gmin - kt * math.log(total)


def build_parser() -> argparse.ArgumentParser:
    """Command-line interface."""
    p = argparse.ArgumentParser(
        description=("Quasi-RRHO thermochemistry from Gaussian frequency logs (teaching "
                     "implementation). 教学版 qRRHO 热化学：Grimme/Truhlar 熵校正、"
                     "Head-Gordon 焓校正、1 M 标准态、频率缩放、高精度单点、Boltzmann 加权。"),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=("examples:\n"
                "  thermo_qrrho.py freq.log\n"
                "  thermo_qrrho.py ts.log --qs grimme --qh --conc 1.0\n"
                "  thermo_qrrho.py r.log --sp r_sp.log --freq-scale 0.986\n"
                "  thermo_qrrho.py c1.log c2.log c3.log --qs grimme --boltz --csv out.csv\n"
                "Energies are in Hartree; dG(1 atm->1 M, 298.15 K) = +1.89 kcal/mol.\n"
                "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License"))
    p.add_argument("logs", nargs="+", help="Gaussian Freq log files / 频率计算输出文件")
    p.add_argument("--qs", choices=["grimme", "truhlar", "none"], default="grimme",
                   help="quasi-harmonic entropy model (default grimme) / 低频熵处理方法")
    p.add_argument("--qh", action="store_true",
                   help="also apply Head-Gordon quasi-harmonic enthalpy / 同时做焓的准谐振插值")
    p.add_argument("--cutoff", type=float, default=100.0,
                   help="nu0 cut-off in cm-1 for qRRHO (default 100) / 截断频率")
    p.add_argument("--alpha", type=float, default=4.0, help="damping exponent (default 4)")
    p.add_argument("-t", "--temp", type=float, default=None,
                   help="temperature in K (default: value in the log) / 温度")
    p.add_argument("-p", "--pressure", type=float, default=None,
                   help="pressure in atm (default: value in the log) / 压力")
    p.add_argument("-c", "--conc", type=float, default=None,
                   help="standard-state concentration in mol/L, e.g. 1.0 / 标准态浓度")
    p.add_argument("-f", "--freq-scale", type=float, default=1.0,
                   help="frequency scale factor for vibrational terms (default 1.0)")
    p.add_argument("--zpe-scale", type=float, default=None,
                   help="separate scale factor for ZPE (default = --freq-scale)")
    p.add_argument("--mult", type=int, default=None,
                   help="override spin multiplicity used for S_elec = R ln(mult)")
    p.add_argument("--sp", nargs="+", default=None,
                   help="high-level single points: one per log, Hartree value or log file")
    p.add_argument("--sp-key", default="auto",
                   help="archive key for SP energies: auto (default), HF, MP2, CCSD(T) ...")
    p.add_argument("--boltz", action="store_true",
                   help="Boltzmann-weight all logs as conformers / 构象 Boltzmann 加权")
    p.add_argument("--csv", default=None, help="write results to CSV")
    p.add_argument("--json", default=None, help="write results to JSON")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="print entropy components and per-mode table / 打印逐模式细节")
    return p


def main(argv: list[str] | None = None) -> int:
    safe_stdout()
    args = build_parser().parse_args(argv)
    if args.sp is not None and len(args.sp) != len(args.logs):
        sys.exit(f"ERROR: --sp needs {len(args.logs)} values (one per log), got {len(args.sp)}")
    if args.conc is not None and args.conc <= 0:
        sys.exit("ERROR: --conc must be > 0")
    if args.cutoff <= 0:
        sys.exit("ERROR: --cutoff must be > 0")

    print(f"# thermo_qrrho.py  qs={args.qs}  qh={'head-gordon' if args.qh else 'off'}  "
          f"cutoff={args.cutoff:g} cm-1  freq-scale={args.freq_scale:g}"
          + (f"  zpe-scale={args.zpe_scale:g}" if args.zpe_scale else "")
          + (f"  conc={args.conc:g} M" if args.conc else "  standard state: pressure of the log"))
    results = []
    for k, path in enumerate(args.logs):
        try:
            data = parse_freq_log(path)
            temp = args.temp if args.temp is not None else data["T_log"]
            press = args.pressure if args.pressure is not None else data["P_log"]
            res = compute_thermo(data, temp, press, args)
            if args.sp is not None:
                e_used, src = energy_token(args.sp[k], args.sp_key)
            else:
                e_used, src = data["E_elec"], "freq log"
        except (OSError, ValueError) as exc:
            sys.exit(f"ERROR: {exc}")
        res.update({"file": path, "E": e_used, "E_source": src, "E_freq_log": data["E_elec"]})
        res["H"] = e_used + res["H_corr"]
        res["qh_H"] = e_used + res["qh_H_corr"]
        res["G"] = e_used + res["G_corr"]
        res["qh_G"] = e_used + res["qh_G_corr"]
        results.append(res)

        name = os.path.basename(path)
        print(f"\n== {name}: T = {temp:.2f} K, P = {press:.5g} atm, mult = {res['mult']}, "
              f"sigma_rot = {res['sigma']}, {res['shape']}, {res['n_vib']} real modes")
        if not data["normal_term"]:
            print("   WARNING: last job step did not end with 'Normal termination'.")
        if res["n_imag"]:
            freqs_txt = ", ".join(f"{v:.1f}" for v in res["imag"])
            print(f"   NOTE: {res['n_imag']} imaginary frequency(ies) ({freqs_txt} cm-1) "
                  "excluded from thermochemistry (1 is expected for a TS only).")
        if res["n_vib"] and res["lowest"] < args.cutoff:
            n_low = sum(1 for v in res["modes"]["nu"] if v < args.cutoff)
            print(f"   {n_low} mode(s) below {args.cutoff:g} cm-1 (lowest {res['lowest']:.1f}).")
        print(f"   E used: {e_used:.6f} Eh ({src})"
              + (f";  E(freq log) = {data['E_elec']:.6f}" if args.sp else ""))
        if args.verbose:
            s = res["S_cal"]
            print(f"   S [cal/(mol K)]: trans {s['trans']:.3f}  rot {s['rot']:.3f}  "
                  f"elec {s['elec']:.3f}  vib(RRHO) {s['vib']:.3f}  vib(used) {s['vib_q']:.3f}")
            print(f"                    total RRHO {s['total']:.3f}  total quasi {s['total_q']:.3f}")
            print_mode_table(res)
        same_settings = (abs(temp - data["T_log"]) < 1e-3 and abs(press - data["P_log"]) < 1e-6
                         and args.freq_scale == 1.0 and args.zpe_scale in (None, 1.0)
                         and (args.mult is None or args.mult == data.get("mult")))
        if same_settings:
            compare_with_gaussian(data, res)
        if args.conc:
            print(f"   dG(standard state -> {args.conc:g} M) = "
                  f"{res['dG_conc'] * HARTREE_KCAL:+.3f} kcal/mol (included in G)")

    # summary table
    print("\n   (Hartree)      E            ZPE        H          T.S      T.qh-S"
          "      G(T)          qh-G(T)")
    for r in results:
        print(f"o {os.path.basename(r['file'])[:24]:24s} {r['E']:14.6f} {r['ZPE']:9.6f} "
              f"{r['H']:14.6f} {r['TS']:9.6f} {r['qh_TS']:9.6f} {r['G']:14.6f} {r['qh_G']:14.6f}")
    if len(results) > 1:
        gmin = min(r["qh_G"] for r in results)
        print("\n   relative qh-G (kcal/mol; meaningful only for isomers/conformers):")
        for r in results:
            print(f"   {os.path.basename(r['file']):30s} {(r['qh_G'] - gmin) * HARTREE_KCAL:8.2f}")
    if args.boltz:
        temps = {round(r["T"], 6) for r in results}
        if len(temps) > 1:
            sys.exit("ERROR: --boltz needs the same temperature for all logs (use --temp).")
        weights, g_ens = boltzmann([r["qh_G"] for r in results], results[0]["T"])
        print(f"\n   Boltzmann weights at {results[0]['T']:.2f} K (from qh-G):")
        for r, w in zip(results, weights):
            r["boltz_weight"] = w
            print(f"   {os.path.basename(r['file']):30s} {w:8.4f}")
        g_avg = sum(w * r["qh_G"] for r, w in zip(results, weights))
        print(f"   weighted average qh-G          = {g_avg:.6f} Eh")
        print(f"   ensemble qh-G (-RT ln sum e^-G/RT) = {g_ens:.6f} Eh "
              f"(includes conformational entropy, {(g_ens - g_avg) * HARTREE_KCAL:+.3f} kcal/mol)")

    cols = ["file", "E_source", "T", "P_atm", "mult", "n_imag", "E", "ZPE", "H_corr",
            "qh_H_corr", "TS", "qh_TS", "dG_conc", "G_corr", "qh_G_corr", "H", "qh_H", "G",
            "qh_G", "boltz_weight"]
    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(cols)
            for r in results:
                writer.writerow([r.get(c, "") for c in cols])
        print(f"\nCSV written: {args.csv}")
    if args.json:
        with open(args.json, "w", encoding="utf-8") as handle:
            json.dump(results, handle, indent=2)
        print(f"JSON written: {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
