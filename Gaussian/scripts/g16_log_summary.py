#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
g16_log_summary.py - summarize Gaussian 16 output files (.log/.out) into a table, CSV or JSON.

批量汇总 Gaussian 16 输出文件：结束状态与出错的 Link、route、最终能量、优化是否收敛、
虚频、ZPE/H/G 热校正、<S**2>（湮灭前/后）、基函数数目、CPU 与墙钟时间。

Usage examples
--------------
  python g16_log_summary.py water.log
  python g16_log_summary.py *.log --csv summary.csv --json summary.json
  python g16_log_summary.py calc_dir/ -r --sort G --relative
  python g16_log_summary.py ts.log --cclib        # optional cross-check with cclib

Notes
-----
* One pass over each file (works for multi-GB logs); multi-step (--Link1--) jobs are
  supported: status counts every job step, energies are taken from the LAST step that
  printed an energy, thermochemistry from the last frequency step.
* "e_final" is the most specific energy of the last step:
  Counterpoise corrected > ONIOM extrapolated > composite (0 K) > CCSD(T) > MP2/double
  hybrid > TD-DFT total > SCF Done. Its type is reported in "e_final_type".
* All energies are in Hartree; relative energies (--relative) in kcal/mol
  (1 Eh = 627.509474 kcal/mol).

Dependencies: Python >= 3.9 standard library; cclib (optional, only for --cclib).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import collections
import csv
import glob
import json
import math
import os
import re
import sys

HARTREE_TO_KCAL = 627.509474

SYMBOLS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
    "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl "
    "Mc Lv Ts Og"
).split()

# --------------------------------------------------------------------------------------
# regular expressions for the Gaussian 16 output format
# --------------------------------------------------------------------------------------
RE_VERSION = re.compile(r"Gaussian (\d\d):\s+(\S+)\s+(\S+)")
RE_REV = re.compile(r"Rev([A-Z]\.\d+)")
RE_CHARGE = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
RE_NBASIS = re.compile(r"^\s*(\d+)\s+basis functions,")
RE_NELEC = re.compile(r"^\s*(\d+)\s+alpha electrons\s+(\d+)\s+beta electrons")
RE_SCF = re.compile(r"SCF Done:\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+)")
RE_MP2 = re.compile(r"EUMP2\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_DH = re.compile(r"E2\(([^)]+)\)\s*=\s*\S+\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_CCSDT = re.compile(r"\bCCSD\(T\)\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_ONIOM = re.compile(r"ONIOM: extrapolated energy\s*=\s*(-?\d+\.\d+)")
RE_CP = re.compile(r"Counterpoise corrected energy\s*=\s*(-?\d+\.\d+)")
RE_BSSE = re.compile(r"BSSE energy\s*=\s*(-?\d+\.\d+)")
RE_TD = re.compile(r"Total Energy, E\(([^)]*)\)\s*=\s*(-?\d+\.\d+)")
RE_COMPOSITE = re.compile(r"^\s*(G4MP2|G4|G3MP2|G3B3|G3|G2MP2|G2|CBS-QB3|CBS-APNO|CBS-4M?|W1BD|W1U|W1RO)"
                          r"\s*\(0 K\)\s*=\s*(-?\d+\.\d+)")
RE_EXCITED = re.compile(r"^ Excited State\s+\d+:")
RE_STEP = re.compile(r"Step number\s+(\d+)\s+out of a maximum of\s+(\d+)")
RE_IRC_PT0 = re.compile(r"Point Number:\s+0\s")
RE_IRC_STEP = re.compile(r"^\s*Pt\s+\d+\s+Step number")
RE_OPT_KW = re.compile(r"(^|[\s#/,])opt(?=$|[\s=(,])", re.IGNORECASE)
RE_IRC_KW = re.compile(r"(^|[\s#/,])irc(?=$|[\s=(,])", re.IGNORECASE)
RE_SCANPT = re.compile(r"on scan point\s+(\d+)\s+out of\s+(\d+)")
RE_CONV_ITEM = re.compile(r"^\s*(Maximum Force|RMS\s+Force|Maximum Displacement|RMS\s+Displacement)"
                          r"\s+(\S+)\s+(\S+)\s+(YES|NO)")
RE_FREQ = re.compile(r"^\s*Frequencies --\s+(.*)$")
RE_TEMP = re.compile(r"Temperature\s+(\d+\.\d+)\s+Kelvin\.\s+Pressure\s+(\d+\.\d+)\s+Atm")
RE_S2_ANN = re.compile(r"S\*\*2 before annihilation\s+(-?\d+\.\d+),\s+after\s+(-?\d+\.\d+)")
RE_S2 = re.compile(r"<S\*\*2>=\s*(-?\d+\.\d+)")
RE_TIME = re.compile(r"(Job cpu time|Elapsed time):\s+(\d+)\s+days\s+(\d+)\s+hours\s+(\d+)\s+minutes\s+"
                     r"(\d+(?:\.\d*)?)\s+seconds")
RE_ERRLINK = re.compile(r"Error termination via Lnk1e\b(?:\s+in\s+(\S+))?")
RE_LINKNUM = re.compile(r"l(\d+)\.exe", re.IGNORECASE)
RE_ERRREQ = re.compile(r"Error termination request processed by link\s+(\d+)")
RE_PG = re.compile(r"Full point group\s+(\S+)")
RE_DASH = re.compile(r"^ -{10,}\s*$")

THERMO_KEYS = (
    ("Zero-point correction=", "zpe"),
    ("Thermal correction to Energy=", "tc_energy"),
    ("Thermal correction to Enthalpy=", "tc_enthalpy"),
    ("Thermal correction to Gibbs Free Energy=", "tc_gibbs"),
    ("Sum of electronic and zero-point Energies=", "e_zpe"),
    ("Sum of electronic and thermal Energies=", "e_thermal"),
    ("Sum of electronic and thermal Enthalpies=", "h_sum"),
    ("Sum of electronic and thermal Free Energies=", "g_sum"),
)

# short list of well-known error lines; g16_error_doctor.py has the full catalogue
KNOWN_ERRORS = (
    "Convergence failure -- run terminated",
    "Number of steps exceeded",
    "corrector steps exceded",
    "is impossible",
    "Atomic number out of range",
    "End of file in ZSymb",
    "End of file in GetChg",
    "Wanted an integer as input",
    "could not allocate memory",
    "Out-of-memory error",
    "Erroneous write",
    "Erroneous read",
    "Wrong number of Negative eigenvalues",
    "SYNTAX ERROR WAS DETECTED",
    "Illegal IType or MSType",
    "FormBX had a problem",
    "Linear angle in",
    "Small interatomic distances",
    "Problem with the distance matrix",
    "Unable to choose the S8 parameter",
    "Inaccurate quadrature",
    "Inv3 failed in PCMMkU",
    "segmentation violation",
    "End of file reading",
)
KNOWN_ERRORS_LOW = tuple(p.lower() for p in KNOWN_ERRORS)

# order of preference for "e_final"
ENERGY_PRIORITY = ("cp", "oniom", "composite", "ccsdt", "mp2", "td", "scf")
ENERGY_LABEL = {"cp": "Counterpoise", "oniom": "ONIOM", "composite": "composite(0K)",
                "ccsdt": "CCSD(T)", "mp2": "MP2/DH", "td": "TD", "scf": "SCF"}

FIELDS = [
    "file", "status", "n_jobs", "n_normal", "error_link", "error_msg", "version", "route",
    "charge", "mult", "natoms", "formula", "nbasis", "point_group",
    "scf_method", "e_scf", "e_final", "e_final_type", "bsse",
    "opt_steps", "opt_status", "max_force", "rms_force", "max_disp", "rms_disp",
    "scan_points", "freq_check", "irc_points", "irc_status", "n_freq", "n_imag", "imag_freqs", "lowest_freq", "stationary",
    "temperature", "pressure", "zpe", "tc_energy", "tc_enthalpy", "tc_gibbs",
    "e_zpe", "e_thermal", "h_sum", "g_sum",
    "s2_before", "s2_after", "s2_ideal", "s2_dev_pct", "n_excited",
    "cpu_time_s", "wall_time_s",
]

DEFAULT_COLUMNS = ["file", "status", "error_link", "e_final", "e_final_type", "opt_status",
                   "n_imag", "g_sum", "s2_before", "wall_time_s"]


def dfloat(text):
    """Convert a Fortran number such as -0.76228D+02 to float."""
    return float(text.replace("D", "E").replace("d", "e"))


def to_float(text):
    try:
        return dfloat(text)
    except ValueError:
        return float("nan")      # e.g. '********' printed by Fortran for overflow


def hill_formula(znums):
    counts = collections.Counter(SYMBOLS[z] if 0 < z < len(SYMBOLS) else "X" for z in znums if z >= 0)
    if not counts:
        return ""
    order = []
    if "C" in counts:
        order.append("C")
        if "H" in counts:
            order.append("H")
    order += sorted(k for k in counts if k not in order)
    return "".join(k + (str(counts[k]) if counts[k] > 1 else "") for k in order)


def read_orientation(lines_iter):
    """Read the atom rows of an orientation block; the header line was already consumed."""
    znums = []
    dashes = 0
    for line in lines_iter:
        if RE_DASH.match(line):
            dashes += 1
            if dashes == 3:
                break
            continue
        if dashes == 2:
            tok = line.split()
            if len(tok) >= 5:
                try:
                    znums.append(int(tok[1]))
                except ValueError:
                    pass
    return znums


def classify_route(route):
    """Return 'opt', 'irc' or 'other' for a route line."""
    if RE_IRC_KW.search(route):
        return "irc"
    if RE_OPT_KW.search(route):
        return "opt"
    return "other"


def empty_record(path):
    rec = {k: None for k in FIELDS}
    rec["file"] = path
    return rec


def parse_log(path):
    """Parse one Gaussian output file and return a dict with the FIELDS keys."""
    rec = empty_record(path)
    n_jobs = 1
    n_normal = 0
    error_link = None
    routes = []
    capturing_route = False
    route_buf = []
    prev = ""
    energies = {}            # type -> (value, job index)
    scf_method = None
    opt_steps = 0
    opt_events = []          # sequence of 'step' / 'done' / 'stopped'
    job_kind = "other"       # 'opt', 'irc' or 'other' (from the route of the current job step)
    freq_check = None        # convergence test printed by a Freq job that follows an Opt
    irc_points = 0
    irc_status = None
    conv = {}
    scan_total = None
    scan_done = 0
    freq_lists = []
    freq_pending_new = True
    thermo = {}
    s2_before = s2_after = s2_plain = None
    n_excited = 0
    cpu = 0.0
    wall = 0.0
    seen_wall = seen_cpu = False
    last_known_err = None
    line_no = 0
    err_line_no = None
    known_err_hits = []
    znums = []

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        it = iter(fh)
        for line in it:
            line_no += 1
            line = line.rstrip("\r\n")

            # ---- route section: ' #...' right after a dashed line, until the next dashed line
            if capturing_route:
                if RE_DASH.match(line):
                    capturing_route = False
                    routes.append("".join(route_buf).strip())
                    job_kind = classify_route(routes[-1])
                    route_buf = []
                else:
                    route_buf.append(line[1:] if line.startswith(" ") else line)
                prev = line
                continue
            if line.startswith(" #") and RE_DASH.match(prev):
                capturing_route = True
                route_buf = [line[1:]]
                prev = line
                continue
            prev = line

            if "Proceeding to internal job step number" in line:
                n_jobs += 1
                freq_pending_new = True
                continue
            if "Normal termination of Gaussian" in line:
                n_normal += 1
                continue
            m = RE_ERRLINK.search(line)
            if m:
                mm = RE_LINKNUM.search(m.group(1) or "")
                error_link = "L" + mm.group(1) if mm else (error_link or "L?")
                err_line_no = line_no
                continue
            m = RE_ERRREQ.search(line)
            if m:
                error_link = error_link or "L" + m.group(1)
                err_line_no = err_line_no or line_no
                continue
            if rec["version"] is None:
                m = RE_VERSION.search(line)
                if m:
                    mm = RE_REV.search(m.group(2))
                    rec["version"] = "G%s %s" % (m.group(1), mm.group(1) if mm else m.group(2))
                    continue
            if ("orientation:" in line) and ("Standard" in line or "Input" in line or "Z-Matrix" in line):
                z = read_orientation(it)
                if z:
                    znums = z
                continue
            m = RE_SCF.search(line)
            if m:
                scf_method = m.group(1)
                energies["scf"] = (float(m.group(2)), n_jobs)
                continue
            if rec["charge"] is None:
                m = RE_CHARGE.search(line)
                if m:
                    rec["charge"], rec["mult"] = int(m.group(1)), int(m.group(2))
                    continue
            m = RE_NBASIS.match(line)
            if m:
                rec["nbasis"] = int(m.group(1))
                continue
            m = RE_PG.search(line)
            if m:
                rec["point_group"] = m.group(1)
                continue
            if "EUMP2" in line:
                m = RE_MP2.search(line)
                if m:
                    energies["mp2"] = (dfloat(m.group(1)), n_jobs)
                continue
            if "E2(" in line:
                m = RE_DH.search(line)
                if m:
                    energies["mp2"] = (dfloat(m.group(3)), n_jobs)
                continue
            if "CCSD(T)" in line:
                m = RE_CCSDT.search(line)
                if m:
                    energies["ccsdt"] = (dfloat(m.group(1)), n_jobs)
                    continue
            if "ONIOM: extrapolated energy" in line:
                m = RE_ONIOM.search(line)
                if m:
                    energies["oniom"] = (float(m.group(1)), n_jobs)
                continue
            if "Counterpoise corrected energy" in line:
                m = RE_CP.search(line)
                if m:
                    energies["cp"] = (float(m.group(1)), n_jobs)
                continue
            if "BSSE energy" in line:
                m = RE_BSSE.search(line)
                if m:
                    rec["bsse"] = float(m.group(1))
                continue
            if "Total Energy, E(" in line:
                m = RE_TD.search(line)
                if m:
                    energies["td"] = (float(m.group(2)), n_jobs)
                continue
            if "(0 K)" in line:
                m = RE_COMPOSITE.search(line)
                if m:
                    energies["composite"] = (float(m.group(2)), n_jobs)
                continue
            if RE_EXCITED.match(line):
                n_excited += 1
                continue
            m = RE_STEP.search(line)
            if m:
                if RE_IRC_STEP.match(line) or job_kind == "irc":
                    continue                    # IRC corrector steps are not optimization steps
                if job_kind != "opt":
                    freq_check = "NOT converged"   # Freq after Opt re-tests convergence
                    continue
                opt_steps += 1
                opt_events.append("step")
                ms = RE_SCANPT.search(line)
                if ms:
                    scan_total = int(ms.group(2))
                continue
            m = RE_CONV_ITEM.match(line)
            if m:
                if job_kind == "opt":
                    key = {"Maximum Force": "max_force", "RMS Force": "rms_force",
                           "Maximum Displacement": "max_disp", "RMS Displacement": "rms_disp"}
                    name = re.sub(r"\s+", " ", m.group(1))
                    conv[key[name]] = to_float(m.group(2))
                continue
            if "Stationary point found" in line:
                if job_kind != "opt":
                    freq_check = "converged"
                    continue
                opt_events.append("done")
                if scan_total:
                    scan_done += 1
                continue
            if "Optimization stopped" in line:
                opt_events.append("stopped")
                continue
            if "Point Number:" in line and "Path Number:" in line:
                if not RE_IRC_PT0.search(line):
                    irc_points += 1
                irc_status = irc_status or "running"
                continue
            if "Reaction path calculation complete" in line:
                irc_status = "complete"
                continue
            if "Harmonic frequencies (cm**-1)" in line:
                freq_pending_new = True
                continue
            m = RE_FREQ.match(line)
            if m:
                if freq_pending_new:
                    freq_lists.append([])
                    freq_pending_new = False
                freq_lists[-1].extend(to_float(x) for x in m.group(1).split())
                continue
            m = RE_TEMP.search(line)
            if m:
                rec["temperature"], rec["pressure"] = float(m.group(1)), float(m.group(2))
                continue
            for label, key in THERMO_KEYS:
                if label in line:
                    try:
                        thermo[key] = float(line.split("=")[1].split()[0])
                    except (IndexError, ValueError):
                        pass
                    break
            if "S**2 before annihilation" in line:
                m = RE_S2_ANN.search(line)
                if m:
                    s2_before, s2_after = float(m.group(1)), float(m.group(2))
                continue
            if "<S**2>=" in line:
                m = RE_S2.search(line)
                if m:
                    s2_plain = float(m.group(1))
                continue
            m = RE_TIME.search(line)
            if m:
                sec = (int(m.group(2)) * 86400 + int(m.group(3)) * 3600 + int(m.group(4)) * 60
                       + float(m.group(5)))
                if m.group(1).startswith("Job"):
                    cpu += sec
                    seen_cpu = True
                else:
                    wall += sec
                    seen_wall = True
                continue
            low = line.lower()
            for pat in KNOWN_ERRORS_LOW:
                if pat in low:
                    known_err_hits.append((line_no, line.strip()))
                    break

    if capturing_route and route_buf:          # truncated file
        routes.append("".join(route_buf).strip())

    # ---- status
    rec["n_jobs"] = n_jobs
    rec["n_normal"] = n_normal
    if error_link:
        rec["status"] = "error"
        rec["error_link"] = error_link
        if known_err_hits:
            before = [h for h in known_err_hits if err_line_no is None or h[0] <= err_line_no]
            last_known_err = (before or known_err_hits)[-1][1]
        rec["error_msg"] = last_known_err
    elif n_normal >= n_jobs and n_normal > 0:
        rec["status"] = "normal"
    else:
        rec["status"] = "incomplete"   # still running, killed (time limit, node failure) or truncated

    user_routes = [r for r in routes if not r.upper().startswith("#N GEOM=ALLCHECK GUESS=TCHECK")]
    rec["route"] = " | ".join(user_routes or routes) or None

    if znums:
        rec["natoms"] = len(znums)
        rec["formula"] = hill_formula(znums)

    # ---- energies (last job step that printed any energy)
    rec["scf_method"] = scf_method
    if "scf" in energies:
        rec["e_scf"] = energies["scf"][0]
    if energies:
        last_job = max(v[1] for v in energies.values())
        for key in ENERGY_PRIORITY:
            if key in energies and energies[key][1] == last_job:
                rec["e_final"] = energies[key][0]
                rec["e_final_type"] = ENERGY_LABEL[key]
                break

    # ---- optimization
    if opt_steps:
        rec["opt_steps"] = opt_steps
        last_evt = opt_events[-1]
        if scan_total:
            rec["scan_points"] = "%d/%d" % (scan_done, scan_total)
        if last_evt == "done":
            rec["opt_status"] = "converged"
        elif last_evt == "stopped" or "stopped" in opt_events:
            rec["opt_status"] = "stopped"
        else:
            rec["opt_status"] = "unfinished"
        for k in ("max_force", "rms_force", "max_disp", "rms_disp"):
            rec[k] = conv.get(k)
    rec["freq_check"] = freq_check
    if irc_points:
        rec["irc_points"] = irc_points
        rec["irc_status"] = irc_status

    # ---- frequencies (last frequency calculation)
    freqs = [f for f in (freq_lists[-1] if freq_lists else []) if not math.isnan(f)]
    if freqs:
        imag = [f for f in freqs if f < 0]
        rec["n_freq"] = len(freqs)
        rec["n_imag"] = len(imag)
        rec["imag_freqs"] = " ".join("%.1f" % f for f in imag) if imag else ""
        rec["lowest_freq"] = min(freqs)
        rec["stationary"] = {0: "minimum", 1: "TS"}.get(len(imag), "saddle(%d)" % len(imag))
    rec.update(thermo)

    # ---- spin
    if s2_before is None and s2_plain is not None:
        s2_before = s2_plain
    rec["s2_before"], rec["s2_after"] = s2_before, s2_after
    if rec["mult"]:
        s = (rec["mult"] - 1) / 2.0
        rec["s2_ideal"] = s * (s + 1)
        if s2_before is not None and rec["s2_ideal"] > 0:
            rec["s2_dev_pct"] = round(100.0 * (s2_before - rec["s2_ideal"]) / rec["s2_ideal"], 1)
    rec["n_excited"] = n_excited or None
    rec["cpu_time_s"] = round(cpu, 1) if seen_cpu else None
    rec["wall_time_s"] = round(wall, 1) if seen_wall else None
    return rec


# --------------------------------------------------------------------------------------
# optional cclib cross-check
# --------------------------------------------------------------------------------------
def cclib_check(path, rec):
    try:
        import cclib  # noqa: F401  (optional dependency)
        from cclib.io import ccread
    except ImportError:
        return "cclib not installed (pip install cclib)"
    try:
        data = ccread(path)
    except Exception as exc:  # cclib raises many different exceptions on odd files
        return "cclib failed: %s" % exc
    if data is None:
        return "cclib could not parse"
    notes = []
    if getattr(data, "scfenergies", None) is not None and len(data.scfenergies) and rec["e_scf"] is not None:
        cc = float(data.scfenergies[-1])
        # cclib <= 1.8 reports eV; later versions may report Hartree - accept both
        if abs(cc - rec["e_scf"]) < 1e-5 or abs(cc / 27.211386245988 - rec["e_scf"]) < 1e-5:
            notes.append("E ok")
        else:
            notes.append("E differs (cclib %.6f)" % cc)
    if getattr(data, "vibfreqs", None) is not None and rec["n_imag"] is not None:
        n_cc = int(sum(1 for f in data.vibfreqs if f < 0))
        notes.append("nimag ok" if n_cc == rec["n_imag"] else "nimag differs (cclib %d)" % n_cc)
    return "; ".join(notes) or "nothing to compare"


# --------------------------------------------------------------------------------------
# input collection and output
# --------------------------------------------------------------------------------------
def collect_files(paths, recursive, exts):
    files = []
    for p in paths:
        if any(ch in p for ch in "*?["):       # Windows cmd does not expand wildcards
            matched = sorted(glob.glob(p))
            if not matched:
                print("WARNING: no file matches %s" % p, file=sys.stderr)
            files.extend(matched)
        elif os.path.isdir(p):
            pattern = "**" if recursive else ""
            for ext in exts:
                files.extend(sorted(glob.glob(os.path.join(p, pattern, "*" + ext), recursive=recursive)))
        elif os.path.isfile(p):
            files.append(p)
        else:
            print("WARNING: %s not found, skipped" % p, file=sys.stderr)
    seen = set()
    unique = []
    for f in files:
        key = os.path.abspath(f)
        if key not in seen:
            seen.add(key)
            unique.append(f)
    return unique


def fmt(value, key):
    if value is None:
        return "-"
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if key in ("e_final", "e_scf", "g_sum", "h_sum", "e_zpe", "e_thermal"):
            return "%.6f" % value
        if key in ("dE_kcal", "dG_kcal"):
            return "%.2f" % value
        if key in ("wall_time_s", "cpu_time_s"):
            return human_time(value)
        if key in ("s2_before", "s2_after"):
            return "%.4f" % value
        return "%.6g" % value
    return str(value)


def human_time(sec):
    sec = int(round(sec))
    h, rem = divmod(sec, 3600)
    m, s = divmod(rem, 60)
    return "%d:%02d:%02d" % (h, m, s)


def print_table(records, columns):
    rows = []
    for r in records:
        row = []
        for c in columns:
            v = r.get(c)
            if c == "file":
                v = os.path.basename(v)
            row.append(fmt(v, c))
        rows.append(row)
    widths = [max(len(c), *(len(row[i]) for row in rows)) if rows else len(c) for i, c in enumerate(columns)]
    widths = [min(w, 48) for w in widths]
    print("  ".join(c.ljust(w) for c, w in zip(columns, widths)))
    print("  ".join("-" * w for w in widths))
    for row in rows:
        print("  ".join((v if len(v) <= w else v[: w - 2] + "..").ljust(w) for v, w in zip(row, widths)))


def add_relative(records):
    for key, out in (("e_final", "dE_kcal"), ("g_sum", "dG_kcal")):
        vals = [r[key] for r in records if r.get(key) is not None]
        ref = min(vals) if vals else None
        for r in records:
            r[out] = (r[key] - ref) * HARTREE_TO_KCAL if (ref is not None and r.get(key) is not None) else None


def json_safe(obj):
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj


def build_parser():
    p = argparse.ArgumentParser(
        description="Summarize Gaussian 16 .log/.out files (status, energies, opt, freq, thermo, <S**2>, "
                    "timing).  批量汇总 Gaussian 输出：结束状态、能量、优化、虚频、热化学、自旋污染、耗时。",
        epilog="Examples:\n"
               "  python g16_log_summary.py job.log\n"
               "  python g16_log_summary.py logs/ -r --csv all.csv --json all.json\n"
               "  python g16_log_summary.py conf*.log --sort G --relative\n"
               "  python g16_log_summary.py a.log --columns file,status,e_final,zpe,g_sum,n_imag\n"
               "Available columns: " + ", ".join(FIELDS + ["dE_kcal", "dG_kcal", "cclib"]),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paths", nargs="+", help="log files, wildcards or directories / 输出文件、通配符或目录")
    p.add_argument("-r", "--recursive", action="store_true", help="search directories recursively / 递归搜索子目录")
    p.add_argument("--ext", default=".log,.out",
                   help="extensions searched in directories (default: .log,.out) / 目录中搜索的扩展名")
    p.add_argument("--csv", metavar="FILE", help="write all fields to a CSV file / 输出 CSV")
    p.add_argument("--json", metavar="FILE", help="write all fields to a JSON file / 输出 JSON")
    p.add_argument("--columns", help="comma-separated columns for the screen table / 屏幕表格显示的列")
    p.add_argument("--sort", choices=["name", "energy", "G", "status"], default="name",
                   help="sort order (default: name) / 排序方式")
    p.add_argument("--relative", action="store_true",
                   help="add dE_kcal and dG_kcal relative to the lowest value / 给出相对能量 (kcal/mol)")
    p.add_argument("--cclib", action="store_true",
                   help="cross-check SCF energy and imaginary-frequency count with cclib (optional) / 用 cclib 交叉验证")
    p.add_argument("-q", "--quiet", action="store_true", help="do not print the table / 不打印表格")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    exts = [e if e.startswith(".") else "." + e for e in args.ext.split(",") if e]
    files = collect_files(args.paths, args.recursive, exts)
    if not files:
        print("ERROR: no Gaussian output files found.", file=sys.stderr)
        return 1

    records = []
    for f in files:
        try:
            rec = parse_log(f)
        except OSError as exc:
            print("ERROR: cannot read %s: %s" % (f, exc), file=sys.stderr)
            continue
        if rec["version"] is None and rec["e_scf"] is None and rec["n_normal"] == 0:
            print("WARNING: %s does not look like a Gaussian output file" % f, file=sys.stderr)
        if args.cclib:
            rec["cclib"] = cclib_check(f, rec)
        records.append(rec)
    if not records:
        return 1

    if args.relative:
        add_relative(records)

    big = float("inf")
    if args.sort == "energy":
        records.sort(key=lambda r: r["e_final"] if r["e_final"] is not None else big)
    elif args.sort == "G":
        records.sort(key=lambda r: r["g_sum"] if r["g_sum"] is not None else big)
    elif args.sort == "status":
        records.sort(key=lambda r: (r["status"] != "error", r["status"] != "incomplete", r["file"]))
    else:
        records.sort(key=lambda r: r["file"])

    fields = list(FIELDS)
    if args.relative:
        fields += ["dE_kcal", "dG_kcal"]
    if args.cclib:
        fields.append("cclib")

    if not args.quiet:
        if args.columns:
            columns = [c.strip() for c in args.columns.split(",") if c.strip()]
            bad = [c for c in columns if c not in fields]
            if bad:
                print("ERROR: unknown column(s): %s" % ", ".join(bad), file=sys.stderr)
                return 2
        else:
            columns = list(DEFAULT_COLUMNS)
            if args.relative:
                columns[4:4] = ["dE_kcal"]
                columns.insert(columns.index("g_sum") + 1, "dG_kcal")
            if args.cclib:
                columns.append("cclib")
        print_table(records, columns)
        n_err = sum(r["status"] == "error" for r in records)
        n_inc = sum(r["status"] == "incomplete" for r in records)
        print("\n%d file(s): %d normal, %d error, %d incomplete" %
              (len(records), len(records) - n_err - n_inc, n_err, n_inc))
        bad_spin = [os.path.basename(r["file"]) for r in records
                    if (r.get("s2_dev_pct") is not None and abs(r["s2_dev_pct"]) > 10)
                    or (r.get("s2_ideal") == 0 and (r.get("s2_before") or 0) > 0.1)]
        if bad_spin:
            print("NOTE: <S**2> deviates > 10%% from S(S+1) (or > 0.1 for singlets: broken symmetry?) in: %s"
                  % ", ".join(bad_spin))
        bad_check = [os.path.basename(r["file"]) for r in records if r.get("freq_check") == "NOT converged"]
        if bad_check:
            print("NOTE: the Freq step re-tested the optimized structure and found it NOT converged in: %s\n"
                  "      (stationary point not fully converged; consider Opt=(ReadFC or CalcFC, Tight) "
                  "and re-run Freq)" % ", ".join(bad_check))
        irc = ["%s (%s points, %s)" % (os.path.basename(r["file"]), r["irc_points"], r["irc_status"])
               for r in records if r.get("irc_points")]
        if irc:
            print("IRC: " + "; ".join(irc))

    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for r in records:
                w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in fields})
        print("CSV written to %s" % args.csv, file=sys.stderr)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump([{k: json_safe(r.get(k)) for k in fields} for r in records], fh, indent=2)
        print("JSON written to %s" % args.json, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
