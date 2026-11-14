#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cp2k_outparse.py - zero-dependency summary of CP2K output files.

一次性汇总 CP2K 输出（.out）中的关键信息：是否正常结束、SCF 是否收敛与步数、总能量、
费米能 / HOMO-LUMO 能隙、Mulliken 电荷（含自旋）、WARNING 汇总、ABORT 信息、内存与耗时；
GEO_OPT / CELL_OPT 还会逐步列出 OPT| 块中的 Max step、RMS step、Max grad、RMS grad、
内压与压力偏差。可输出 JSON 或 CSV，适合批量检查。

Usage examples
--------------
  python cp2k_outparse.py h2o_energy.out                 # human-readable summary
  python cp2k_outparse.py h2o_geo_opt_bfgs.out --opt     # + table of every OPT| step
  python cp2k_outparse.py fe_afm.out --mulliken          # + Mulliken charges / spin moments
  python cp2k_outparse.py run.out --warnings 20          # list up to 20 distinct warnings
  python cp2k_outparse.py */output.out --brief           # one line per file (batch check)
  python cp2k_outparse.py si_cell_opt.out --json > si.json
  python cp2k_outparse.py si_cell_opt.out --opt-csv opt_steps.csv

Exit status: 0 if every file ended normally ("PROGRAM ENDED AT") without an ABORT,
1 otherwise (handy in shell loops: `python cp2k_outparse.py x.out --brief || echo check x.out`).

What is parsed (strings taken from the CP2K 2026.2 sources; earlier 2023.x-2025.x
versions print the same lines unless noted)
-------------------------------------------------------------------------
* "CP2K| version string", "GLOBAL| Project name", "GLOBAL| Run type",
  "GLOBAL| Total number of message passing processes", "GLOBAL| Number of threads ..."
* "- Atoms:" in the TOTAL NUMBERS AND MAXIMUM NUMBERS block
* "*** SCF run converged in N steps ***", "outer SCF loop converged in ...",
  "outer SCF loop FAILED to converge", "SCF run NOT converged" (warning or ABORT)
* "ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]" (2025.1+) and
  "... energy [a.u.]:" (2024.x and older); the unit can be changed by
  &FORCE_EVAL/&PRINT/&PROGRAM_RUN_INFO/ENERGY_UNIT and is read from the line
* "Fermi energy:" (Hartree, printed with smearing), "MO| E(Fermi):", "MO| Band gap:",
  "HOMO - LUMO gap [eV] :"
* "Mulliken Population Analysis" (restricted and spin-polarized layouts)
* "*** WARNING in <file>:<line> :: <message> ***" (multi-line messages are joined),
  "The number of warnings for this run is"
* the [ABORT] box, "MEMORY| Estimated peak process memory [MiB]",
  the timing table (row "CP2K", total time, maximum over ranks), "PROGRAM ENDED AT"
* OPT| blocks of GEO_OPT / CELL_OPT (CP2K 2023.1 and newer); the older
  "Informations at step" layout is parsed on a best-effort basis
* "GEOMETRY OPTIMIZATION COMPLETED", "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED",
  the L-BFGS specific convergence message, "STRESS| 1/3 Trace", "CELL| Volume"

Dependencies: Python >= 3.9, standard library only.

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import csv
import json
import os
import re
import sys
from collections import OrderedDict

HARTREE_EV = 27.211386245988  # CODATA 2018

FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"

RE_VERSION = re.compile(r"CP2K\|\s+version string:\s+(.+?)\s*$")
RE_PROJECT = re.compile(r"GLOBAL\|\s+Project name\s+(\S+)")
RE_RUNTYPE = re.compile(r"GLOBAL\|\s+Run type\s+(\S+)")
RE_NPROC = re.compile(r"GLOBAL\|\s+Total number of message passing processes\s+(\d+)")
RE_NTHREAD = re.compile(r"GLOBAL\|\s+Number of threads for this process\s+(\d+)")
RE_NATOMS = re.compile(r"^\s*-\s+Atoms:\s+(\d+)\s*$")
RE_SCF_OK = re.compile(r"\*\*\* SCF run converged in\s+(\d+)\s+steps")
RE_SCF_LEAVE = re.compile(r"Leaving inner SCF loop after reaching\s+(\d+)\s+steps")
RE_OUTER_OK = re.compile(r"outer SCF loop converged in\s+(\d+)\s+iterations or\s+(\d+)\s+steps")
RE_OUTER_FAIL = re.compile(r"outer SCF loop FAILED to converge after\s+(\d+)\s+iterations or\s+(\d+)\s+steps")
RE_SCF_NOT = re.compile(r"SCF run NOT converged")
RE_ENERGY = re.compile(r"ENERGY\|\s+Total FORCE_EVAL\s+\(\s*(\S+)\s*\)\s+energy\s+\[([^\]]+)\]:?\s+(" + FLOAT + r")")
RE_FREE = re.compile(r"ENERGY\|\s+Total FORCE_EVAL\s+\(\s*(\S+)\s*\)\s+free energy\s+\[([^\]]+)\]:?\s+(" + FLOAT + r")")
RE_FERMI = re.compile(r"^\s*Fermi energy:\s+(" + FLOAT + r")")
RE_MO_FERMI = re.compile(r"MO\|\s+E\(Fermi\):\s+(" + FLOAT + r")\s+a\.u\.\s+(" + FLOAT + r")\s+eV")
RE_MO_GAP = re.compile(r"MO\|\s+Band gap:\s+(" + FLOAT + r")\s+a\.u\.\s+(" + FLOAT + r")\s+eV")
RE_HL_GAP = re.compile(r"HOMO - LUMO gap \[eV\] :\s*(" + FLOAT + r")")
RE_WARN_START = re.compile(r"^\s*\*\*\* WARNING in (\S+) :: (.*?)\s*\*\*\*\s*$")
RE_STAR_ROW = re.compile(r"^\s*\*\*\* (.*?)\s*\*\*\*\s*$")
RE_NWARN = re.compile(r"number of warnings for this run is\s*:?\s*(\d+)")
RE_MEMORY = re.compile(r"MEMORY\|\s+Estimated peak process memory \[MiB\]\s+(\d+)")
RE_TIMING_CP2K = re.compile(r"^\s*CP2K\s+\d+\s+" + FLOAT + r"\s+(" + FLOAT + r")\s+(" + FLOAT + r")\s+("
                            + FLOAT + r")\s+(" + FLOAT + r")")
RE_STRESS_TRACE = re.compile(r"STRESS\|\s+1/3 Trace\s+(" + FLOAT + r")")
RE_STRESS_UNIT = re.compile(r"STRESS\|\s+(Analytical|Numerical) stress tensor \[([^\]]+)\]")
RE_CELL_VOL = re.compile(r"CELL\|\s+Volume \[([^\]]+)\]:?\s+(" + FLOAT + r")")

# OPT| block (CP2K >= 2023.1), labels from gopt_f_methods.F
OPT_FLOAT_KEYS = [
    (re.compile(r"OPT\|\s+Total energy \[([^\]]+)\]\s+(" + FLOAT + r")"), "energy"),
    (re.compile(r"OPT\|\s+Internal pressure \[([^\]]+)\]\s+(" + FLOAT + r")"), "pressure"),
    (re.compile(r"OPT\|\s+Effective energy change \[([^\]]+)\]\s+(" + FLOAT + r")"), "energy_change"),
    (re.compile(r"OPT\|\s+Predicted energy change \[([^\]]+)\]\s+(" + FLOAT + r")"), "predicted_change"),
    (re.compile(r"OPT\|\s+Pressure deviation(?: with constraint| without constraint)? \[([^\]]+)\]\s+(" + FLOAT + r")"),
     "pressure_deviation"),
    (re.compile(r"OPT\|\s+Pressure tolerance \[([^\]]+)\]\s+(" + FLOAT + r")"), "pressure_tolerance"),
]
OPT_PLAIN_KEYS = [
    (re.compile(r"OPT\|\s+Step size\s+(" + FLOAT + r")"), "step_size"),
    (re.compile(r"OPT\|\s+Trust radius\s+(" + FLOAT + r")"), "trust_radius"),
    (re.compile(r"OPT\|\s+Used time \[s\]\s+(" + FLOAT + r")"), "used_time_s"),
    (re.compile(r"OPT\|\s+Maximum step size\s+(" + FLOAT + r")"), "max_step"),
    (re.compile(r"OPT\|\s+Convergence limit for maximum step size\s+(" + FLOAT + r")"), "max_step_limit"),
    (re.compile(r"OPT\|\s+RMS step size\s+(" + FLOAT + r")"), "rms_step"),
    (re.compile(r"OPT\|\s+Convergence limit for RMS step size\s+(" + FLOAT + r")"), "rms_step_limit"),
    (re.compile(r"OPT\|\s+Maximum gradient\s+(" + FLOAT + r")"), "max_grad"),
    (re.compile(r"OPT\|\s+Convergence limit for maximum gradient\s+(" + FLOAT + r")"), "max_grad_limit"),
    (re.compile(r"OPT\|\s+RMS gradient\s+(" + FLOAT + r")"), "rms_grad"),
    (re.compile(r"OPT\|\s+Convergence limit for RMS gradient\s+(" + FLOAT + r")"), "rms_grad_limit"),
]
OPT_FLAG_KEYS = [
    (re.compile(r"OPT\|\s+Maximum step size is converged\s+(YES|NO)"), "max_step_ok"),
    (re.compile(r"OPT\|\s+RMS step size is converged\s+(YES|NO)"), "rms_step_ok"),
    (re.compile(r"OPT\|\s+Maximum gradient is converged\s+(YES|NO)"), "max_grad_ok"),
    (re.compile(r"OPT\|\s+RMS gradient is converged\s+(YES|NO)"), "rms_grad_ok"),
    (re.compile(r"OPT\|\s+Pressure is converged\s+(YES|NO)"), "pressure_ok"),
    (re.compile(r"OPT\|\s+Decrease in energy\s+(YES|NO)"), "energy_decreased"),
]
RE_OPT_STEP = re.compile(r"OPT\|\s+Step number\s+(\d+)")
RE_OPT_METHOD = re.compile(r"OPT\|\s+Optimization method\s+(\S+)")

# legacy layout (CP2K <= 2022.x), best effort
RE_LEG_STEP = re.compile(r"Informations at step\s*=\s*(\d+)")
LEGACY_KEYS = [
    (re.compile(r"^\s*Total Energy\s*=\s*(" + FLOAT + r")"), "energy"),
    (re.compile(r"^\s*Internal Pressure \[bar\]\s*=\s*(" + FLOAT + r")"), "pressure"),
    (re.compile(r"^\s*Real energy change\s*=\s*(" + FLOAT + r")"), "energy_change"),
    (re.compile(r"^\s*Max\. step size\s*=\s*(" + FLOAT + r")"), "max_step"),
    (re.compile(r"^\s*RMS step size\s*=\s*(" + FLOAT + r")"), "rms_step"),
    (re.compile(r"^\s*Max\. gradient\s*=\s*(" + FLOAT + r")"), "max_grad"),
    (re.compile(r"^\s*RMS gradient\s*=\s*(" + FLOAT + r")"), "rms_grad"),
    (re.compile(r"^\s*Pressure Deviation \[bar\]\s*=\s*(" + FLOAT + r")"), "pressure_deviation"),
]

STEP_COLUMNS = ["step", "method", "energy", "energy_change", "max_step", "rms_step",
                "max_grad", "rms_grad", "pressure", "pressure_deviation",
                "max_step_ok", "rms_step_ok", "max_grad_ok", "rms_grad_ok", "pressure_ok",
                "used_time_s"]


def to_float(text):
    return float(text.replace("D", "E").replace("d", "e"))


def energy_to_ev(value, unit):
    u = unit.strip().lower()
    if u in ("hartree", "a.u.", "au", "ha"):
        return value * HARTREE_EV
    if u == "ev":
        return value
    if u in ("kjmol", "kj/mol"):
        return value / 96.48533212
    if u in ("kcalmol", "kcal/mol"):
        return value / 23.06054783
    if u == "ry":
        return value * HARTREE_EV / 2.0
    return None


def read_lines(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"file not found: {path}")
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        return fh.read().splitlines()


def parse_abort(lines, i):
    """Collect the text of the [ABORT] box that contains line i."""
    texts = []
    j = i
    while j < len(lines):
        raw = lines[j].rstrip()
        stripped = raw.strip()
        if j > i and (stripped.startswith("****") or not stripped.startswith("*")):
            break
        body = stripped.strip("*").strip()
        # drop the little ASCII-art figure printed at the left side of the box
        body = re.sub(r"^(\[ABORT\]|\\___/|O/\||/\| \||/ \\|\||___|/   \\)\s*", "", body)
        if body:
            texts.append(body)
        j += 1
    return " ".join(texts).strip()


def parse_output(path):
    lines = read_lines(path)
    res = OrderedDict()
    res["file"] = path
    res["cp2k_version"] = None
    res["project"] = None
    res["run_type"] = None
    res["mpi_processes"] = None
    res["omp_threads"] = None
    res["natoms"] = None
    res["program_started_count"] = 0
    res["ended_normally"] = False
    res["abort"] = None
    scf = OrderedDict(n_scf_runs=0, converged_runs=0, not_converged_runs=0, steps_per_run=[],
                      outer_loops=0, outer_failed=0, last_converged=None)
    energies = []
    free_energies = []
    fermi_ha = []
    mo_fermi_ev = []
    gaps_ev = []
    hl_gaps_ev = []
    warnings = OrderedDict()
    res["n_warnings_reported"] = None
    res["peak_memory_mib"] = None
    res["walltime_s"] = None
    stress = []
    stress_unit = None
    volume = None
    mulliken = None
    opt_steps = []
    opt_status = None
    cur = None
    last_nc = [-100]  # line index of the last counted "SCF run NOT converged" event

    def count_not_converged(idx):
        # the same failure can appear as a plain line, a WARNING and an ABORT box
        if idx - last_nc[0] > 12:
            scf["n_scf_runs"] += 1
            scf["not_converged_runs"] += 1
        scf["last_converged"] = False
        last_nc[0] = idx

    n = len(lines)
    i = 0
    while i < n:
        line = lines[i]

        if "PROGRAM STARTED AT" in line:
            res["program_started_count"] += 1
        elif "PROGRAM ENDED AT" in line:
            res["ended_normally"] = True

        m = RE_VERSION.search(line)
        if m and res["cp2k_version"] is None:
            res["cp2k_version"] = m.group(1)
        m = RE_PROJECT.search(line)
        if m and res["project"] is None:
            res["project"] = m.group(1)
        m = RE_RUNTYPE.search(line)
        if m and res["run_type"] is None:
            res["run_type"] = m.group(1)
        m = RE_NPROC.search(line)
        if m and res["mpi_processes"] is None:
            res["mpi_processes"] = int(m.group(1))
        m = RE_NTHREAD.search(line)
        if m and res["omp_threads"] is None:
            res["omp_threads"] = int(m.group(1))
        m = RE_NATOMS.match(line)
        if m and res["natoms"] is None:
            res["natoms"] = int(m.group(1))

        # ---------------- SCF ----------------
        m = RE_SCF_OK.search(line)
        if m:
            scf["n_scf_runs"] += 1
            scf["converged_runs"] += 1
            scf["steps_per_run"].append(int(m.group(1)))
            scf["last_converged"] = True
        m = RE_SCF_LEAVE.search(line)
        if m:
            scf["steps_per_run"].append(int(m.group(1)))
        m = RE_OUTER_OK.search(line)
        if m:
            scf["outer_loops"] += 1
            scf["last_converged"] = True
        m = RE_OUTER_FAIL.search(line)
        if m:
            scf["outer_loops"] += 1
            scf["outer_failed"] += 1
            scf["last_converged"] = False
        if RE_SCF_NOT.search(line) and "WARNING" not in line and "ABORT" not in line \
                and "To continue" not in line:
            # plain "*** SCF run NOT converged ***" line (older versions)
            count_not_converged(i)

        # ---------------- energies ----------------
        m = RE_ENERGY.search(line)
        if m:
            energies.append((to_float(m.group(3)), m.group(2).strip()))
        m = RE_FREE.search(line)
        if m:
            free_energies.append((to_float(m.group(3)), m.group(2).strip()))
        m = RE_FERMI.match(line)
        if m:
            fermi_ha.append(to_float(m.group(1)))
        m = RE_MO_FERMI.search(line)
        if m:
            mo_fermi_ev.append(to_float(m.group(2)))
        m = RE_MO_GAP.search(line)
        if m:
            gaps_ev.append(to_float(m.group(2)))
        m = RE_HL_GAP.search(line)
        if m:
            hl_gaps_ev.append(to_float(m.group(1)))

        # ---------------- warnings ----------------
        m = RE_WARN_START.match(line)
        if m:
            loc, msg = m.group(1), m.group(2)
            j = i + 1
            while j < n:
                mm = RE_STAR_ROW.match(lines[j])
                if not mm or RE_WARN_START.match(lines[j]):
                    break
                msg += " " + mm.group(1)
                j += 1
            msg = re.sub(r"\s+", " ", msg).strip()
            key = msg
            if key in warnings:
                warnings[key]["count"] += 1
            else:
                warnings[key] = {"count": 1, "location": loc}
            if "SCF run NOT converged" in msg:
                count_not_converged(i)
            i = j
            continue
        m = RE_NWARN.search(line)
        if m:
            res["n_warnings_reported"] = int(m.group(1))

        if "[ABORT]" in line:
            res["abort"] = parse_abort(lines, i)
            if res["abort"] and "SCF run NOT converged" in res["abort"]:
                count_not_converged(i)

        m = RE_MEMORY.search(line)
        if m:
            res["peak_memory_mib"] = int(m.group(1))
        m = RE_TIMING_CP2K.match(line)
        if m:
            res["walltime_s"] = to_float(m.group(4))  # total time, maximum over ranks

        m = RE_STRESS_UNIT.search(line)
        if m:
            stress_unit = m.group(2)
        m = RE_STRESS_TRACE.search(line)
        if m:
            stress.append(to_float(m.group(1)))
        m = RE_CELL_VOL.search(line)
        if m:
            volume = (to_float(m.group(2)), m.group(1))

        # ---------------- Mulliken (keep the last block) ----------------
        if "Mulliken Population Analysis" in line:
            blk, j = parse_mulliken(lines, i + 1)
            if blk:
                mulliken = blk
            i = j
            continue

        # ---------------- OPT| blocks ----------------
        m = RE_OPT_STEP.search(line)
        if m:
            cur = OrderedDict(step=int(m.group(1)))
            opt_steps.append(cur)
        elif cur is not None and "OPT|" in line:
            m = RE_OPT_METHOD.search(line)
            if m:
                cur["method"] = m.group(1)
            for rx, key in OPT_FLOAT_KEYS:
                mm = rx.search(line)
                if mm:
                    cur[key] = to_float(mm.group(2))
                    cur[key + "_unit"] = mm.group(1)
                    break
            else:
                for rx, key in OPT_PLAIN_KEYS:
                    mm = rx.search(line)
                    if mm:
                        cur[key] = to_float(mm.group(1))
                        break
                else:
                    for rx, key in OPT_FLAG_KEYS:
                        mm = rx.search(line)
                        if mm:
                            cur[key] = (mm.group(1) == "YES")
                            break
        m = RE_LEG_STEP.search(line)
        if m:
            cur = OrderedDict(step=int(m.group(1)), legacy=True)
            opt_steps.append(cur)
        elif cur is not None and cur.get("legacy"):
            for rx, key in LEGACY_KEYS:
                mm = rx.match(line)
                if mm:
                    cur[key] = to_float(mm.group(1))
                    break

        if "GEOMETRY OPTIMIZATION COMPLETED" in line:
            opt_status = "converged"
        elif "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in line:
            opt_status = "max_iter_reached"
        elif "Specific L-BFGS convergence criteria" in line and opt_status is None:
            opt_status = "lbfgs_specific_criteria"
        i += 1

    res["scf"] = scf
    if energies:
        e, u = energies[-1]
        res["final_energy"] = e
        res["energy_unit"] = u
        res["final_energy_eV"] = energy_to_ev(e, u)
        res["n_energy_evaluations"] = len(energies)
    else:
        res["final_energy"] = None
        res["energy_unit"] = None
        res["final_energy_eV"] = None
        res["n_energy_evaluations"] = 0
    if free_energies:
        res["final_free_energy"] = free_energies[-1][0]
    res["fermi_energy_eV"] = (fermi_ha[-1] * HARTREE_EV if fermi_ha
                              else (mo_fermi_ev[-1] if mo_fermi_ev else None))
    if hl_gaps_ev:
        nspin = 2 if mulliken and mulliken.get("spin_polarized") else 1
        res["homo_lumo_gap_eV"] = hl_gaps_ev[-nspin:]
    elif gaps_ev:
        res["homo_lumo_gap_eV"] = gaps_ev[-2:] if len(gaps_ev) >= 2 else gaps_ev
    else:
        res["homo_lumo_gap_eV"] = None
    res["warnings"] = [{"message": k, "count": v["count"], "location": v["location"]}
                       for k, v in warnings.items()]
    res["n_warnings_parsed"] = sum(v["count"] for v in warnings.values())
    res["pressure"] = stress[-1] if stress else None
    res["pressure_unit"] = stress_unit if stress else None
    res["cell_volume"] = volume[0] if volume else None
    res["cell_volume_unit"] = volume[1] if volume else None
    res["mulliken"] = mulliken
    res["opt_steps"] = opt_steps
    res["opt_status"] = opt_status
    return res


def parse_mulliken(lines, start):
    """Parse the atom table following 'Mulliken Population Analysis'."""
    n = len(lines)
    j = start
    header = None
    while j < n and j < start + 8:
        if lines[j].lstrip().startswith("#  Atom") or lines[j].lstrip().startswith("# Orbital"):
            header = lines[j]
            break
        j += 1
    if header is None:
        return None, start
    spin = "alpha,beta" in header
    atoms = []
    total = None
    j += 1
    while j < n:
        s = lines[j].strip()
        if s.startswith("# Total charge"):
            nums = re.findall(FLOAT, s.replace("# Total charge and spin", "").replace("# Total charge", ""))
            total = [to_float(x) for x in nums]
            j += 1
            break
        tok = s.split()
        if len(tok) >= 5 and tok[0].isdigit() and tok[2].isdigit():
            try:
                rec = OrderedDict(index=int(tok[0]), element=tok[1], kind=int(tok[2]))
                vals = [to_float(x) for x in tok[3:]]
            except ValueError:
                j += 1
                continue
            if spin and len(vals) >= 4:
                rec["pop_alpha"], rec["pop_beta"], rec["charge"], rec["spin_moment"] = vals[:4]
            else:
                rec["population"], rec["charge"] = vals[0], vals[-1]
            atoms.append(rec)
        elif s and not s.startswith("#") and atoms and not tok[0].isdigit():
            break
        j += 1
    return OrderedDict(spin_polarized=spin, atoms=atoms, total=total), j


def fmt(value, spec=".6f", none="n/a"):
    if value is None:
        return none
    if isinstance(value, bool):
        return "YES" if value else "NO"
    if isinstance(value, float):
        return format(value, spec)
    return str(value)


def scf_word(scf):
    if scf["not_converged_runs"] or scf["outer_failed"]:
        return "NOT converged ({} of {} SCF runs)".format(
            scf["not_converged_runs"] + scf["outer_failed"], max(scf["n_scf_runs"], 1))
    if scf["converged_runs"] or scf["outer_loops"]:
        return "converged"
    return "unknown"


def print_summary(r, args):
    print("=" * 78)
    print(f"File            : {r['file']}")
    print(f"CP2K version    : {fmt(r['cp2k_version'])}")
    print(f"Project / run   : {fmt(r['project'])} / {fmt(r['run_type'])}")
    print(f"MPI x OpenMP    : {fmt(r['mpi_processes'])} x {fmt(r['omp_threads'])}")
    print(f"Atoms           : {fmt(r['natoms'])}")
    status = "ended normally" if r["ended_normally"] else "NOT finished (no 'PROGRAM ENDED AT')"
    print(f"Status          : {status}")
    if r["program_started_count"] > 1:
        print(f"  !! 'PROGRAM STARTED AT' appears {r['program_started_count']} times: several independent "
              "serial runs wrote into one file (MPI launcher / library mismatch) or the file was appended.")
    if r["abort"]:
        print(f"ABORT           : {r['abort']}")
    scf = r["scf"]
    steps = scf["steps_per_run"]
    print(f"SCF             : {scf_word(scf)}; {len(steps)} SCF cycle(s) with step counts, "
          f"total {sum(steps)} steps, last {steps[-1] if steps else 'n/a'}")
    if scf["outer_loops"]:
        print(f"Outer SCF loops : {scf['outer_loops']} (failed: {scf['outer_failed']})")
    if r["final_energy"] is not None:
        ev = r["final_energy_eV"]
        ev_txt = f" = {ev:.6f} eV" if ev is not None and r["energy_unit"].lower() != "ev" else ""
        print(f"Final energy    : {r['final_energy']:.10f} [{r['energy_unit']}]{ev_txt}"
              f"   ({r['n_energy_evaluations']} ENERGY| line(s))")
    else:
        print("Final energy    : n/a (no 'ENERGY| Total FORCE_EVAL' line)")
    if r.get("final_free_energy") is not None:
        print(f"Free energy     : {r['final_free_energy']:.10f} [{r['energy_unit']}]")
    print(f"Fermi energy    : {fmt(r['fermi_energy_eV'], '.4f')} eV")
    gap = r["homo_lumo_gap_eV"]
    print(f"HOMO-LUMO gap   : {', '.join(f'{g:.4f}' for g in gap) + ' eV' if gap else 'n/a'}")
    if r["pressure"] is not None:
        print(f"Pressure        : {r['pressure']:.4f} [{r['pressure_unit'] or '?'}] (STRESS| 1/3 Trace)")
    if r["cell_volume"] is not None:
        print(f"Cell volume     : {r['cell_volume']:.4f} [{r['cell_volume_unit']}] (last CELL| Volume)")
    print(f"Peak memory     : {fmt(r['peak_memory_mib'])} MiB")
    print(f"Wall time       : {fmt(r['walltime_s'], '.2f')} s (timing table, row CP2K)")
    nrep = r["n_warnings_reported"]
    print(f"Warnings        : {r['n_warnings_parsed']} parsed"
          + (f", CP2K reports {nrep}" if nrep is not None else "")
          + f", {len(r['warnings'])} distinct")
    for w in r["warnings"][:args.warnings]:
        print(f"   [{w['count']:>4}x] {w['message']}  ({w['location']})")
    if len(r["warnings"]) > args.warnings:
        print(f"   ... {len(r['warnings']) - args.warnings} more (use --warnings N)")

    if r["opt_steps"]:
        last = r["opt_steps"][-1]
        status = {"converged": "GEOMETRY OPTIMIZATION COMPLETED",
                  "max_iter_reached": "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED",
                  "lbfgs_specific_criteria": "L-BFGS specific criteria met (check the last OPT| block!)"
                  }.get(r["opt_status"], "running / interrupted")
        print(f"Optimization    : {len(r['opt_steps'])} OPT| block(s), status: {status}")
        print(f"  last step {last['step']}: max_grad {fmt(last.get('max_grad'), '.2e')}"
              f" (limit {fmt(last.get('max_grad_limit'), '.1e')}), rms_grad {fmt(last.get('rms_grad'), '.2e')},"
              f" max_step {fmt(last.get('max_step'), '.2e')}, rms_step {fmt(last.get('rms_step'), '.2e')}")
        if last.get("pressure_deviation") is not None:
            print(f"  pressure deviation {last['pressure_deviation']:.3f} "
                  f"(tolerance {fmt(last.get('pressure_tolerance'), '.1f')})")
    if args.opt and r["opt_steps"]:
        print_opt_table(r["opt_steps"])
    if args.mulliken:
        print_mulliken(r["mulliken"])


def print_opt_table(steps):
    print("-" * 78)
    print(f"{'step':>5} {'energy':>18} {'dE':>11} {'max_step':>10} {'rms_step':>10} "
          f"{'max_grad':>10} {'rms_grad':>10} {'P_dev':>10} conv")
    for s in steps:
        flags = "".join("Y" if s.get(k) else ("n" if k in s else "-")
                        for k in ("max_step_ok", "rms_step_ok", "max_grad_ok", "rms_grad_ok", "pressure_ok"))
        print(f"{s['step']:>5} {fmt(s.get('energy'), '.10f'):>18} {fmt(s.get('energy_change'), '.2e'):>11} "
              f"{fmt(s.get('max_step'), '.2e'):>10} {fmt(s.get('rms_step'), '.2e'):>10} "
              f"{fmt(s.get('max_grad'), '.2e'):>10} {fmt(s.get('rms_grad'), '.2e'):>10} "
              f"{fmt(s.get('pressure_deviation'), '.1f'):>10} {flags}")
    print("conv flags: max_step rms_step max_grad rms_grad pressure (Y = converged, n = not, - = not printed)")


def print_mulliken(mul):
    print("-" * 78)
    if not mul:
        print("Mulliken: no 'Mulliken Population Analysis' block found")
        return
    if mul["spin_polarized"]:
        print(f"{'#':>5} {'El':>3} {'kind':>4} {'pop_alpha':>11} {'pop_beta':>11} {'charge':>10} {'spin':>10}")
        for a in mul["atoms"]:
            print(f"{a['index']:>5} {a['element']:>3} {a['kind']:>4} {a['pop_alpha']:>11.6f} {a['pop_beta']:>11.6f} "
                  f"{a['charge']:>10.6f} {a['spin_moment']:>10.6f}")
    else:
        print(f"{'#':>5} {'El':>3} {'kind':>4} {'population':>12} {'charge':>10}")
        for a in mul["atoms"]:
            print(f"{a['index']:>5} {a['element']:>3} {a['kind']:>4} {a['population']:>12.6f} {a['charge']:>10.6f}")
    if mul["total"]:
        print("total: " + "  ".join(f"{x:.6f}" for x in mul["total"]))
    print("Mulliken charges depend strongly on the basis set: compare trends only.")


def print_brief(results):
    print(f"{'file':<40} {'ended':>5} {'scf':>6} {'energy':>20} {'unit':>8} {'opt':>9} {'warn':>5} {'time_s':>9}")
    for r in results:
        scf = r["scf"]
        scf_txt = "NO" if (scf["not_converged_runs"] or scf["outer_failed"]) else (
            "yes" if (scf["converged_runs"] or scf["outer_loops"]) else "?")
        opt = {"converged": "done", "max_iter_reached": "maxiter",
               "lbfgs_specific_criteria": "lbfgs?"}.get(r["opt_status"], "-" if not r["opt_steps"] else "running")
        name = r["file"] if len(r["file"]) <= 40 else "..." + r["file"][-37:]
        print(f"{name:<40} {'yes' if r['ended_normally'] else 'NO':>5} {scf_txt:>6} "
              f"{fmt(r['final_energy'], '.10f'):>20} {fmt(r['energy_unit']):>8} {opt:>9} "
              f"{r['n_warnings_parsed']:>5} {fmt(r['walltime_s'], '.1f'):>9}")


def write_opt_csv(results, path):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["file"] + STEP_COLUMNS)
        for r in results:
            for s in r["opt_steps"]:
                w.writerow([r["file"]] + ["" if s.get(c) is None else s.get(c) for c in STEP_COLUMNS])
    print(f"written: {path}", file=sys.stderr)


def build_parser():
    p = argparse.ArgumentParser(
        description="Summarise CP2K output files (SCF, energies, Fermi level/gap, Mulliken, "
                    "warnings, timing, GEO_OPT/CELL_OPT steps). Standard library only.",
        epilog="Examples:\n"
               "  python cp2k_outparse.py h2o.out\n"
               "  python cp2k_outparse.py opt.out --opt --opt-csv steps.csv\n"
               "  python cp2k_outparse.py */output.out --brief\n"
               "  python cp2k_outparse.py run.out --json > run.json",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("outputs", nargs="+", help="CP2K output file(s) written with -o (or redirected stdout)")
    p.add_argument("--json", action="store_true", help="print everything as JSON (one object, or a list)")
    p.add_argument("--brief", action="store_true", help="one summary line per file")
    p.add_argument("--opt", action="store_true", help="print the table of all OPT| steps")
    p.add_argument("--opt-csv", metavar="FILE", help="write all OPT| steps of all files to a CSV file")
    p.add_argument("--mulliken", action="store_true", help="print the last Mulliken table")
    p.add_argument("--warnings", type=int, default=10, metavar="N",
                   help="number of distinct warnings to list (default 10)")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    results = []
    status = 0
    for path in args.outputs:
        try:
            r = parse_output(path)
        except (OSError, UnicodeError) as exc:
            print(f"ERROR: {path}: {exc}", file=sys.stderr)
            status = 1
            continue
        if r["cp2k_version"] is None and r["final_energy"] is None and not r["opt_steps"]:
            print(f"WARNING: {path} does not look like a CP2K output file", file=sys.stderr)
        if not r["ended_normally"] or r["abort"]:
            status = 1
        results.append(r)
    if not results:
        return 1
    if args.json:
        json.dump(results[0] if len(results) == 1 else results, sys.stdout, indent=2)
        print()
    elif args.brief:
        print_brief(results)
    else:
        for r in results:
            print_summary(r, args)
    if args.opt_csv:
        write_opt_csv(results, args.opt_csv)
    return status


if __name__ == "__main__":
    sys.exit(main())
