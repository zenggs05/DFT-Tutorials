#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
job_status.py -- scan calculation directories and summarize the state of
                 VASP / Gaussian 16 / CP2K / OpenMX jobs in one table.
                 批量扫描子目录，识别作业状态（正常结束 / 收敛 / 未收敛 / 报错 /
                 运行中 / 中断），汇总能量、离子步（优化步）数和耗时。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

How each code is judged (all files are streamed line by line, so multi-GB
OUTCARs are fine)
--------------------------------------------------------------------------
VASP 6.x (every directory with an OUTCAR)
    finished   : "General timing and accounting" in OUTCAR; time = "Elapsed time (sec):"
                 (unfinished runs: sum of "LOOP+ ... real time")
    energy     : last "energy(sigma->0)" (fallback: "free  energy   TOTEN")
    steps      : number of ionic steps ("free  energy   TOTEN" lines)
    converged  : relaxations (IBRION 1/2/3, NSW > 0) need "reached required accuracy";
                 ionic steps whose electronic loop ran NELM iterations are reported
    errors     : ZBRENT, EDDDAV/EDDRMM (ZHEGV), Sub-Space-Matrix, BRMIX, VERY BAD NEWS,
                 Tetrahedron method fails, ZPOTRF, SIGSEGV, Slurm TIME LIMIT / OOM ...
                 searched in OUTCAR and in the stdout files (vasp.out, *.out, *.log,
                 slurm-*.out, *.o<jobid>, *.err); VASP 6 boxed messages that are
                 split over several lines are joined before matching
Gaussian 16 (*.log, *.out with the Gaussian banner)
    "Normal termination" for every Link1 step / "Error termination" (+ link number),
    last "SCF Done", optimization steps ("Step number N out of a maximum of M"),
    "Optimization completed" / "Number of steps exceeded", imaginary frequencies,
    "Elapsed time:" (summed over Link1 steps)
CP2K (*.out, *.log with "PROGRAM STARTED AT" / "CP2K|")
    "PROGRAM ENDED AT", [ABORT] boxes, "SCF run NOT converged",
    "ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]" (2025.1+) or "[a.u.]" (older),
    "OPT| Step number" (2023.1+) or "Informations at step" (older), "MD| Step number",
    "GEOMETRY OPTIMIZATION COMPLETED" / "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED",
    time from the "CP2K" row of the timing table (or the START/END time stamps)
OpenMX 4.0 / 3.9 (System.Name.out)
    "Elapsed.Time." present = finished, last "Utot." (Hartree), "scf_conv=0/1"
    flags in System.Name.md (SCF convergence of every MD/optimization step),
    number of steps vs. MD.maxIter of the input .dat (optimizations stopped at
    MD.maxIter are flagged)

Status values
    DONE           finished normally and converged (or nothing to converge)
    NOT_CONVERGED  finished, but SCF or geometry optimization did not converge
    ERROR          an error message was found / error termination
    RUNNING        unfinished and the output was modified < --running-min minutes ago
    INCOMPLETE     unfinished and stale (killed, time limit, node failure ...)

Usage examples
--------------
    python job_status.py                       # scan the current directory tree
    python job_status.py relax/ scf/ band/ --depth 2
    python job_status.py ~/projects/TiO2 --only-bad
    python job_status.py . --csv status.csv --json status.json
    python job_status.py . --unit Ha --sort energy
    python job_status.py . --xlsx status.xlsx  # needs pandas + openpyxl

Exit status: 0 = all jobs DONE/RUNNING, 1 = at least one ERROR/NOT_CONVERGED/
INCOMPLETE, 2 = no job found.

Dependencies: Python >= 3.9, standard library; pandas (+ openpyxl) only for --xlsx.
"""

import argparse
import csv
import datetime
import json
import os
import re
import sys
import time

HA2EV = 27.211386245988          # CODATA 2018
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
STATUS_ORDER = {"ERROR": 0, "INCOMPLETE": 1, "NOT_CONVERGED": 2, "RUNNING": 3, "DONE": 4}

# --------------------------------------------------------------------------
# VASP
# --------------------------------------------------------------------------
VASP_ERRORS = [
    # (pattern searched in whitespace-normalised text, short message)
    ("ZBRENT: fatal error", "ZBRENT bracketing failed (try IBRION=1/POTIM smaller, restart from CONTCAR)"),
    ("Error EDDDAV", "EDDDAV/ZHEGV failed (ALGO=Normal/All, check structure)"),
    ("EDDRMM: call to ZHEGV failed", "EDDRMM/ZHEGV failed (ALGO=Normal)"),
    ("Sub-Space-Matrix is not hermitian", "Sub-space matrix not hermitian (ALGO, NCORE, structure)"),
    ("BRMIX: very serious problems", "BRMIX mixing problem (restart without WAVECAR/CHGCAR)"),
    ("VERY BAD NEWS", "VERY BAD NEWS internal error (often symmetry: ISYM/SYMPREC)"),
    ("Tetrahedron method fails", "Tetrahedron method fails (too few k-points for ISMEAR=-5)"),
    ("ZPOTRF", "ZPOTRF failed (overlapping atoms or bad WAVECAR)"),
    ("Error reading item", "INCAR/KPOINTS read error"),
    ("M_divide: can not subdivide", "MPI ranks not divisible (NCORE/KPAR/IMAGES)"),
    ("internal error in subroutine", "internal error (see stdout)"),
]
GENERIC_ERRORS = [
    ("DUE TO TIME LIMIT", "killed by the scheduler: TIME LIMIT"),
    ("oom-kill", "killed: out of memory (oom-kill)"),
    ("Out Of Memory", "killed: out of memory"),
    ("exceeded memory limit", "killed: memory limit"),
    ("SIGSEGV, segmentation fault", "segmentation fault (ulimit -s unlimited, MPI/OMP settings)"),
    ("BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES", "MPI process crashed"),
    ("NODE FAILURE", "node failure"),
]
RE_VASP_E0 = re.compile(r"energy\(sigma->0\)\s*=\s*(" + FLOAT + ")")
RE_VASP_TOTEN = re.compile(r"free  energy   TOTEN\s*=\s*(" + FLOAT + r")\s*eV")
RE_VASP_ITER = re.compile(r"-{3,}\s*Iteration\s+(\d+)\(\s*(\d+)\)")
RE_VASP_PAR = re.compile(r"^\s*(NSW|IBRION|NELM)\s*=\s*(-?\d+)")
RE_VASP_ELAPSED = re.compile(r"Elapsed time \(sec\):\s*(" + FLOAT + ")")
RE_VASP_LOOPP = re.compile(r"LOOP\+:.*?real time\s*(" + FLOAT + ")")
RE_VASP_VERSION = re.compile(r"^\s*vasp\.(\d+\.\d+\.\d+)")


def norm_ws(text):
    """Join VASP 6 boxed messages: drop '|' borders and collapse whitespace."""
    return re.sub(r"\s+", " ", text.replace("|", " "))


def tail_text(path, nbytes=300000):
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            size = fh.tell()
            fh.seek(max(0, size - nbytes))
            return fh.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def stdout_candidates(d, exclude=()):
    out = []
    try:
        names = sorted(os.listdir(d))
    except OSError:
        return out
    for n in names:
        p = os.path.join(d, n)
        if p in exclude or not os.path.isfile(p):
            continue
        if (n in ("vasp.out", "stdout", "out", "log") or n.endswith((".out", ".log", ".err", ".std"))
                or n.startswith("slurm-") or re.search(r"\.[oe]\d+$", n)):
            out.append(p)
    return out


def find_errors(texts, table):
    found = []
    for t in texts:
        nt = norm_ws(t)
        for key, msg in table:
            if key in nt and msg not in found:
                found.append(msg)
    return found


def parse_vasp(d):
    outcar = os.path.join(d, "OUTCAR")
    r = new_row(outcar, "VASP")
    par = {}
    e0 = toten = None
    nion = 0
    finished = reached = False
    loop_real = 0.0
    cur_ion, cur_elec, nelm_hits = 0, 0, 0
    ion_max = {}
    with open(outcar, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "Iteration" in line:
                m = RE_VASP_ITER.search(line)
                if m:
                    cur_ion, cur_elec = int(m.group(1)), int(m.group(2))
                    ion_max[cur_ion] = max(ion_max.get(cur_ion, 0), cur_elec)
                    continue
            if "TOTEN" in line:
                m = RE_VASP_TOTEN.search(line)
                if m:
                    toten = float(m.group(1))
                    nion += 1
                    continue
            if "sigma->0" in line:
                m = RE_VASP_E0.search(line)
                if m:
                    e0 = float(m.group(1))
                continue
            if "LOOP+" in line:
                m = RE_VASP_LOOPP.search(line)
                if m:
                    loop_real += float(m.group(1))
                continue
            if "reached required accuracy" in line:
                reached = True
            elif "General timing and accounting" in line:
                finished = True
            elif "Elapsed time (sec)" in line:
                m = RE_VASP_ELAPSED.search(line)
                if m:
                    r["time_s"] = float(m.group(1))
            elif len(par) < 3 and "=" in line:
                m = RE_VASP_PAR.match(line)
                if m and m.group(1) not in par:
                    par[m.group(1)] = int(m.group(2))
            elif r.get("version") is None and line.lstrip().startswith("vasp."):
                m = RE_VASP_VERSION.match(line)
                if m:
                    r["version"] = m.group(1)
    energy = e0 if e0 is not None else toten
    if energy is not None:
        r["energy_eV"] = energy
    r["steps"] = nion
    if r["time_s"] is None and loop_real > 0:
        r["time_s"] = loop_real
        r["notes"].append("time = sum of LOOP+ (unfinished)")
    nelm = par.get("NELM")
    if nelm:
        nelm_hits = sum(1 for v in ion_max.values() if v >= nelm)
    nsw, ibrion = par.get("NSW", 0), par.get("IBRION", -1)
    relax = nsw > 0 and ibrion in (1, 2, 3)
    texts = [tail_text(outcar)] + [tail_text(p, 100000) for p in stdout_candidates(d, (outcar,))]
    errs = find_errors(texts, VASP_ERRORS + GENERIC_ERRORS)
    if finished:
        conv_ok = True
        if relax and not reached:
            conv_ok = False
            r["notes"].append("NSW=%d reached without EDIFFG convergence" % nsw)
        elif relax:
            r["notes"].append("relaxation converged")
        if nelm_hits:
            conv_ok = False if (nsw == 0 or ion_max.get(max(ion_max), 0) >= nelm) else conv_ok
            r["notes"].append("%d ionic step(s) hit NELM=%d (SCF not converged)" % (nelm_hits, nelm))
        r["status"] = "DONE" if conv_ok else "NOT_CONVERGED"
        r["notes"] += ["warning: " + e for e in errs]
    else:
        if errs:
            r["status"] = "ERROR"
            r["notes"] += errs
        else:
            r["status"] = running_or_incomplete(outcar)
    r["task"] = "relax" if relax else ("md" if nsw > 0 and ibrion == 0 else
                                       ("static" if nsw <= 0 or ibrion == -1 else "IBRION=%d" % ibrion))
    return r


# --------------------------------------------------------------------------
# Gaussian
# --------------------------------------------------------------------------
RE_G_SCF = re.compile(r"SCF Done:\s+E\(([^)]+)\)\s*=\s*(" + FLOAT + ")")
RE_G_STEP = re.compile(r"Step number\s+(\d+)\s+out of a maximum of\s+(\d+)")
RE_G_IMAG = re.compile(r"(\d+)\s+imaginary frequencies \(negative Signs\)")
RE_G_TIME = re.compile(r"Elapsed time:\s+(\d+)\s+days\s+(\d+)\s+hours\s+(\d+)\s+minutes\s+(" + FLOAT + r")\s+seconds")
RE_G_ERRLINK = re.compile(r"Error termination via Lnk1e in \S*?(l\d+)\.exe")
RE_G_ERRREQ = re.compile(r"Error termination request processed by link\s+(\d+)")
GAUSSIAN_HINTS = [
    ("Convergence failure -- run terminated", "SCF convergence failure (l502): SCF=XQC/QC, better guess"),
    ("Number of steps exceeded", "optimization steps exceeded: Opt=(MaxCycles=N), restart from chk"),
    ("could not allocate memory", "memory allocation failed: lower %Mem or request more"),
    ("Erroneous write", "disk full/unwritable scratch (GAUSS_SCRDIR)"),
    ("Bend failed for angle", "linear angle in internal coordinates: Opt=Cartesian or GIC"),
    ("FormBX had a problem", "internal coordinate problem: Opt=Cartesian"),
    ("The combination of multiplicity", "charge/multiplicity inconsistent with the electrons"),
    ("Atoms too close", "atoms too close: check geometry/units"),
    ("End of file in ZSymb", "input ended early: missing blank line after the coordinates?"),
    ("Small interatomic distances encountered", "atoms too close in the geometry"),
]


def parse_gaussian(path):
    r = new_row(path, "Gaussian")
    nlink1 = nnormal = 0
    error = False
    errlink = None
    energy = method = None
    nsteps = maxsteps = 0
    opt_done = opt_fail = False
    imag = None
    tsum = 0.0
    tfound = False
    hints = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "Entering Link 1" in line:
                nlink1 += 1
                continue
            if "SCF Done" in line:
                m = RE_G_SCF.search(line)
                if m:
                    method, energy = m.group(1), float(m.group(2))
                continue
            if "Step number" in line:
                m = RE_G_STEP.search(line)
                if m:
                    nsteps, maxsteps = int(m.group(1)), int(m.group(2))
                continue
            if "Optimization completed" in line or "Stationary point found" in line:
                opt_done = True
            elif "Optimization stopped" in line:
                opt_fail = True
            elif "imaginary frequencies" in line:
                m = RE_G_IMAG.search(line)
                if m:
                    imag = int(m.group(1))
            elif "Normal termination of Gaussian" in line:
                nnormal += 1
            elif "Error termination" in line:
                error = True
                m = RE_G_ERRLINK.search(line) or RE_G_ERRREQ.search(line)
                if m:
                    errlink = m.group(1) if m.group(1).startswith("l") else "l" + m.group(1)
            elif "Elapsed time:" in line:
                m = RE_G_TIME.search(line)
                if m:
                    tfound = True
                    tsum += (int(m.group(1)) * 86400 + int(m.group(2)) * 3600
                             + int(m.group(3)) * 60 + float(m.group(4)))
            for key, msg in GAUSSIAN_HINTS:
                if key in line and msg not in hints:
                    hints.append(msg)
    if energy is not None:
        r["energy_eV"] = energy * HA2EV
        r["task"] = method
    r["steps"] = nsteps if nsteps else None
    if tfound:
        r["time_s"] = tsum
    nlink1 = max(nlink1, 1)
    if error:
        r["status"] = "ERROR"
        r["notes"].append("Error termination" + (" in %s" % errlink if errlink else ""))
    elif nnormal >= nlink1:
        r["status"] = "DONE"
        if opt_fail and not opt_done:
            r["status"] = "NOT_CONVERGED"
    else:
        r["status"] = running_or_incomplete(path)
        if nnormal:
            r["notes"].append("%d of %d Link1 steps finished" % (nnormal, nlink1))
    if nsteps:
        r["notes"].append("opt %s (%d/%d steps)" % ("converged" if opt_done else
                                                     ("NOT converged" if opt_fail else "running"),
                                                     nsteps, maxsteps))
    if imag is not None:
        r["notes"].append("%d imaginary freq" % imag)
    r["notes"] += hints
    return r


# --------------------------------------------------------------------------
# CP2K
# --------------------------------------------------------------------------
RE_C_ENERGY = re.compile(r"ENERGY\|\s+Total FORCE_EVAL\s+\(\s*\S+\s*\)\s+energy\s+\[([^\]]+)\]:?\s+(" + FLOAT + ")")
RE_C_RUNTYPE = re.compile(r"GLOBAL\|\s+Run type\s+(\S+)")
RE_C_VERSION = re.compile(r"CP2K\|\s+version string:\s+(.+?)\s*$")
RE_C_OPTSTEP = re.compile(r"OPT\|\s+Step number\s+(\d+)|Informations at step\s*=\s*(\d+)")
RE_C_MDSTEP = re.compile(r"MD\|\s+Step number\s+(\d+)")
RE_C_TIMING = re.compile(r"^\s*CP2K\s+\d+\s+" + FLOAT + r"\s+" + FLOAT + r"\s+" + FLOAT
                         + r"\s+" + FLOAT + r"\s+(" + FLOAT + r")\s*$")
RE_C_STAMP = re.compile(r"PROGRAM (STARTED|ENDED) AT\s+(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})")


def parse_cp2k(path):
    r = new_row(path, "CP2K")
    ended = abort = False
    abort_msg = []
    n_notconv = 0
    energy = None
    runtype = None
    steps = 0
    md_steps = 0
    opt_done = opt_max = False
    timing = None
    stamps = {}
    grab = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if grab:
                txt = line.strip(" *\n")
                if txt and not set(txt) <= set("_-/\\|*() "):
                    abort_msg.append(txt)
                grab -= 1
            if "ENERGY|" in line:
                m = RE_C_ENERGY.search(line)
                if m:
                    unit = m.group(1).lower()
                    energy = float(m.group(2)) * (HA2EV if unit in ("hartree", "a.u.") else 1.0)
                continue
            if "Step number" in line or "Informations at step" in line:
                m = RE_C_OPTSTEP.search(line)
                if m:
                    steps = int(m.group(1) or m.group(2))
                    continue
                m = RE_C_MDSTEP.search(line)
                if m:
                    md_steps = int(m.group(1))
                    continue
            if "SCF run NOT converged" in line:
                n_notconv += 1
            elif "PROGRAM ENDED AT" in line or "PROGRAM STARTED AT" in line:
                if "ENDED" in line:
                    ended = True
                m = RE_C_STAMP.search(line)
                if m:
                    stamps[m.group(1)] = m.group(2)
            elif "ABORT" in line:
                abort = True
                grab = 6
            elif "GEOMETRY OPTIMIZATION COMPLETED" in line:
                opt_done = True
            elif "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in line:
                opt_max = True
            elif runtype is None and "Run type" in line:
                m = RE_C_RUNTYPE.search(line)
                if m:
                    runtype = m.group(1)
            elif r.get("version") is None and "version string" in line:
                m = RE_C_VERSION.search(line)
                if m:
                    r["version"] = m.group(1)
            elif line.lstrip().startswith("CP2K "):
                m = RE_C_TIMING.match(line)
                if m:
                    timing = float(m.group(1))
    r["task"] = runtype
    if energy is not None:
        r["energy_eV"] = energy
    r["steps"] = steps or md_steps or None
    if timing is not None:
        r["time_s"] = timing
    elif "STARTED" in stamps and "ENDED" in stamps:
        fmt = "%Y-%m-%d %H:%M:%S"
        try:
            dt = (datetime.datetime.strptime(stamps["ENDED"], fmt)
                  - datetime.datetime.strptime(stamps["STARTED"], fmt))
            r["time_s"] = dt.total_seconds()
        except ValueError:
            pass
    if abort:
        r["status"] = "ERROR"
        r["notes"].append("ABORT: " + " ".join(abort_msg)[:160] if abort_msg else "ABORT")
    elif ended:
        bad = n_notconv > 0 or opt_max
        r["status"] = "NOT_CONVERGED" if bad else "DONE"
    else:
        r["status"] = running_or_incomplete(path)
    if n_notconv:
        r["notes"].append("SCF run NOT converged x%d" % n_notconv)
    if opt_done:
        r["notes"].append("optimization completed")
    if opt_max:
        r["notes"].append("MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED")
    if not ended and not abort:
        errs = find_errors([tail_text(p, 50000) for p in stdout_candidates(os.path.dirname(path) or ".",
                                                                              (path,))], GENERIC_ERRORS)
        if errs:
            r["status"] = "ERROR"
            r["notes"] += errs
    return r


# --------------------------------------------------------------------------
# OpenMX
# --------------------------------------------------------------------------
RE_O_UTOT = re.compile(r"^\s*Utot\.\s+(" + FLOAT + ")")
RE_O_TIME = re.compile(r"Elapsed\.Time\.\s+(" + FLOAT + ")")
OPT_TYPES = {"opt", "diis", "bfgs", "rf", "ef", "rfc5", "optc1", "optc2", "optc3", "optc4",
             "optc5", "optc6", "optc7", "rfc6", "rfc7"}


def omx_input_for(outpath):
    """Find the .dat whose System.Name equals the stem of System.Name.out."""
    d = os.path.dirname(outpath) or "."
    stem = os.path.basename(outpath)[:-4]
    for n in sorted(os.listdir(d)):
        if not n.endswith(".dat"):
            continue
        p = os.path.join(d, n)
        kw = {}
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                w = line.split("#", 1)[0].split()
                if len(w) >= 2:
                    kw[w[0].lower()] = w[1]
        if kw.get("system.name") == stem:
            return p, kw
    return None, {}


def parse_openmx(path):
    r = new_row(path, "OpenMX")
    utot = None
    finished = False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "Utot." in line:
                m = RE_O_UTOT.match(line)
                if m:
                    utot = float(m.group(1))
            elif "Elapsed.Time." in line:
                m = RE_O_TIME.search(line)
                finished = True
                if m:
                    r["time_s"] = float(m.group(1))
    if utot is not None:
        r["energy_eV"] = utot * HA2EV
    dat, kw = omx_input_for(path)
    mdtype = kw.get("md.type", "nomd").lower()
    r["task"] = kw.get("md.type", "NOMD")
    mdfile = path[:-4] + ".md"
    nconv0 = nblocks = 0
    if os.path.isfile(mdfile):
        with open(mdfile, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if "scf_conv=" in line:
                    nblocks += 1
                    if "scf_conv=0" in line:
                        nconv0 += 1
        r["steps"] = nblocks or None
    if finished:
        bad = nconv0 > 0
        if nconv0:
            r["notes"].append("%d step(s) with scf_conv=0" % nconv0)
        if mdtype in OPT_TYPES and nblocks:
            try:
                maxiter = int(kw.get("md.maxiter", "1"))
            except ValueError:
                maxiter = 0
            if maxiter and nblocks >= maxiter:
                bad = True
                r["notes"].append("optimization stopped at MD.maxIter=%d (check History of ... table)"
                                  % maxiter)
        r["status"] = "NOT_CONVERGED" if bad else "DONE"
        if not os.path.isfile(mdfile):
            r["notes"].append("no %s: SCF convergence not checked" % os.path.basename(mdfile))
    else:
        r["status"] = running_or_incomplete(path)
        d = os.path.dirname(path) or "."
        errs = find_errors([tail_text(p, 50000) for p in stdout_candidates(d, (path,))], GENERIC_ERRORS)
        if errs:
            r["status"] = "ERROR"
            r["notes"] += errs
    if dat is None:
        r["notes"].append("input .dat not found (MD.Type/MD.maxIter unknown)")
    return r


# --------------------------------------------------------------------------
# common
# --------------------------------------------------------------------------
RUNNING_MIN = 10.0


def new_row(path, code):
    return {"path": path, "code": code, "status": "?", "task": None, "energy_eV": None,
            "steps": None, "time_s": None, "version": None, "notes": []}


def running_or_incomplete(path):
    try:
        age = (time.time() - os.path.getmtime(path)) / 60.0
    except OSError:
        return "INCOMPLETE"
    return "RUNNING" if age < RUNNING_MIN else "INCOMPLETE"


def head_text(path, nbytes=40000):
    try:
        with open(path, "rb") as fh:
            return fh.read(nbytes).decode("utf-8", errors="replace")
    except OSError:
        return ""


def classify(path):
    """Return 'gaussian', 'cp2k', 'openmx' or None for a .log/.out file."""
    h = head_text(path)
    if "Entering Gaussian System" in h or "Gaussian, Inc." in h:
        return "gaussian"
    if "PROGRAM STARTED AT" in h or "CP2K|" in h:
        return "cp2k"
    if path.endswith(".out"):
        if "OpenMX" in h:
            return "openmx"
        try:
            if omx_input_for(path)[0]:      # a .dat with System.Name = stem of this file
                return "openmx"
        except OSError:
            pass
    return None


def scan(roots, depth, codes):
    rows = []
    seen = set()
    for root in roots:
        if os.path.isfile(root):
            entries = [(os.path.dirname(root) or ".", [os.path.basename(root)])]
        elif os.path.isdir(root):
            entries = []
            for d, dirs, files in os.walk(root):
                dirs[:] = sorted(x for x in dirs if not x.startswith(".") and not x.endswith("_rst"))
                rel_d = os.path.relpath(d, root)
                level = 0 if rel_d == "." else len(rel_d.split(os.sep))
                if depth is not None and level >= depth:
                    dirs[:] = []
                entries.append((d, sorted(files)))
        else:
            print("warning: %s does not exist" % root, file=sys.stderr)
            continue
        for d, files in entries:
            fs = set(files)
            if "OUTCAR" in fs and "vasp" in codes:
                p = os.path.join(d, "OUTCAR")
                if p not in seen:
                    seen.add(p)
                    rows.append(safe(parse_vasp, d, p))
            for f in files:
                if not f.endswith((".log", ".out")) or f.startswith("slurm-") or f == "vasp.out":
                    continue
                p = os.path.join(d, f)
                if p in seen:
                    continue
                kind = classify(p)
                if kind is None or kind not in codes:
                    continue
                seen.add(p)
                func = {"gaussian": parse_gaussian, "cp2k": parse_cp2k, "openmx": parse_openmx}[kind]
                rows.append(safe(func, p, p))
    return rows


def safe(func, arg, path):
    try:
        return func(arg)
    except (OSError, UnicodeError, ValueError) as exc:
        r = new_row(path, func.__name__.replace("parse_", ""))
        r["status"] = "ERROR"
        r["notes"].append("could not parse: %s" % exc)
        return r


def fmt_time(sec):
    if sec is None:
        return ""
    sec = int(round(sec))
    d, rem = divmod(sec, 86400)
    h, rem = divmod(rem, 3600)
    m, s = divmod(rem, 60)
    return ("%dd%02d:%02d:%02d" % (d, h, m, s)) if d else ("%02d:%02d:%02d" % (h, m, s))


def build_parser():
    p = argparse.ArgumentParser(
        prog="job_status.py",
        description="Summarize VASP / Gaussian / CP2K / OpenMX job states below one or more directories.",
        epilog=("Examples:\n"
                "  python job_status.py\n"
                "  python job_status.py relax/ scf/ band/ --depth 2\n"
                "  python job_status.py ~/projects/TiO2 --only-bad\n"
                "  python job_status.py . --csv status.csv --json status.json --unit Ha\n"
                "Status: DONE, NOT_CONVERGED, ERROR, RUNNING, INCOMPLETE (see the script header)."),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("roots", nargs="*", default=["."], help="directories (or single output files) to scan")
    p.add_argument("--depth", type=int, help="maximum directory depth below each root (default: unlimited)")
    p.add_argument("--codes", default="vasp,gaussian,cp2k,openmx",
                   help="comma separated subset of vasp,gaussian,cp2k,openmx")
    p.add_argument("--unit", choices=("eV", "Ha"), default="eV", help="energy unit of the table (default eV)")
    p.add_argument("--sort", choices=("path", "status", "energy", "time"), default="path",
                   help="row order (default path)")
    p.add_argument("--only-bad", action="store_true", help="show only ERROR / NOT_CONVERGED / INCOMPLETE rows")
    p.add_argument("--running-min", type=float, default=10.0,
                   help="unfinished output modified within this many minutes counts as RUNNING (default 10)")
    p.add_argument("--csv", metavar="FILE", help="write the table as CSV")
    p.add_argument("--json", metavar="FILE", help="write all fields as JSON")
    p.add_argument("--xlsx", metavar="FILE", help="write an Excel file (needs pandas + openpyxl)")
    p.add_argument("--no-notes", action="store_true", help="hide the notes column in the terminal")
    return p


def main(argv=None):
    global RUNNING_MIN
    a = build_parser().parse_args(argv)
    RUNNING_MIN = a.running_min
    codes = {c.strip().lower() for c in a.codes.split(",") if c.strip()}
    unknown = codes - {"vasp", "gaussian", "cp2k", "openmx"}
    if unknown:
        print("job_status.py: error: unknown code(s) in --codes: %s" % ", ".join(sorted(unknown)),
              file=sys.stderr)
        return 2
    rows = scan(a.roots, a.depth, codes)
    if not rows:
        print("no VASP/Gaussian/CP2K/OpenMX output found below: %s" % " ".join(a.roots))
        return 2
    fac = 1.0 / HA2EV if a.unit == "Ha" else 1.0
    for r in rows:
        r["energy"] = None if r["energy_eV"] is None else r["energy_eV"] * fac
    key = {"path": lambda r: r["path"],
           "status": lambda r: (STATUS_ORDER.get(r["status"], 9), r["path"]),
           "energy": lambda r: (r["energy"] is None, r["energy"] or 0.0),
           "time": lambda r: (r["time_s"] is None, -(r["time_s"] or 0.0))}[a.sort]
    rows.sort(key=key)
    shown = [r for r in rows if not a.only_bad or r["status"] in ("ERROR", "NOT_CONVERGED", "INCOMPLETE")]

    cwd = os.getcwd()

    def rel(pth):
        try:
            return os.path.relpath(pth, cwd)
        except ValueError:      # different drive on Windows
            return pth

    w = max([len(rel(r["path"])) for r in shown] + [4])
    ecol = "E(%s)" % a.unit
    print("%-*s  %-8s %-13s %-10s %18s %6s %11s  %s"
          % (w, "path", "code", "status", "task", ecol, "steps", "time", "" if a.no_notes else "notes"))
    for r in shown:
        e = "" if r["energy"] is None else ("%18.6f" % r["energy"])
        print("%-*s  %-8s %-13s %-10s %18s %6s %11s  %s"
              % (w, rel(r["path"]), r["code"], r["status"], (r["task"] or "")[:10], e,
                 "" if r["steps"] is None else r["steps"], fmt_time(r["time_s"]),
                 "" if a.no_notes else "; ".join(r["notes"])))
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    print("\n%d job(s): %s" % (len(rows), ", ".join("%s=%d" % (k, counts[k])
                                                   for k in sorted(counts, key=lambda s: STATUS_ORDER.get(s, 9)))))

    fields = ["path", "code", "status", "task", "energy_eV", "steps", "time_s", "version", "notes"]
    if a.csv:
        with open(a.csv, "w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh)
            wr.writerow(fields)
            for r in rows:
                wr.writerow([("; ".join(r[f]) if f == "notes" else ("" if r[f] is None else r[f]))
                             for f in fields])
        print("written %s" % a.csv)
    if a.json:
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump([{f: r[f] for f in fields} for r in rows], fh, indent=1)
        print("written %s" % a.json)
    if a.xlsx:
        try:
            import pandas as pd
        except ImportError:
            print("--xlsx needs pandas (and openpyxl): pip install pandas openpyxl", file=sys.stderr)
            return 2
        df = pd.DataFrame([{f: ("; ".join(r[f]) if f == "notes" else r[f]) for f in fields} for r in rows])
        try:
            df.to_excel(a.xlsx, index=False)
        except ImportError:
            print("--xlsx needs openpyxl: pip install openpyxl", file=sys.stderr)
            return 2
        print("written %s" % a.xlsx)
    bad = any(r["status"] in ("ERROR", "NOT_CONVERGED", "INCOMPLETE") for r in rows)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
