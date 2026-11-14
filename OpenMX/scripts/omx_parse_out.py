#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_parse_out.py -- parse OpenMX 4.0 / 3.9 output (System.Name.out, .std, .md).
                    解析 OpenMX 主输出：总能 Utot 及能量分解、化学势、SCF 历史（NormRD/Uele、
                    dUele）、是否收敛、Mulliken 布居/电荷/自旋磁矩、原子力、结构/晶胞优化历史、
                    网格与耗时，输出文字摘要、JSON 或 CSV，供其它脚本复用。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Sections recognised in System.Name.out (formats of OpenMX 4.0, also valid for 3.9)
------------------------------------------------------------------------------------
  echoed input            keywords before "Required cutoff energy"
  grid                    "Required cutoff energy (Ryd) ...", "Num. of grids of a-, b-, and c-axes = ..."
  SCF history at MD= n    "   SCF=   1  NormRD=  1.000000000000  Uele= -3.523169099737"
  Total energy (Hartree) at MD = n   lines "  Utot.   -8.033540953061" (Uele, Ukin, ..., Enpy)
  Chemical potential (Hartree)       value
  Eigenvalues (Hartree) for SCF KS-eq.   HOMO and eigenvalues (Cluster calculations only)
  Mulliken populations    per atom Up/Down/Sum/Diff (+ theta/phi for nc), total spin moment
  Dipole moment (Debye)   "Absolute D"
  <coordinates.forces ... coordinates.forces>   x y z (Ang), Fx Fy Fz (Hartree/Bohr)
  History of geometry optimization / History of cell optimization (table)
  Computational Time (second)   "Elapsed.Time." and "Total Computational Time"

System.Name.std is used for the termination message ("The calculation was
normally finished.") and for runs that crashed before the .out was complete
("<DFT>  Uele   = ... dUele = ..." / "<DFT>  NormRD = ... Criterion = ...").
System.Name.md (manual Sec. 18.6) supplies scf_conv=0/1 flags per MD step.

SCF convergence
---------------
OpenMX leaves the SCF loop either when the criterion is met or after
scf.maxIter steps.  A block is flagged "not converged" when its number of
steps reaches scf.maxIter (read from the echoed input, default 40) or when the
.md file reports scf_conv=0 for that step.

Charges and magnetic moments
----------------------------
Mulliken charge = valence electrons of the VPS - Mulliken population (Sum).
The valence charge is read from DATA.PATH/VPS/<vps>.vps when available
(--data-path, or DATA.PATH from the echoed input), otherwise from the
DFT_DATA19 table in omx_make_input.py (same folder).  Spin moments are the
"Diff" column in the OpenMX convention (opposite in sign to the physical
magnetic moment of electrons, see chapter 11).

Usage examples
--------------
  python omx_parse_out.py met.out
  python omx_parse_out.py si_ecut_*.out --summary-csv ecut_summary.csv
  python omx_parse_out.py fe.out --json fe.json --csv-dir fe_tables
  python omx_parse_out.py crashed.std                  # .std only (crashed runs)
  python omx_parse_out.py nio.out --data-path /opt/openmx4.0/DFT_DATA19 --strict

  # as a module (other scripts of this folder do this)
  from omx_parse_out import parse_out
  res = parse_out("met.out"); print(res["utot_Ha"], res["scf_converged"])

Dependencies: python >= 3.9, standard library only.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys

HA2EV = 27.211386245988          # CODATA 2018
BOHR2ANG = 0.529177210903
HA_BOHR2EV_ANG = HA2EV / BOHR2ANG

_NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?"
FINISHED_MSG = "The calculation was normally finished."
ERROR_PATTERNS = re.compile(
    r"Could not find|can not open|Found an undefined species|ERROR key=|errors in the inputfile|"
    r"Format error|Segmentation fault|signal 11|BAD TERMINATION|Killed|out of memory|oom-kill|"
    r"error while loading shared libraries|Invalid communicator|was incorrect on entry")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _read(path: str) -> str:
    if not os.path.isfile(path):
        raise FileNotFoundError("file not found: %s" % path)
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read()


def _f(x: str) -> float:
    return float(x)


def parse_keywords(text: str) -> tuple[dict, dict]:
    """Keywords and blocks of an OpenMX input (or of the input echoed in a .out).

    Returns (keywords, blocks); keys are lower case, values are the raw strings
    after the keyword (comments removed); blocks map to lists of token lists.
    """
    kw, blocks = {}, {}
    current = None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if current is not None:
            if line.lower() == current + ">" or line.lower().endswith(current + ">"):
                current = None
                continue
            blocks[current].append(line.split())
            continue
        if line.startswith("<"):
            current = line[1:].split()[0].lower()
            blocks[current] = []
            continue
        tok = line.split(None, 1)
        kw[tok[0].lower()] = tok[1].strip() if len(tok) > 1 else ""
    return kw, blocks


def _valence_lookup():
    """Import valence_of() from omx_make_input.py of the same folder (optional)."""
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    try:
        from omx_make_input import valence_of  # type: ignore
        return valence_of
    except Exception:
        return None


# ---------------------------------------------------------------------------
# .md file
# ---------------------------------------------------------------------------
def parse_md(path: str) -> list[dict]:
    """Blocks of System.Name.md: natoms / comment / atom lines (manual Sec. 18.6)."""
    lines = _read(path).splitlines()
    steps, i = [], 0
    while i < len(lines):
        tok = lines[i].split()
        if len(tok) == 1 and tok[0].isdigit():
            nat = int(tok[0])
            if i + 2 + nat > len(lines):
                break                                   # incomplete last block
            comment = lines[i + 1]
            step = {"natoms": nat}
            for key in ("time", "Energy", "Temperature"):
                m = re.search(r"\b%s=\s*(%s)" % (key, _NUM), comment)
                if m:
                    step[key.lower()] = float(m.group(1))
            m = re.search(r"scf_conv=\s*(\d)", comment)
            if m:
                step["scf_conv"] = int(m.group(1))
            m = re.search(r"Cell_Vectors=\s*((?:%s\s+){8}%s)" % (_NUM, _NUM), comment)
            if m:
                v = [float(x) for x in m.group(1).split()]
                step["cell"] = [v[0:3], v[3:6], v[6:9]]
            fmax = 0.0
            ok = True
            for j in range(nat):
                if i + 2 + j >= len(lines):
                    ok = False
                    break
                t = lines[i + 2 + j].split()
                if len(t) >= 7:
                    try:
                        fx, fy, fz = float(t[4]), float(t[5]), float(t[6])
                        fmax = max(fmax, (fx * fx + fy * fy + fz * fz) ** 0.5)
                    except ValueError:
                        pass
            if not ok:
                break
            step["max_force_Ha_Bohr"] = fmax
            steps.append(step)
            i += 2 + nat
        else:
            i += 1
    return steps


# ---------------------------------------------------------------------------
# .std file
# ---------------------------------------------------------------------------
def parse_std(path: str) -> dict:
    text = _read(path)
    res = {"file": path, "finished": FINISHED_MSG in text, "scf_std": [], "errors": []}
    md, scf = None, None
    cur = None
    for line in text.splitlines():
        m = re.search(r"\*+\s*MD=\s*(\d+)\s+SCF=\s*(\d+)\s*\*+", line)
        if m:
            md, scf = int(m.group(1)), int(m.group(2))
            cur = {"md": md, "scf": scf}
            res["scf_std"].append(cur)
            continue
        m = re.match(r"\s*<DFT>\s+Uele\s*=\s*(%s)\s+dUele\s*=\s*(%s)" % (_NUM, _NUM), line)
        if m and cur is not None:
            cur["Uele"], cur["dUele"] = _f(m.group(1)), _f(m.group(2))
            continue
        m = re.match(r"\s*<DFT>\s+NormRD\s*=\s*(%s)\s+Criterion\s*=\s*(%s)" % (_NUM, _NUM), line)
        if m and cur is not None:
            cur["NormRD"], cur["criterion"] = _f(m.group(1)), _f(m.group(2))
            continue
        m = re.match(r"\s*<DFT>\s+MD=\s*(\d+)\s+SCF=\s*(\d+)\s+NormRD\s*=\s*(%s)\s+Criterion\s*=\s*(%s)"
                     % (_NUM, _NUM), line)                       # level.of.stdout 0
        if m:
            res["scf_std"].append({"md": int(m.group(1)), "scf": int(m.group(2)),
                                   "NormRD": _f(m.group(3)), "criterion": _f(m.group(4))})
            continue
        if ERROR_PATTERNS.search(line) and len(res["errors"]) < 20:
            res["errors"].append(line.strip())
    return res


# ---------------------------------------------------------------------------
# .out file
# ---------------------------------------------------------------------------
def parse_out(path: str, data_path: str | None = None, std: str | None = None,
              md: str | None = None) -> dict:
    """Parse System.Name.out; .std/.md next to it are used automatically if present."""
    text = _read(path)
    res: dict = {"file": path}
    stem = path[:-4] if path.endswith(".out") else path
    m = re.search(r"performed by OpenMX Ver\.\s*(\S+)\s+using\s+(\d+)\s+MPI processes and\s+(\d+)\s+OpenMP",
                  text[:3000])
    if m:
        res["openmx_version"] = m.group(1)
        res["mpi_processes"], res["omp_threads"] = int(m.group(2)), int(m.group(3))

    # ---- echoed input
    cut = len(text)
    for marker in ("Required cutoff energy", "SCF history at MD"):
        k = text.find(marker)
        if 0 <= k < cut:
            cut = k
    kw, blocks = parse_keywords(text[:cut])
    res["system_name"] = kw.get("system.name", os.path.basename(stem))
    try:
        maxiter = int(float(kw.get("scf.maxiter", "40").split()[0]))
    except ValueError:
        maxiter = 40
    res["scf_maxIter"] = maxiter
    for key in ("scf.criterion", "scf.energycutoff"):
        if key in kw:
            try:
                res[key.replace("scf.", "input_")] = float(kw[key].split()[0])
            except ValueError:
                pass
    for key in ("scf.spinpolarization", "scf.eigenvaluesolver", "scf.xctype", "md.type",
                "scf.kgrid", "scf.ngrid", "data.path"):
        if key in kw:
            res["input_" + key.split(".", 1)[1].lower()] = kw[key]
    species_vps = {}
    for row in blocks.get("definition.of.atomic.species", []):
        if len(row) >= 3:
            species_vps[row[0]] = row[2]
    res["species"] = {k: {"vps": v} for k, v in species_vps.items()}
    try:
        res["natoms"] = int(kw.get("atoms.number", "0").split()[0])
    except ValueError:
        res["natoms"] = 0

    # ---- grid
    m = re.search(r"Required cutoff energy \(Ryd\) for 3D-grids\s*=\s*(%s)" % _NUM, text)
    if m:
        res["required_cutoff_Ry"] = float(m.group(1))
    m = re.search(r"Used cutoff energy \(Ryd\) for 3D-grids\s*=\s*([^\n]+)", text)
    if m:
        res["used_cutoff_Ry"] = [float(x) for x in re.findall(_NUM, m.group(1))]
    m = re.search(r"Num\. of grids of a-, b-, and c-axes\s*=\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)", text)
    if m:
        res["ngrid"] = [int(m.group(i)) for i in (1, 2, 3)]

    # ---- SCF histories
    hist = []
    parts = re.split(r"SCF history at MD=\s*(\d+)", text)
    for j in range(1, len(parts), 2):
        mdstep = int(parts[j])
        rows = re.findall(r"^\s*SCF=\s*(\d+)\s+NormRD=\s*(%s)\s+Uele=\s*(%s)" % (_NUM, _NUM),
                          parts[j + 1], re.M)
        steps = [{"scf": int(a), "NormRD": float(b), "Uele": float(c)} for a, b, c in rows]
        for k in range(1, len(steps)):
            steps[k]["dUele"] = abs(steps[k]["Uele"] - steps[k - 1]["Uele"])
        hist.append({"md": mdstep, "steps": steps})
    res["scf_history"] = hist

    # ---- energies
    eblocks = []
    parts = re.split(r"Total energy \(Hartree\) at MD\s*=\s*(\d+)", text)
    for j in range(1, len(parts), 2):
        body = parts[j + 1].split("Note:", 1)[0]
        terms = {name: float(val) for name, val in
                 re.findall(r"^\s*(U[A-Za-z0-9]+|Enpy)\.\s+(%s)\s*$" % _NUM, body, re.M)}
        eblocks.append({"md": int(parts[j]), "terms": terms})
    res["energy_blocks"] = eblocks
    if eblocks:
        last = eblocks[-1]["terms"]
        res["energies_Ha"] = last
        if "Utot" in last:
            res["utot_Ha"] = last["Utot"]
            res["utot_eV"] = last["Utot"] * HA2EV
        res["utot_history_Ha"] = [b["terms"].get("Utot") for b in eblocks]
    mu = re.findall(r"Chemical potential \(Hartree\)\s+(%s)" % _NUM, text)
    if mu:
        res["chemical_potential_Ha"] = float(mu[-1])
        res["chemical_potential_eV"] = float(mu[-1]) * HA2EV

    # ---- eigenvalues
    #   Cluster:      "Eigenvalues (Hartree) for SCF KS-eq."          (+ "HOMO = n")
    #   Cluster, nc:  "Eigenvalues (Hartree) for non-collinear KS-eq." (one column)
    #   Band:         "Eigenvalues (Hartree) of SCF KS-eq."           (kloop=... blocks)
    em = list(re.finditer(r"Eigenvalues \(Hartree\) (?:for|of) [^\n]*KS-eq\.", text))
    if em:
        block = text[em[-1].end():]
        end = re.search(r"^\s*Mulliken populations\s*$", block, re.M)
        block = block[:end.start()] if end else block[:2000000]
        rows = []
        for line in block.splitlines():
            t = line.split()
            if len(t) in (2, 3) and t[0].isdigit():
                try:
                    rows.append((int(t[0]), [float(x) for x in t[1:]]))
                except ValueError:
                    continue
        mh = re.search(r"HOMO\s*=\s*(\d+)", block)
        if mh and "kloop" not in block:
            eig = dict(rows)
            homo = int(mh.group(1))
            res["homo_index"] = homo
            if homo in eig:
                res["homo_Ha"] = eig[homo]
                if homo + 1 in eig:
                    res["lumo_Ha"] = eig[homo + 1]
                    res["gap_eV"] = (min(eig[homo + 1]) - max(eig[homo])) * HA2EV
        elif "kloop" in block and mu:
            # band edges on the SCF k-grid relative to the chemical potential (estimate only)
            ef = float(mu[-1])
            vals = [e for _, es in rows for e in es]
            occ = [e for e in vals if e <= ef]
            emp = [e for e in vals if e > ef]
            if occ and emp:
                res["vbm_scf_kgrid_eV"] = max(occ) * HA2EV
                res["cbm_scf_kgrid_eV"] = min(emp) * HA2EV
                res["gap_scf_kgrid_eV"] = (min(emp) - max(occ)) * HA2EV

    # ---- Mulliken
    mm = list(re.finditer(r"^\s*Mulliken populations\s*$", text, re.M))
    if mm:
        block = text[mm[-1].end():]
        end = re.search(r"Decomposed Mulliken populations|Dipole moment \(Debye\)|<coordinates\.forces",
                        block)
        block = block[:end.start()] if end else block[:2000000]
        atoms = []
        for line in block.splitlines():
            t = line.split()
            if len(t) in (6, 8) and t[0].isdigit():
                try:
                    nums = [float(x) for x in t[2:]]
                except ValueError:
                    continue
                a = {"index": int(t[0]), "species": t[1], "up": nums[0], "down": nums[1],
                     "sum": nums[2], "moment": nums[3]}
                if len(nums) == 6:
                    a["theta_deg"], a["phi_deg"] = nums[4], nums[5]
                atoms.append(a)
        m = re.search(r"Total spin moment \(muB\)\s+(%s)" % _NUM, block)
        if m:
            res["total_spin_moment_muB"] = float(m.group(1))
        else:
            m = re.search(r"Total spin S\s*=\s*(%s)" % _NUM, block)
            if m:
                res["total_spin_moment_muB"] = 2.0 * float(m.group(1))
        m = re.search(r"Angles \(Deg\)\s+(%s)\s+(%s)" % (_NUM, _NUM), block)
        if m:
            res["total_spin_angles_deg"] = [float(m.group(1)), float(m.group(2))]
        m = re.search(r"total=\s*(%s)\s+ideal\(neutral\)=\s*(%s)" % (_NUM, _NUM), block)
        if m:
            res["mulliken_total"] = float(m.group(1))
            res["electrons_neutral"] = float(m.group(2))
        # charges
        dp = data_path or kw.get("data.path")
        if dp and not os.path.isabs(dp):
            dp = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)), dp))
        valence_of = _valence_lookup()
        for a in atoms:
            vps = species_vps.get(a["species"])
            val = valence_of(vps, dp) if (valence_of and vps) else None
            if val is not None:
                a["valence"] = val
                a["charge"] = val - a["sum"]
        res["mulliken"] = atoms

    # ---- dipole
    m = list(re.finditer(r"Absolute D\s+(%s)" % _NUM, text))
    if m:
        res["dipole_abs_Debye"] = float(m[-1].group(1))
        mt = re.findall(r"^\s*Total\s+(%s)\s+(%s)\s+(%s)\s*$" % (_NUM, _NUM, _NUM), text[m[-1].start():], re.M)
        if mt:
            res["dipole_Debye"] = [float(x) for x in mt[0]]

    # ---- forces (last block)
    blocks_f = re.findall(r"<coordinates\.forces(.*?)coordinates\.forces>", text, re.S)
    if blocks_f:
        forces = []
        for line in blocks_f[-1].splitlines():
            t = line.split()
            if len(t) >= 8 and t[0].isdigit():
                try:
                    v = [float(x) for x in t[2:8]]
                except ValueError:
                    continue
                f = (v[3] ** 2 + v[4] ** 2 + v[5] ** 2) ** 0.5
                forces.append({"index": int(t[0]), "species": t[1], "xyz_Ang": v[:3],
                               "force_Ha_Bohr": v[3:], "abs_force_Ha_Bohr": f})
        res["forces"] = forces
        if forces:
            fm = max(x["abs_force_Ha_Bohr"] for x in forces)
            res["max_force_Ha_Bohr"] = fm
            res["max_force_eV_Ang"] = fm * HA_BOHR2EV_ANG

    # ---- optimization history
    m = re.search(r"History of (geometry|cell) optimization", text)
    if m:
        body = text[m.end():]
        hdr = re.search(r"^\s*MD_iter\s+(.+)$", body, re.M)
        cols = ["MD_iter"] + (re.split(r"\s{2,}", hdr.group(1).strip()) if hdr else [])
        rows = []
        started = False
        for line in body.splitlines():
            t = line.split()
            if t and t[0].isdigit() and len(t) >= 5:
                try:
                    rows.append([int(t[0])] + [float(x) for x in t[1:]])
                    started = True
                    continue
                except ValueError:
                    pass
            if started and line.strip() and not line.strip()[0].isdigit():
                break
        res["opt_kind"] = m.group(1)
        res["opt_columns"] = cols
        res["opt_history"] = rows

    # ---- timing / termination
    m = re.findall(r"Elapsed\.Time\.\s+(%s)" % _NUM, text)
    if m:
        res["elapsed_s"] = float(m[-1])
    m = re.search(r"Total Computational Time\s*=\s*\d+\s+(%s)\s+\d+\s+(%s)" % (_NUM, _NUM), text)
    if m:
        res["total_time_max_s"] = float(m.group(2))
    res["out_complete"] = "Computational Time (second)" in text

    # ---- .std and .md companions
    std = std or (stem + ".std")
    if os.path.isfile(std):
        s = parse_std(std)
        res["std_file"] = std
        res["finished"] = s["finished"]
        res["errors"] = s["errors"]
    else:
        res["finished"] = res["out_complete"]
    md = md or (stem + ".md")
    md_steps = []
    if os.path.isfile(md):
        try:
            md_steps = parse_md(md)
        except (OSError, ValueError):
            md_steps = []
        res["md_file"] = md
        res["md_nsteps"] = len(md_steps)

    # ---- convergence verdict
    flags = []
    for i, h in enumerate(hist):
        n = len(h["steps"])
        conv = n < maxiter
        if i < len(md_steps) and "scf_conv" in md_steps[i]:
            conv = conv and md_steps[i]["scf_conv"] == 1
        flags.append(conv)
        h["converged"] = conv
        h["nsteps"] = n
        if h["steps"]:
            h["last_NormRD"] = h["steps"][-1]["NormRD"]
            if "dUele" in h["steps"][-1]:
                h["last_dUele"] = h["steps"][-1]["dUele"]
    res["scf_converged"] = bool(flags) and flags[-1]
    res["scf_all_converged"] = bool(flags) and all(flags)
    res["scf_nsteps_last"] = len(hist[-1]["steps"]) if hist else 0
    return res


def parse_any(path: str, data_path: str | None = None) -> dict:
    """parse_out() for .out files, parse_std() for .std files."""
    if path.endswith(".std"):
        out = path[:-4] + ".out"
        if os.path.isfile(out):
            return parse_out(out, data_path=data_path, std=path)
        s = parse_std(path)
        s["system_name"] = os.path.basename(path[:-4])
        last = s["scf_std"][-1] if s["scf_std"] else {}
        s["scf_nsteps_last"] = sum(1 for x in s["scf_std"] if x.get("md") == last.get("md"))
        return s
    return parse_out(path, data_path=data_path)


def quick_utot(path: str) -> float | None:
    """Last 'Utot.' value of a .out file in Hartree (None if missing)."""
    try:
        text = _read(path)
    except OSError:
        return None
    m = re.findall(r"^\s*Utot\.\s+(%s)\s*$" % _NUM, text, re.M)
    return float(m[-1]) if m else None


# ---------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------
def summary_text(r: dict) -> str:
    L = ["==== %s  (System.Name = %s)" % (r.get("file"), r.get("system_name"))]
    if "scf_std" in r and "scf_history" not in r:
        L.append("  .std only: finished = %s, SCF lines = %d" % (r["finished"], len(r["scf_std"])))
        if r["scf_std"]:
            last = r["scf_std"][-1]
            L.append("  last: MD=%s SCF=%s NormRD=%s dUele=%s" % (last.get("md"), last.get("scf"),
                     last.get("NormRD"), last.get("dUele")))
        for e in r.get("errors", []):
            L.append("  error/warning: %s" % e)
        return "\n".join(L)
    L.append("  finished normally : %s%s" % (r.get("finished"), "" if "std_file" in r else
                                            " (no .std found; judged from the .out)"))
    if "utot_Ha" in r:
        L.append("  Utot              : %.10f Ha = %.6f eV" % (r["utot_Ha"], r["utot_eV"]))
    if "chemical_potential_Ha" in r:
        L.append("  chemical potential: %.8f Ha = %.4f eV" % (r["chemical_potential_Ha"],
                                                           r["chemical_potential_eV"]))
    if "gap_eV" in r:
        L.append("  HOMO-LUMO gap     : %.4f eV (HOMO index %d, Cluster)" % (r["gap_eV"], r["homo_index"]))
    if "gap_scf_kgrid_eV" in r:
        L.append("  gap on SCF k-grid : %.4f eV (estimate; use a band/DOS calculation for the real gap)"
                 % r["gap_scf_kgrid_eV"])
    if r.get("scf_history"):
        h = r["scf_history"][-1]
        L.append("  SCF (last MD=%d)   : %d steps (scf.maxIter %d), last NormRD %.3e, |dUele| %s -> %s"
                 % (h["md"], h["nsteps"], r["scf_maxIter"], h.get("last_NormRD", float("nan")),
                    ("%.3e" % h["last_dUele"]) if "last_dUele" in h else "n/a",
                    "converged" if h["converged"] else "NOT CONVERGED"))
        if len(r["scf_history"]) > 1:
            bad = [x["md"] for x in r["scf_history"] if not x["converged"]]
            L.append("  MD/opt steps      : %d, SCF not converged at MD = %s"
                     % (len(r["scf_history"]), bad if bad else "none"))
    if "ngrid" in r:
        L.append("  grid              : %s (required cutoff %s Ry)" % ("x".join(map(str, r["ngrid"])),
                                                                     r.get("required_cutoff_Ry")))
    if "total_spin_moment_muB" in r:
        L.append("  total spin moment : %.4f muB (OpenMX convention)" % r["total_spin_moment_muB"])
    if r.get("mulliken"):
        L.append("  Mulliken          : idx species      up       down      sum     moment   charge")
        for a in r["mulliken"][:50]:
            L.append("                      %3d %-6s %9.4f %9.4f %9.4f %9.4f  %s"
                     % (a["index"], a["species"], a["up"], a["down"], a["sum"], a["moment"],
                        ("%8.4f" % a["charge"]) if "charge" in a else "   n/a"))
        if len(r["mulliken"]) > 50:
            L.append("                      ... (%d atoms, see --json/--csv-dir)" % len(r["mulliken"]))
    if "dipole_abs_Debye" in r:
        L.append("  dipole            : %.4f Debye" % r["dipole_abs_Debye"])
    if "max_force_Ha_Bohr" in r:
        L.append("  max |F|           : %.6f Ha/Bohr = %.4f eV/Ang" % (r["max_force_Ha_Bohr"], r["max_force_eV_Ang"]))
    if r.get("opt_history"):
        last = r["opt_history"][-1]
        L.append("  %s optimization  : %d steps, last |Fmax| %.6f Ha/Bohr, last Utot %.8f Ha"
                 % (r["opt_kind"], len(r["opt_history"]), last[2], last[4]))
    if "elapsed_s" in r:
        L.append("  elapsed time      : %.1f s" % r["elapsed_s"])
    for e in r.get("errors", []):
        L.append("  error/warning     : %s" % e)
    return "\n".join(L)


def write_csv(path: str, header: list, rows: list) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def write_tables(r: dict, outdir: str) -> list:
    os.makedirs(outdir, exist_ok=True)
    name = r.get("system_name", "omx")
    written = []
    if r.get("scf_history"):
        rows = [[h["md"], s["scf"], s["NormRD"], s["Uele"], s.get("dUele", "")]
                for h in r["scf_history"] for s in h["steps"]]
        p = os.path.join(outdir, name + "_scf.csv")
        write_csv(p, ["md", "scf", "NormRD", "Uele_Ha", "dUele_Ha"], rows)
        written.append(p)
    if r.get("energy_blocks"):
        keys = sorted({k for b in r["energy_blocks"] for k in b["terms"]})
        rows = [[b["md"]] + [b["terms"].get(k, "") for k in keys] for b in r["energy_blocks"]]
        p = os.path.join(outdir, name + "_energies.csv")
        write_csv(p, ["md"] + [k + "_Ha" for k in keys], rows)
        written.append(p)
    if r.get("mulliken"):
        rows = [[a["index"], a["species"], a["up"], a["down"], a["sum"], a["moment"],
                 a.get("theta_deg", ""), a.get("phi_deg", ""), a.get("valence", ""), a.get("charge", "")]
                for a in r["mulliken"]]
        p = os.path.join(outdir, name + "_mulliken.csv")
        write_csv(p, ["index", "species", "up", "down", "sum", "moment_muB", "theta_deg", "phi_deg",
                      "valence", "charge"], rows)
        written.append(p)
    if r.get("forces"):
        rows = [[f["index"], f["species"]] + f["xyz_Ang"] + f["force_Ha_Bohr"] + [f["abs_force_Ha_Bohr"]]
                for f in r["forces"]]
        p = os.path.join(outdir, name + "_forces.csv")
        write_csv(p, ["index", "species", "x_Ang", "y_Ang", "z_Ang", "Fx_Ha_Bohr", "Fy_Ha_Bohr",
                      "Fz_Ha_Bohr", "absF_Ha_Bohr"], rows)
        written.append(p)
    if r.get("opt_history"):
        p = os.path.join(outdir, name + "_opt.csv")
        ncol = max(len(x) for x in r["opt_history"])
        cols = r.get("opt_columns", [])
        header = [(cols[i] if i < len(cols) else "col%d" % i) for i in range(ncol)]
        write_csv(p, header, r["opt_history"])
        written.append(p)
    return written


SUMMARY_FIELDS = ["file", "system_name", "finished", "scf_converged", "scf_nsteps_last", "scf_maxIter",
                  "utot_Ha", "utot_eV", "chemical_potential_eV", "total_spin_moment_muB", "gap_eV",
                  "max_force_eV_Ang", "ngrid", "input_energycutoff", "input_kgrid", "elapsed_s"]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Parse OpenMX output files (System.Name.out / .std / .md) and print a summary, "
                    "JSON or CSV tables (Utot, chemical potential, SCF history, Mulliken, forces, "
                    "optimization history).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python omx_parse_out.py met.out\n"
               "  python omx_parse_out.py si_ecut_*.out --summary-csv ecut.csv\n"
               "  python omx_parse_out.py fe.out --json fe.json --csv-dir tables\n"
               "  python omx_parse_out.py crashed.std")
    p.add_argument("files", nargs="+", help="System.Name.out (or .std) files")
    p.add_argument("--data-path", default=None,
                   help="DFT_DATA19 directory for valence charges (default: DATA.PATH of the input echo)")
    p.add_argument("--json", default=None, metavar="FILE", help="write all parsed data as JSON ('-' = stdout)")
    p.add_argument("--csv-dir", default=None, metavar="DIR",
                   help="write per-file tables <name>_scf/_energies/_mulliken/_forces/_opt.csv")
    p.add_argument("--summary-csv", default=None, metavar="FILE", help="one row per input file")
    p.add_argument("-q", "--quiet", action="store_true", help="no text summary")
    p.add_argument("--strict", action="store_true",
                   help="exit status 2 if any run did not finish normally or its last SCF did not converge")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    results, bad = [], False
    for f in args.files:
        try:
            r = parse_any(f, data_path=args.data_path)
        except FileNotFoundError as exc:
            print("ERROR: %s" % exc, file=sys.stderr)
            bad = True
            continue
        results.append(r)
        if not r.get("finished") or ("scf_converged" in r and not r["scf_converged"]):
            bad = True
        if not args.quiet:
            print(summary_text(r))
        if args.csv_dir:
            for p in write_tables(r, args.csv_dir):
                print("  wrote %s" % p, file=sys.stderr)
    if not results:
        return 1
    if args.json:
        data = results[0] if len(results) == 1 else results
        if args.json == "-":
            json.dump(data, sys.stdout, indent=1)
            print()
        else:
            with open(args.json, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=1)
            print("wrote %s" % args.json, file=sys.stderr)
    if args.summary_csv:
        rows = []
        for r in results:
            row = []
            for k in SUMMARY_FIELDS:
                v = r.get(k, "")
                row.append("x".join(map(str, v)) if isinstance(v, list) else v)
            rows.append(row)
        write_csv(args.summary_csv, SUMMARY_FIELDS, rows)
        print("wrote %s" % args.summary_csv, file=sys.stderr)
    return 2 if (args.strict and bad) else 0


if __name__ == "__main__":
    sys.exit(main())
