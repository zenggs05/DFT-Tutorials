#!/usr/bin/env python3
"""Scan calculation directories and report the status of VASP / Gaussian /
CP2K / OpenMX jobs.  Standard library only (Python >= 3.9).

Markers used (search the files yourself with grep to double-check):
  VASP     OUTCAR "General timing and accounting"   -> run finished
           OUTCAR "reached required accuracy"       -> ionic relaxation converged
           stdout "ZBRENT: fatal error", "Error EDDDAV", "VERY BAD NEWS" ... -> error
           OSZICAR: #SCF steps of last ionic step == NELM -> SCF not converged
  Gaussian "Normal termination" / "Error termination"
  CP2K     "PROGRAM ENDED AT", "SCF run NOT converged", "ABORT"
  OpenMX   stdout "The SCF was achieved" / "The SCF was not fully achieved",
           .out "Utot." and "Elapsed.Time."

Usage:
  python job_status_scanner.py ROOT [--csv status.csv]
"""
import argparse
import csv
import glob
import os
import re

VASP_ERRORS = {
    "ZBRENT: fatal error": "ZBRENT (ionic step bracketing failed)",
    "Error EDDDAV": "EDDDAV (Davidson/ZHEGV failed)",
    "Sub-Space-Matrix is not hermitian": "Sub-space matrix not hermitian",
    "VERY BAD NEWS": "VERY BAD NEWS (see message)",
    "BRMIX: very serious problems": "BRMIX mixing problem",
    "WARNING: Sub-Space-Matrix": "Sub-space matrix warning",
    "Error reading item": "INCAR/KPOINTS read error",
    "M_divide: can not subdivide": "MPI ranks not divisible (NEB IMAGES/KPAR/NCORE)",
}
HA2EV = 27.211386245988


def tail(path, nbytes=200000):
    """Return the last nbytes of a text file (fast for huge OUTCARs)."""
    with open(path, "rb") as f:
        f.seek(0, os.SEEK_END)
        size = f.tell()
        f.seek(max(0, size - nbytes))
        return f.read().decode("utf-8", errors="ignore")


def read(path):
    with open(path, errors="ignore") as f:
        return f.read()


def vasp_status(d):
    outcar = os.path.join(d, "OUTCAR")
    txt = read(outcar)
    end = tail(outcar)
    info = {"code": "VASP", "status": "RUNNING/KILLED", "energy_eV": "", "steps": "", "note": ""}
    m = re.findall(r"free  energy   TOTEN\s+=\s+(-?\d+\.\d+)", end)
    if m:
        info["energy_eV"] = m[-1]
    if "General timing and accounting" in end:
        info["status"] = "FINISHED"
    nsw = re.search(r"NSW\s*=\s*(\d+)", txt)
    if nsw and int(nsw.group(1)) > 0 and info["status"] == "FINISHED":
        conv = "reached required accuracy" in end
        info["status"] = "CONVERGED" if conv else "NOT_CONVERGED(NSW reached)"
    # electronic convergence from OSZICAR
    osz = os.path.join(d, "OSZICAR")
    nelm = re.search(r"NELM\s*=\s*(\d+)", txt)
    if os.path.exists(osz) and nelm:
        blocks = re.split(r"^\s*\d+ F=", read(osz), flags=re.M)
        last = blocks[-2] if len(blocks) > 1 else blocks[-1]
        nscf = len(re.findall(r"^(DAV|RMM|CG|SDA|ESD|DIA):", last, flags=re.M))
        info["steps"] = str(len(blocks) - 1)
        if nscf >= int(nelm.group(1)):
            info["note"] += "last ionic step hit NELM (SCF not converged); "
    # errors usually go to stdout: vasp.out, slurm-*.out, *.log
    for f in [outcar] + glob.glob(os.path.join(d, "*.out")) + glob.glob(os.path.join(d, "slurm-*")):
        t = tail(f, 50000)
        for key, msg in VASP_ERRORS.items():
            if key in t and msg not in info["note"]:
                info["note"] += msg + "; "
                if not key.startswith("WARNING"):
                    info["status"] = "ERROR"
    return info


def gaussian_status(path):
    t = tail(path)
    info = {"code": "Gaussian", "status": "RUNNING/KILLED", "energy_eV": "", "steps": "", "note": ""}
    m = re.findall(r"SCF Done:\s+E\(\S+\)\s+=\s+(-?\d+\.\d+)", read(path))
    if m:
        info["energy_eV"] = f"{float(m[-1]) * HA2EV:.6f}"
        info["steps"] = str(len(m))
    if "Normal termination" in t[-3000:]:
        info["status"] = "FINISHED"
    elif "Error termination" in t:
        info["status"] = "ERROR"
        link = re.findall(r"Error termination via Lnk1e in \S+/(l\d+)\.exe", t)
        info["note"] = f"failed in {link[-1]}" if link else ""
    if "-- Stationary point found." in t:
        info["note"] += "opt converged; "
    nimag = re.findall(r"(\d+) imaginary frequencies \(negative Signs\)", t)
    if nimag:
        info["note"] += f"{nimag[-1]} imaginary freq; "
    return info


def cp2k_status(path):
    t = read(path)
    info = {"code": "CP2K", "status": "RUNNING/KILLED", "energy_eV": "", "steps": "", "note": ""}
    # 2025.1+: "... energy [hartree]"; 2024.x and older: "... energy [a.u.]:"
    m = re.findall(r"ENERGY\| Total FORCE_EVAL \( \w+ \) energy \[(?:hartree|a\.u\.)\]:?\s+(-?\d+\.\d+)", t)
    if m:
        info["energy_eV"] = f"{float(m[-1]) * HA2EV:.6f}"
        info["steps"] = str(len(m))
    if "PROGRAM ENDED AT" in t:
        info["status"] = "FINISHED"
    if "SCF run NOT converged" in t:
        info["note"] += "SCF run NOT converged; "
    if "ABORT" in t[-5000:]:
        info["status"] = "ERROR"
    if "GEOMETRY OPTIMIZATION COMPLETED" in t:
        info["note"] += "geo-opt converged; "
    return info


def openmx_status(outfile, d):
    t = read(outfile)
    info = {"code": "OpenMX", "status": "RUNNING/KILLED", "energy_eV": "", "steps": "", "note": ""}
    m = re.findall(r"Utot\.\s+(-?\d+\.\d+)", t)
    if m:
        info["energy_eV"] = f"{float(m[-1]) * HA2EV:.6f}"
    if "Elapsed.Time." in t:
        info["status"] = "FINISHED"
    for std in glob.glob(os.path.join(d, "*.std")) + glob.glob(os.path.join(d, "*.log")):
        s = tail(std)
        if "The SCF was not fully achieved" in s:
            info["note"] += "SCF not fully achieved; "
        if "The geometry optimization was achieved" in s:
            info["note"] += "geo-opt converged; "
    return info


def scan(root):
    rows = []
    for d, _dirs, files in sorted(os.walk(root)):
        fs = set(files)
        if "OUTCAR" in fs:
            rows.append((d, vasp_status(d)))
        for f in sorted(files):
            p = os.path.join(d, f)
            if not (f.endswith(".log") or f.endswith(".out")) or f.startswith("slurm"):
                continue
            head = read(p)[:20000]
            if "Gaussian, Inc." in head or "Entering Gaussian System" in head:
                rows.append((p, gaussian_status(p)))
            elif "CP2K|" in head or "**** **** ******  **  PROGRAM STARTED AT" in head:
                rows.append((p, cp2k_status(p)))
            elif "OpenMX" in head and f.endswith(".out"):
                rows.append((p, openmx_status(p, d)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?", default=".")
    ap.add_argument("--csv")
    a = ap.parse_args()
    rows = scan(a.root)
    w = max([len(os.path.relpath(p, a.root)) for p, _ in rows] + [4])
    print(f"{'path':<{w}}  {'code':<8} {'status':<26} {'energy(eV)':>16} {'steps':>5}  note")
    for p, i in rows:
        print(f"{os.path.relpath(p, a.root):<{w}}  {i['code']:<8} {i['status']:<26} "
              f"{i['energy_eV']:>16} {i['steps']:>5}  {i['note']}")
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["path", "code", "status", "energy_eV", "steps", "note"])
            for p, i in rows:
                wr.writerow([p, i["code"], i["status"], i["energy_eV"], i["steps"], i["note"]])
        print("written", a.csv)


if __name__ == "__main__":
    main()
