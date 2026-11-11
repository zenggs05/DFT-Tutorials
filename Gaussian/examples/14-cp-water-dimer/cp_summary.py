#!/usr/bin/env python3
"""Summarize Counterpoise (BSSE) results from a Gaussian 16 log.

For every job step that contains a counterpoise summary, the LAST summary of
that step is used (a CP optimization prints one per geometry step):

    Counterpoise corrected energy =   ...
                      BSSE energy =   ...
                 sum of fragments =   ...

Derived quantities (all in kcal/mol):
    dE_int(raw) = E_raw - sum_of_fragments,  E_raw = E_CP - BSSE
    dE_int(CP)  = E_CP  - sum_of_fragments
    (sum_of_fragments = monomers at the complex geometry in their own basis)
If the relaxed monomer energies are given with --mono (Hartree, one value per
fragment, same method/basis), the deformation and binding energies are added:
    E_def  = sum_of_fragments - sum(E_mono)
    dE_bind(CP) = dE_int(CP) + E_def

Usage:
    python cp_summary.py water_dimer_counterpoise.log
    python cp_summary.py water_dimer_counterpoise.log --mono -76.4xxxx -76.4xxxx
Only the Python standard library is used (Python >= 3.8).
"""
import argparse
import re

HARTREE_TO_KCAL = 627.509474

RE_ROUTE = re.compile(r"^\s*#\S*\s+(.*)$")
RE_CP = re.compile(r"Counterpoise corrected energy\s*=\s*(-?\d+\.\d+)")
RE_BSSE = re.compile(r"BSSE energy\s*=\s*(-?\d+\.\d+)")
RE_SUM = re.compile(r"sum of (?:fragments|monomers)\s*=\s*(-?\d+\.\d+)")
RE_END = re.compile(r"(Normal|Error) termination of Gaussian")


def parse(path):
    steps, cur, route = [], {}, None
    with open(path, errors="replace") as fh:
        for line in fh:
            if route is None:
                m = RE_ROUTE.match(line)
                if m and line.lstrip().startswith("#"):
                    route = line.strip()
            for key, rx in (("cp", RE_CP), ("bsse", RE_BSSE), ("sum", RE_SUM)):
                m = rx.search(line)
                if m:
                    cur[key] = float(m.group(1))
            if RE_END.search(line):
                cur["route"] = route
                cur["status"] = "normal" if "Normal" in line else "ERROR"
                steps.append(cur)
                cur, route = {}, None
    if cur:
        cur["route"] = route
        cur["status"] = "unfinished"
        steps.append(cur)
    return steps


def main():
    ap = argparse.ArgumentParser(description="Counterpoise/BSSE summary of a Gaussian log")
    ap.add_argument("log")
    ap.add_argument("--mono", type=float, nargs="+",
                    help="relaxed monomer energies in Hartree (same level), one per fragment")
    a = ap.parse_args()
    k = HARTREE_TO_KCAL
    found = False
    for i, s in enumerate(parse(a.log), 1):
        if "cp" not in s:
            continue
        found = True
        e_cp, bsse = s["cp"], s.get("bsse", float("nan"))
        e_raw = e_cp - bsse
        print(f"--- job step {i} ({s['status']}): {s['route']}")
        print(f"  E(complex, CP corrected)   = {e_cp:.8f} Eh")
        print(f"  BSSE                       = {bsse:.8f} Eh = {bsse * k:.2f} kcal/mol")
        if "sum" in s:
            e_sum = s["sum"]
            print(f"  dE_int raw                 = {(e_raw - e_sum) * k:8.2f} kcal/mol")
            print(f"  dE_int CP-corrected        = {(e_cp - e_sum) * k:8.2f} kcal/mol")
            if a.mono:
                e_def = e_sum - sum(a.mono)
                print(f"  E_def (deformation)        = {e_def * k:8.2f} kcal/mol")
                print(f"  dE_bind CP-corrected       = {(e_cp - e_sum + e_def) * k:8.2f} kcal/mol")
        else:
            print("  (no 'sum of fragments' line found in this step)")
    if not found:
        print("No counterpoise summary found. Was Counterpoise=N in the route?")


if __name__ == "__main__":
    main()
