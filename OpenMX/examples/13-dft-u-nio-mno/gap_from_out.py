#!/usr/bin/env python3
"""Summarize OpenMX .out files of DFT(+U) band calculations.

For every .out file this script reports
  * the largest Hubbard U found in the echoed <Hubbard.U.values block (0 if +U is off),
  * Utot and Uhub (Hartree),
  * VBM, CBM and band gap estimated from the eigenvalues printed on the SCF k-mesh
    (overall and per spin channel), relative to the chemical potential,
  * Mulliken spin moments (column "Diff", muB) of every atom.

Usage:
  python3 gap_from_out.py nio_gga.out nio_dudarev_u5.out
  python3 gap_from_out.py --csv scan_U_summary.csv nio_u0.out nio_u2.out ...

Notes:
  * Only the k-points of the SCF mesh are used, so the gap is an estimate.
    Confirm with a band structure (Band.dispersion) or DOS calculation.
  * Written for collinear (scf.SpinPolarization on/off) Band calculations of
    OpenMX 3.9/4.0. Requires only the Python standard library.
"""
import argparse
import csv
import re
import sys

HARTREE2EV = 27.211386245988
FLOAT = r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?"


def max_hubbard_u(lines):
    """Largest U (eV) in the echoed input; 0.0 if scf.Hubbard.U is not on."""
    hub_on = False
    for ln in lines:
        m = re.match(r"^\s*scf\.Hubbard\.U\s+(\S+)", ln, re.IGNORECASE)
        if m:
            hub_on = m.group(1).lower() == "on"
            break
    if not hub_on:
        return 0.0
    umax, inside = 0.0, False
    for ln in lines:
        s = ln.split("#", 1)[0].strip()
        if s.startswith("<Hubbard.U.values"):
            inside = True
            continue
        if s.startswith("Hubbard.U.values>"):
            break
        if inside and s:
            tokens = s.split()[1:]          # drop species name
            for i in range(1, len(tokens), 2):
                try:
                    umax = max(umax, float(tokens[i]))
                except ValueError:
                    pass
    return umax


def parse_out(fname):
    with open(fname, errors="replace") as f:
        lines = f.read().splitlines()
    res = {"file": fname, "U": max_hubbard_u(lines)}

    for ln in lines:
        m = re.match(r"^\s*Utot\.\s+(" + FLOAT + r")", ln)
        if m:
            res["Utot_Ha"] = float(m.group(1))
        m = re.match(r"^\s*Uhub\.\s+(" + FLOAT + r")", ln)
        if m:
            res["Uhub_Ha"] = float(m.group(1))

    # eigenvalue block (last occurrence)
    start = None
    for i, ln in enumerate(lines):
        if "Eigenvalues (Hartree) of SCF KS-eq." in ln:
            start = i
    mu, eigs = None, [[], []]
    if start is not None:
        for ln in lines[start + 1:]:
            if "Mulliken populations" in ln or "Dipole moment" in ln:
                break
            m = re.search(r"Chemical Potential \(Hartree\)\s*=\s*(" + FLOAT + r")", ln)
            if m:
                mu = float(m.group(1))
                continue
            m = re.match(r"^\s*\d+\s+(" + FLOAT + r")(?:\s+(" + FLOAT + r"))?\s*$", ln)
            if m:
                eigs[0].append(float(m.group(1)))
                if m.group(2) is not None:
                    eigs[1].append(float(m.group(2)))
    res["mu_Ha"] = mu

    def edges(values):
        occ = [e for e in values if e < mu]
        emp = [e for e in values if e > mu]
        if not occ or not emp:
            return None, None
        return max(occ), min(emp)

    if mu is not None and eigs[0]:
        chans = [eigs[0]] + ([eigs[1]] if eigs[1] else [])
        allv = eigs[0] + eigs[1]
        vbm, cbm = edges(allv)
        if vbm is not None:
            res["vbm_eV"] = (vbm - mu) * HARTREE2EV
            res["cbm_eV"] = (cbm - mu) * HARTREE2EV
            res["gap_eV"] = (cbm - vbm) * HARTREE2EV
        for tag, ch in zip(("up", "dn"), chans):
            v, c = edges(ch)
            if v is not None:
                res["gap_%s_eV" % tag] = (c - v) * HARTREE2EV

    # Mulliken spin moments (first Mulliken block)
    moments = []
    for i, ln in enumerate(lines):
        if "Mulliken populations" in ln:
            for ln2 in lines[i + 1:]:
                if "Sum of MulP" in ln2:
                    break
                m = re.match(r"^\s*(\d+)\s+([A-Za-z]\S*)\s+(" + FLOAT + r")\s+(" + FLOAT +
                             r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")\s*$", ln2)
                if m:
                    moments.append((int(m.group(1)), m.group(2), float(m.group(6))))
            break
    res["moments"] = moments
    return res


def main():
    ap = argparse.ArgumentParser(description="Gap / moments / energies from OpenMX .out files")
    ap.add_argument("files", nargs="+", help="OpenMX .out files")
    ap.add_argument("--csv", help="also write a CSV table to this file")
    args = ap.parse_args()

    results = [parse_out(f) for f in args.files]
    nat = max((len(r["moments"]) for r in results), default=0)

    print("%-28s %6s %16s %10s %8s %8s %8s  moments(muB)" %
          ("file", "U(eV)", "Utot(Ha)", "Uhub(Ha)", "gap", "gap_up", "gap_dn"))
    for r in results:
        def fmt(key, f="%8.3f"):
            return f % r[key] if r.get(key) is not None else " " * (len(f % 0.0) - 3) + "n/a"
        mom = " ".join("%s%d:%+.3f" % (sp, n, d) for n, sp, d in r["moments"][:8])
        if len(r["moments"]) > 8:
            mom += " ... (%d atoms, see --csv)" % len(r["moments"])
        print("%-28s %6.2f %16s %10s %8s %8s %8s  %s" % (
            r["file"][-28:], r["U"], fmt("Utot_Ha", "%16.8f"), fmt("Uhub_Ha", "%10.6f"),
            fmt("gap_eV"), fmt("gap_up_eV"), fmt("gap_dn_eV"), mom))
        if r.get("gap_eV") is not None and r["gap_eV"] < 0.05:
            print("   warning: gap < 0.05 eV on the SCF mesh -> metallic or not converged?")
        if r.get("mu_Ha") is None:
            print("   warning: no eigenvalue block found (cluster run or level.of.fileout too low?)")

    if args.csv:
        keys = ["file", "U", "Utot_Ha", "Uhub_Ha", "mu_Ha", "vbm_eV", "cbm_eV",
                "gap_eV", "gap_up_eV", "gap_dn_eV"] + ["moment_atom%d" % (i + 1) for i in range(nat)]
        with open(args.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in results:
                row = {k: r.get(k) for k in keys if not k.startswith("moment_atom")}
                for n, sp, d in r["moments"]:
                    row["moment_atom%d" % n] = d
                w.writerow(row)
        print("CSV written to", args.csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
