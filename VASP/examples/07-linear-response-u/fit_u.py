#!/usr/bin/env python3
"""Fit the linear-response U from VASP LDAUTYPE=3 runs.

Reads the d occupation of atom 1 from the LAST "total charge" table of each
OUTCAR (needs LORBIT = 11), fits N_d(V) for the non-self-consistent (nscf_V)
and self-consistent (scf_V) runs, and reports

    chi0 = dN_nscf/dV ,  chi = dN_scf/dV ,  U = 1/chi - 1/chi0

Directory names must look like  nscf_<V>  and  scf_<V>  (V in eV).
The ground state (V = 0) is added to both data sets if --gs is given.

Usage:
    python fit_u.py --gs gs/OUTCAR --dirs nscf_* scf_*  [--atom 1] [--plot]
"""
import argparse
import os
import re
import sys

import numpy as np


def d_occupation(outcar, atom=1, column="d"):
    """Return the occupation in `column` of `atom` from the last 'total charge' table."""
    try:
        lines = open(outcar).read().splitlines()
    except OSError as err:
        sys.exit(f"cannot read {outcar}: {err}")
    starts = [i for i, l in enumerate(lines) if l.strip() == "total charge"]
    if not starts:
        sys.exit(f"{outcar}: no 'total charge' table (LORBIT = 11 set? run finished?)")
    i = starts[-1]
    # find the header line '# of ion  s  p  d  (f)  tot'
    while not lines[i].lstrip().startswith("# of ion"):
        i += 1
    header = lines[i].split()[3:]          # drop '#', 'of', 'ion'
    col = header.index(column) + 1         # +1 because data lines start with ion index
    for line in lines[i + 2:]:
        parts = line.split()
        if parts and parts[0] == str(atom):
            return float(parts[col])
    sys.exit(f"{outcar}: atom {atom} not found in the 'total charge' table")


def collect(dirs, prefix, atom):
    data = []
    for d in dirs:
        name = os.path.basename(os.path.normpath(d))
        m = re.fullmatch(prefix + r"_(-?\d+(?:\.\d+)?)", name)
        if m:
            data.append((float(m.group(1)), d_occupation(os.path.join(d, "OUTCAR"), atom)))
    return data


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dirs", nargs="+", required=True, help="nscf_V and scf_V directories")
    p.add_argument("--gs", help="OUTCAR of the unperturbed ground state (V = 0)")
    p.add_argument("--atom", type=int, default=1, help="index of the perturbed atom (default 1)")
    p.add_argument("--plot", action="store_true", help="save Nd_vs_V.png (needs matplotlib)")
    a = p.parse_args()

    nscf = collect(a.dirs, "nscf", a.atom)
    scf = collect(a.dirs, "scf", a.atom)
    if a.gs:
        n0 = d_occupation(a.gs, a.atom)
        nscf.append((0.0, n0))
        scf.append((0.0, n0))
    if len(nscf) < 2 or len(scf) < 2:
        sys.exit("need at least two points for each of nscf and scf")

    out = {}
    for label, data in (("nscf", nscf), ("scf", scf)):
        v, n = np.array(sorted(data)).T
        slope, icpt = np.polyfit(v, n, 1)
        resid = n - (slope * v + icpt)
        out[label] = (v, n, slope, icpt)
        print(f"{label:5s}: points = {len(v)}, dN/dV = {slope: .5f} 1/eV, "
              f"max |residual| = {np.abs(resid).max():.2e}")
        for vi, ni in zip(v, n):
            print(f"        V = {vi: .3f} eV   N_d = {ni:.4f}")

    chi0 = out["nscf"][2]
    chi = out["scf"][2]
    u = 1.0 / chi - 1.0 / chi0
    print(f"\nchi0 = {chi0:.5f} 1/eV, chi = {chi:.5f} 1/eV")
    print(f"U = 1/chi - 1/chi0 = {u:.2f} eV")

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.5, 3.5))
        for label, marker in (("nscf", "o"), ("scf", "s")):
            v, n, s, c = out[label]
            ax.plot(v, n, marker, label=f"{label}: slope {s:.3f}")
            vv = np.linspace(v.min(), v.max(), 50)
            ax.plot(vv, s * vv + c, "-", lw=1)
        ax.set_xlabel("V (eV)")
        ax.set_ylabel(f"d occupation of atom {a.atom}")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig("Nd_vs_V.png", dpi=200)
        print("saved Nd_vs_V.png")


if __name__ == "__main__":
    main()
