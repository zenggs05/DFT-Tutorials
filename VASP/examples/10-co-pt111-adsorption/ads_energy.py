#!/usr/bin/env python3
"""Adsorption energy from three VASP runs.

    E_ads = E(slab+ads) - E(slab) - n * E(gas)

The energy of each run is the last "energy(sigma->0)" in its OUTCAR
(i.e. the extrapolated zero-smearing energy of the final ionic step).
Negative E_ads means the adsorbed state is more stable than slab + gas.

Examples
    python ads_energy.py --slab 1_slab --ads 3_co_top --gas 2_co_gas
    python ads_energy.py --slab 1_slab --ads 3_co_top 4_co_fcc --gas 2_co_gas
    # dissociative adsorption of H2 -> 2H*: E_ads per H uses n = 0.5
    python ads_energy.py --slab slab --ads H_fcc --gas H2_gas -n 0.5

Each argument may be a directory (OUTCAR is read from it) or an OUTCAR file.
Only the Python standard library is needed.
"""
import argparse
import os
import re
import sys

PAT = re.compile(r"energy\(sigma->0\)\s*=\s*([-+0-9.Ee]+)")


def last_energy(path):
    f = os.path.join(path, "OUTCAR") if os.path.isdir(path) else path
    if not os.path.isfile(f):
        sys.exit("ERROR: %s not found" % f)
    e = None
    with open(f, errors="replace") as fh:
        for line in fh:
            m = PAT.search(line)
            if m:
                e = float(m.group(1))
    if e is None:
        sys.exit("ERROR: no 'energy(sigma->0)' line in %s (run not finished?)" % f)
    with open(f, errors="replace") as fh:
        txt = fh.read()
    ib = re.search(r"IBRION\s*=\s*(-?\d+)", txt)
    nsw = re.search(r"NSW\s*=\s*(\d+)", txt)
    if (ib and int(ib.group(1)) == -1) or (nsw and int(nsw.group(1)) == 0):
        conv = "single point"
    elif "reached required accuracy" in txt:
        conv = "relaxed"
    else:
        conv = "NOT CONVERGED"
    return e, conv


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slab", required=True, help="clean slab run")
    ap.add_argument("--ads", required=True, nargs="+", help="slab+adsorbate run(s)")
    ap.add_argument("--gas", required=True, help="gas-phase reference run")
    ap.add_argument("-n", type=float, default=1.0,
                    help="number of gas molecules per adsorbate (default 1)")
    args = ap.parse_args()

    es, cs = last_energy(args.slab)
    eg, cg = last_energy(args.gas)
    print("%-28s %16.6f eV  (%s)" % ("E(slab)  " + args.slab, es, cs))
    print("%-28s %16.6f eV  (%s)" % ("E(gas)   " + args.gas, eg, cg))
    print("-" * 66)
    for a in args.ads:
        ea, ca = last_energy(a)
        eads = ea - es - args.n * eg
        print("%-28s %16.6f eV  (%s)" % ("E(ads)   " + a, ea, ca))
        print("%-28s %16.3f eV" % ("  -> E_ads", eads))
    if any("NOT" in c for c in (cs, cg)):
        print("WARNING: a reference run did not reach the force criterion.")


if __name__ == "__main__":
    main()
