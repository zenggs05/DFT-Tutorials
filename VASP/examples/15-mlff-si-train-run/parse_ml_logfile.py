#!/usr/bin/env python3
"""Summarise a VASP ML_LOGFILE (machine-learned force fields, VASP >= 6.3).

    python parse_ml_logfile.py                    # ./ML_LOGFILE
    python parse_ml_logfile.py -f train/ML_LOGFILE --plot

Every data line in ML_LOGFILE starts with a keyword (STATUS, ERR, BEEF, LCONF,
NORME, SFF, ...).  The column names are documented in header lines such as
"# ERR   nstep   rmse_energy   rmse_force   rmse_stress"; this script reads
those names from the file itself, so it adapts to the VASP version.

Prints
  * how many MD steps were 'accurate', 'threshold', 'learning', 'critical'
    or 'predict' (STATUS lines) and how many ab initio calculations were done;
  * the last ERR line (RMSE of energy/force/stress on the training set);
  * the last BEEF line (Bayesian error estimate and current ML_CTIFOR);
  * the last LCONF line (local reference configurations per species);
  * spilling-factor statistics if SFF lines are present (ML_MODE = run).
Writes <keyword>.dat for every keyword found; --plot writes ml_errors.png.
Only numpy (and matplotlib for --plot) is required.
"""
import argparse
import collections
import sys

import numpy as np


def parse(path):
    headers = {}
    data = collections.defaultdict(list)
    with open(path, errors="replace") as fh:
        for line in fh:
            tok = line.split()
            if not tok:
                continue
            if tok[0] == "#" and len(tok) > 2 and tok[2] == "nstep":
                headers[tok[1]] = tok[2:]          # column names after the keyword
            elif tok[0].isupper() and tok[0].isalpha() and len(tok) > 1 and tok[1].isdigit():
                data[tok[0]].append(tok[1:])
    return headers, data


def numeric(rows):
    out = []
    for r in rows:
        try:
            out.append([float(x) for x in r])
        except ValueError:
            return None
    return np.array(out) if out else None


def col(headers, key, name, fallback):
    names = headers.get(key, [])
    return names.index(name) if name in names else fallback


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--file", default="ML_LOGFILE")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    headers, data = parse(args.file)
    if not data:
        sys.exit(f"no data lines found in {args.file}")
    print(f"keywords found: {', '.join(sorted(data))}")

    if "STATUS" in data:
        ist = col(headers, "STATUS", "state", 1)
        iab = col(headers, "STATUS", "doabin", 3)
        states = collections.Counter(r[ist] for r in data["STATUS"])
        nab = sum(1 for r in data["STATUS"] if r[iab] == "T")
        print(f"\nMD steps: {len(data['STATUS'])}   ab initio calculations: {nab}")
        for s, n in states.most_common():
            print(f"  {s:10s} {n:8d}")

    for key in ("ERR", "BEEF", "NORME", "LCONF"):
        if key in data:
            names = headers.get(key, [])
            print(f"\nlast {key} line" + (f" ({' '.join(names)})" if names else "") + ":")
            print("  " + "  ".join(data[key][-1]))

    for key, rows in data.items():
        arr = numeric(rows)
        if arr is not None:
            hdr = " ".join(headers.get(key, []))
            np.savetxt(f"{key}.dat", arr, fmt="%.8e", header=hdr)

    sff = numeric(data.get("SFF", []))
    if sff is not None and sff.shape[1] > 1:
        print(f"\nspilling factor (SFF, column '{(headers.get('SFF') or ['', '?'])[1]}'): "
              f"max {sff[:, 1].max():.3f}, mean {sff[:, 1].mean():.3f} over {len(sff)} samples")
        if sff[:, 1].max() > 0.5:
            print("  -> some structures are poorly described: consider continuing the training there")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4))
        err = numeric(data.get("ERR", []))
        if err is not None:
            i = col(headers, "ERR", "rmse_force", 2)
            ax.plot(err[:, 0], err[:, i], label="ERR rmse_force")
        beef = numeric(data.get("BEEF", []))
        if beef is not None:
            i = col(headers, "BEEF", "bee_max_force", 2)
            j = col(headers, "BEEF", "threshold", 4)
            ax.plot(beef[:, 0], beef[:, i], lw=0.8, label="BEEF bee_max_force")
            ax.plot(beef[:, 0], beef[:, j], label="threshold (ML_CTIFOR)")
        ax.set_xlabel("MD step")
        ax.set_ylabel("force error (eV/A)")
        ax.legend()
        fig.tight_layout()
        fig.savefig("ml_errors.png", dpi=150)
        print("\nwrote ml_errors.png")


if __name__ == "__main__":
    main()
