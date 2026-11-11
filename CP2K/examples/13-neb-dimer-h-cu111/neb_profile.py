#!/usr/bin/env python3
"""Print (and optionally plot) the last NEB energy profile from a CP2K BAND output.

Reads the compact summary blocks written to the main output by
&BAND/&PROGRAM_RUN_INFO (format of CP2K 2023.x - 2026.2 with PLOT_REL_ENERGY F):

    DISTANCES REP =   ...   (NUMBER_OF_REPLICA - 1 values, bohr, 3N-dim distance)
    ENERGIES [au] =   ...   (NUMBER_OF_REPLICA values, hartree)

Usage:
    python3 neb_profile.py h_diffusion_cineb.out            # table only
    python3 neb_profile.py h_diffusion_cineb.out --png mep.png
Only the Python standard library is needed; --png also needs matplotlib.
"""
import argparse
import re
import sys

HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903
NUM = r"[-+]?\d+\.\d*(?:[EeDd][-+]?\d+)?"


def read_blocks(path):
    """Return the (distances, energies) lists of the last complete summary block."""
    blocks = []
    cur_key, cur_vals, dist = None, [], None
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            m = re.match(r"\s*(DISTANCES REP|ENERGIES \[au\])\s*=(.*)", line)
            if m:
                if cur_key == "ENERGIES [au]":
                    blocks.append((dist, cur_vals))
                elif cur_key == "DISTANCES REP":
                    dist = cur_vals
                cur_key, cur_vals = m.group(1), [float(x.replace("D", "E")) for x in re.findall(NUM, m.group(2))]
                continue
            if cur_key and re.match(r"^\s+" + NUM + r"(\s+" + NUM + r")*\s*$", line):
                cur_vals += [float(x.replace("D", "E")) for x in re.findall(NUM, line)]
                continue
            if cur_key == "DISTANCES REP":
                dist = cur_vals
            elif cur_key == "ENERGIES [au]":
                blocks.append((dist, cur_vals))
            cur_key, cur_vals = None, []
    if cur_key == "ENERGIES [au]":
        blocks.append((dist, cur_vals))
    blocks = [b for b in blocks if b[0] is not None and len(b[0]) == len(b[1]) - 1]
    if not blocks:
        sys.exit("no 'DISTANCES REP' / 'ENERGIES [au]' blocks found "
                 "(PRINT_LEVEL too low, or PLOT_REL_ENERGY T format?)")
    return blocks


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("output", help="main CP2K output of a RUN_TYPE BAND job")
    ap.add_argument("--png", help="write a plot of the profile to this file (needs matplotlib)")
    args = ap.parse_args()

    blocks = read_blocks(args.output)
    dist, ener = blocks[-1]
    print(f"# {len(blocks)} band steps found, showing the last one")
    s = [0.0]
    for d in dist:
        s.append(s[-1] + d * BOHR_ANG)
    e0 = ener[0]
    rel = [(e - e0) * HARTREE_EV for e in ener]
    imax = max(range(1, len(rel) - 1), key=lambda i: rel[i]) if len(rel) > 2 else 0
    print(f"{'replica':>7s} {'s (A, 3N-dim)':>14s} {'E (Ha)':>18s} {'E-E1 (eV)':>11s}")
    for i, (si, e, r) in enumerate(zip(s, ener, rel), 1):
        mark = "  <- highest image" if i - 1 == imax else ""
        print(f"{i:7d} {si:14.4f} {e:18.8f} {r:11.4f}{mark}")
    print(f"# forward barrier  E(max) - E(1)  = {rel[imax]:.4f} eV")
    print(f"# reverse barrier  E(max) - E(N)  = {rel[imax] - rel[-1]:.4f} eV")
    print(f"# reaction energy  E(N)   - E(1)  = {rel[-1]:.4f} eV")

    if args.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.5, 3.2))
        ax.plot(s, rel, "o-", color="#1f5f99")
        ax.plot([s[imax]], [rel[imax]], "o", color="#c0392b")
        ax.set_xlabel("reaction coordinate s (A)")
        ax.set_ylabel("E - E(initial) (eV)")
        fig.tight_layout()
        fig.savefig(args.png, dpi=200)
        print(f"# plot written to {args.png}")


if __name__ == "__main__":
    main()
