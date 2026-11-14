#!/usr/bin/env python3
"""List on-site energies and the largest hoppings stored in an OpenMX .HWR file.

The .HWR file (MLWF run, manual 47.2.5) holds <m,0|H|n,R> on the lattice
vectors R of the Wigner-Seitz supercell:
    line 1   : comment
    line 2   : "Number of Wannier Function  N"
    line 3   : "Number of Wigner-Seitz supercell  NR"
    line 4-7 : "Lattice vector (in Bohr)" + 3 vectors
    line 8   : collinear / non-collinear, spinsize
    line 9   : "Fermi level  E_F" (Hartree)
    then     : for each spin and R: "R ( n1 n2 n3 )  degeneracy" followed by
               N*N lines "m n Re Im" (Hartree)
Energies are printed in eV; on-site energies relative to E_F.
The degeneracy is printed as stored; R = (0,0,0) always has degeneracy 1.

Usage:
    python3 hwr_hopping.py graphene_pz.HWR            # 15 largest |H| with R != 0 or m != n
    python3 hwr_hopping.py graphene_pz.HWR --top 30
Standard library only.
"""
import argparse
import math
import re
import sys

HA2EV = 27.2113845
NUM = r"[-+]?\d+\.?\d*(?:[EeDd][-+]?\d+)?"


def read_hwr(fname):
    with open(fname, errors="replace") as f:
        lines = f.read().splitlines()
    nwf = int(re.findall(r"\d+", lines[1])[-1])
    cell = [[float(v) for v in lines[i].split()[:3]] for i in (4, 5, 6)]
    efermi = float(re.findall(NUM, lines[8])[-1])
    nr = int(re.findall(r"\d+", lines[2])[-1])
    blocks, cur = [], None
    for line in lines[9:]:
        if line.lstrip().startswith("R"):
            m = re.match(r"\s*R\s*\(\s*(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*\)\s*(\d+)", line)
            if not m:
                continue
            r = tuple(int(m.group(i)) for i in (1, 2, 3))
            # spin is the outermost loop: the first nr blocks are spin 0
            cur = dict(spin=len(blocks) // nr, R=r, deg=int(m.group(4)), H={})
            blocks.append(cur)
        elif cur is not None:
            p = line.split()
            if len(p) >= 4:
                cur["H"][(int(p[0]), int(p[1]))] = complex(float(p[2]), float(p[3]))
    return nwf, cell, efermi, blocks


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("hwr", help="System.Name.HWR")
    ap.add_argument("--top", type=int, default=15, help="number of hoppings to list")
    a = ap.parse_args()

    nwf, cell, ef, blocks = read_hwr(a.hwr)
    if not blocks:
        sys.exit("no R blocks found in %s" % a.hwr)
    print("Wannier functions: %d   R vectors: %d   E_F = %.5f Ha" % (nwf, len(blocks), ef))
    for b in blocks:
        if b["R"] == (0, 0, 0):
            print("\nspin %d on-site energies (eV, relative to E_F):" % b["spin"])
            for m in range(1, nwf + 1):
                h = b["H"].get((m, m), 0j)
                print("  WF %3d : %9.4f" % (m, (h.real - ef) * HA2EV))
    hops = []
    for b in blocks:
        rlen = math.sqrt(sum((b["R"][0] * cell[0][i] + b["R"][1] * cell[1][i]
                              + b["R"][2] * cell[2][i]) ** 2 for i in range(3)))
        for (m, n), h in b["H"].items():
            if b["R"] == (0, 0, 0) and m == n:
                continue
            hops.append((abs(h), b["spin"], b["R"], b["deg"], m, n, h, rlen))
    hops.sort(key=lambda t: -t[0])
    print("\nlargest hoppings (eV):")
    print(" spin       R          deg   m   n      Re(H)      Im(H)   |R| (Bohr)")
    for t in hops[:a.top]:
        _, s, r, deg, m, n, h, rlen = t
        print("  %d  (%3d %3d %3d)  %3d  %3d %3d  %9.4f  %9.4f  %8.3f"
              % (s, r[0], r[1], r[2], deg, m, n, h.real * HA2EV, h.imag * HA2EV, rlen))


if __name__ == "__main__":
    main()
