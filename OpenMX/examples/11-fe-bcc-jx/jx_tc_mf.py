#!/usr/bin/env python3
"""Group the exchange couplings printed by OpenMX 'jx' into neighbour shells
and estimate the mean-field Curie temperature from the listed pairs.

Usage:
    python jx_tc_mf.py jx.log Fe_bcc_jx.dat

jx.log      : standard output of  jx System.Name.scfout jx.config
              (bulk mode, Flag.PeriodicSum off: lines "i j c1 c2 c3 J[meV] J[mRy] ...")
input .dat  : the OpenMX input of the SCF run (unit vectors and coordinates are read
              from <Atoms.UnitVectors and <Atoms.SpeciesAndCoordinates)

Mean-field estimate (OpenMX manual, Eq. 48 of the jx section):
    k_B T_C = (2/3) * lambda_max( Jt ),  Jt_ij = sum_R J(i0,jR)  (R = 0 excluded for i = j)
Only the pairs present in jx.log enter the sum, i.e. the result is truncated at the
listed shells. Long-ranged (RKKY-like) couplings in metals make this a rough estimate.
"""
import re
import sys
from collections import defaultdict

import numpy as np

KB_MEV = 0.08617333262   # Boltzmann constant in meV/K
BOHR = 0.529177210903     # Angstrom


def read_block(lines, name):
    start = None
    for k, line in enumerate(lines):
        s = line.strip()
        if start is None and s.lower().startswith("<" + name.lower()):
            start = k + 1
        elif start is not None and s.lower().startswith(name.lower() + ">"):
            return lines[start:k]
    raise ValueError(f"block <{name} not found")


def read_keyword(lines, key, default=None):
    for line in lines:
        tok = line.split()
        if tok and tok[0].lower() == key.lower() and len(tok) > 1:
            return tok[1]
    return default


def read_structure(datfile):
    lines = [l.split("#")[0] for l in open(datfile, encoding="ascii", errors="replace")]
    cell = np.array([[float(x) for x in l.split()[:3]]
                     for l in read_block(lines, "Atoms.UnitVectors") if l.strip()])
    if read_keyword(lines, "Atoms.UnitVectors.Unit", "Ang").lower() == "au":
        cell *= BOHR
    unit = read_keyword(lines, "Atoms.SpeciesAndCoordinates.Unit", "Ang").lower()
    pos = []
    for l in read_block(lines, "Atoms.SpeciesAndCoordinates"):
        tok = l.split()
        if len(tok) >= 5:
            pos.append([float(x) for x in tok[2:5]])
    pos = np.array(pos)
    if unit == "frac":
        pos = pos @ cell
    elif unit == "au":
        pos = pos * BOHR
    return cell, pos


def read_jx(logfile):
    pat = re.compile(r"^\s*(\d+)\s+(\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s+")
    out = []
    for line in open(logfile, errors="replace"):
        m = pat.match(line)
        if m:
            i, j, c1, c2, c3 = (int(m.group(k)) for k in range(1, 6))
            out.append((i, j, (c1, c2, c3), float(m.group(6))))
    return out


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    cell, pos = read_structure(sys.argv[2])
    pairs = read_jx(sys.argv[1])
    if not pairs:
        print("no bulk-mode J lines (i j c1 c2 c3 J[meV] ...) found in", sys.argv[1])
        sys.exit(1)

    shells = defaultdict(list)
    print(f"{'i':>3} {'j':>3} {'R (l1 l2 l3)':>14} {'dist(Ang)':>10} {'J(meV)':>12}")
    for i, j, R, J in pairs:
        d = np.linalg.norm(pos[j - 1] + np.array(R) @ cell - pos[i - 1])
        shells[(i, j, round(d, 3))].append(J)
        print(f"{i:3d} {j:3d} {R[0]:4d}{R[1]:5d}{R[2]:5d} {d:10.4f} {J:12.5f}")

    print("\nshell summary (same i, j and distance):")
    print(f"{'i':>3} {'j':>3} {'dist(Ang)':>10} {'n':>4} {'<J>(meV)':>12} {'spread(meV)':>12}")
    for (i, j, d), Js in sorted(shells.items(), key=lambda x: (x[0][2], x[0][0], x[0][1])):
        Js = np.array(Js)
        print(f"{i:3d} {j:3d} {d:10.4f} {len(Js):4d} {Js.mean():12.5f} {Js.max() - Js.min():12.5f}")

    atoms = sorted({p[0] for p in pairs} | {p[1] for p in pairs})
    idx = {a: k for k, a in enumerate(atoms)}
    M = np.zeros((len(atoms), len(atoms)))
    for i, j, R, J in pairs:
        if i == j and R == (0, 0, 0):
            continue
        M[idx[i], idx[j]] += J
    # if only (i, j) and not (j, i) was listed, use J(j0,iR') = J(i0,j(-R)) to fill the gap
    for a in range(len(atoms)):
        for b in range(len(atoms)):
            if a != b and M[a, b] == 0.0 and M[b, a] != 0.0:
                M[a, b] = M[b, a]
    M = 0.5 * (M + M.T)
    lam = np.linalg.eigvalsh(M).max()
    print(f"\nmagnetic atoms: {atoms}")
    print(f"sum of listed J per sublattice matrix, largest eigenvalue = {lam:.4f} meV")
    print(f"mean-field Tc (truncated at the listed shells) = {2.0 * lam / (3.0 * KB_MEV):.0f} K")


if __name__ == "__main__":
    main()
