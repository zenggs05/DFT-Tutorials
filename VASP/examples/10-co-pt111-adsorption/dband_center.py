#!/usr/bin/env python3
"""d-band centre, width and filling from a VASP DOSCAR (LORBIT = 11).

    eps_d = int E n_d(E) dE / int n_d(E) dE           (first moment)
    W_d   = sqrt( int (E-eps_d)^2 n_d dE / int n_d dE ) (square root of 2nd moment)
    f_d   = int_{E<E_F} n_d dE / int n_d dE             (filling)

Energies are referenced to the Fermi level written in DOSCAR (E_F = 0).
Works for ISPIN = 1 and ISPIN = 2 (non-collinear DOSCAR is not supported).

Examples
    python dband_center.py DOSCAR -a 9-16            # top two Pt layers
    python dband_center.py DOSCAR -a 13-16 -w -10 2  # surface layer, window [-10, 2] eV
    python dband_center.py DOSCAR -a 13-16 --dat pdos_d.dat

Only numpy is required.
"""
import argparse
import sys

import numpy as np


def parse_indices(tokens):
    out = []
    for tok in tokens:
        for part in tok.split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                a, b = part.split("-", 1)
                out += list(range(int(a), int(b) + 1))
            else:
                out.append(int(part))
    return sorted(set(out))


def read_doscar(path):
    with open(path) as fh:
        lines = fh.readlines()
    natoms = int(lines[0].split()[0])
    head = lines[5].split()
    nedos = int(head[2])
    efermi = float(head[3])
    # total DOS block: E, DOS(up) [, DOS(dn)], integrated ...
    tot = np.array([[float(x) for x in l.split()] for l in lines[6:6 + nedos]])
    energies = tot[:, 0] - efermi
    pdos = []
    start = 6 + nedos
    for ia in range(natoms):
        blk_start = start + ia * (nedos + 1) + 1
        if blk_start + nedos > len(lines):
            sys.exit("ERROR: DOSCAR has no projected DOS for all atoms - "
                     "set LORBIT = 11 (and do not use LORBIT = 0).")
        blk = np.array([[float(x) for x in l.split()]
                        for l in lines[blk_start:blk_start + nedos]])
        pdos.append(blk)
    return energies, efermi, pdos


def d_columns(ncol):
    """Return list of column groups (one list per spin) holding the 5 d orbitals.

    LORBIT=11 column layout (first column is energy):
      ISPIN=1:  s py pz px dxy dyz dz2 dxz dx2 [f ...]      -> 10 or 17 columns
      ISPIN=2:  s_u s_d py_u py_d ... dx2_u dx2_d [f ...]   -> 19 or 33 columns
    """
    if ncol in (10, 17):
        return [list(range(5, 10))]
    if ncol in (19, 33):
        up = [9, 11, 13, 15, 17]
        dn = [10, 12, 14, 16, 18]
        return [up, dn]
    sys.exit("ERROR: unexpected number of PDOS columns (%d). LORBIT=11 needed; "
             "non-collinear DOSCAR is not supported." % ncol)


_trap = getattr(np, "trapezoid", None) or np.trapz   # numpy >= 2.0 renamed trapz


def moments(e, n, window):
    m = (e >= window[0]) & (e <= window[1])
    e, n = e[m], n[m]
    norm = _trap(n, e)
    if norm <= 0:
        return float("nan"), float("nan"), float("nan"), 0.0
    c = _trap(n * e, e) / norm
    w = np.sqrt(_trap(n * (e - c) ** 2, e) / norm)
    occ = e <= 0.0
    fill = _trap(n[occ], e[occ]) / norm if occ.sum() > 1 else 0.0
    return c, w, fill, norm


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("doscar", nargs="?", default="DOSCAR")
    ap.add_argument("-a", "--atoms", nargs="+", required=True,
                    help="1-based atom indices, ranges allowed (e.g. 9-16 or 13,14)")
    ap.add_argument("-w", "--window", nargs=2, type=float, default=None,
                    metavar=("EMIN", "EMAX"),
                    help="integration window relative to E_F in eV "
                         "(default: whole DOSCAR energy range)")
    ap.add_argument("--dat", default=None,
                    help="write summed d-PDOS (E-E_F, up[, down]) to this file")
    args = ap.parse_args()

    energies, efermi, pdos = read_doscar(args.doscar)
    atoms = parse_indices(args.atoms)
    if atoms[0] < 1 or atoms[-1] > len(pdos):
        sys.exit("ERROR: atom index out of range 1..%d" % len(pdos))
    groups = d_columns(pdos[0].shape[1])
    window = args.window or (energies.min(), energies.max())

    spin_dos = []
    for cols in groups:
        s = np.zeros_like(energies)
        for ia in atoms:
            s += pdos[ia - 1][:, cols].sum(axis=1)
        spin_dos.append(s)

    print("DOSCAR: %s   E_F = %.4f eV   atoms: %s" % (args.doscar, efermi,
                                                      " ".join(map(str, atoms))))
    print("window (E - E_F): [%.2f, %.2f] eV" % tuple(window))
    labels = ["total"] if len(spin_dos) == 1 else ["up", "down"]
    for lab, n in zip(labels, spin_dos):
        c, w, f, norm = moments(energies, n, window)
        print("  %-5s  eps_d = %8.3f eV   W_d = %6.3f eV   filling = %5.3f   "
              "N_d(window) = %7.3f" % (lab, c, w, f, norm))
    if len(spin_dos) == 2:
        c, w, f, norm = moments(energies, spin_dos[0] + spin_dos[1], window)
        print("  %-5s  eps_d = %8.3f eV   W_d = %6.3f eV   filling = %5.3f   "
              "N_d(window) = %7.3f" % ("sum", c, w, f, norm))

    if args.dat:
        data = np.column_stack([energies] + spin_dos)
        hdr = "E-E_F(eV) " + " ".join("d_" + l for l in labels)
        np.savetxt(args.dat, data, fmt="%12.5f", header=hdr)
        print("summed d-PDOS written to", args.dat)


if __name__ == "__main__":
    main()
