#!/usr/bin/env python3
"""Effective masses from a VASP band calculation on k-points around an extremum.

The k-points (e.g. KPOINTS_fd) contain the extremum k0 plus a few points along
several straight lines starting at k0. For each direction the band energy is fitted to

    E(k0 + dk) = E(k0) + a * |dk|^2  (+ b * |dk|^4 with --quartic)
    m*/m0 = hbar^2 / (2 m0 a) = 3.80998 eV*A^2 / a

with |dk| in 1/Angstrom (2*pi included). Negative m* = hole-like (band maximum).

Input : EIGENVAL (or EIGENVAL_OPT of a KPOINTS_OPT run) + POSCAR of the same calculation.
EIGENVAL k-points are in fractional reciprocal coordinates; POSCAR gives the lattice.
Requires numpy only.

Examples
--------
    python fit_effmass.py                       # CBM and VBM at Gamma
    python fit_effmass.py --band 8 --band 9     # explicit band indices (1-based)
    python fit_effmass.py --npoints 3 --quartic # fewer points, include non-parabolicity
"""
import argparse
import sys

import numpy as np

HB2_2M = 3.80998212   # hbar^2 / (2 m_e) in eV * Angstrom^2


def read_lattice(poscar):
    with open(poscar) as fh:
        lines = fh.readlines()
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
    if scale < 0:
        return lat * (abs(scale) / abs(np.linalg.det(lat))) ** (1 / 3)
    return lat * scale


def read_eigenval(path):
    """Return kfrac (nk,3), energies (nspin,nk,nb), occupations (nspin,nk,nb) or None, nelect."""
    with open(path) as fh:
        lines = fh.readlines()
    ispin = int(lines[0].split()[3])
    nelect, nk, nb = lines[5].split()[:3]
    nelect, nk, nb = float(nelect), int(nk), int(nb)
    kf = np.zeros((nk, 3))
    en = np.zeros((ispin, nk, nb))
    occ = np.zeros((ispin, nk, nb))
    has_occ = True
    i = 6
    for ik in range(nk):
        while not lines[i].split():
            i += 1
        kf[ik] = [float(x) for x in lines[i].split()[:3]]
        i += 1
        for ib in range(nb):
            t = lines[i].split()
            i += 1
            vals = [float(x) for x in t[1:]]
            if ispin == 1:
                en[0, ik, ib] = vals[0]
                if len(vals) >= 2:
                    occ[0, ik, ib] = vals[1]
                else:
                    has_occ = False
            else:
                en[0, ik, ib], en[1, ik, ib] = vals[0], vals[1]
                if len(vals) >= 4:
                    occ[0, ik, ib], occ[1, ik, ib] = vals[2], vals[3]
                else:
                    has_occ = False
    return kf, en, (occ if has_occ else None), nelect


def pick_band(spec, en, occ, nelect, ik0, spin):
    nb = en.shape[2]
    if spec not in ("vbm", "cbm"):
        b = int(spec) - 1
        if not 0 <= b < nb:
            sys.exit("band %s out of range 1..%d" % (spec, nb))
        return b
    if occ is not None:
        nocc = int(np.sum(occ[spin, ik0] > 0.5))   # EIGENVAL occupations are normalised to 1
    else:
        nocc = int(round(nelect / 2))          # no occupations in old EIGENVAL: assume non-magnetic
    b = nocc - 1 if spec == "vbm" else nocc
    if not 0 <= b < nb:
        sys.exit("cannot locate %s (increase NBANDS?)" % spec)
    return b


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--eigenval", default="EIGENVAL")
    p.add_argument("--poscar", default="POSCAR")
    p.add_argument("--k0", nargs=3, type=float, default=[0.0, 0.0, 0.0],
                   help="extremum in fractional reciprocal coordinates (default Gamma)")
    p.add_argument("--band", action="append", default=None,
                   help="'vbm', 'cbm' or a 1-based band index; can be repeated (default: cbm and vbm)")
    p.add_argument("--spin", type=int, default=1, choices=[1, 2])
    p.add_argument("--npoints", type=int, default=0, help="use only the N closest points per direction (0 = all)")
    p.add_argument("--quartic", action="store_true", help="fit a*dk^2 + b*dk^4 instead of a*dk^2")
    args = p.parse_args()

    lat = read_lattice(args.poscar)
    brec = 2 * np.pi * np.linalg.inv(lat).T                 # rows = b1, b2, b3 (1/A, 2*pi included)
    kf, en, occ, nelect = read_eigenval(args.eigenval)
    spin = args.spin - 1
    if spin >= en.shape[0]:
        sys.exit("EIGENVAL is not spin polarized")
    kc = kf @ brec
    k0 = np.array(args.k0) @ brec
    d0 = np.linalg.norm(kc - k0, axis=1)
    ik0 = int(np.argmin(d0))
    if d0[ik0] > 1e-6:
        sys.exit("k0 is not in the k-point list (closest point %.2e 1/A away)" % d0[ik0])

    # group the other k-points by direction from k0
    groups = {}
    for ik in range(len(kc)):
        dk = kc[ik] - k0
        n = np.linalg.norm(dk)
        if n < 1e-8:
            continue
        key = tuple(np.round(dk / n, 3) + 0.0)
        groups.setdefault(key, []).append((n, ik))

    bands = args.band or ["cbm", "vbm"]
    for spec in bands:
        b = pick_band(spec, en, occ, nelect, ik0, spin)
        e0 = en[spin, ik0, b]
        degenerate = [j + 1 for j in range(en.shape[2])
                      if j != b and abs(en[spin, ik0, j] - e0) < 2e-3]
        print("\nband %d (%s), E(k0) = %.4f eV, spin %d" % (b + 1, spec, e0, args.spin))
        if degenerate:
            print("  NOTE: degenerate at k0 with band(s) %s -> bands may swap along a line;"
                  " masses refer to energy-ordered bands" % degenerate)
        print("  %-26s %6s %12s %12s %10s" % ("direction (cartesian)", "npts", "m*_fit/m0", "m*_2pt/m0", "max|dk|"))
        for key, pts in sorted(groups.items()):
            pts.sort()
            if args.npoints:
                pts = pts[:args.npoints]
            x = np.array([q[0] for q in pts])
            y = np.array([en[spin, q[1], b] for q in pts]) - e0
            if args.quartic and len(x) >= 2:
                A = np.column_stack([x ** 2, x ** 4])
                a = np.linalg.lstsq(A, y, rcond=None)[0][0]
            else:
                a = np.sum(x ** 2 * y) / np.sum(x ** 4)
            m_fit = HB2_2M / a if abs(a) > 1e-12 else float("inf")
            m_2pt = HB2_2M * x[0] ** 2 / y[0] if abs(y[0]) > 1e-12 else float("inf")
            print("  (%7.4f %7.4f %7.4f)    %6d %12.4f %12.4f %10.4f"
                  % (key[0], key[1], key[2], len(x), m_fit, m_2pt, x.max()))
    print("\nm* > 0: electron-like (minimum); m* < 0: hole-like (maximum). |dk| in 1/A (2*pi included).")


if __name__ == "__main__":
    main()
