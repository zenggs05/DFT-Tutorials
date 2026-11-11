#!/usr/bin/env python3
"""Fit the energy curve written by OpenMX with MD.Type EvsLC (file System.Name.EvsLC).

Each data line of System.Name.EvsLC contains
    ax ay az  bx by bz  cx cy cz  Etot
with the lattice vectors in Angstrom and the total energy in Hartree.
Lines starting with '#' or containing fewer than 10 numbers are ignored, so several
files (e.g. a restarted scan) can be given at once.

Modes
  volume (default)  third-order Birch-Murnaghan fit of E(V): V0, E0, B0, B0'
                    (linear least squares of E against V^(-2/3), cubic polynomial)
  c                 polynomial fit of E against |c| for c-axis scans
                    (e.g. MD.EvsLC.flag 0 0 1 for layered materials): c0, E_min

Examples
  python fit_evslc.py cdia_evslc.EvsLC --lattice fcc --natoms 2
  python fit_evslc.py graphite_d3.EvsLC --mode c --natoms 4 --plot graphite_c.png

Only numpy is required; matplotlib is needed for --plot.
"""
import argparse
import sys

import numpy as np

HA2EV = 27.211386245988
EVA3_TO_GPA = 160.21766208


def read_evslc(paths):
    rows = []
    for path in paths:
        with open(path) as fh:
            for line in fh:
                tok = line.split("#", 1)[0].split()
                if len(tok) < 10:
                    continue
                try:
                    rows.append([float(t) for t in tok[:10]])
                except ValueError:
                    continue
    if len(rows) < 4:
        sys.exit("ERROR: need at least 4 data points, found %d" % len(rows))
    data = np.array(rows)
    cells = data[:, :9].reshape(-1, 3, 3)
    energy = data[:, 9] * HA2EV
    return cells, energy


def fit_bm3(volume, energy):
    """E(V) = sum_k c_k x^k with x = V^(-2/3) (cubic) is exactly the 3rd-order BM form."""
    x = volume ** (-2.0 / 3.0)
    poly = np.poly1d(np.polyfit(x, energy, 3))
    d1, d2, d3 = poly.deriv(1), poly.deriv(2), poly.deriv(3)
    cands = [r.real for r in d1.r if abs(r.imag) < 1e-10 and r.real > 0 and d2(r.real) > 0]
    if not cands:
        raise RuntimeError("no energy minimum found in the fitted curve")
    x0 = min(cands, key=lambda r: abs(r - x.mean()))
    v0 = x0 ** -1.5
    e0 = poly(x0)
    # chain rule for x(V) = V^(-2/3)
    xv1 = -2.0 / 3.0 * v0 ** (-5.0 / 3.0)
    xv2 = 10.0 / 9.0 * v0 ** (-8.0 / 3.0)
    xv3 = -80.0 / 27.0 * v0 ** (-11.0 / 3.0)
    f1, f2, f3 = d1(x0), d2(x0), d3(x0)
    ev2 = f2 * xv1 ** 2 + f1 * xv2
    ev3 = f3 * xv1 ** 3 + 3.0 * f2 * xv1 * xv2 + f1 * xv3
    b0 = v0 * ev2                       # eV/A^3
    b0p = -1.0 - v0 * ev3 / ev2
    resid = energy - poly(x)
    return v0, e0, b0, b0p, resid, (lambda v: poly(v ** (-2.0 / 3.0)))


def lattice_constant(v0, lattice):
    factor = {"fcc": 4.0, "bcc": 2.0, "sc": 1.0}
    return (factor[lattice] * v0) ** (1.0 / 3.0) if lattice in factor else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="System.Name.EvsLC file(s)")
    ap.add_argument("--mode", choices=["volume", "c"], default="volume")
    ap.add_argument("--lattice", choices=["fcc", "bcc", "sc", "none"], default="none",
                    help="primitive-cell type used to convert V0 into a cubic lattice constant")
    ap.add_argument("--natoms", type=int, default=1, help="atoms per cell (for per-atom energies)")
    ap.add_argument("--plot", metavar="PNG", help="save a plot of data and fit")
    args = ap.parse_args()

    cells, energy = read_evslc(args.files)
    order = None

    if args.mode == "volume":
        volume = np.abs(np.linalg.det(cells))
        order = np.argsort(volume)
        volume, energy = volume[order], energy[order]
        v0, e0, b0, b0p, resid, model = fit_bm3(volume, energy)
        print("# 3rd-order Birch-Murnaghan fit of %d points" % len(volume))
        print("V0  = %.4f A^3/cell   (%.4f A^3/atom)" % (v0, v0 / args.natoms))
        a0 = lattice_constant(v0, args.lattice)
        if a0 is not None:
            print("a0  = %.4f A          (%s, from V0)" % (a0, args.lattice))
        print("E0  = %.6f eV/cell   (%.6f eV/atom)" % (e0, e0 / args.natoms))
        print("B0  = %.2f GPa" % (b0 * EVA3_TO_GPA))
        print("B0' = %.2f" % b0p)
        print("RMS residual = %.3f meV/cell" % (1000.0 * np.sqrt(np.mean(resid ** 2))))
        print("#     V(A^3)    E-E0(meV/cell)  fit-residual(meV)")
        for v, e, r in zip(volume, energy, resid):
            print("  %10.4f  %14.3f  %14.3f" % (v, 1000.0 * (e - e0), 1000.0 * r))
        if not volume.min() < v0 < volume.max():
            print("WARNING: V0 lies outside the scanned range; extend the scan", file=sys.stderr)
        elif v0 < volume[1] or v0 > volume[-2]:
            print("WARNING: V0 is close to the edge of the scanned range", file=sys.stderr)
        xdat, xlab = volume, "V (A^3/cell)"
        xfit = np.linspace(volume.min(), volume.max(), 200)
        yfit = model(xfit)
        ymin = e0
    else:
        clen = np.linalg.norm(cells[:, 2, :], axis=1)
        order = np.argsort(clen)
        clen, energy = clen[order], energy[order]
        deg = 4 if len(clen) >= 7 else 3
        poly = np.poly1d(np.polyfit(clen, energy, deg))
        grid = np.linspace(clen.min(), clen.max(), 4001)
        i = int(np.argmin(poly(grid)))
        c0, emin = grid[i], poly(grid[i])
        print("# polynomial (degree %d) fit of E(|c|), %d points" % (deg, len(clen)))
        print("c0    = %.4f A" % c0)
        print("E_min = %.6f eV/cell" % emin)
        print("E(c_max) - E_min = %.2f meV/cell (%.2f meV/atom); a lower bound of the binding energy"
              % (1000 * (energy[-1] - emin), 1000 * (energy[-1] - emin) / args.natoms))
        print("#     |c|(A)   E-E_min(meV/cell)")
        for c, e in zip(clen, energy):
            print("  %10.4f  %14.3f" % (c, 1000.0 * (e - emin)))
        j = int(np.argmin(energy))
        if i in (0, len(grid) - 1) or j in (0, len(energy) - 1):
            print("WARNING: the lowest energy is at the edge of the scan; no bound minimum inside "
                  "the range (extend the scan, or the interaction is not binding)", file=sys.stderr)
        xdat, xlab = clen, "|c| (A)"
        xfit = grid
        yfit = poly(grid)
        ymin = emin

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.5, 3.4))
        ax.plot(xfit, 1000.0 * (yfit - ymin), "-", lw=1.5, label="fit")
        ax.plot(xdat, 1000.0 * (energy - ymin), "o", ms=4, label="OpenMX")
        ax.set_xlabel(xlab)
        ax.set_ylabel("E - E_min (meV/cell)")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200)
        print("plot written to %s" % args.plot)


if __name__ == "__main__":
    main()
