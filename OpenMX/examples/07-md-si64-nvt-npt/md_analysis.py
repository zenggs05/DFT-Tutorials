#!/usr/bin/env python3
"""Analyse an OpenMX MD run: temperature, conserved quantities, g(r), MSD, NPT volume.

Input files (System.Name = prefix):
    <prefix>.ene      one line per MD step (manual Sec. 20): col 2 time (fs),
                      14 Ukc, 15 Utot, 16 Utot+Ukc, 22 Nose-Hoover Hamiltonian (Hartree),
                      18 target T, 19 calculated T (K)
    <prefix>.md       trajectory in the OpenMX 4.0 md format (manual Sec. 18.6)
    <prefix>.npt.ene  (optional, NPT runs) col 3 volume (Bohr^3), 4 T (K), 5 pressure (GPa)

Usage examples:
    python3 md_analysis.py Si64_NVT_VS --skip 400
    python3 md_analysis.py Si64_NVT_NH --skip 250 --rmax 5.0 --plot
    python3 md_analysis.py Si8_NPT_VS_PR --skip 200
    python3 md_analysis.py Si64_NVT_NH --skip 250 --dpdata deepmd_si64 --stride 5

Requires numpy; matplotlib for --plot; dpdata for --dpdata.
"""
import argparse
import os

import numpy as np

HA2MEV = 27211.386
BOHR = 0.529177210903


def read_md(fname):
    """Return lists: symbols (first frame), positions [nfr,N,3] (Ang), cells [nfr,3,3], times."""
    pos, cells, times, sym = [], [], [], None
    with open(fname) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.strip():
                continue
            nat = int(line.split()[0])
            head = f.readline()
            cell = None
            if "Cell_Vectors=" in head:
                cell = np.array(head.split("Cell_Vectors=")[1].split()[:9], float).reshape(3, 3)
            t = np.nan
            tok = head.split()
            if "time=" in tok:
                t = float(tok[tok.index("time=") + 1])
            s, p = [], []
            for _ in range(nat):
                q = f.readline().split()
                s.append(q[0])
                p.append([float(x) for x in q[1:4]])
            sym = sym or s
            pos.append(p)
            cells.append(cell)
            times.append(t)
    return sym, np.array(pos), cells, np.array(times)


def perpendicular_widths(cell):
    vol = abs(np.linalg.det(cell))
    return np.array([vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3]))
                     for i in range(3)])


def rdf(pos, cells, sym, pair, rmax, nbins):
    a, b = pair
    ia = np.array([i for i, s in enumerate(sym) if s == a])
    ib = np.array([i for i, s in enumerate(sym) if s == b])
    edges = np.linspace(0.0, rmax, nbins + 1)
    hist = np.zeros(nbins)
    vol_sum = 0.0
    for x, cell in zip(pos, cells):
        inv = np.linalg.inv(cell)
        fa, fb = x[ia] @ inv, x[ib] @ inv
        d = fa[:, None, :] - fb[None, :, :]
        d -= np.round(d)
        r = np.linalg.norm(d @ cell, axis=2)
        if a == b:
            r = r[np.triu_indices(len(ia), k=1)]
        else:
            r = r.ravel()
        h, _ = np.histogram(r, bins=edges)
        hist += h
        vol_sum += abs(np.linalg.det(cell))
    nfr = len(pos)
    vol = vol_sum / nfr
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    npairs = len(ia) * (len(ia) - 1) / 2.0 if a == b else len(ia) * len(ib)
    g = hist / nfr / (npairs / vol) / shell
    rc = 0.5 * (edges[1:] + edges[:-1])
    # running coordination number: average number of b neighbours around an a atom
    cn = np.cumsum(hist / nfr) * (2.0 if a == b else 1.0) / len(ia)
    return rc, g, cn


def unwrap(pos, cells):
    """Unwrap coordinates using the minimum-image displacement between successive frames."""
    out = np.empty_like(pos)
    out[0] = pos[0]
    for k in range(1, len(pos)):
        inv = np.linalg.inv(cells[k])
        df = (pos[k] - pos[k - 1]) @ inv
        df -= np.round(df)
        out[k] = out[k - 1] + df @ cells[k]
    return out


def msd(xu, maxlag, stride=1):
    lags = np.arange(1, maxlag + 1)
    m = np.empty(len(lags))
    for i, lag in enumerate(lags):
        d = xu[lag::stride] - xu[:-lag:stride][: len(xu[lag::stride])]
        m[i] = np.mean(np.sum(d ** 2, axis=2))
    return lags, m


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="System.Name of the MD run")
    ap.add_argument("--skip", type=int, default=0, help="steps/frames to discard (equilibration)")
    ap.add_argument("--pair", default="", help="element pair for g(r), e.g. Si-Si (default: first element)")
    ap.add_argument("--rmax", type=float, default=0.0, help="g(r) cutoff in Ang (default: half cell width)")
    ap.add_argument("--nbins", type=int, default=200)
    ap.add_argument("--plot", action="store_true", help="save <prefix>_analysis.png")
    ap.add_argument("--dpdata", default="", help="export frames to deepmd/npy in this directory")
    ap.add_argument("--stride", type=int, default=5, help="frame stride for --dpdata")
    args = ap.parse_args()
    p = args.prefix

    # ---------- .ene ----------
    if os.path.exists(p + ".ene"):
        ene = np.loadtxt(p + ".ene", ndmin=2)[args.skip:]
        t, temp = ene[:, 1], ene[:, 18]
        print("ene: %d steps used (skip %d)" % (len(t), args.skip))
        print("<T> = %.1f +- %.1f K (target %.1f K at last step)" % (temp.mean(), temp.std(), ene[-1, 17]))
        if len(t) > 2:
            for col, name in ((15, "Utot+Ukc"), (21, "H_NoseHoover")):
                slope = np.polyfit(t / 1000.0, ene[:, col], 1)[0] * HA2MEV
                print("drift of %-12s = %10.3f meV/ps (meaningful for NVE / NVT_NH respectively)"
                      % (name, slope))

    # ---------- .npt.ene ----------
    if os.path.exists(p + ".npt.ene"):
        npt = np.loadtxt(p + ".npt.ene", ndmin=2)[args.skip:]
        v = npt[:, 2] * BOHR ** 3
        print("npt: <V> = %.3f +- %.3f A^3, <T> = %.1f K, <P_int> = %.2f +- %.2f GPa"
              % (v.mean(), v.std(), npt[:, 3].mean(), npt[:, 4].mean(), npt[:, 4].std()))
        print("     cubic-cell equivalent a = <V>^(1/3) = %.4f A" % (v.mean() ** (1.0 / 3.0)))

    # ---------- .md ----------
    if not os.path.exists(p + ".md"):
        return
    sym, pos, cells, times = read_md(p + ".md")
    pos, cells, times = pos[args.skip:], cells[args.skip:], times[args.skip:]
    if any(c is None for c in cells):
        print("md: no Cell_Vectors= in the md file (cluster run?); g(r)/MSD skipped")
        return
    dt = np.median(np.diff(times)) if len(times) > 1 else np.nan
    print("md: %d frames used, dt = %.3f fs, %d atoms" % (len(pos), dt, len(sym)))

    first = sym[0]
    pair = tuple(args.pair.split("-")) if args.pair else (first, first)
    wmin = min(perpendicular_widths(c).min() for c in cells)
    rmax = args.rmax if args.rmax > 0 else 0.5 * wmin
    if rmax > 0.5 * wmin:
        print("warning: rmax %.2f > half cell width %.2f, minimum image not valid" % (rmax, 0.5 * wmin))
    rc, g, cn = rdf(pos, cells, sym, pair, rmax, args.nbins)
    # first peak / first minimum on a lightly smoothed g(r)
    gs = np.convolve(g, np.ones(5) / 5.0, mode="same")
    ipk = int(np.argmax(gs > 0.3 * gs.max()))
    while ipk + 1 < len(gs) and gs[ipk + 1] >= gs[ipk]:
        ipk += 1
    imin = ipk
    while imin + 1 < len(gs) and gs[imin + 1] <= gs[imin]:
        imin += 1
    print("g(r) %s-%s: first peak %.3f A, first minimum %.3f A, coordination %.2f"
          % (pair[0], pair[1], rc[ipk], rc[imin], cn[imin]))
    np.savetxt(p + "_rdf.dat", np.column_stack([rc, g, cn]), header="r(A) g(r) n(r)")

    xu = unwrap(pos, cells)
    maxlag = max(1, len(xu) // 2)
    lags, m = msd(xu, maxlag)
    tl = lags * dt
    np.savetxt(p + "_msd.dat", np.column_stack([tl, m]), header="t(fs) MSD(A^2)")
    half = len(tl) // 2
    if half > 2:
        slope = np.polyfit(tl[half:], m[half:], 1)[0]          # A^2 / fs
        print("MSD at %.0f fs = %.4f A^2; D from linear fit (second half) = %.3e cm^2/s"
              % (tl[-1], m[-1], slope / 6.0 * 0.1))
        print("(a plateau with D ~ 0 means a solid; a linear MSD means diffusion)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
        if os.path.exists(p + ".ene"):
            ax[0].plot(t, temp, lw=0.8, label="T calc")
            ax[0].plot(t, ene[:, 17], lw=1.2, label="T target")
            ax[0].set_xlabel("time (fs)")
            ax[0].set_ylabel("T (K)")
            ax[0].legend()
        ax[1].plot(rc, g)
        ax[1].set_xlabel("r (A)")
        ax[1].set_ylabel("g(r) %s-%s" % pair)
        ax[2].plot(tl, m)
        ax[2].set_xlabel("t (fs)")
        ax[2].set_ylabel("MSD (A$^2$)")
        fig.tight_layout()
        fig.savefig(p + "_analysis.png", dpi=150)
        print("figure written to %s_analysis.png" % p)

    if args.dpdata:
        import dpdata
        ls = dpdata.LabeledSystem(p, fmt="openmx/md")
        sub = ls[args.skip::args.stride]
        sub.to("deepmd/npy", args.dpdata)
        print("dpdata: %d frames written to %s/" % (len(sub), args.dpdata))


if __name__ == "__main__":
    main()
