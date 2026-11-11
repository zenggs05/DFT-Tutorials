#!/usr/bin/env python3
"""Check an OpenMX NVE run: energy conservation and vibrational frequencies.

Reads <prefix>.ene (one line per MD step) and <prefix>.md (OpenMX 4.0 md format).
Column numbers of .ene follow the OpenMX manual (Sec. 20.1-20.3):
    2: time (fs), 14: Ukc, 15: Utot, 16: Utot+Ukc (Hartree), 19: temperature (K)

Usage:
    python3 nve_check.py H2O_NVE
    python3 nve_check.py H2O_NVE --bonds 1-2,1-3 --angles 2-1-3 --skip 0 --plot

Requires numpy; matplotlib only for --plot.
"""
import argparse

import numpy as np

HA2MEV = 27211.386
FS2CM = 1.0e15 / 2.99792458e10      # (1/fs) -> cm^-1


def read_md(fname):
    """Return list of (symbols, positions[Ang], cell or None, header line)."""
    frames = []
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
                vals = head.split("Cell_Vectors=")[1].split()[:9]
                cell = np.array(vals, float).reshape(3, 3)
            sym, pos = [], []
            for _ in range(nat):
                p = f.readline().split()
                sym.append(p[0])
                pos.append([float(x) for x in p[1:4]])
            frames.append((sym, np.array(pos), cell, head))
    return frames


def spectrum(signal, dt_fs):
    """Power spectrum of a real signal; returns wavenumbers (cm^-1) and power."""
    x = np.asarray(signal) - np.mean(signal)
    x = x * np.hanning(len(x))
    pw = np.abs(np.fft.rfft(x)) ** 2
    freq = np.fft.rfftfreq(len(x), d=dt_fs) * FS2CM
    return freq, pw


def top_peaks(freq, pw, fmin=300.0, fmax=5000.0, n=3, rel=0.05):
    """Local maxima between fmin and fmax stronger than rel * strongest peak."""
    idx = [i for i in range(1, len(pw) - 1)
           if pw[i] > pw[i - 1] and pw[i] >= pw[i + 1] and fmin < freq[i] < fmax]
    if not idx:
        return []
    pmax = max(pw[i] for i in idx)
    idx = [i for i in idx if pw[i] >= rel * pmax]
    idx.sort(key=lambda i: pw[i], reverse=True)
    return sorted(freq[i] for i in idx[:n])


def pairs(text, size):
    out = []
    for item in text.split(","):
        item = item.strip()
        if item:
            ids = [int(x) - 1 for x in item.split("-")]
            if len(ids) != size:
                raise SystemExit("bad atom list: %s" % item)
            out.append(ids)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="System.Name of the run")
    ap.add_argument("--skip", type=int, default=0, help="number of initial steps to discard")
    ap.add_argument("--bonds", default="1-2,1-3", help="atom pairs (1-based), e.g. 1-2,1-3")
    ap.add_argument("--angles", default="2-1-3", help="atom triples (1-based), e.g. 2-1-3")
    ap.add_argument("--plot", action="store_true", help="save nve_check.png")
    args = ap.parse_args()

    ene = np.loadtxt(args.prefix + ".ene", ndmin=2)[args.skip:]
    t = ene[:, 1]
    dt = np.median(np.diff(t)) if len(t) > 1 else float("nan")
    ukc, econs, temp = ene[:, 13], ene[:, 15], ene[:, 18]
    print("steps = %d, dt = %.3f fs" % (len(t), dt))
    print("<T> = %.1f K" % temp.mean())
    if len(t) > 2:
        slope = np.polyfit(t / 1000.0, econs, 1)[0]          # Hartree / ps
        print("E_cons drift = %.3f meV/ps, peak-to-peak = %.3f meV"
              % (slope * HA2MEV, np.ptp(econs) * HA2MEV))
        print("Ukc peak-to-peak = %.1f meV" % (np.ptp(ukc) * HA2MEV))

    frames = read_md(args.prefix + ".md")[args.skip:]
    if not frames:
        return
    pos = np.array([f[1] for f in frames])
    nfr = min(len(pos), len(t))
    series = []
    for i, j in pairs(args.bonds, 2):
        d = np.linalg.norm(pos[:nfr, i] - pos[:nfr, j], axis=1)
        series.append(("bond %d-%d" % (i + 1, j + 1), d))
    for i, j, k in pairs(args.angles, 3):
        u = pos[:nfr, i] - pos[:nfr, j]
        w = pos[:nfr, k] - pos[:nfr, j]
        c = np.sum(u * w, axis=1) / np.linalg.norm(u, axis=1) / np.linalg.norm(w, axis=1)
        series.append(("angle %d-%d-%d" % (i + 1, j + 1, k + 1),
                       np.degrees(np.arccos(np.clip(c, -1, 1)))))
    for name, s in series:
        f, p = spectrum(s, dt)
        peaks = ", ".join("%.0f" % x for x in top_peaks(f, p))
        print("%-14s mean = %8.4f  spectrum peaks near %s cm^-1" % (name, s.mean(), peaks))

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].plot(t, (econs - econs[0]) * HA2MEV, label="Utot+Ukc")
        ax[0].plot(t, (ukc - ukc[0]) * HA2MEV, label="Ukc", alpha=0.6)
        ax[0].set_xlabel("time (fs)")
        ax[0].set_ylabel("energy change (meV)")
        ax[0].legend()
        for name, s in series:
            f, p = spectrum(s, dt)
            ax[1].plot(f, p / p.max(), label=name)
        ax[1].set_xlim(0, 4500)
        ax[1].set_xlabel("wavenumber (cm$^{-1}$)")
        ax[1].set_ylabel("normalized power")
        ax[1].legend()
        fig.tight_layout()
        fig.savefig("nve_check.png", dpi=150)
        print("figure written to nve_check.png")


if __name__ == "__main__":
    main()
