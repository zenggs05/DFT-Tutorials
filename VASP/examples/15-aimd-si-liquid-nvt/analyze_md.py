#!/usr/bin/env python3
"""Minimal analysis of a VASP MD run (XDATCAR + OSZICAR).

    python analyze_md.py --potim 2.0                    # thermo + RDF + MSD
    python analyze_md.py --potim 2.0 --skip 500 --pair Si-Si --vdos
    python analyze_md.py --potim 0.5 --nblock 1 --species O

Outputs
  thermo.png / thermo.dat  temperature and energies from OSZICAR (T=, E=, F=)
  rdf.dat / rdf.png        g(r) and running coordination number n(r)
  msd.dat / msd.png        mean-squared displacement per species, D from 6Dt
  vdos.dat / vdos.png      vibrational DOS from the velocity autocorrelation
                           (velocities by finite differences, needs NBLOCK=1)

Notes
  * --potim is the MD time step in fs, --nblock the XDATCAR output interval
    (NBLOCK, or ML_OUTBLOCK in MLFF production runs).
  * Works for fixed-cell and variable-cell (ISIF=3) XDATCAR files.
  * Atoms are unwrapped with the minimum-image convention between successive
    frames, so atoms must move less than half a box length between frames.
  * RDF uses the fractional minimum image: use rmax below half of the
    smallest cell height (default) and near-orthogonal cells.
Requires numpy and matplotlib.
"""
import argparse
import re
import sys

import numpy as np


# ----------------------------------------------------------------- reading
def read_xdatcar(path):
    """Return species list, counts, lattices [F,3,3] (A), frac [F,N,3]."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    i = 0
    lat = None
    species = counts = None
    lats, frames = [], []
    while i < len(lines):
        line = lines[i].strip()
        if line.lower().startswith(("direct configuration", "cartesian configuration")):
            if line.lower().startswith("cartesian"):
                sys.exit("Cartesian XDATCAR not supported")
            n = sum(counts)
            frac = np.array([[float(x) for x in lines[i + 1 + j].split()[:3]] for j in range(n)])
            frames.append(frac)
            lats.append(lat)
            i += n + 1
            continue
        if not line:
            i += 1
            continue
        # header block: comment, scale, 3 lattice lines, species, counts
        scale = float(lines[i + 1].split()[0])
        lat = np.array([[float(x) for x in lines[i + k].split()[:3]] for k in (2, 3, 4)])
        lat = lat * scale if scale > 0 else lat * (abs(scale) / abs(np.linalg.det(lat))) ** (1 / 3)
        species = lines[i + 5].split()
        counts = [int(x) for x in lines[i + 6].split()]
        i += 7
    if not frames:
        sys.exit(f"no configurations found in {path}")
    return species, counts, np.array(lats), np.array(frames)


def read_oszicar(path):
    pat = re.compile(r"^\s*(\d+)\s+T=\s*([-\d.]+)\s+E=\s*(\S+)\s+F=\s*(\S+)")
    rows = []
    try:
        with open(path) as fh:
            for line in fh:
                m = pat.match(line)
                if m:
                    rows.append([float(m.group(k)) for k in (1, 2, 3, 4)])
    except FileNotFoundError:
        return None
    return np.array(rows) if rows else None


# ----------------------------------------------------------------- analysis
def unwrap(frac, lats):
    d = np.diff(frac, axis=0)
    d -= np.round(d)
    ufrac = np.concatenate([frac[:1], frac[:1] + np.cumsum(d, axis=0)])
    return np.einsum("fni,fij->fnj", ufrac, lats)  # Cartesian, A


def rdf(frac, lats, idx_a, idx_b, rmax, nbins):
    edges = np.linspace(0.0, rmax, nbins + 1)
    hist = np.zeros(nbins)
    norm = 0.0
    same = np.array_equal(idx_a, idx_b)
    for f in range(len(frac)):
        lat = lats[f]
        vol = abs(np.linalg.det(lat))
        fa, fb = frac[f, idx_a], frac[f, idx_b]
        d = fb[None, :, :] - fa[:, None, :]
        d -= np.round(d)
        r = np.linalg.norm(d @ lat, axis=-1)
        if same:
            r = r[~np.eye(len(idx_a), dtype=bool)]
        hist += np.histogram(r.ravel(), bins=edges)[0]
        norm += len(idx_a) * len(idx_b) / vol if not same else len(idx_a) * (len(idx_a) - 1) / vol
    rc = 0.5 * (edges[1:] + edges[:-1])
    shell = 4.0 / 3.0 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = hist / (shell * norm)
    rho_b = np.mean([len(idx_b) / abs(np.linalg.det(l)) for l in lats])
    coord = np.cumsum(g * rho_b * shell)
    return rc, g, coord


def msd(cart, idx, max_lag, origin_step):
    x = cart[:, idx, :]
    lags = np.arange(max_lag + 1)
    out = np.zeros(len(lags))
    for k, lag in enumerate(lags):
        origins = np.arange(0, len(x) - lag, origin_step)
        out[k] = np.mean(np.sum((x[origins + lag] - x[origins]) ** 2, axis=-1))
    return lags, out


def vdos(cart, dt_fs):
    v = (cart[2:] - cart[:-2]) / (2 * dt_fs)          # A/fs, central differences
    nt = len(v)
    vf = np.fft.rfft(v, n=2 * nt, axis=0)
    acf = np.fft.irfft(np.abs(vf) ** 2, axis=0)[:nt].sum(axis=(1, 2))
    acf /= np.arange(nt, 0, -1)
    acf /= acf[0]
    window = np.hanning(2 * nt)[nt:]
    spec = np.abs(np.fft.rfft(acf * window, n=4 * nt))
    freq = np.fft.rfftfreq(4 * nt, d=dt_fs) * 1000.0  # 1/fs -> THz
    return freq, spec / (spec.sum() * (freq[1] - freq[0]))


# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xdatcar", default="XDATCAR")
    ap.add_argument("--oszicar", default="OSZICAR")
    ap.add_argument("--potim", type=float, required=True, help="MD time step (fs)")
    ap.add_argument("--nblock", type=int, default=1, help="steps between XDATCAR frames")
    ap.add_argument("--skip", type=int, default=0, help="frames discarded as equilibration")
    ap.add_argument("--pair", help="RDF pair, e.g. Si-Si or O-H (default: all-all)")
    ap.add_argument("--rmax", type=float, help="RDF cutoff (A)")
    ap.add_argument("--nbins", type=int, default=200)
    ap.add_argument("--species", help="species for MSD/VDOS (default: every species)")
    ap.add_argument("--fit", nargs=2, type=float, default=(0.2, 0.6),
                    help="fraction of the MSD lag range used for the linear fit")
    ap.add_argument("--vdos", action="store_true")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    species, counts, lats, frac = read_xdatcar(args.xdatcar)
    labels = np.repeat(species, counts)
    print(f"{len(frac)} frames, {sum(counts)} atoms ({' '.join(f'{s}{n}' for s, n in zip(species, counts))})")
    lats, frac = lats[args.skip:], frac[args.skip:]
    if len(frac) < 10:
        sys.exit("fewer than 10 frames left after --skip")
    dt = args.potim * args.nblock
    vols = np.abs(np.linalg.det(lats))
    print(f"analysed frames: {len(frac)}  (dt between frames = {dt:g} fs, "
          f"total {dt * (len(frac) - 1) / 1000:.2f} ps); <V> = {vols.mean():.2f} A^3")

    # ---- thermodynamics from OSZICAR
    th = read_oszicar(args.oszicar)
    if th is not None:
        np.savetxt("thermo.dat", th, fmt="%14.6f", header="step T(K) E(eV) F(eV)")
        start = min(args.skip * args.nblock, len(th) - 1)
        print(f"<T> = {th[start:, 1].mean():.1f} +- {th[start:, 1].std():.1f} K (after skip)")
        fig, ax = plt.subplots(2, 1, figsize=(6, 5), sharex=True)
        ax[0].plot(th[:, 0] * args.potim / 1000, th[:, 1])
        ax[0].set_ylabel("T (K)")
        ax[1].plot(th[:, 0] * args.potim / 1000, th[:, 2], label="E (conserved)")
        ax[1].plot(th[:, 0] * args.potim / 1000, th[:, 3], label="F (potential)")
        ax[1].set_ylabel("energy (eV)")
        ax[1].set_xlabel("time (ps)")
        ax[1].legend()
        fig.tight_layout()
        fig.savefig("thermo.png", dpi=150)

    # ---- RDF
    if args.pair:
        a, b = args.pair.split("-")
        ia, ib = np.where(labels == a)[0], np.where(labels == b)[0]
        if len(ia) == 0 or len(ib) == 0:
            sys.exit(f"species in --pair not found; available: {species}")
    else:
        ia = ib = np.arange(len(labels))
    heights = [abs(np.linalg.det(l)) / np.linalg.norm(np.cross(l[(k + 1) % 3], l[(k + 2) % 3]))
               for l in lats for k in range(3)]
    rmax = args.rmax or 0.5 * min(heights)
    rc, g, coord = rdf(frac, lats, ia, ib, rmax, args.nbins)
    np.savetxt("rdf.dat", np.column_stack([rc, g, coord]), fmt="%12.5f", header="r(A) g(r) n(r)")
    # first peak = first local maximum with g > 1; first minimum = next local minimum
    gs = np.convolve(g, np.ones(3) / 3, mode="same")
    ipk = next((i for i in range(1, len(gs) - 1)
                if gs[i] > 1.0 and gs[i] >= gs[i - 1] and gs[i] >= gs[i + 1]), int(np.argmax(gs)))
    imin = next((i for i in range(ipk + 1, len(gs) - 1)
                 if gs[i] <= gs[i - 1] and gs[i] <= gs[i + 1]), len(gs) - 1)
    print(f"RDF {args.pair or 'all-all'}: first peak at {rc[ipk]:.2f} A, "
          f"first minimum ~{rc[imin]:.2f} A, coordination n(r_min) = {coord[imin]:.2f}")
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(rc, g, label="g(r)")
    ax.plot(rc, coord, "--", label="n(r)")
    ax.set_xlabel("r (A)")
    ax.legend()
    fig.tight_layout()
    fig.savefig("rdf.png", dpi=150)

    # ---- MSD and diffusion coefficients
    cart = unwrap(frac, lats)
    sel = [args.species] if args.species else species
    max_lag = len(cart) // 2
    origin_step = max(1, len(cart) // 200)
    cols, header = [], []
    fig, ax = plt.subplots(figsize=(5, 4))
    for s in sel:
        idx = np.where(labels == s)[0]
        lags, m = msd(cart, idx, max_lag, origin_step)
        t = lags * dt
        lo, hi = int(args.fit[0] * max_lag), int(args.fit[1] * max_lag)
        slope = np.polyfit(t[lo:hi + 1], m[lo:hi + 1], 1)[0]   # A^2/fs
        d_cm2_s = slope / 6.0 * 0.1                              # 1 A^2/fs = 0.1 cm^2/s
        print(f"MSD {s}: D = {d_cm2_s:.3e} cm^2/s  (fit {t[lo] / 1000:.2f}-{t[hi] / 1000:.2f} ps)")
        ax.plot(t / 1000, m, label=f"{s}: D={d_cm2_s:.2e} cm$^2$/s")
        cols.append(m)
        header.append(f"MSD_{s}(A^2)")
    np.savetxt("msd.dat", np.column_stack([t] + cols), fmt="%14.6f",
               header="t(fs) " + " ".join(header))
    ax.set_xlabel("time (ps)")
    ax.set_ylabel(r"MSD ($\AA^2$)")
    ax.legend()
    fig.tight_layout()
    fig.savefig("msd.png", dpi=150)

    # ---- vibrational DOS
    if args.vdos:
        if args.nblock != 1:
            print("warning: VDOS from positions needs NBLOCK = 1 (every MD step stored)")
        idx = np.where(labels == args.species)[0] if args.species else np.arange(len(labels))
        freq, spec = vdos(cart[:, idx, :], dt)
        keep = freq < 1000.0 / (2 * dt)
        np.savetxt("vdos.dat", np.column_stack([freq[keep], spec[keep]]), fmt="%14.6f",
                   header="frequency(THz) VDOS(normalised)")
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.plot(freq[keep], spec[keep])
        ax.set_xlabel("frequency (THz)")
        ax.set_ylabel("VDOS (arb. units)")
        fig.tight_layout()
        fig.savefig("vdos.png", dpi=150)
        print("wrote vdos.dat / vdos.png (1 THz = 33.356 cm^-1)")
    print("wrote thermo/rdf/msd .dat and .png files")


if __name__ == "__main__":
    main()
