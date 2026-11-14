#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_band_plot.py -- plot OpenMX band structures and report VBM / CBM / gap.
                    读取 OpenMX 的 System.Name.Band（Hartree，未平移）或 bandgnu13 输出的
                    BANDDAT1/2（eV，已以化学势为零点），画自旋分辨能带，给出 VBM/CBM 与带隙。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Supported inputs (format checked against OpenMX 4.0 Band_DFT_kpath.c / bandgnu13.c)
------------------------------------------------------------------------------------
1. System.Name.Band (written when Band.dispersion = on), plain text:
       n_eig  spin_flag  mu                  mu = chemical potential (Hartree)
       b11 b12 ... b33                       reciprocal vectors of the SCF cell (Bohr^-1)
       Nkpath
       n k1 k2 k3 k1' k2' k3' lab1 lab2      one line per path segment
       n k1 k2 k3                            then for every segment, k point and spin:
       e1 e2 ... en                          eigenvalues (Hartree)
   spin_flag = 1 for collinear spin (up block, then down block per k point);
   spin_flag = 0 for non-spin AND for non-collinear / SOC runs (then n = 2 x basis size).
   Energies are converted with 27.2113845 eV/Hartree (same factor as bandgnu13)
   and shifted so that mu = 0.  The x axis is the cumulative path length in
   Bohr^-1 (2*pi included), with zero-length jumps between segments, exactly as
   bandgnu13 computes it.
2. System.Name.BANDDAT1 (+ BANDDAT2 for collinear spin) written by bandgnu13:
   blocks separated by blank lines; the first Nkpath-1 blocks are vertical
   separator lines (constant x) and are skipped.  Tick labels are taken from
   System.Name.GNUBAND ("set xtics (...)") when it exists, else from --dat.

Labels
------
* .Band files already contain the labels of <Band.kpath.
* --dat FILE reads Band.Nkpath / <Band.kpath from an OpenMX input file and
  uses its labels (useful for BANDDAT files without GNUBAND, or to relabel).
* --labels "L G X W K G" overrides everything: give Nkpath+1 labels for a
  continuous path or 2*Nkpath labels (start/end of each segment).
  g / G / Gamma / GM are drawn as the Greek Gamma.

Gap analysis
------------
States with E - mu <= 0 are counted as occupied (or use --nocc N: number of
occupied bands per spin channel; for non-collinear runs N = number of
electrons).  If a band has states on both sides of E_F the system is reported
as metallic along the path.  VBM/CBM are searched ONLY on the k points of the
path - a true gap needs a dense k mesh / DOS as well.

Usage examples
--------------
    python omx_band_plot.py si.Band
    python omx_band_plot.py si.Band --emin -13 --emax 6 -o si_band.png --csv si_band.csv
    python omx_band_plot.py fe.Band --spin both --emin -8 --emax 4
    python omx_band_plot.py gaas_soc.Band --compare gaas.Band --compare-label "no SOC"
    python omx_band_plot.py cdia.BANDDAT1                  # labels from cdia.GNUBAND
    python omx_band_plot.py cdia.BANDDAT1 --dat Cdia_band.dat
    python omx_band_plot.py si.Band --labels "L G X W K G" --ref vbm

Dependencies: Python >= 3.9, numpy, matplotlib.
"""

import argparse
import os
import re
import sys

HA2EV = 27.2113845  # identical to bandgnu13.c, so numbers agree with BANDDAT


def die(msg):
    sys.exit("omx_band_plot.py: error: " + msg)


# ----------------------------------------------------------------------------
# readers
# ----------------------------------------------------------------------------
def read_band_file(fname):
    """Parse System.Name.Band. Returns a dict with x, eig[spin,k,band] (eV, mu=0)."""
    import numpy as np

    try:
        with open(fname) as f:
            tok = f.read().split()
    except OSError as exc:
        die("cannot read %s (%s)" % (fname, exc))
    pos = [0]

    def nxt():
        if pos[0] >= len(tok):
            die("%s ends unexpectedly (truncated file or not an OpenMX .Band file)" % fname)
        pos[0] += 1
        return tok[pos[0] - 1]

    try:
        neig, spin_flag, mu = int(nxt()), int(nxt()), float(nxt())
        rtv = np.array([float(nxt()) for _ in range(9)]).reshape(3, 3)
        nseg = int(nxt())
        segs = []
        for _ in range(nseg):
            nk = int(nxt())
            k0 = [float(nxt()) for _ in range(3)]
            k1 = [float(nxt()) for _ in range(3)]
            segs.append([nk, k0, k1, nxt(), nxt()])
        nspin = 2 if spin_flag == 1 else 1
        kfrac, eig = [], [[] for _ in range(nspin)]
        for nk, _, _, _, _ in segs:
            for _ in range(nk):
                for s in range(nspin):
                    n = int(nxt())
                    k = [float(nxt()) for _ in range(3)]
                    e = [float(nxt()) for _ in range(n)]
                    if s == 0:
                        kfrac.append(k)
                    eig[s].append(e)
    except ValueError as exc:
        die("%s is malformed near token %d (%s)" % (fname, pos[0], exc))
    if pos[0] != len(tok):
        print("warning: %d trailing tokens in %s ignored" % (len(tok) - pos[0], fname))
    # bandgnu13 uses the smallest number of eigenvalues over all k points
    nmin = min(len(e) for es in eig for e in es)
    if nmin < neig:
        print("note: some k points have fewer than %d eigenvalues; using %d bands" % (neig, nmin))
    eig = np.array([[e[:nmin] for e in es] for es in eig])
    eig = (eig - mu) * HA2EV
    kfrac = np.array(kfrac)
    kcart = kfrac @ rtv
    x = np.zeros(len(kfrac))
    bounds = []
    i = 0
    for nk, _, _, _, _ in segs:
        for j in range(nk):
            if j == 0:
                x[i] = x[i - 1] if i > 0 else 0.0
            else:
                x[i] = x[i - 1] + np.linalg.norm(kcart[i] - kcart[i - 1])
            i += 1
        bounds.append((i - nk, i - 1))
    labels = [(s[3], s[4]) for s in segs]
    return dict(kind="Band", mu_ev=mu * HA2EV, x=x, eig=eig, kfrac=kfrac,
                bounds=bounds, seg_labels=labels, noncollinear_possible=(spin_flag == 0))


def _read_blocks(fname):
    import numpy as np

    blocks, cur = [], []
    try:
        with open(fname) as f:
            for line in f:
                p = line.split()
                if len(p) >= 2:
                    try:
                        cur.append((float(p[0]), float(p[1])))
                    except ValueError:
                        die("non-numeric line in %s: %s" % (fname, line.strip()))
                elif cur:
                    blocks.append(np.array(cur))
                    cur = []
    except OSError as exc:
        die("cannot read %s (%s)" % (fname, exc))
    if cur:
        blocks.append(np.array(cur))
    if not blocks:
        die("%s contains no data blocks" % fname)
    return blocks


def read_banddat(fname):
    """Parse one bandgnu13 BANDDAT file -> (x, eig[k,band], bounds)."""
    import numpy as np

    blocks = _read_blocks(fname)
    nvert = 0
    while nvert < len(blocks) and np.ptp(blocks[nvert][:, 0]) < 1e-10 and len(blocks[nvert]) > 2:
        nvert += 1
    nseg = nvert + 1
    data = blocks[nvert:]
    if len(data) % nseg != 0:
        die("%s: %d band blocks is not a multiple of %d segments" % (fname, len(data), nseg))
    nband = len(data) // nseg
    seglen = [len(b) for b in data[:nseg]]
    x = np.concatenate([b[:, 0] for b in data[:nseg]])
    eig = np.zeros((len(x), nband))
    for ib in range(nband):
        part = data[ib * nseg:(ib + 1) * nseg]
        if [len(b) for b in part] != seglen:
            die("%s: band %d has a different number of k points" % (fname, ib + 1))
        eig[:, ib] = np.concatenate([b[:, 1] for b in part])
    bounds, start = [], 0
    for n in seglen:
        bounds.append((start, start + n - 1))
        start += n
    return x, eig, bounds


def read_gnuband_tics(fname):
    """Return list of (label, x) from 'set xtics ("L" 0.0, "G" 0.5, ...)'."""
    if not os.path.exists(fname):
        return []
    with open(fname) as f:
        for line in f:
            if line.strip().lower().startswith("set xtics"):
                return [(lab, float(xv)) for lab, xv in
                        re.findall(r'"([^"]*)"\s+([-+0-9.eE]+)', line)]
    return []


def read_banddat_set(fname):
    import numpy as np

    base = re.sub(r"\.BANDDAT[12]$", "", fname)
    files = [fname]
    if fname.endswith(".BANDDAT1") and os.path.exists(base + ".BANDDAT2"):
        files.append(base + ".BANDDAT2")
    xs, eigs, bounds = [], [], None
    for fn in files:
        x, e, b = read_banddat(fn)
        xs.append(x)
        eigs.append(e)
        bounds = bounds or b
    if len(eigs) == 2:
        if eigs[0].shape != eigs[1].shape:
            n = min(eigs[0].shape[1], eigs[1].shape[1])
            eigs = [eigs[0][:, :n], eigs[1][:, :n]]
            if eigs[0].shape != eigs[1].shape:
                die("BANDDAT1 and BANDDAT2 have different k grids")
    eig = np.array(eigs)
    tics = read_gnuband_tics(base + ".GNUBAND")
    print("read %s (%d spin file(s)); GNUBAND labels: %s"
          % (fname, len(files), "yes" if tics else "no"))
    return dict(kind="BANDDAT", mu_ev=None, x=xs[0], eig=eig, kfrac=None,
                bounds=bounds, seg_labels=None, gnuband_tics=tics,
                noncollinear_possible=False)


def read_any(fname):
    if not os.path.exists(fname):
        die("file not found: %s" % fname)
    if re.search(r"\.BANDDAT[12]$", fname):
        return read_banddat_set(fname)
    with open(fname) as f:
        first = f.readline().split()
    if len(first) == 3:
        return read_band_file(fname)
    if len(first) == 2:
        return read_banddat_set(fname)
    die("%s is neither a .Band file nor a BANDDAT file" % fname)


def read_dat_kpath(datfile):
    """Labels of <Band.kpath in an OpenMX input file -> list of (lab1, lab2)."""
    try:
        lines = open(datfile).read().splitlines()
    except OSError as exc:
        die("cannot read %s (%s)" % (datfile, exc))
    out, inside = [], False
    for raw in lines:
        s = raw.split("#", 1)[0].strip()
        low = s.lower()
        if low.startswith("<band.kpath") and not low.startswith("<band.kpath.unitcell"):
            inside = True
            continue
        if inside and low.startswith("band.kpath>"):
            break
        if inside and s:
            p = s.split()
            if len(p) < 9:
                die("%s: <Band.kpath line needs 9 items: %s" % (datfile, s))
            out.append((p[7], p[8]))
    if not out:
        die("no <Band.kpath block found in %s" % datfile)
    return out


# ----------------------------------------------------------------------------
# labels / ticks
# ----------------------------------------------------------------------------
def pretty(label):
    return "G" if label.lower() in ("g", "gamma", "gm", "gam") else label


def tex(label):
    return "|".join(r"$\Gamma$" if p == "G" else p for p in label.split("|"))


def make_ticks(data, args):
    x, bounds = data["x"], data["bounds"]
    nseg = len(bounds)
    seg_labels = data.get("seg_labels")
    if args.dat:
        lab = read_dat_kpath(args.dat)
        if len(lab) != nseg:
            print("warning: %s has %d k-path segments but the band file has %d; labels ignored"
                  % (args.dat, len(lab), nseg))
        else:
            seg_labels = lab
    if args.labels:
        toks = args.labels.split()
        if len(toks) == nseg + 1:
            seg_labels = [(toks[i], toks[i + 1]) for i in range(nseg)]
        elif len(toks) == 2 * nseg:
            seg_labels = [(toks[2 * i], toks[2 * i + 1]) for i in range(nseg)]
        else:
            die("--labels needs %d (continuous) or %d labels, got %d" % (nseg + 1, 2 * nseg, len(toks)))
    ticks = []
    if seg_labels:
        for (i0, i1), (l0, l1) in zip(bounds, seg_labels):
            ticks.append((x[i0], pretty(l0)))
            ticks.append((x[i1], pretty(l1)))
    elif data.get("gnuband_tics"):
        ticks = [(xv, "|".join(pretty(p) for p in lab.split(","))) for lab, xv in data["gnuband_tics"]]
    else:
        for i0, i1 in bounds:
            ticks.append((x[i0], ""))
            ticks.append((x[i1], ""))
    merged = []
    for xv, lab in ticks:
        if merged and abs(merged[-1][0] - xv) < 1e-6:
            if lab and lab != merged[-1][1]:
                merged[-1] = (xv, (merged[-1][1] + "|" + lab) if merged[-1][1] else lab)
        else:
            merged.append((xv, lab))
    # label of every k index that is a segment end point (for the gap report)
    point_label = {}
    if seg_labels:
        for (i0, i1), (l0, l1) in zip(bounds, seg_labels):
            point_label[i0] = pretty(l0)
            point_label[i1] = pretty(l1)
    return merged, point_label


# ----------------------------------------------------------------------------
# analysis
# ----------------------------------------------------------------------------
def describe_k(data, idx, point_label):
    s = "x=%.4f" % data["x"][idx]
    if data.get("kfrac") is not None:
        s += " k=(%.4f %.4f %.4f)" % tuple(data["kfrac"][idx])
    if idx in point_label:
        s = point_label[idx] + " " + s
    return s


def gap_report(data, point_label, nocc=None):
    import numpy as np

    eig = data["eig"]
    nspin, nk, nb = eig.shape
    if nocc is not None:
        if not 0 < nocc < nb:
            die("--nocc must be between 1 and %d" % (nb - 1))
        occ_mask = np.zeros_like(eig, dtype=bool)
        occ_mask[:, :, :nocc] = True
    else:
        occ_mask = eig <= 0.0
    crossing = []
    for s in range(nspin):
        for b in range(nb):
            col = occ_mask[s, :, b]
            if col.any() and (~col).any():
                crossing.append((s, b))
    res = {}
    if data.get("mu_ev") is not None:
        print("chemical potential mu = %.4f eV (absolute; set to 0 in the plot)" % data["mu_ev"])
    if crossing and nocc is None:
        txt = ", ".join("spin%d/band%d" % (s, b + 1) for s, b in crossing[:8])
        print("bands crossing E_F: %s%s -> metallic (or semimetal) along this path"
              % (txt, " ..." if len(crossing) > 8 else ""))
    occ = np.where(occ_mask, eig, -np.inf)
    emp = np.where(~occ_mask, eig, np.inf)
    if not np.isfinite(occ).any() or not np.isfinite(emp).any():
        print("no occupied or no empty state in the energy list; gap not evaluated")
        return res
    channels = [("all", slice(None))]
    if nspin == 2:
        channels = [("spin up", slice(0, 1)), ("spin down", slice(1, 2)), ("overall", slice(None))]
    for name, sl in channels:
        o, e = occ[sl], emp[sl]
        sv, kv, bv = np.unravel_index(np.argmax(o), o.shape)
        sc, kc, bc = np.unravel_index(np.argmin(e), e.shape)
        vbm, cbm = o[sv, kv, bv], e[sc, kc, bc]
        sv += sl.start or 0
        sc += sl.start or 0
        print("[%s]" % name)
        print("  VBM = %9.4f eV at %s (spin %d, band %d)" % (vbm, describe_k(data, kv, point_label), sv, bv + 1))
        print("  CBM = %9.4f eV at %s (spin %d, band %d)" % (cbm, describe_k(data, kc, point_label), sc, bc + 1))
        if crossing and nocc is None:
            continue
        kind = "direct" if kv == kc else "indirect"
        direct = np.min(np.min(e, axis=2).min(axis=0) - np.max(o, axis=2).max(axis=0))
        print("  gap along the path = %.4f eV (%s); smallest direct gap = %.4f eV" % (cbm - vbm, kind, direct))
        res[name] = (vbm, cbm)
    print("note: VBM/CBM are searched only on the k points of the path")
    return res


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        prog="omx_band_plot.py",
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bandfile", help="System.Name.Band or System.Name.BANDDAT1")
    ap.add_argument("--dat", default="", help="OpenMX input (.dat) to take <Band.kpath labels from")
    ap.add_argument("--labels", default="", help='override labels, e.g. "L G X W K G"')
    ap.add_argument("--emin", type=float, default=-6.0, help="lower energy limit in eV (default -6)")
    ap.add_argument("--emax", type=float, default=6.0, help="upper energy limit in eV (default 6)")
    ap.add_argument("--ref", choices=["mu", "vbm"], default="mu",
                    help="energy zero: chemical potential (default) or VBM (insulators only)")
    ap.add_argument("--shift", type=float, default=0.0, help="extra rigid shift in eV added to all bands")
    ap.add_argument("--nocc", type=int, default=None,
                    help="number of occupied bands per spin channel (instead of E<=mu)")
    ap.add_argument("--spin", choices=["both", "up", "down"], default="both",
                    help="which spin channel to draw for collinear runs (default both)")
    ap.add_argument("--compare", default="", help="second .Band/BANDDAT1 drawn as dashed lines (same path)")
    ap.add_argument("--compare-label", default="compare", help="legend label for --compare")
    ap.add_argument("--label", default="", help="legend label for the main file")
    ap.add_argument("--title", default="", help="figure title")
    ap.add_argument("--color", default="", help="line color for the main bands (non-spin)")
    ap.add_argument("--figsize", type=float, nargs=2, default=[5.0, 4.5], metavar=("W", "H"))
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("-o", "--out", default="", help="output image (default <stem>_band.png)")
    ap.add_argument("--csv", default="", help="write x and band energies (eV) to this CSV file")
    ap.add_argument("--no-plot", action="store_true", help="only print the analysis")
    ap.add_argument("--show", action="store_true", help="open an interactive window instead of Agg")
    args = ap.parse_args()

    import numpy as np

    data = read_any(args.bandfile)
    nspin, nk, nb = data["eig"].shape
    print("%s: %d segment(s), %d k points, %d spin channel(s), %d bands"
          % (args.bandfile, len(data["bounds"]), nk, nspin, nb))
    if data["kind"] == "Band" and data["noncollinear_possible"]:
        print("note: spin_flag = 0 -> non-spin-polarized OR non-collinear/SOC run "
              "(in the latter case every band holds one electron)")
    ticks, point_label = make_ticks(data, args)
    print("ticks: " + "  ".join("%s=%.4f" % (lab or "?", xv) for xv, lab in ticks))
    res = gap_report(data, point_label, args.nocc)

    zero = 0.0
    if args.ref == "vbm":
        key = "overall" if "overall" in res else "all"
        if key not in res:
            die("--ref vbm needs a gapped system (or use --nocc)")
        zero = res[key][0]
        print("energy zero moved to the VBM (%.4f eV relative to mu)" % zero)
    eig = data["eig"] - zero + args.shift

    comp = None
    if args.compare:
        comp = read_any(args.compare)
        ce = comp["eig"] - zero + args.shift
        cx = comp["x"]
        if abs(cx[-1] - data["x"][-1]) > 1e-3 * max(abs(data["x"][-1]), 1e-12):
            print("warning: path lengths differ (%.4f vs %.4f); x of --compare rescaled"
                  % (cx[-1], data["x"][-1]))
            cx = cx / cx[-1] * data["x"][-1] if cx[-1] else cx
        comp = (cx, ce)

    if args.csv:
        cols, head = [data["x"]], ["x_Bohr^-1"]
        for s in range(nspin):
            for b in range(nb):
                cols.append(eig[s, :, b])
                head.append(("up_" if nspin == 2 and s == 0 else "dn_" if nspin == 2 else "") + "b%d" % (b + 1))
        np.savetxt(args.csv, np.column_stack(cols), delimiter=",", fmt="%.6f",
                   header=",".join(head), comments="")
        print("band energies written to %s (eV, zero = %s)" % (args.csv, args.ref))

    if args.no_plot:
        return
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    x = data["x"]
    spins = range(nspin)
    if nspin == 2 and args.spin != "both":
        spins = [0] if args.spin == "up" else [1]
    colors = {0: "#1f4e79", 1: "#c0504d"} if nspin == 2 else {0: args.color or "#1f4e79"}
    names = {0: "spin up", 1: "spin down"} if nspin == 2 else {0: args.label or None}
    for s in spins:
        for b in range(nb):
            ax.plot(x, eig[s, :, b], color=colors[s], lw=1.1,
                    label=(names[s] if b == 0 else None))
    if comp is not None:
        cx, ce = comp
        for s in range(ce.shape[0]):
            for b in range(ce.shape[2]):
                ax.plot(cx, ce[s, :, b], color="0.45", lw=0.9, ls="--",
                        label=(args.compare_label if (s == 0 and b == 0) else None))
    for xv, _ in ticks:
        ax.axvline(xv, color="k", lw=0.6)
    ax.axhline(0.0, color="0.5", lw=0.8, ls=":")
    ax.set_xticks([xv for xv, _ in ticks])
    ax.set_xticklabels([tex(lab) for _, lab in ticks])
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(args.emin, args.emax)
    ax.set_ylabel(r"$E - E_\mathrm{F}$ (eV)" if args.ref == "mu" else r"$E - E_\mathrm{VBM}$ (eV)")
    if args.title:
        ax.set_title(args.title)
    if nspin == 2 or comp is not None or args.label:
        ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    if args.show:
        plt.show()
    else:
        stem = re.sub(r"\.(Band|BANDDAT[12])$", "", os.path.basename(args.bandfile))
        out = args.out or stem + "_band.png"
        fig.savefig(out, dpi=args.dpi)
        print("figure written to %s" % out)


if __name__ == "__main__":
    main()
