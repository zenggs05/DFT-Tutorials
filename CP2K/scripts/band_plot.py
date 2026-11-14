#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
band_plot.py -- plot a CP2K band structure (&DFT/&PRINT/&BAND_STRUCTURE).
                解析 CP2K 能带文件（.bs 的 Set/Point/Band 块），按 KPOINT_SET 生成横轴和
                高对称点标签，以 VBM 或费米能为零点作图，判断直接/间接带隙并导出 CSV。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Input
-----
* The file named by &BAND_STRUCTURE/FILE_NAME (written in append mode!):
      # Set 1: 4 special points, 91 k-points, 20 bands
      #  Special point 1   0.5 0.5 0.5  L
      #  Point 1  Spin 1:  0.5 0.5 0.5  0.0109
      #   Band    Energy [eV]     Occupation
             1     -5.8     2.0
  If the same file contains several runs (a repeated "# Set 1"), only the last
  run is used.  Energies are absolute and in eV.
* Or the main CP2K output when FILE_NAME was not given: then the "#  Point"
  blocks are printed there and the special points are taken from the
  "KPOINTS| Special point" lines.

k-path (x axis)
---------------
k-points are given in fractions of the reciprocal lattice vectors (B_VECTOR).
To get a true Cartesian path length (1/Angstrom) the cell is needed:
  --cell-from run.out   reads the "CELL| Vector a/b/c [angstrom]" lines
  --lattice 9 numbers   a1 a2 a3 b1 b2 b3 c1 c2 c3 in Angstrom
If the .bs data come from a main output, the cell is read from it.
Without a cell the path length is computed in fractional units (a warning
is printed; the shape of the plot is then distorted for non-cubic cells).
Consecutive KPOINT_SETs are joined; if the end point of one set differs
from the start of the next, the labels are merged as "U|K" (path break).

Energy zero (--zero)
--------------------
  vbm   : highest occupied band energy (occupation > 0.5 * max occupation;
          use --nocc N to fix the number of occupied bands per spin instead)
  fermi : --efermi X (eV) or the last "Fermi energy:" line of --fermi-from OUT
  none  : absolute energies

Usage examples
--------------
    python band_plot.py si_band.bs --cell-from si_band_structure.out
    python band_plot.py si_band.bs --lattice 0 2.7155 2.7155 2.7155 0 2.7155 2.7155 2.7155 0 \
                         --erange -13 6 --csv si_band.csv
    python band_plot.py graphene.bs --cell-from graphene.out --zero fermi --fermi-from graphene.out
    python band_plot.py si_band_structure.out --nocc 4 -o si_band.png

Dependencies: python >= 3.9, numpy, matplotlib.
"""
import argparse
import os
import re
import sys

HA2EV = 27.211386245988
NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
COLORS = ["#2a78d6", "#eb6834"]   # spin 1, spin 2


def die(msg):
    sys.exit("band_plot.py: error: " + msg)


def need_numpy():
    try:
        import numpy as np
    except ImportError:
        die("numpy is required (pip install numpy)")
    return np


# --------------------------------------------------------------------------- parsing
def new_set(idx):
    return {"index": idx, "special": [], "points": []}


def read_band_file(path):
    """Return a list of k-point sets. Each set: special=[(label, k)], points=[{k, w, spin, e, occ}]."""
    try:
        with open(path, "r", errors="replace") as fh:
            lines = fh.readlines()
    except OSError as exc:
        die(f"cannot read {path}: {exc}")

    sets, cur, pt = [], None, None
    runs = 0
    main_output = False
    re_set = re.compile(r"#\s*Set\s+(\d+):\s*(\d+)\s+special points,\s*(\d+)\s+k-points,\s*(\d+)\s+bands")
    re_sp = re.compile(r"#\s*Special point\s+(\d+)\s+(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*(.*)$")
    re_pt = re.compile(r"#\s*Point\s+(\d+)\s+Spin\s+(\d+):\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")")
    re_kpset = re.compile(r"KPOINTS\|\s*Number of k-points in set\s+(\d+)\s+(\d+)")
    re_kpsp = re.compile(r"KPOINTS\|\s*Special point\s+(\d+)\s+(.*?)\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*$")

    for raw in lines:
        line = raw.rstrip("\n")
        s = line.strip()
        m = re_set.match(s)
        if m:
            idx = int(m.group(1))
            if idx == 1 and sets:
                runs += 1
                sets = []
            cur = new_set(idx)
            sets.append(cur)
            pt = None
            continue
        m = re_kpset.search(line)
        if m:
            main_output = True
            idx = int(m.group(1))
            if idx == 1 and sets:
                runs += 1
                sets = []
            cur = new_set(idx)
            sets.append(cur)
            pt = None
            continue
        m = re_sp.match(s)
        if m and cur is not None:
            k = [float(m.group(i)) for i in (2, 3, 4)]
            cur["special"].append((m.group(5).strip(), k))
            continue
        m = re_kpsp.search(line)
        if m and cur is not None and main_output:
            k = [float(m.group(i)) for i in (3, 4, 5)]
            cur["special"].append((m.group(2).strip(), k))
            continue
        m = re_pt.match(s)
        if m:
            if cur is None:            # data without any set header
                cur = new_set(1)
                sets.append(cur)
            pt = {"ik": int(m.group(1)), "spin": int(m.group(2)),
                  "k": [float(m.group(i)) for i in (3, 4, 5)], "w": float(m.group(6)), "e": [], "occ": []}
            cur["points"].append(pt)
            continue
        if s.startswith("#"):
            continue
        if pt is not None:
            toks = s.split()
            if len(toks) == 3 and re.fullmatch(r"\d+", toks[0]):
                try:
                    pt["e"].append(float(toks[1]))
                    pt["occ"].append(float(toks[2]))
                except ValueError:
                    pt = None
            else:
                pt = None          # end of a band block (e.g. other text in the main output)
    if runs:
        print(f"note: {path} contains {runs + 1} appended runs; using the last one "
              "(delete the .bs file before re-running CP2K)")
    sets = [st for st in sets if st["points"]]
    if not sets:
        die(f"{path}: no '#  Point ... Spin ...' blocks found - is this a CP2K band-structure file?")
    return sets, main_output


def read_cell(path):
    """Cell vectors (rows, Angstrom) from the first 'CELL| Vector x [angstrom]:' lines."""
    np = need_numpy()
    vec = {}
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                m = re.search(r"\bCELL\|\s*Vector\s+([abc])\s+\[angstrom\]:\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")", line)
                if m and m.group(1) not in vec:
                    vec[m.group(1)] = [float(m.group(i)) for i in (2, 3, 4)]
                if len(vec) == 3:
                    break
    except OSError as exc:
        die(f"cannot read {path}: {exc}")
    if len(vec) < 3:
        return None
    return np.array([vec["a"], vec["b"], vec["c"]])


def read_fermi(path):
    ef = None
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                m = re.search(r"Fermi energy:\s*(" + NUM + r")", line)
                if m:
                    ef = float(m.group(1)) * HA2EV
    except OSError as exc:
        die(f"cannot read {path}: {exc}")
    if ef is None:
        die(f"{path}: no 'Fermi energy:' line (it is printed only with smearing); use --efermi")
    return ef


# --------------------------------------------------------------------------- main
def build_parser():
    ap = argparse.ArgumentParser(
        prog="band_plot.py",
        description="Plot a CP2K band structure from a &BAND_STRUCTURE .bs file (or the main output). "
                    "解析 .bs 文件作能带图、判断带隙、导出 CSV。",
        epilog="Examples:\n"
               "  python band_plot.py si_band.bs --cell-from si_band_structure.out --csv si_band.csv\n"
               "  python band_plot.py graphene.bs --cell-from graphene.out --zero fermi --fermi-from graphene.out\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bsfile", help=".bs file of &BAND_STRUCTURE (or a main output containing the band data)")
    ap.add_argument("--cell-from", metavar="OUT", help="CP2K main output to read the cell vectors from")
    ap.add_argument("--lattice", nargs=9, type=float, metavar="X",
                    help="cell vectors a b c as 9 numbers in Angstrom (row-wise)")
    ap.add_argument("--zero", choices=["vbm", "fermi", "none"], default="vbm", help="energy zero (default vbm)")
    ap.add_argument("--efermi", type=float, metavar="EV", help="Fermi energy in eV (absolute) for --zero fermi")
    ap.add_argument("--fermi-from", metavar="OUT", help="read the last 'Fermi energy:' (hartree) from this output")
    ap.add_argument("--nocc", type=int, metavar="N",
                    help="number of occupied bands per spin channel (instead of using the occupations)")
    ap.add_argument("--erange", nargs=2, type=float, default=[-6.0, 6.0], metavar=("EMIN", "EMAX"),
                    help="energy window relative to the zero, eV (default -6 6)")
    ap.add_argument("--labels", nargs="+", metavar="L",
                    help="override the special-point labels along the whole path (use 'U|K' for breaks)")
    ap.add_argument("--csv", metavar="FILE",
                    help="write x, kx, ky, kz and all bands (columns per spin) to CSV; ticks go to FILE_ticks.csv")
    ap.add_argument("-o", "--output", default="band.png", help="figure file (default band.png; 'none' = no figure)")
    ap.add_argument("--title", default=None, help="figure title")
    ap.add_argument("--dpi", type=int, default=200, help="figure resolution (default 200)")
    ap.add_argument("--show", action="store_true", help="open an interactive window instead of only saving")
    return ap


def pretty(label):
    lab = label.strip() or "?"
    rep = {"GAMMA": r"$\Gamma$", "G": r"$\Gamma$", "GM": r"$\Gamma$", "SIGMA": r"$\Sigma$", "DELTA": r"$\Delta$",
           "LAMBDA": r"$\Lambda$"}
    return "|".join(rep.get(p.upper(), p) for p in lab.split("|"))


def main():
    args = build_parser().parse_args()
    np = need_numpy()
    if args.erange[0] >= args.erange[1]:
        die("--erange: EMIN must be smaller than EMAX")
    sets, main_output = read_band_file(args.bsfile)

    # reciprocal lattice ---------------------------------------------------
    cell = None
    if args.lattice:
        cell = np.array(args.lattice).reshape(3, 3)
    elif args.cell_from:
        cell = read_cell(args.cell_from)
        if cell is None:
            die(f"{args.cell_from}: no 'CELL| Vector a/b/c [angstrom]' lines found")
    elif main_output:
        cell = read_cell(args.bsfile)
    if cell is not None:
        if abs(np.linalg.det(cell)) < 1e-8:
            die("cell vectors are linearly dependent")
        bmat = 2.0 * np.pi * np.linalg.inv(cell).T     # rows b1 b2 b3, 1/Angstrom
        xunit = r"k-path (1/$\mathrm{\AA}$)"
    else:
        print("warning: no cell given (--cell-from / --lattice); path length in fractional units")
        bmat = np.eye(3)
        xunit = "k-path (fractional units)"

    # arrays per set ------------------------------------------------------
    nspin = max(p["spin"] for st in sets for p in st["points"])
    nband = min(len(p["e"]) for st in sets for p in st["points"])
    if nband == 0:
        die("empty band blocks")
    if any(len(p["e"]) != nband for st in sets for p in st["points"]):
        print(f"warning: different band counts per k-point; truncating to {nband}")

    xs, ks, es, occs, setid = [], [], [], [], []
    ticks = []          # (x, label)
    x0 = 0.0
    prev_end_k = None
    for st in sets:
        kp = {}
        for p in st["points"]:
            kp.setdefault(p["ik"], {})[p["spin"]] = p
        iks = sorted(kp)
        if any(len(kp[i]) != nspin for i in iks):
            die(f"set {st['index']}: missing spin blocks for some k-points")
        kfrac = np.array([kp[i][1]["k"] for i in iks])
        e = np.array([[kp[i][s]["e"][:nband] for i in iks] for s in range(1, nspin + 1)])
        o = np.array([[kp[i][s]["occ"][:nband] for i in iks] for s in range(1, nspin + 1)])
        kc = kfrac @ bmat
        seg = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kc, axis=0), axis=1))]) if len(kc) > 1 else np.zeros(1)
        x = x0 + seg
        # ticks: locate each special point along the set
        spec = st["special"]
        jstart = 0
        for j, (lab, ksp) in enumerate(spec):
            d = np.linalg.norm(kfrac[jstart:] - np.array(ksp), axis=1)
            if d.size == 0:
                break
            jj = jstart + int(np.argmin(d))
            if d.min() > 1e-4:
                print(f"warning: special point {lab or j + 1} of set {st['index']} not found on the path")
                continue
            ticks.append([float(x[jj]), lab or f"P{j + 1}", st["index"], j == 0])
            jstart = jj
        if prev_end_k is not None and spec:
            # merge first tick of this set with the last tick of the previous one
            first = [t for t in ticks if t[2] == st["index"] and t[3]]
            prev = [t for t in ticks if t[2] != st["index"]]
            if first and prev:
                same = np.linalg.norm(np.array(spec[0][1]) - prev_end_k) < 1e-6
                if not same:
                    prev[-1][1] = prev[-1][1] + "|" + first[0][1]
                ticks.remove(first[0])
        prev_end_k = kfrac[-1]
        xs.append(x)
        ks.append(kfrac)
        es.append(e)
        occs.append(o)
        setid.append(np.full(len(x), st["index"]))
        x0 = x[-1]

    X = np.concatenate(xs)
    K = np.concatenate(ks)
    E = np.concatenate(es, axis=1)        # (nspin, nk, nband)
    O = np.concatenate(occs, axis=1)
    S = np.concatenate(setid)
    print(f"{len(sets)} k-point set(s), {X.size} k-points, {nband} bands, {nspin} spin channel(s)")

    if args.labels:
        if len(args.labels) != len(ticks):
            die(f"--labels: {len(args.labels)} labels given but the path has {len(ticks)} ticks")
        for t, lab in zip(ticks, args.labels):
            t[1] = lab

    # occupied / empty -----------------------------------------------------
    if args.nocc is not None:
        if not 0 < args.nocc < nband:
            die(f"--nocc must be between 1 and {nband - 1}")
        occmask = np.zeros_like(E, dtype=bool)
        occmask[:, :, :args.nocc] = True
    else:
        omax = O.max()
        if omax <= 0:
            die("all occupations are zero; use --nocc N")
        occmask = O > 0.5 * omax
    vbm = E[occmask].max() if occmask.any() else None
    cbm = E[~occmask].min() if (~occmask).any() else None

    def where(val, mask):
        idx = np.argwhere(np.isclose(E, val, atol=1e-6) & mask)[0]
        s, ik, ib = idx
        return int(s), int(ik), int(ib)

    if vbm is not None and cbm is not None:
        sv, kv, bv = where(vbm, occmask)
        sc, kc_, bc = where(cbm, ~occmask)
        gap = cbm - vbm
        print(f"VBM = {vbm:.4f} eV  band {bv + 1}, spin {sv + 1}, k = ({K[kv][0]:.4f} {K[kv][1]:.4f} {K[kv][2]:.4f})")
        print(f"CBM = {cbm:.4f} eV  band {bc + 1}, spin {sc + 1}, k = ({K[kc_][0]:.4f} {K[kc_][1]:.4f} {K[kc_][2]:.4f})")
        if gap <= 0:
            print(f"bands overlap (CBM - VBM = {gap:.4f} eV): metal / semimetal on this path")
        else:
            direct = np.linalg.norm(K[kv] - K[kc_]) < 1e-6
            dmin = (np.where(~occmask, E, np.inf).min(axis=2) - np.where(occmask, E, -np.inf).max(axis=2))
            print(f"band gap on this path = {gap:.4f} eV ({'direct' if direct else 'indirect'}); "
                  f"smallest direct gap = {dmin.min():.4f} eV")
            print("note: a gap read from a k-path is only as good as the path sampling (NPOINTS)")

    # zero ---------------------------------------------------------------
    if args.zero == "vbm":
        if vbm is None:
            die("no occupied band found; use --nocc or --zero fermi/none")
        ref, zl = vbm, r"E - E$_\mathrm{VBM}$ (eV)"
    elif args.zero == "fermi":
        if args.efermi is not None:
            ref = args.efermi
        elif args.fermi_from:
            ref = read_fermi(args.fermi_from)
        elif main_output:
            ref = read_fermi(args.bsfile)
        else:
            die("--zero fermi needs --efermi EV or --fermi-from OUT")
        zl = r"E - E$_\mathrm{F}$ (eV)"
    else:
        ref, zl = 0.0, "E (eV)"
    print(f"energy zero = {ref:.4f} eV ({args.zero})")
    Er = E - ref

    # CSV -------------------------------------------------------------------
    if args.csv:
        cols = [X, K[:, 0], K[:, 1], K[:, 2], S]
        heads = ["x", "kx_frac", "ky_frac", "kz_frac", "set"]
        for s in range(nspin):
            for b in range(nband):
                cols.append(Er[s, :, b])
                heads.append(f"s{s + 1}_b{b + 1}")
        np.savetxt(args.csv, np.column_stack(cols), delimiter=",", header=",".join(heads), comments="", fmt="%.6f")
        root, ext = os.path.splitext(args.csv)
        tfile = root + "_ticks" + (ext or ".csv")
        with open(tfile, "w") as fh:
            fh.write("x,label\n")
            for t in ticks:
                fh.write(f"{t[0]:.6f},{t[1]}\n")
        print(f"wrote {args.csv} and {tfile} (energies relative to the zero, eV)")

    if args.output.lower() == "none" and not args.show:
        return
    try:
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        die("matplotlib is required for plotting (pip install matplotlib); use -o none --csv FILE otherwise")

    fig, ax = plt.subplots(figsize=(4.8, 4.2))
    names = ["spin up", "spin down"] if nspin == 2 else [None]
    for s in range(nspin):
        for b in range(nband):
            for sid in np.unique(S):            # do not connect separate sets
                m = S == sid
                ax.plot(X[m], Er[s, m, b], color=COLORS[s], lw=1.2, ls="-" if s == 0 else "--",
                        label=names[s] if (b == 0 and sid == S[0]) else None)
    for t in ticks:
        ax.axvline(t[0], color="#b5b4ad", lw=0.6)
    if args.zero != "none":
        ax.axhline(0.0, color="#8a8a86", lw=0.8, ls=":")
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([pretty(t[1]) for t in ticks])
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(args.erange)
    ax.set_ylabel(zl)
    ax.set_xlabel(xunit, fontsize=8, color="#5f5e5a")
    if nspin == 2:
        ax.legend(frameon=False, fontsize=8, loc="upper right")
    if args.title:
        ax.set_title(args.title)
    fig.tight_layout()
    if args.output.lower() != "none":
        fig.savefig(args.output, dpi=args.dpi)
        print(f"wrote {args.output}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
