#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pdos_plot.py -- plot CP2K PDOS files (legacy &PRINT/&PDOS and 2026.2 &PRINT/&DOS).
                读取 CP2K 的 .pdos 文件（旧版逐态格式 + 2026.2 的 &CURVE 展宽曲线），
                展宽、对齐能量零点、按元素/角动量求和，分自旋作图并导出 CSV。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Supported input files
---------------------
1. Raw (unbroadened) per-state files, one per atomic kind (and per spin):
     <PROJECT>-k1-1.pdos, <PROJECT>-ALPHA_k1-1.pdos, <PROJECT>-list1-1.pdos ...
   * CP2K <= 2026.1 (&PRINT/&PDOS):
       # Projected DOS for atomic kind Ti at iteration step i = 0, E(Fermi) = -0.1 a.u.
       #     MO Eigenvalue [a.u.]      Occupation     s     p     d
   * CP2K 2026.2 (&PRINT/&DOS/&PDOS): E(Fermi) and E(HOCO) on their own
     lines and an extra "Total" column after "Occupation".
   Columns are located by the header names, never by fixed position.
   Energies are absolute and in hartree; this script converts them to eV and
   broadens them (Gaussian or Lorentzian, FWHM given by --fwhm).
2. Broadened curve files of CP2K 2026.2 (&DOS/&CURVE):
     <PROJECT>-k1_curve-1.pdos, <PROJECT>-ALPHA_k1_curve-1.pdos ...
       # Energy zero: HOCO
       # Energy[eV]  total  s  py  pz  px ...
   They are already broadened and shifted by CP2K (ENERGY_ZERO); they are
   plotted as they are unless --zero asks for another reference.
Raw and curve files cannot be mixed in one plot (different zero/broadening).
If a file contains several appended blocks (e.g. &EACH in GEO_OPT), the last
block is used.

Energy zero (--zero)
--------------------
  fermi : E(Fermi) from the file headers (max over spins, as CP2K does)
  hoco  : highest occupied state; header E(HOCO) (2026.2) or, for old files,
          the highest state with occupation > 0.5 * max occupation
  none  : absolute energies
  file  : (curve files only, default for them) keep the zero CP2K used
  --efermi X overrides the reference with an absolute value X in eV.

Usage examples
--------------
    python pdos_plot.py                              # all *.pdos in the current dir
    python pdos_plot.py tio2_pdos-k*-1.pdos --group element-l --erange -8 6
    python pdos_plot.py *-k*-1.pdos --zero hoco --fwhm 0.3 --lorentzian -o h2o_pdos.png
    python pdos_plot.py h2o_pdos-k*_curve-1.pdos --group orbital --elements O
    python pdos_plot.py fe-ALPHA_k1-1.pdos fe-BETA_k1-1.pdos --group l --csv fe_pdos.csv
    python pdos_plot.py *.pdos --group total --total --show

Dependencies: python >= 3.9, numpy, matplotlib.
"""
import argparse
import glob
import os
import re
import sys

HA2EV = 27.211386245988
NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"

# categorical colours in a fixed order (never cycled; extra series get new line styles)
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
STYLES = ["-", "--", ":", "-."]
L_ORDER = ["s", "p", "d", "f", "g", "h"]


def die(msg):
    sys.exit("pdos_plot.py: error: " + msg)


def need_numpy():
    try:
        import numpy as np  # noqa: F401
    except ImportError:
        die("numpy is required (pip install numpy)")
    return np


# --------------------------------------------------------------------------- parsing
class PdosFile:
    """Content of one .pdos file (last block)."""

    def __init__(self, path):
        self.path = path
        self.kind = "?"           # atomic kind name, 'list1', 'r_list1'
        self.spin = 0             # 0 = closed shell, 1 = ALPHA, 2 = BETA
        self.curve = False
        self.efermi = None        # eV, absolute
        self.hoco = None          # eV, absolute
        self.zero_label = None    # curve files: FERMI / HOCO / ABSOLUTE
        self.energy = None        # eV (absolute for raw files, as written for curves)
        self.occ = None           # raw files only
        self.total = None         # curve: total column; raw: Total column if present
        self.channels = {}        # label -> weights (raw) or density (curve)


def spin_from_name(path):
    base = os.path.basename(path).upper()
    if re.search(r"(^|[-_])ALPHA_", base):
        return 1
    if re.search(r"(^|[-_])BETA_", base):
        return 2
    return 0


def split_blocks(lines):
    """Split a file into blocks, each starting with a '# ... DOS ...' title line."""
    blocks, cur = [], []
    for ln in lines:
        if ln.startswith("#") and re.match(r"#\s*(Projected DOS|DOS at)", ln):
            if cur:
                blocks.append(cur)
            cur = [ln]
        else:
            cur.append(ln)
    if cur:
        blocks.append(cur)
    return blocks


def parse_float_row(line):
    vals = []
    for tok in line.split():
        try:
            vals.append(float(tok.replace("D", "E").replace("d", "e")))
        except ValueError:
            # Fortran overflow like 0.12345-100 (exponent without E)
            m = re.fullmatch(r"([-+]?\d*\.\d+)([-+]\d{3})", tok)
            if m:
                vals.append(float(m.group(1) + "E" + m.group(2)))
            else:
                return None
    return vals


def read_pdos(path):
    np = need_numpy()
    try:
        with open(path, "r", errors="replace") as fh:
            lines = fh.readlines()
    except OSError as exc:
        die(f"cannot read {path}: {exc}")
    if not lines:
        die(f"{path} is empty")
    blocks = split_blocks(lines)
    blk = blocks[-1]
    if len(blocks) > 1:
        print(f"note: {os.path.basename(path)} contains {len(blocks)} appended blocks, using the last one")

    pf = PdosFile(path)
    pf.spin = spin_from_name(path)
    header = [ln for ln in blk if ln.startswith("#")]
    text = "".join(header)
    m = re.search(r"atomic kind\s+(\S+)", text)
    if m:
        pf.kind = m.group(1)
    else:
        m = re.search(r"for list\s+(\d+)", text)
        if m:
            pf.kind = "list" + m.group(1)
        elif "in real space" in text:
            pf.kind = "r_ldos"
        else:
            pf.kind = os.path.basename(path)
    m = re.search(r"E\(Fermi\)\s*=\s*(" + NUM + r")\s*a\.u\.", text)
    if m:
        pf.efermi = float(m.group(1)) * HA2EV
    m = re.search(r"E\(HOCO\)\s*=\s*(" + NUM + r")\s*a\.u\.", text)
    if m:
        pf.hoco = float(m.group(1)) * HA2EV
    m = re.search(r"Energy zero:\s*(.+)", text)
    if m:
        pf.zero_label = m.group(1).split()[-1].upper()   # 'AUTO -> FERMI' -> FERMI

    # column header = last comment line containing 'Occupation' or starting with '# Energy['
    colhead = None
    for ln in header:
        if "Occupation" in ln or re.match(r"#\s*Energy\[", ln):
            colhead = ln
    if colhead is None:
        die(f"{path}: no column header line found - is this a CP2K .pdos file?")
    toks = colhead.lstrip("#").split()

    rows = []
    for ln in blk:
        if ln.startswith("#") or not ln.strip():
            continue
        r = parse_float_row(ln)
        if r is None:
            print(f"warning: {os.path.basename(path)}: skipping malformed line: {ln.strip()[:60]}")
            continue
        rows.append(r)
    if not rows:
        die(f"{path}: no data lines")
    ncol = min(len(r) for r in rows)
    data = np.array([r[:ncol] for r in rows], dtype=float)

    if "Occupation" in toks and toks[0] == "MO":
        # raw per-state file: MO, energy (a.u.), occupation, [Total], channels...
        labels = toks[toks.index("Occupation") + 1:]
        pf.energy = data[:, 1] * HA2EV
        pf.occ = data[:, 2]
        for i, lab in enumerate(labels):
            col = 3 + i
            if col >= ncol:
                break
            if lab.lower() == "total":
                pf.total = data[:, col]
            else:
                pf.channels[lab] = data[:, col]
        if not pf.channels and pf.total is not None:
            pf.channels["total"] = pf.total
    else:
        # broadened curve: Energy[unit], total, components...
        pf.curve = True
        unit_ev = toks[0].lower().startswith("energy[ev")
        efac = 1.0 if unit_ev else HA2EV
        dfac = 1.0 if unit_ev else 1.0 / HA2EV
        pf.energy = data[:, 0] * efac
        labels = toks[1:]
        for i, lab in enumerate(labels):
            col = 1 + i
            if col >= ncol:
                break
            if lab.lower() == "total":
                pf.total = data[:, col] * dfac
            else:
                pf.channels[lab] = data[:, col] * dfac
        if pf.total is None and pf.channels:
            pf.total = sum(pf.channels.values())
        if pf.zero_label is None:
            pf.zero_label = "ABSOLUTE"
    if pf.efermi is None and not pf.curve:
        print(f"note: {os.path.basename(path)}: no E(Fermi) in the header")
    return pf


# --------------------------------------------------------------------------- helpers
def element_of(kind, by_kind):
    if by_kind:
        return kind
    m = re.match(r"([A-Z][a-z]?)", kind)
    if m and not kind.startswith("list") and not kind.startswith("r_"):
        return m.group(1)
    return kind


def l_of(label):
    lab = label.strip().lower()
    if lab in ("total", "ldos"):
        return "total"
    return lab[0] if lab and lab[0] in L_ORDER else lab


def broaden(np, grid, centers, weights, fwhm, shape):
    """Sum of normalised peaks (states/eV) on grid, evaluated in chunks."""
    out = np.zeros_like(grid)
    if centers.size == 0:
        return out
    chunk = max(1, int(2.0e6 // max(1, grid.size)))
    if shape == "lorentzian":
        gam = 0.5 * fwhm
        for i in range(0, centers.size, chunk):
            x = grid[:, None] - centers[None, i:i + chunk]
            out += (gam / np.pi / (x * x + gam * gam)) @ weights[i:i + chunk]
    else:
        sig = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        norm = 1.0 / (sig * np.sqrt(2.0 * np.pi))
        for i in range(0, centers.size, chunk):
            x = (grid[:, None] - centers[None, i:i + chunk]) / sig
            out += (norm * np.exp(-0.5 * x * x)) @ weights[i:i + chunk]
    return out


def collect_files(args):
    files = []
    pats = args.files or [os.path.join(args.dir, "*.pdos")]
    for p in pats:
        hits = sorted(glob.glob(p))
        if not hits and os.path.isfile(p):
            hits = [p]
        if not hits:
            if args.files:
                die(f"no file matches {p}")
        files += hits
    files = list(dict.fromkeys(files))
    if not files:
        die("no .pdos files found (give file names or use -d DIR)")
    return files


def reference_energy(np, pfs, how, user_ef):
    """Absolute reference energy (eV) for raw files."""
    if user_ef is not None:
        return user_ef, f"user value {user_ef:.4f} eV"
    if how == "none":
        return 0.0, "absolute energies"
    if how == "fermi":
        vals = [p.efermi for p in pfs if p.efermi is not None]
        if not vals:
            die("--zero fermi: no E(Fermi) in the file headers; use --efermi or --zero hoco/none")
        return max(vals), "E(Fermi)"
    if how == "hoco":
        vals = [p.hoco for p in pfs if p.hoco is not None]
        if vals:
            return max(vals), "E(HOCO) from header"
        est = []
        for p in pfs:
            if p.occ is None or p.occ.size == 0:
                continue
            occmax = p.occ.max()
            sel = p.occ > 0.5 * occmax if occmax > 0 else np.zeros_like(p.occ, bool)
            if sel.any():
                est.append(p.energy[sel].max())
        if not est:
            die("--zero hoco: cannot determine the highest occupied state")
        return max(est), "highest state with occ > 0.5*max (old format)"
    die(f"--zero {how} is not valid for raw files")
    return 0.0, ""


def curve_offset(pfs, how, user_ef):
    """Shift to add to curve energies so that they refer to the requested zero."""
    if how == "file" and user_ef is None:
        return 0.0, f"as written by CP2K ({pfs[0].zero_label})"
    labels = {p.zero_label for p in pfs}
    if len(labels) > 1:
        die(f"curve files use different energy zeros {sorted(labels)}; use --zero file")
    lab = labels.pop()
    ef = max([p.efermi for p in pfs if p.efermi is not None], default=None)
    ho = max([p.hoco for p in pfs if p.hoco is not None], default=None)
    old = {"ABSOLUTE": 0.0, "FERMI": ef, "HOCO": ho}.get(lab)
    if old is None:
        die(f"cannot undo energy zero '{lab}' (missing header value); use --zero file")
    if user_ef is not None:
        new, txt = user_ef, f"user value {user_ef:.4f} eV"
    elif how == "none":
        new, txt = 0.0, "absolute energies"
    elif how == "fermi":
        if ef is None:
            die("no E(Fermi) in the curve headers")
        new, txt = ef, "E(Fermi)"
    elif how == "hoco":
        if ho is None:
            die("no E(HOCO) in the curve headers")
        new, txt = ho, "E(HOCO)"
    else:
        new, txt = old, lab
    return old - new, txt


def gap_report(np, pfs, ref):
    """Print HOCO/LUMO estimate from raw per-state files."""
    occ_e, emp_e = [], []
    for p in pfs:
        if p.occ is None or p.occ.size == 0:
            continue
        occmax = p.occ.max()
        if occmax <= 0:
            continue
        sel = p.occ > 0.5 * occmax
        occ_e += list(p.energy[sel])
        emp_e += list(p.energy[~sel])
    if occ_e and emp_e:
        homo = max(occ_e)
        above = [e for e in emp_e if e > homo - 1e-6]
        if above:
            lumo = min(above)
            print(f"highest occupied state  = {homo:10.4f} eV (absolute), {homo - ref:+.4f} eV vs zero")
            print(f"lowest unoccupied state = {lumo:10.4f} eV (absolute), {lumo - ref:+.4f} eV vs zero")
            print(f"gap estimate (states in these files) = {lumo - homo:.4f} eV")


# --------------------------------------------------------------------------- main
def build_parser():
    ap = argparse.ArgumentParser(
        prog="pdos_plot.py",
        description="Plot CP2K PDOS: legacy/2026.2 per-state .pdos files (broadened here) "
                    "or 2026.2 &CURVE files (already broadened). 支持旧版 .pdos 与 2026.2 曲线文件。",
        epilog="Examples:\n"
               "  python pdos_plot.py tio2_pdos-k*-1.pdos --group element-l --erange -8 6\n"
               "  python pdos_plot.py h2o_pdos-k*_curve-1.pdos --group orbital --elements O\n"
               "  python pdos_plot.py *ALPHA_k*-1.pdos *BETA_k*-1.pdos --group l --csv pdos.csv\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", help="PDOS files or glob patterns (default: all *.pdos in -d DIR)")
    ap.add_argument("-d", "--dir", default=".", help="directory searched when no files are given (default: .)")
    ap.add_argument("--mode", choices=["auto", "raw", "curve"], default="auto",
                    help="which file type to use when both are present (auto: raw if any raw file, else curve)")
    ap.add_argument("--group", choices=["total", "element", "l", "element-l", "orbital"], default="element",
                    help="how to sum: total, per element, per angular momentum l, per element and l, "
                         "or every column (default: element)")
    ap.add_argument("--by-kind", action="store_true",
                    help="keep CP2K kind names (Fe1, Fe2 ...) instead of merging them into elements")
    ap.add_argument("--elements", nargs="+", metavar="EL", help="only plot these elements/kinds")
    ap.add_argument("--zero", choices=["fermi", "hoco", "none", "file"], default=None,
                    help="energy zero (default: fermi for raw files, file for curve files)")
    ap.add_argument("--efermi", type=float, metavar="EV", help="absolute reference energy in eV (overrides --zero)")
    ap.add_argument("--fwhm", type=float, default=None, help="broadening FWHM in eV for raw files (default 0.2)")
    ap.add_argument("--lorentzian", action="store_true", help="Lorentzian instead of Gaussian broadening (raw files)")
    ap.add_argument("--de", type=float, default=0.01, help="energy grid step in eV for raw files (default 0.01)")
    ap.add_argument("--erange", nargs=2, type=float, default=[-10.0, 5.0], metavar=("EMIN", "EMAX"),
                    help="energy window relative to the zero, eV (default -10 5)")
    ap.add_argument("--spin", choices=["both", "up", "down", "sum"], default="both",
                    help="spin-polarised files: both (beta mirrored below zero), up, down or sum (default both)")
    ap.add_argument("--total", action="store_true", help="also draw the sum of all loaded files")
    ap.add_argument("--scale", type=float, default=1.0, help="multiply all densities by this factor")
    ap.add_argument("--fill", action="store_true", help="fill the area under each curve (light alpha)")
    ap.add_argument("-o", "--output", default="pdos.png", help="figure file (default pdos.png; 'none' = no figure)")
    ap.add_argument("--csv", metavar="FILE", help="write energy + all plotted series to a CSV file")
    ap.add_argument("--title", default=None, help="figure title")
    ap.add_argument("--ymax", type=float, default=None, help="upper limit of the DOS axis")
    ap.add_argument("--dpi", type=int, default=200, help="figure resolution (default 200)")
    ap.add_argument("--show", action="store_true", help="open an interactive window instead of only saving")
    return ap


def main():
    args = build_parser().parse_args()
    np = need_numpy()
    if args.de <= 0:
        die("--de must be positive")
    if args.erange[0] >= args.erange[1]:
        die("--erange: EMIN must be smaller than EMAX")

    files = collect_files(args)
    pfs = [read_pdos(f) for f in files]
    raw = [p for p in pfs if not p.curve]
    cur = [p for p in pfs if p.curve]
    if args.mode == "raw" or (args.mode == "auto" and raw):
        if not raw:
            die("--mode raw: no raw per-state .pdos files among the inputs")
        if cur:
            print(f"note: ignoring {len(cur)} curve file(s); use --mode curve to plot them")
        use, is_curve = raw, False
    else:
        if not cur:
            die("--mode curve: no &CURVE files among the inputs")
        if raw and args.mode == "auto":
            print(f"note: ignoring {len(raw)} raw file(s)")
        use, is_curve = cur, True
    if args.elements:
        want = set(args.elements)
        use = [p for p in use if element_of(p.kind, args.by_kind) in want or p.kind in want]
        if not use:
            die("no file left after --elements filter")

    spins = sorted({p.spin for p in use})
    if 0 in spins and len(spins) > 1:
        die("mixing closed-shell and ALPHA/BETA files from different runs")
    polarized = spins != [0]

    # energy grid and zero -------------------------------------------------
    if is_curve:
        zero = args.zero or "file"
        shift, ztxt = curve_offset(use, zero, args.efermi)
        if args.fwhm is not None or args.lorentzian:
            print("note: curve files are already broadened by CP2K; --fwhm/--lorentzian ignored")
        grid = use[0].energy + shift
        for p in use[1:]:
            if p.energy.size != grid.size or np.abs(p.energy + shift - grid).max() > 1e-6:
                die(f"{p.path}: energy grid differs from {use[0].path}")
    else:
        zero = args.zero or "fermi"
        if zero == "file":
            die("--zero file only applies to curve files")
        ref, ztxt = reference_energy(np, use, zero, args.efermi)
        fwhm = 0.2 if args.fwhm is None else args.fwhm
        if fwhm <= 0:
            die("--fwhm must be positive")
        if args.de > fwhm / 3:
            print("warning: --de is coarse compared with --fwhm; peaks will look jagged")
        pad = 5.0 * fwhm
        grid = np.arange(args.erange[0] - pad, args.erange[1] + pad + 0.5 * args.de, args.de)
        shape = "lorentzian" if args.lorentzian else "gaussian"
        gap_report(np, use, ref)

    print(f"files used: {len(use)}  ({'curve' if is_curve else 'raw per-state'}, "
          f"{'spin-polarised' if polarized else 'closed shell'})")
    print(f"energy zero: {ztxt}")
    for p in use:
        efs = f"{p.efermi:.4f}" if p.efermi is not None else "n/a"
        print(f"  {os.path.basename(p.path):40s} kind={p.kind:8s} spin={p.spin} "
              f"E_F={efs} eV channels={','.join(p.channels)}")

    # series ---------------------------------------------------------------
    series = {}   # (name, spin) -> density on grid

    def add(name, spin, dens):
        key = (name, spin)
        series[key] = series.get(key, 0.0) + dens * args.scale

    for p in use:
        el = element_of(p.kind, args.by_kind)
        if is_curve:
            chans = dict(p.channels)
            tot = p.total
        else:
            e_rel = p.energy - ref
            win = (e_rel > grid[0] - 10 * fwhm) & (e_rel < grid[-1] + 10 * fwhm)
            chans = {lab: broaden(np, grid, e_rel[win], w[win], fwhm, shape) for lab, w in p.channels.items()}
            tot = sum(chans.values()) if chans else np.zeros_like(grid)
        for lab, dens in chans.items():
            if args.group == "orbital":
                add(f"{el} {lab}", p.spin, dens)
            elif args.group == "l":
                add(l_of(lab), p.spin, dens)
            elif args.group == "element-l":
                add(f"{el} {l_of(lab)}", p.spin, dens)
        if args.group == "element":
            add(el, p.spin, tot)
        if args.group == "total" or args.total:
            add("total", p.spin, tot)

    if args.spin == "sum" and polarized:
        merged = {}
        for (name, s), d in series.items():
            merged[(name, 0)] = merged.get((name, 0), 0.0) + d
        series = merged
    elif args.spin in ("up", "down") and polarized:
        keep = 1 if args.spin == "up" else 2
        series = {k: v for k, v in series.items() if k[1] == keep}

    def sort_key(k):
        name, s = k
        parts = name.split()
        lpos = L_ORDER.index(parts[-1]) if parts[-1] in L_ORDER else 9
        return (name == "total", parts[0] if len(parts) > 1 else "", lpos, name, s)

    keys = sorted(series, key=sort_key)
    names = list(dict.fromkeys(k[0] for k in keys))

    sel = (grid >= args.erange[0]) & (grid <= args.erange[1])
    x = grid[sel]

    # CSV --------------------------------------------------------------------
    if args.csv:
        cols, heads = [x], ["E_eV"]
        for k in keys:
            cols.append(series[k][sel])
            heads.append(k[0].replace(" ", "_") + ("" if k[1] == 0 else ("_up" if k[1] == 1 else "_down")))
        np.savetxt(args.csv, np.column_stack(cols), delimiter=",", header=",".join(heads), comments="", fmt="%.6f")
        print(f"wrote {args.csv} ({len(heads) - 1} series, states/eV)")

    if args.output.lower() == "none" and not args.show:
        return

    try:
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        die("matplotlib is required for plotting (pip install matplotlib); use -o none --csv FILE otherwise")

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for i, name in enumerate(names):
        color = "#4a4a48" if name == "total" else COLORS[i % len(COLORS)]
        style = "-" if name == "total" else STYLES[(i // len(COLORS)) % len(STYLES)]
        for s in (0, 1, 2):
            if (name, s) not in series:
                continue
            y = series[(name, s)][sel]
            if s == 2 and args.spin == "both":
                y = -y
            lab = name if s != 2 or (name, 1) not in series else None
            ax.plot(x, y, style, color=color, lw=1.5, label=lab)
            if args.fill:
                ax.fill_between(x, 0, y, color=color, alpha=0.15, lw=0)
    if polarized and args.spin == "both":
        ax.axhline(0.0, color="#8a8a86", lw=0.8)
        ax.text(0.01, 0.97, "spin up", transform=ax.transAxes, va="top", fontsize=8, color="#4a4a48")
        ax.text(0.01, 0.03, "spin down", transform=ax.transAxes, va="bottom", fontsize=8, color="#4a4a48")
    absolute = zero == "none" or (zero == "file" and use[0].zero_label == "ABSOLUTE")
    if not absolute or args.efermi is not None:
        ax.axvline(0.0, color="#8a8a86", lw=0.8, ls="--")
    ax.set_xlim(args.erange)
    if args.ymax:
        if polarized and args.spin == "both":
            ax.set_ylim(-args.ymax, args.ymax)
        else:
            ax.set_ylim(0, args.ymax)
    elif not (polarized and args.spin == "both"):
        ax.set_ylim(bottom=0)
    zl = {"fermi": r"E - E$_\mathrm{F}$", "hoco": r"E - E$_\mathrm{HOCO}$", "none": "E",
          "file": "E (CP2K zero: %s)" % use[0].zero_label}.get(zero, "E")
    if args.efermi is not None:
        zl = "E - E$_\\mathrm{ref}$"
    ax.set_xlabel(zl + " (eV)")
    ax.set_ylabel("PDOS (states/eV)")
    ax.grid(alpha=0.25, lw=0.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if len(names) > 1:
        ax.legend(frameon=False, fontsize=8, ncol=2 if len(names) > 6 else 1)
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
