#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_eos_fit.py -- equation of state for OpenMX: scale the lattice, run, fit Birch-Murnaghan.
                  按晶格常数缩放生成一组 OpenMX 输入（固定 scf.Ngrid，避免实空间网格数突变导致
                  能量曲线跳变），依次运行，拟合三阶 Birch-Murnaghan（或 Murnaghan / Vinet）状态
                  方程，得到 V0、a0、B0、B0'；也可直接拟合 MD.Type EvsLC 生成的 .EvsLC 文件。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Workflow
--------
1. The template (a converged single-point input with <Atoms.UnitVectors) is
   scaled homogeneously: every unit vector and, for Ang/AU coordinates, every
   Cartesian coordinate is multiplied by s (FRAC coordinates are unchanged).
   Scales: --scales 0.97 0.98 ... or --range -3 3 7 (linear strain in %).
2. The real-space grid must be identical for all volumes (manual Sec. 14.2):
     * scf.Ngrid of the template is kept if present, or
     * --ngrid N1 N2 N3, or
     * --ngrid-from FILE.out ("Num. of grids of a-, b-, and c-axes = ..."), or
     * automatically: the LARGEST cell is run once with scf.maxIter 1 and its
       grid is taken, so the requested cutoff is met (or exceeded) everywhere.
3. Each input <prefix>_s<scale>.dat is run as  <mpirun> <openmx> in.dat -nt <nt>.
4. E(V) is fitted; B0 is converted with 1 eV/Ang^3 = 160.21766 GPa.
   The 3rd-order Birch-Murnaghan fit is a cubic polynomial in V^(-2/3) (exact,
   numpy only); Murnaghan and Vinet use scipy.optimize.curve_fit.

Usage examples
--------------
  python omx_eos_fit.py Si_scf.dat --range -3 3 7 --a-ref 5.43 --mpirun "mpirun -np 8"
  python omx_eos_fit.py Fe_bcc.dat --scales 0.96 0.98 1.00 1.02 1.04 --ngrid 32 32 32
  python omx_eos_fit.py Si_scf.dat --range -3 3 7 --dry-run --ngrid 24 24 24   # inputs only
  python omx_eos_fit.py Si_scf.dat --range -3 3 7 --collect --a-ref 5.43       # fit existing runs
  python omx_eos_fit.py --evslc cdia_evslc.EvsLC --natoms 2 --a-ref-factor 4    # fit an EvsLC file

--a-ref A gives the lattice constant (Ang) that corresponds to scale 1.0, so that
a0 = A * s0 with s0 = (V0/V_ref)^(1/3).  For an EvsLC file use --a-ref-factor F
to get a0 = (F * V0)^(1/3) (F = 4 for an fcc primitive cell, 2 for bcc, 1 for a
conventional cubic cell).

Environment: MPIRUN (default "mpirun -np 4"), OPENMX (default "openmx").
Dependencies: python >= 3.9, numpy, matplotlib (plot), scipy (Murnaghan/Vinet only).
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shlex
import subprocess
import sys

HA2EV = 27.211386245988
BOHR2ANG = 0.529177210903
EVA3_TO_GPA = 160.21766208
FINISHED_MSG = "The calculation was normally finished."
HERE = os.path.dirname(os.path.abspath(__file__))


def _np():
    try:
        import numpy as np
    except ImportError:
        sys.exit("ERROR: numpy is required (pip install numpy)")
    return np


# ---------------------------------------------------------------------------
# input editing
# ---------------------------------------------------------------------------
def kw_get(lines, key):
    low = key.lower()
    for line in lines:
        tok = line.split("#", 1)[0].split()
        if tok and tok[0].lower() == low:
            return tok[1:]
    return None


def kw_set(lines, key, value):
    low = key.lower()
    for i, line in enumerate(lines):
        tok = line.split("#", 1)[0].split()
        if tok and tok[0].lower() == low:
            lines[i] = "%-27s %s\n" % (key, value)
            return
    lines.append("%-27s %s\n" % (key, value))


def block_range(lines, name):
    low = name.lower()
    start = None
    for i, line in enumerate(lines):
        s = line.split("#", 1)[0].strip().lower()
        if start is None and s.startswith("<") and s[1:].split()[:1] == [low]:
            start = i
        elif start is not None and s.endswith(low + ">"):
            return start, i
    return None


def cell_of(lines):
    rng = block_range(lines, "Atoms.UnitVectors")
    if rng is None:
        sys.exit("ERROR: <Atoms.UnitVectors block not found; the EOS needs a periodic cell")
    unit = (kw_get(lines, "Atoms.UnitVectors.Unit") or ["Ang"])[0].lower()
    rows = []
    for line in lines[rng[0] + 1: rng[1]]:
        tok = line.split("#", 1)[0].split()
        if len(tok) >= 3:
            rows.append([float(x) for x in tok[:3]])
    if len(rows) != 3:
        sys.exit("ERROR: Atoms.UnitVectors must contain three vectors")
    fac = BOHR2ANG if unit == "au" else 1.0
    return [[x * fac for x in r] for r in rows]


def scaled_lines(lines, s):
    out = list(lines)
    rng = block_range(out, "Atoms.UnitVectors")
    for i in range(rng[0] + 1, rng[1]):
        tok = out[i].split("#", 1)[0].split()
        if len(tok) >= 3:
            out[i] = "  %18.12f %18.12f %18.12f\n" % tuple(float(x) * s for x in tok[:3])
    unit = (kw_get(out, "Atoms.SpeciesAndCoordinates.Unit") or ["Ang"])[0].lower()
    if unit in ("ang", "au"):
        rng = block_range(out, "Atoms.SpeciesAndCoordinates")
        if rng is None:
            sys.exit("ERROR: <Atoms.SpeciesAndCoordinates block not found")
        for i in range(rng[0] + 1, rng[1]):
            tok = out[i].split("#", 1)[0].split()
            if len(tok) >= 5:
                xyz = ["%16.10f" % (float(x) * s) for x in tok[2:5]]
                out[i] = "  " + "  ".join(tok[:2] + xyz + tok[5:]) + "\n"
    return out


# ---------------------------------------------------------------------------
# EOS fits (E in eV, V in Ang^3)
# ---------------------------------------------------------------------------
def fit_bm3_poly(V, E):
    """3rd-order Birch-Murnaghan as a cubic polynomial in x = V^(-2/3)."""
    np = _np()
    x = V ** (-2.0 / 3.0)
    c = np.polyfit(x, E, 3)
    p = np.poly1d(c)
    d1, d2, d3 = p.deriv(1), p.deriv(2), p.deriv(3)
    roots = [r.real for r in d1.r if abs(r.imag) < 1e-12 and r.real > 0 and d2(r.real) > 0]
    if not roots:
        raise RuntimeError("the fitted curve has no minimum; widen or recenter the volume range")
    x0 = min(roots, key=lambda r: abs(r - x.mean()))
    V0 = x0 ** -1.5
    g1 = -2.0 / 3.0 * V0 ** (-5.0 / 3.0)
    g2 = 10.0 / 9.0 * V0 ** (-8.0 / 3.0)
    g3 = -80.0 / 27.0 * V0 ** (-11.0 / 3.0)
    Ex, Exx, Exxx = d1(x0), d2(x0), d3(x0)
    EVV = Exx * g1 ** 2 + Ex * g2
    EVVV = Exxx * g1 ** 3 + 3.0 * Exx * g1 * g2 + Ex * g3
    B0 = V0 * EVV
    Bp = -(1.0 + V0 * EVVV / EVV)
    return {"E0": float(p(x0)), "V0": float(V0), "B0": float(B0), "Bp": float(Bp),
            "model": lambda v: p(v ** (-2.0 / 3.0))}


def e_murnaghan(V, E0, V0, B0, Bp):
    return E0 + B0 * V / Bp * ((V0 / V) ** Bp / (Bp - 1.0) + 1.0) - B0 * V0 / (Bp - 1.0)


def e_vinet(V, E0, V0, B0, Bp):
    np = _np()
    eta = (V / V0) ** (1.0 / 3.0)
    return E0 + 2.0 * B0 * V0 / (Bp - 1.0) ** 2 * (
        2.0 - (5.0 + 3.0 * Bp * (eta - 1.0) - 3.0 * eta) * np.exp(-1.5 * (Bp - 1.0) * (eta - 1.0)))


def fit_eos(V, E, kind):
    start = fit_bm3_poly(V, E)
    if kind == "bm3":
        return start
    try:
        from scipy.optimize import curve_fit
    except ImportError:
        sys.exit("ERROR: --eos %s needs scipy (pip install scipy); --eos bm3 works with numpy only" % kind)
    f = {"murnaghan": e_murnaghan, "vinet": e_vinet}[kind]
    p0 = [start["E0"], start["V0"], start["B0"], start["Bp"] if 1.5 < start["Bp"] < 10 else 4.0]
    popt, _ = curve_fit(f, V, E, p0=p0, maxfev=20000)
    E0, V0, B0, Bp = (float(x) for x in popt)
    return {"E0": E0, "V0": V0, "B0": B0, "Bp": Bp, "model": lambda v: f(v, E0, V0, B0, Bp)}


# ---------------------------------------------------------------------------
# running / parsing
# ---------------------------------------------------------------------------
def run(workdir, name, args):
    cmd = shlex.split(args.mpirun) + [args.openmx, name + ".dat", "-nt", str(args.nt)]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(args.nt)
    print(">>> (%s) %s" % (workdir, " ".join(cmd)), flush=True)
    try:
        with open(os.path.join(workdir, name + ".std"), "w") as fh:
            return subprocess.run(cmd, cwd=workdir, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode
    except FileNotFoundError as exc:
        sys.exit("ERROR: cannot start '%s': %s" % (cmd[0], exc))


def is_finished(workdir, name):
    p = os.path.join(workdir, name + ".std")
    if not os.path.isfile(p):
        return False
    with open(p, errors="replace") as fh:
        return FINISHED_MSG in fh.read()


def read_out(path):
    """(Utot in Ha or None, SCF converged or None, grid list or None)."""
    if not os.path.isfile(path):
        return None, None, None
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    try:
        import omx_parse_out  # type: ignore
        r = omx_parse_out.parse_out(path)
        return r.get("utot_Ha"), r.get("scf_converged"), r.get("ngrid")
    except ImportError:
        pass
    except Exception as exc:
        print("WARNING: could not parse %s: %s" % (path, exc), file=sys.stderr)
        return None, None, None
    with open(path, errors="replace") as fh:
        text = fh.read()
    m = re.findall(r"^\s*Utot\.\s+([-+0-9.eE]+)\s*$", text, re.M)
    g = re.search(r"Num\. of grids of a-, b-, and c-axes\s*=\s*(\d+),\s*(\d+),\s*(\d+)", text)
    return (float(m[-1]) if m else None), None, ([int(x) for x in g.groups()] if g else None)


def read_evslc(path):
    """Rows 'ax ay az bx by bz cx cy cz Etot' (Ang, Hartree) of System.Name.EvsLC."""
    np = _np()
    if not os.path.isfile(path):
        sys.exit("ERROR: file not found: %s" % path)
    V, E = [], []
    with open(path, errors="replace") as fh:
        for line in fh:
            tok = line.split("#", 1)[0].split()
            if len(tok) < 10:
                continue
            try:
                v = [float(t) for t in tok[:10]]
            except ValueError:
                continue
            V.append(abs(np.linalg.det(np.array(v[:9]).reshape(3, 3))))
            E.append(v[9] * HA2EV)
    if len(V) < 4:
        sys.exit("ERROR: %s contains fewer than 4 data lines (ax ... cz Etot)" % path)
    return np.array(V), np.array(E)


# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="E(V) equation of state with OpenMX: scaled inputs with a fixed scf.Ngrid, runs, "
                    "and a Birch-Murnaghan (or Murnaghan/Vinet) fit.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python omx_eos_fit.py Si_scf.dat --range -3 3 7 --a-ref 5.43\n"
               "  python omx_eos_fit.py Fe_bcc.dat --scales 0.96 0.98 1.00 1.02 1.04 --ngrid 32 32 32\n"
               "  python omx_eos_fit.py Si_scf.dat --range -3 3 7 --collect --a-ref 5.43\n"
               "  python omx_eos_fit.py --evslc cdia_evslc.EvsLC --natoms 2 --a-ref-factor 4")
    p.add_argument("template", nargs="?", help="converged single-point OpenMX input (.dat)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--scales", nargs="+", type=float, help="lattice scale factors, e.g. 0.97 0.98 1.00 1.02")
    g.add_argument("--range", nargs=3, type=float, metavar=("MIN%", "MAX%", "N"),
                   help="linear strain range in percent and number of points, e.g. -3 3 7")
    p.add_argument("--evslc", default=None, help="fit an existing System.Name.EvsLC file instead (no runs)")
    p.add_argument("--ngrid", nargs=3, type=int, default=None, metavar=("N1", "N2", "N3"),
                   help="fixed scf.Ngrid for all volumes")
    p.add_argument("--ngrid-from", default=None, metavar="OUT", help="take the grid from this .out file")
    p.add_argument("--workdir", default="eos", help="directory for inputs/outputs (default eos)")
    p.add_argument("--prefix", default=None, help="System.Name prefix (default: System.Name of the template)")
    p.add_argument("--mpirun", default=os.environ.get("MPIRUN", "mpirun -np 4"),
                   help="MPI launcher with options (default $MPIRUN or 'mpirun -np 4'); '' for serial")
    p.add_argument("--openmx", default=os.environ.get("OPENMX", "openmx"),
                   help="OpenMX executable (default $OPENMX or 'openmx')")
    p.add_argument("--nt", type=int, default=1, help="OpenMP threads per MPI process (-nt), default 1")
    p.add_argument("--data-path", default=None, help="override DATA.PATH in the generated inputs")
    p.add_argument("--keep-fileout", action="store_true", help="keep level.of.fileout (default: set 0)")
    p.add_argument("--dry-run", action="store_true", help="only write the inputs (needs a known grid)")
    p.add_argument("--collect", action="store_true", help="only fit existing outputs, no runs")
    p.add_argument("--rerun", action="store_true", help="repeat runs that already finished")
    p.add_argument("--eos", choices=["bm3", "murnaghan", "vinet"], default="bm3",
                   help="EOS form (default bm3 = 3rd-order Birch-Murnaghan)")
    p.add_argument("--a-ref", type=float, default=None,
                   help="lattice constant (Ang) corresponding to scale 1.0; reports a0 = a_ref * s0")
    p.add_argument("--a-ref-factor", type=float, default=None,
                   help="EvsLC mode: a0 = (factor * V0)^(1/3), e.g. 4 for an fcc primitive cell")
    p.add_argument("--natoms", type=int, default=None, help="atoms per cell (default: Atoms.Number)")
    p.add_argument("--csv", default=None, help="CSV file (default <workdir>/eos.csv)")
    p.add_argument("--plot", default=None, help="PNG file (default <workdir>/eos.png)")
    p.add_argument("--no-plot", action="store_true", help="no plot")
    p.add_argument("--show", action="store_true", help="show the plot window")
    return p


def report_and_plot(V, E, labels, args, natoms, ref_volume, title):
    np = _np()
    order = np.argsort(V)
    V, E = V[order], E[order]
    labels = [labels[i] for i in order]
    if len(V) < 4:
        sys.exit("ERROR: need at least 4 E(V) points, have %d" % len(V))
    try:
        fit = fit_eos(V, E, args.eos)
    except (RuntimeError, ValueError) as exc:
        sys.exit("ERROR: EOS fit failed: %s" % exc)
    except Exception as exc:          # scipy convergence problems
        sys.exit("ERROR: EOS fit failed: %s" % exc)
    resid = E - fit["model"](V)
    rms = float(np.sqrt(np.mean(resid ** 2))) * 1000.0
    print("\nEOS fit (%s), %d points" % (args.eos, len(V)))
    print("  E0  = %.6f eV/cell" % fit["E0"])
    print("  V0  = %.4f Ang^3/cell%s" % (fit["V0"], ("  (%.4f Ang^3/atom)" % (fit["V0"] / natoms)) if natoms else ""))
    print("  B0  = %.2f GPa" % (fit["B0"] * EVA3_TO_GPA))
    print("  B0' = %.3f" % fit["Bp"])
    if ref_volume:
        s0 = (fit["V0"] / ref_volume) ** (1.0 / 3.0)
        print("  s0  = %.5f (lattice scale at V0)" % s0)
        if args.a_ref:
            print("  a0  = %.4f Ang (a_ref %.4f x s0)" % (args.a_ref * s0, args.a_ref))
    if args.a_ref_factor:
        print("  a0  = %.4f Ang ((%.3g x V0)^(1/3))" % ((args.a_ref_factor * fit["V0"]) ** (1.0 / 3.0), args.a_ref_factor))
    print("  RMS residual = %.3f meV/cell" % rms)
    if not (V.min() < fit["V0"] < V.max()):
        print("WARNING: V0 lies outside the sampled volumes; shift/widen the scale range.")
    elif fit["V0"] < V[1] or fit["V0"] > V[-2]:
        print("WARNING: V0 is close to the edge of the sampled range; add points on that side.")
    if rms > 1.0:
        print("WARNING: RMS residual > 1 meV/cell; check for grid jumps (scf.Ngrid), unconverged SCF "
              "or a too wide strain range.")

    csv_path = args.csv or os.path.join(args.workdir, "eos.csv")
    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["label", "volume_Ang3", "energy_eV", "fit_eV", "residual_meV"])
        for lab, v, e, r in zip(labels, V, E, resid):
            w.writerow([lab, "%.6f" % v, "%.8f" % e, "%.8f" % (e - r), "%.4f" % (r * 1000)])
        w.writerow([])
        w.writerow(["# eos", args.eos, "E0_eV", fit["E0"], "V0_Ang3", fit["V0"],
                    "B0_GPa", fit["B0"] * EVA3_TO_GPA, "B0p", fit["Bp"]])
    print("wrote %s" % csv_path)

    if args.no_plot:
        return
    try:
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not installed; no plot", file=sys.stderr)
        return
    vv = np.linspace(V.min() - 0.02 * (V.max() - V.min()), V.max() + 0.02 * (V.max() - V.min()), 300)
    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.plot(V, E - fit["E0"], "o", color="#1f77b4", label="OpenMX")
    ax.plot(vv, fit["model"](vv) - fit["E0"], "-", color="#d62728",
            label="%s: V0=%.3f, B0=%.1f GPa, B0'=%.2f" % (args.eos, fit["V0"], fit["B0"] * EVA3_TO_GPA, fit["Bp"]))
    ax.axvline(fit["V0"], color="0.6", ls=":", lw=1)
    ax.set_xlabel("Volume (Ang$^3$/cell)")
    ax.set_ylabel("E - E0 (eV/cell)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    png = args.plot or os.path.join(args.workdir, "eos.png")
    fig.savefig(png, dpi=150)
    print("wrote %s" % png)
    if args.show:
        plt.show()


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    np = _np()

    if args.evslc:
        V, E = read_evslc(args.evslc)
        report_and_plot(V, E, ["pt%d" % (i + 1) for i in range(len(V))], args, args.natoms, None,
                        os.path.basename(args.evslc))
        return 0

    if not args.template:
        sys.exit("ERROR: give a template .dat (or --evslc FILE)")
    if not os.path.isfile(args.template):
        sys.exit("ERROR: template not found: %s" % args.template)
    with open(args.template, encoding="utf-8-sig", errors="replace") as fh:
        tlines = fh.read().replace("\r\n", "\n").splitlines(True)
    if tlines and not tlines[-1].endswith("\n"):
        tlines[-1] += "\n"

    mdtype = (kw_get(tlines, "MD.Type") or ["nomd"])[0].lower()
    if mdtype.startswith("optc") or mdtype.startswith("rfc") or mdtype in ("evslc", "df"):
        sys.exit("ERROR: MD.Type %s changes/scans the cell; use nomd (or an atomic relaxation such as RF)" % mdtype)
    if mdtype not in ("nomd",):
        print("NOTE: MD.Type %s: atoms are relaxed at every volume (fine for internal coordinates); "
              "the last Utot of each run is used" % mdtype)

    if args.scales:
        scales = sorted(set(args.scales))
    elif args.range:
        lo, hi, n = args.range
        n = int(n)
        if n < 4:
            sys.exit("ERROR: --range needs N >= 4")
        scales = [1.0 + (lo + (hi - lo) * i / (n - 1)) / 100.0 for i in range(n)]
    else:
        sys.exit("ERROR: give --scales or --range")
    if len(scales) < 4:
        sys.exit("ERROR: need at least 4 scale factors")

    cell = np.array(cell_of(tlines))
    vref = abs(np.linalg.det(cell))
    natoms = args.natoms
    if natoms is None:
        v = kw_get(tlines, "Atoms.Number")
        natoms = int(float(v[0])) if v else None
    prefix = args.prefix or (kw_get(tlines, "System.Name") or ["omx"])[0]
    os.makedirs(args.workdir, exist_ok=True)

    dp = args.data_path or ((kw_get(tlines, "DATA.PATH") or ["../DFT_DATA19"])[0])
    if not os.path.isabs(dp):
        dp = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.template)), dp))

    def common(lines, name):
        kw_set(lines, "System.Name", name)
        kw_set(lines, "DATA.PATH", dp)
        if not args.keep_fileout:
            kw_set(lines, "level.of.fileout", "0")
        return lines

    # ---- fixed grid
    ngrid = args.ngrid
    if ngrid is None and args.ngrid_from:
        _, _, ngrid = read_out(args.ngrid_from)
        if not ngrid:
            sys.exit("ERROR: no 'Num. of grids' line in %s" % args.ngrid_from)
    if ngrid is None and kw_get(tlines, "scf.Ngrid"):
        ngrid = [int(x) for x in kw_get(tlines, "scf.Ngrid")[:3]]
    if ngrid is None:
        probe = "%s_gridprobe" % prefix
        if args.dry_run or args.collect:
            pout = os.path.join(args.workdir, probe + ".out")
            if os.path.isfile(pout):
                _, _, ngrid = read_out(pout)
            if not ngrid:
                sys.exit("ERROR: the grid is unknown; give --ngrid N1 N2 N3 or --ngrid-from FILE.out")
        else:
            smax = max(scales)
            lines = common(scaled_lines(tlines, smax), probe)
            kw_set(lines, "scf.maxIter", "1")
            kw_set(lines, "MD.Type", "nomd")
            with open(os.path.join(args.workdir, probe + ".dat"), "w", encoding="utf-8", newline="\n") as fh:
                fh.write("".join(lines))
            print("grid probe: largest cell (scale %.4f), scf.maxIter 1" % smax)
            run(args.workdir, probe, args)
            _, _, ngrid = read_out(os.path.join(args.workdir, probe + ".out"))
            if not ngrid:
                sys.exit("ERROR: grid probe failed; see %s/%s.std" % (args.workdir, probe))
    print("fixed scf.Ngrid = %d %d %d" % tuple(ngrid))

    # ---- inputs
    names = []
    for s in scales:
        name = "%s_s%.4f" % (prefix, s)
        lines = common(scaled_lines(tlines, s), name)
        kw_set(lines, "scf.Ngrid", "%d %d %d" % tuple(ngrid))
        with open(os.path.join(args.workdir, name + ".dat"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("".join(lines))
        names.append(name)
    print("wrote %d inputs to %s/" % (len(names), args.workdir))
    if args.dry_run:
        return 0

    if not args.collect:
        for name in names:
            if is_finished(args.workdir, name) and not args.rerun:
                print("skip %s (finished)" % name)
                continue
            run(args.workdir, name, args)

    V, E, labels = [], [], []
    print("\n%10s %14s %20s %8s" % ("scale", "V (Ang^3)", "Utot (Ha)", "SCF"))
    for s, name in zip(scales, names):
        utot, conv, grid = read_out(os.path.join(args.workdir, name + ".out"))
        vol = vref * s ** 3
        print("%10.4f %14.4f %20s %8s" % (s, vol, ("%.10f" % utot) if utot is not None else "missing",
                                          {True: "ok", False: "NOCONV", None: "?"}[conv]))
        if grid and list(grid) != list(ngrid):
            print("WARNING: %s used grid %s instead of %s" % (name, grid, ngrid))
        if utot is None:
            continue
        if conv is False:
            print("WARNING: %s did not converge; point excluded" % name)
            continue
        V.append(vol)
        E.append(utot * HA2EV)
        labels.append("%.4f" % s)
    report_and_plot(np.array(V), np.array(E), labels, args, natoms, vref,
                    "%s (fixed grid %s)" % (prefix, "x".join(map(str, ngrid))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
