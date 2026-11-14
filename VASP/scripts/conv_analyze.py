#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conv_analyze.py -- analyse VASP convergence tests (ENCUT / k-points / KSPACING /
                   SIGMA / any INCAR tag) and plot the convergence curves.
                   解析收敛测试目录（OUTCAR/OSZICAR），输出 dE/atom、最大受力、T*S 表并作图。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Input
-----
Directories created by conv_test.sh (each holding conv_info.txt), or any set of
VASP run directories; the tested parameter is taken from conv_info.txt, from
the directory name (``encut_400``, ``kpoints_8x8x8``, ``sigma_0.05`` ...) or
from OUTCAR (``--param ENCUT``).

For every run: energy(sigma->0) (or TOTEN with ``--energy toten``), energy per
atom, deviation from the reference run (default: the most accurate one =
largest ENCUT / densest k-mesh / smallest KSPACING / smallest SIGMA) in
meV/atom, max force and max force deviation, external pressure, T*S per atom,
number of irreducible k-points and elapsed time.  The first value from which on
|dE| stays below ``--tol`` (default 1 meV/atom) is reported as converged; for
SIGMA tests the run must in addition have |T*S| below ``--tstol`` (default
1 meV/atom).

Reference for k-point tests: the x value of a k-mesh run is the PRODUCT
N1*N2*N3 of the mesh.  ``--ref`` accepts a mesh ("8x8x8"), the product ("512")
or a single N ("8" -> the run whose mesh is 8x8x8, or 8x8x1 for slabs, i.e. all
divisions > 1 equal to N); a single N is tried first, then the product.

Usage examples
--------------
    python conv_analyze.py conv_encut
    python conv_analyze.py conv_kpoints --tol 0.5 --ftol 0.01
    python conv_analyze.py encut_* --param ENCUT -o encut_test.png
    python conv_analyze.py conv_sigma --energy toten --tstol 0.5
    python conv_analyze.py conv_kpoints --ref 8x8x8   # same as --ref 8 or --ref 512
    python conv_analyze.py conv_encut --show       # interactive window

Outputs: a table on screen, <prefix>.dat (whitespace table) and <prefix>.png.

Dependencies: python >= 3.9, numpy, matplotlib.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

__version__ = "1.1.0"

FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"
KNOWN = {"encut": "ENCUT", "kpoints": "KPOINTS", "kspacing": "KSPACING", "sigma": "SIGMA"}
# light-mode categorical slots + text/grid colours (validated reference palette)
C1, C2 = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def fnum(s):
    try:
        return float(str(s).replace("D", "E"))
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------
def read_conv_info(d: Path) -> dict:
    f = d / "conv_info.txt"
    out = {}
    if f.is_file():
        for line in f.read_text(errors="replace").splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                out[k.strip().lower()] = v.strip()
    return out


def read_poscar_flags(path: Path):
    """(lattice, flags) of a POSCAR with selective dynamics, else None."""
    try:
        lines = path.read_text(errors="replace").splitlines()
        s = [float(x) for x in lines[1].split()[:3]]
        lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
        lat = lat * (s[0] if len(s) == 1 and s[0] > 0 else (np.array(s) if len(s) == 3 else 1.0))
        i = 5 if all(t.isdigit() for t in lines[5].split()) else 6
        nat = sum(int(x) for x in lines[i].split())
        i += 1
        if lines[i].strip()[:1] not in ("s", "S"):
            return None
        i += 2
        flags = []
        for j in range(nat):
            t = lines[i + j].split()
            flags.append([x.upper().startswith("T") for x in t[3:6]] if len(t) >= 6 else [True] * 3)
        return lat, np.array(flags)
    except (OSError, IndexError, ValueError):
        return None


def project(forces, lat, flags):
    out = np.zeros_like(forces)
    for k, (f, fl) in enumerate(zip(forces, flags)):
        free = lat[fl]
        if len(free) == 3:
            out[k] = f
        elif len(free):
            q, _ = np.linalg.qr(free.T)
            out[k] = q @ (q.T @ f)
    return out


def parse_outcar(path: Path) -> dict:
    r = {"e0": None, "toten": None, "ts": None, "nions": None, "nkpts": None, "forces": None,
         "pressure": None, "encut": None, "sigma": None, "ismear": None, "kspacing": None,
         "elapsed": None, "finished": False}
    nions = None
    with open(path, errors="replace") as fh:
        it = iter(fh)
        for line in it:
            if nions is None and "NIONS" in line:
                m = re.search(r"NIONS\s*=\s*(\d+)", line)
                if m:
                    nions = r["nions"] = int(m.group(1))
            if r["nkpts"] is None and "NKPTS" in line:
                m = re.search(r"NKPTS\s*=\s*(\d+)", line)
                if m:
                    r["nkpts"] = int(m.group(1))
            if r["encut"] is None and re.search(r"^\s*ENCUT\s*=", line):
                r["encut"] = fnum(re.search(r"ENCUT\s*=\s*(" + FLOAT + ")", line).group(1))
            if r["ismear"] is None and "ISMEAR" in line and "SIGMA" in line:
                m1 = re.search(r"ISMEAR\s*=\s*(-?\d+)", line)
                m2 = re.search(r"SIGMA\s*=\s*(" + FLOAT + ")", line)
                if m1 and m2:
                    r["ismear"], r["sigma"] = int(m1.group(1)), fnum(m2.group(1))
            if r["kspacing"] is None and re.search(r"^\s*KSPACING\s*=", line):
                m = re.search(r"KSPACING\s*=\s*(" + FLOAT + ")", line)
                if m:
                    r["kspacing"] = fnum(m.group(1))
            if "EENTRO" in line:
                m = re.search(r"EENTRO\s*=\s*(" + FLOAT + ")", line)
                if m:
                    r["ts"] = fnum(m.group(1))
            elif "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in line:
                for _ in range(6):
                    nxt = next(it, "")
                    m = re.search(r"TOTEN\s*=\s*(" + FLOAT + ")", nxt)
                    if m:
                        r["toten"] = fnum(m.group(1))
                    m = re.search(r"energy\(sigma->0\)\s*=\s*(" + FLOAT + ")", nxt)
                    if m:
                        r["e0"] = fnum(m.group(1))
                        break
            elif "TOTAL-FORCE" in line and "POSITION" in line and nions:
                next(it, "")
                rows = []
                for _ in range(nions):
                    t = next(it, "").split()
                    if len(t) >= 6:
                        rows.append([fnum(x) for x in t[3:6]])
                if len(rows) == nions and all(None not in x for x in rows):
                    r["forces"] = np.array(rows)
            elif "external pressure" in line:
                m = re.search(r"external pressure\s*=\s*(" + FLOAT + ")", line)
                if m:
                    r["pressure"] = fnum(m.group(1))
            elif "Elapsed time (sec)" in line:
                r["elapsed"] = fnum(line.split(":")[-1])
            elif "General timing and accounting" in line:
                r["finished"] = True
    return r


def parse_oszicar_e0(path: Path):
    e0 = None
    with open(path, errors="replace") as fh:
        for line in fh:
            m = re.search(r"\bE0=\s*(" + FLOAT + ")", line)
            if m:
                e0 = fnum(m.group(1))
    return e0


def natoms_from_poscar(path: Path):
    try:
        lines = path.read_text(errors="replace").splitlines()
        i = 5 if all(t.isdigit() for t in lines[5].split()) else 6
        return sum(int(x) for x in lines[i].split())
    except (OSError, IndexError, ValueError):
        return None


def mesh_from_kpoints(path: Path):
    try:
        lines = path.read_text(errors="replace").splitlines()
        if int(lines[1].split()[0]) == 0 and lines[2].strip()[:1].lower() in ("g", "m"):
            return [int(x) for x in lines[3].split()[:3]]
    except (OSError, IndexError, ValueError):
        pass
    return None


def collect(paths, args):
    """Find run directories and parse them."""
    runs = []
    cand = []
    for p in paths:
        p = Path(p)
        if p.is_file():
            p = p.parent
        if not p.is_dir():
            print(f"WARNING: '{p}' is not a directory, skipped", file=sys.stderr)
            continue
        subs = sorted(x for x in p.iterdir() if x.is_dir() and ((x / "OUTCAR").is_file() or (x / "conv_info.txt").is_file()))
        if (p / "OUTCAR").is_file() or (p / "conv_info.txt").is_file():
            cand.append(p)
        elif subs:
            cand += subs
        else:
            print(f"WARNING: no OUTCAR/conv_info.txt in '{p}'", file=sys.stderr)
    seen = set()
    for d in cand:
        if d.resolve() in seen:
            continue
        seen.add(d.resolve())
        info = read_conv_info(d)
        param = (args.param or info.get("param") or "").upper()
        value = info.get("value")
        if not param or value is None:
            m = re.match(r"^([A-Za-z]+)_(.+)$", d.name)
            if m:
                param = param or KNOWN.get(m.group(1).lower(), m.group(1).upper())
                value = value if value is not None else m.group(2)
        run = {"dir": d, "param": param, "label": value, "x": None}
        oc = d / "OUTCAR"
        if oc.is_file():
            run.update(parse_outcar(oc))
        else:
            run.update({"e0": None, "toten": None, "nions": None, "forces": None, "ts": None,
                        "nkpts": None, "pressure": None, "elapsed": None, "finished": False})
        if run["e0"] is None and (d / "OSZICAR").is_file():
            run["e0"] = parse_oszicar_e0(d / "OSZICAR")
        if run["nions"] is None:
            run["nions"] = natoms_from_poscar(d / "POSCAR")
        # numeric x value
        if param == "KPOINTS":
            mesh = None
            if info.get("mesh"):
                mesh = [int(x) for x in info["mesh"].split()]
            elif value and re.fullmatch(r"\d+([xX]\d+[xX]\d+)?", str(value)):
                mesh = [int(x) for x in re.split(r"[xX]", value)] if "x" in value.lower() else [int(value)] * 3
            else:
                mesh = mesh_from_kpoints(d / "KPOINTS")
            run["mesh"] = mesh
            run["x"] = float(np.prod(mesh)) if mesh else (float(run["nkpts"]) if run["nkpts"] else None)
            run["label"] = "x".join(str(v) for v in mesh) if mesh else str(value)
        elif param:
            v = fnum(value) if value is not None else None
            if v is None:
                v = {"ENCUT": run.get("encut"), "SIGMA": run.get("sigma"), "KSPACING": run.get("kspacing")}.get(param)
            run["x"] = v
            if run["label"] is None and v is not None:
                run["label"] = f"{v:g}"
        if run["forces"] is not None and not args.no_constraints:
            pf = read_poscar_flags(d / "POSCAR")
            if pf is not None and len(pf[1]) == len(run["forces"]):
                run["forces"] = project(run["forces"], pf[0], pf[1])
        runs.append(run)
    return runs


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------
def pick_reference(valid, param, text):
    """Run selected by --ref.  For KPOINTS the x value is N1*N2*N3, so a mesh
    ('8x8x8'), a single N ('8' -> 8x8x8 or 8x8x1) or the product ('512') is
    accepted; for other tags the nearest numeric value is used."""
    text = str(text).strip()
    if param == "KPOINTS":
        if re.fullmatch(r"\d+[xX*]\d+[xX*]\d+", text):
            mesh = [int(v) for v in re.split(r"[xX*]", text)]
            exact = [r for r in valid if r.get("mesh") and list(r["mesh"]) == mesh]
            if exact:
                return exact[0]
            target = float(np.prod(mesh))
        elif re.fullmatch(r"\d+", text):
            n = int(text)
            same = [r for r in valid if r.get("mesh")
                    and [m for m in r["mesh"] if m > 1]
                    and all(m == n for m in r["mesh"] if m > 1)]
            if same:
                return max(same, key=lambda r: r["x"])
            target = float(n)          # not a uniform N: interpret as the product N1*N2*N3
        else:
            raise SystemExit(f"ERROR: --ref '{text}' is not a k-mesh (e.g. 8x8x8), N or N1*N2*N3")
        hit = [r for r in valid if abs(r["x"] - target) < 1e-6]
        if not hit:
            print(f"WARNING: no k-mesh with N1*N2*N3 = {target:g}; using the nearest one", file=sys.stderr)
        return hit[0] if hit else min(valid, key=lambda r: abs(r["x"] - target))
    v = fnum(text)
    if v is None:
        raise SystemExit(f"ERROR: --ref '{text}' is not a number")
    return min(valid, key=lambda r: abs(r["x"] - v))


def analyse(runs, args):
    params = {r["param"] for r in runs if r["param"]}
    if len(params) > 1:
        print(f"WARNING: mixed parameters {sorted(params)}; use --param or analyse them separately", file=sys.stderr)
    param = sorted(params)[0] if params else "VALUE"
    good = [r for r in runs if r["x"] is not None]
    numeric = bool(good)
    if not good:
        # text values (e.g. LREAL = Auto / .FALSE.): keep the given order
        good = [r for r in runs if r["label"] is not None]
        for i, r in enumerate(good):
            r["x"] = float(i)
    else:
        for r in runs:
            if r["x"] is None:
                print(f"WARNING: cannot determine the {param} value of {r['dir']} - skipped", file=sys.stderr)
        good.sort(key=lambda r: r["x"])
    if not good:
        raise SystemExit("ERROR: could not determine the tested values (use --param)")
    key = "toten" if args.energy == "toten" else "e0"
    for r in good:
        r["E"] = r.get(key)
        r["Ea"] = r["E"] / r["nions"] if (r["E"] is not None and r["nions"]) else None
        r["fmax"] = float(np.max(np.linalg.norm(r["forces"], axis=1))) if r["forces"] is not None else None
        r["ts_a"] = 1000.0 * r["ts"] / r["nions"] if (r["ts"] is not None and r["nions"]) else None
    # reference = most accurate
    valid = [r for r in good if r["Ea"] is not None]
    if not valid:
        raise SystemExit("ERROR: no finished run with an energy found (check OUTCAR files)")
    if args.ref is not None:
        ref = pick_reference(valid, param, args.ref)
    elif param in ("KSPACING", "SIGMA"):
        ref = valid[0]          # smallest spacing / smallest smearing
    else:
        ref = valid[-1]         # largest ENCUT / densest mesh / last value
    for r in good:
        r["dE"] = 1000.0 * (r["Ea"] - ref["Ea"]) if r["Ea"] is not None else None
        if r["forces"] is not None and ref["forces"] is not None and r["forces"].shape == ref["forces"].shape:
            r["dF"] = float(np.max(np.abs(r["forces"] - ref["forces"])))
        else:
            r["dF"] = None
        r["dP"] = (r["pressure"] - ref["pressure"]) if (r.get("pressure") is not None and ref.get("pressure") is not None) else None
    # converged value: walking from the coarse side, first value after which
    # every value (towards the reference) satisfies the tolerances
    order = list(good) if ref is good[-1] else list(reversed(good))
    conv = None
    for i, r in enumerate(order):
        tail = [t for t in order[i:] if t is not ref]
        ok = all(t["dE"] is not None and abs(t["dE"]) < args.tol for t in tail)
        if args.ftol is not None:
            ok = ok and all(t["dF"] is not None and t["dF"] < args.ftol for t in tail)
        if param == "SIGMA":
            # the smearing entropy must be small as well (unknown T*S cannot be checked)
            ok = ok and all(t["ts_a"] is None or abs(t["ts_a"]) < args.tstol for t in tail)
        if ok and r is not ref:
            conv = r
            break
    return param, good, ref, conv, numeric


def print_table(param, runs, ref, conv, args):
    unit = {"ENCUT": "eV", "SIGMA": "eV", "KSPACING": "1/A"}.get(param, "")
    head = (f"{param + (' (' + unit + ')' if unit else ''):>16} {'NKPTS':>6} {'E (eV)':>15} {'E/atom (eV)':>13} "
            f"{'dE (meV/at)':>12} {'Fmax':>8} {'dF':>8} {'P (kB)':>8} {'T*S (meV/at)':>13} {'time':>8}")
    print(head)
    print("-" * len(head))
    rows = []
    for r in runs:
        def f(v, fmt):
            return format(v, fmt) if isinstance(v, (int, float)) and v is not None else "-"
        mark = " <- ref" if r is ref else (" <- converged" if r is conv else "")
        t = r.get("elapsed")
        print(f"{str(r['label']):>16} {f(r.get('nkpts'), 'd'):>6} {f(r['E'], '.6f'):>15} {f(r['Ea'], '.6f'):>13} "
              f"{f(r['dE'], '.3f'):>12} {f(r['fmax'], '.4f'):>8} {f(r['dF'], '.4f'):>8} {f(r.get('pressure'), '.2f'):>8} "
              f"{f(r['ts_a'], '.3f'):>13} {f(t, '.0f'):>8}{mark}")
        if not r.get("finished", True):
            print(f"{'':>16} (warning: {r['dir']} did not finish normally)")
        rows.append(r)
    print()
    crit = f"|dE| < {args.tol:g} meV/atom" + (f" and dF < {args.ftol:g} eV/A" if args.ftol is not None else "")
    if param == "SIGMA":
        crit += f" and |T*S| < {args.tstol:g} meV/atom"
    if conv is not None:
        print(f"converged ({crit}) from {param} = {conv['label']} on (reference {ref['label']}).")
    else:
        print(f"no value satisfies {crit} relative to the reference {ref['label']} - extend the test range.")
    if param == "SIGMA":
        bad = [r for r in runs if r["ts_a"] is not None and abs(r["ts_a"]) >= args.tstol]
        if bad:
            print(f"|T*S| >= {args.tstol:g} meV/atom for SIGMA = " + ", ".join(str(r["label"]) for r in bad)
                  + " -> too large for total energies (compare energy(sigma->0)).")
    if param == "ENCUT":
        print("note: for variable-cell relaxations use ENCUT >= 1.3 x max(ENMAX) regardless (Pulay stress).")


def write_dat(path, param, runs):
    with open(path, "w", encoding="ascii", errors="replace") as fh:
        fh.write(f"# {param}  x  NKPTS  E(eV)  E/atom(eV)  dE(meV/atom)  Fmax(eV/A)  dF(eV/A)  P(kB)  TS(meV/atom)\n")
        for r in runs:
            vals = [r["x"], r.get("nkpts"), r["E"], r["Ea"], r["dE"], r["fmax"], r["dF"], r.get("pressure"), r["ts_a"]]
            fh.write(f"{str(r['label']):>14} " + " ".join("nan" if v is None else f"{v:.8g}" for v in vals) + "\n")


def plot(path, param, runs, ref, conv, args, numeric):
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [str(r["label"]) for r in runs]
    categorical = param == "KPOINTS" or not numeric
    # k-meshes: evenly spaced positions ordered by density (labels carry the mesh)
    x = np.arange(len(runs), dtype=float) if categorical else np.array([r["x"] for r in runs], dtype=float)
    xlab = {"ENCUT": "ENCUT (eV)", "KPOINTS": "k-mesh (ordered by number of k-points)",
            "KSPACING": "KSPACING (1/A)", "SIGMA": "SIGMA (eV)"}.get(param, param)
    have_f = any(r["dF"] is not None and r["dF"] > 1e-6 for r in runs) and len(runs) > 1
    have_ts = param == "SIGMA" and any(r["ts_a"] is not None for r in runs)
    npan = 2 + int(have_f) + int(have_ts)
    fig, axes = plt.subplots(npan, 1, figsize=(6.4, 2.3 * npan + 0.6), sharex=True)
    plt.rcParams.update({"font.size": 10})

    def style(ax, ylabel):
        ax.set_ylabel(ylabel, color=INK)
        ax.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(INK2)
        ax.tick_params(colors=INK2)

    def line(ax, y, color):
        y = np.array([np.nan if v is None else v for v in y], dtype=float)
        ax.plot(x, y, "-o", color=color, lw=1.6, ms=6, mec="white", mew=1.0)
        return y

    ax = axes[0]
    line(ax, [r["Ea"] for r in runs], C1)
    style(ax, "E per atom (eV)")
    ax.set_title(f"{param} convergence  (reference {ref['label']})", color=INK, fontsize=11, loc="left")

    ax = axes[1]
    dE = np.array([np.nan if r["dE"] is None else abs(r["dE"]) for r in runs], dtype=float)
    mask = np.array([r is not ref for r in runs])
    ax.semilogy(x[mask], np.clip(dE[mask], 1e-4, None), "-o", color=C1, lw=1.6, ms=6, mec="white", mew=1.0)
    ax.axhline(args.tol, color=INK2, lw=1.0, ls="--")
    ax.text(x.min(), args.tol * 1.15, f"{args.tol:g} meV/atom", color=INK2, fontsize=8, va="bottom")
    style(ax, "|dE| (meV/atom)")
    if conv is not None:
        xc = x[runs.index(conv)]
        for a in axes:
            a.axvline(xc, color=C2, lw=1.2, ls=":")
        axes[0].annotate(f"converged: {conv['label']}", (xc, 1.0), xycoords=("data", "axes fraction"),
                         xytext=(4, -12), textcoords="offset points", color=INK2, fontsize=8)
    k = 2
    if have_f:
        ax = axes[k]
        y = np.array([np.nan if r["dF"] is None else r["dF"] for r in runs], dtype=float)
        ax.semilogy(x[mask], np.clip(y[mask], 1e-5, None), "-o", color=C1, lw=1.6, ms=6, mec="white", mew=1.0)
        if args.ftol is not None:
            ax.axhline(args.ftol, color=INK2, lw=1.0, ls="--")
        style(ax, "max |dF| (eV/A)")
        k += 1
    if have_ts:
        ax = axes[k]
        line(ax, [r["ts_a"] for r in runs], C1)
        ax.axhline(args.tstol, color=INK2, lw=1.0, ls="--")
        ax.axhline(-args.tstol, color=INK2, lw=1.0, ls="--")
        style(ax, "T*S (meV/atom)")
    axes[-1].set_xlabel(xlab, color=INK)
    if categorical:
        axes[-1].set_xticks(x)
        axes[-1].set_xticklabels(labels, rotation=30 if len(labels) > 5 else 0, fontsize=8)
    if param == "KSPACING":
        axes[-1].invert_xaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=args.dpi)
    print(f"figure written to {path}")
    if args.show:
        plt.show()
    plt.close(fig)


def build_parser():
    ap = argparse.ArgumentParser(
        prog="conv_analyze.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Analyse VASP convergence tests created by conv_test.sh (or any run directories).\n"
                     "分析收敛测试：dE/atom、受力、压强、T*S，并给出收敛参数、画图。"),
        epilog=("Examples / 示例:\n  python conv_analyze.py conv_encut\n"
                "  python conv_analyze.py conv_kpoints --tol 0.5 --ftol 0.005\n"
                "  python conv_analyze.py encut_* --param ENCUT -o encut\n"
                "  python conv_analyze.py conv_kpoints --ref 8x8x8   (same as --ref 8 or --ref 512)\n"
                "  python conv_analyze.py conv_sigma --tstol 0.5\n"
                "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License"))
    ap.add_argument("paths", nargs="+", help="test parent directory or run directories")
    ap.add_argument("--param", help="tested INCAR tag if not in conv_info.txt/dir names (ENCUT, SIGMA, ...)")
    ap.add_argument("--energy", choices=("e0", "toten"), default="e0",
                    help="energy(sigma->0) (default) or free energy TOTEN")
    ap.add_argument("--tol", type=float, default=1.0, help="energy tolerance in meV/atom (default 1)")
    ap.add_argument("--ftol", type=float, help="optional force tolerance in eV/A (e.g. 0.01)")
    ap.add_argument("--ref", metavar="VALUE",
                    help="reference run (default: most accurate run). Numeric value of the tag "
                         "(e.g. 600 for ENCUT). For k-point tests the x value is the PRODUCT "
                         "N1*N2*N3: give the mesh '8x8x8', a single N '8' (matches 8x8x8, or "
                         "8x8x1 for slabs) or the product '512'")
    ap.add_argument("--tstol", type=float, default=1.0,
                    help="SIGMA tests: |T*S| per atom must also be below this value (meV/atom) "
                         "for the '<- converged' mark (default 1)")
    ap.add_argument("--no-constraints", action="store_true", help="ignore Selective dynamics flags for forces")
    ap.add_argument("-o", "--output", help="output prefix (default conv_<param>)")
    ap.add_argument("--no-plot", action="store_true", help="table only")
    ap.add_argument("--show", action="store_true", help="show the figure interactively")
    ap.add_argument("--dpi", type=int, default=200, help="figure resolution (default 200)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return ap


def main(argv=None) -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except AttributeError:
        pass
    args = build_parser().parse_args(argv)
    runs = collect(args.paths, args)
    if not runs:
        print("ERROR: no run directories found", file=sys.stderr)
        return 2
    param, good, ref, conv, numeric = analyse(runs, args)
    print_table(param, good, ref, conv, args)
    prefix = args.output or f"conv_{param.lower()}"
    if prefix.endswith(".png"):
        prefix = prefix[:-4]
    write_dat(prefix + ".dat", param, good)
    print(f"table written to {prefix}.dat")
    if not args.no_plot and len(good) >= 2:
        plot(prefix + ".png", param, good, ref, conv, args, numeric)
    return 0


if __name__ == "__main__":
    sys.exit(main())
