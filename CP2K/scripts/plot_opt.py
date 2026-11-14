#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_opt.py - convergence plots for CP2K GEO_OPT / CELL_OPT runs.

读取 CP2K 输出中每一步的 OPT| 收敛块，画出能量、最大/RMS 梯度（力）、最大/RMS 位移、
压力偏差随优化步数的变化，并用虚线标出收敛阈值（实心点 = 该判据已满足）；
如果存在 <PROJECT>-1.cell（&MOTION/&PRINT/&CELL ON），再加一幅晶胞体积随步数的曲线。
可以在优化运行中途使用。

Usage examples
--------------
  python plot_opt.py h2o_geo_opt_bfgs.out                     # -> h2o_geo_opt_bfgs_opt.png
  python plot_opt.py si_cell_opt.out                          # finds si_cell-1.cell automatically
  python plot_opt.py si_cell_opt.out --cell si_cell-1.cell -o si_cell.png --dpi 200
  python plot_opt.py opt.out --energy-unit kcal/mol --csv opt_steps.csv
  python plot_opt.py opt.out --show                           # interactive window

Notes
-----
* The OPT| block layout (CP2K 2023.1 and newer, including 2026.2):
    OPT| Step number / Total energy [hartree] / Internal pressure [bar] (CELL_OPT) /
    Maximum step size + Convergence limit for maximum step size + ... is converged YES|NO /
    RMS step size / Maximum gradient / RMS gradient (same triplets) /
    Pressure deviation [bar] / Pressure tolerance [bar] / Pressure is converged (CELL_OPT)
  Step 0 only carries the starting energy.  Gradients and step sizes are in atomic units
  (hartree/bohr and bohr).  The older "Informations at step" layout (<= 2022) is read on a
  best-effort basis (no convergence limits are shown then).
* The .cell file has the columns  Step  Time[fs]  Ax Ay Az Bx By Bz Cx Cy Cz  Volume[Angstrom^3].
* Energies are plotted relative to the lowest value.

Dependencies: Python >= 3.9, numpy, matplotlib.

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import csv
import glob
import os
import re
import sys

HARTREE_EV = 27.211386245988
ENERGY_FACTORS = {"hartree": 1.0, "ev": HARTREE_EV, "mev": HARTREE_EV * 1000.0,
                  "kcal/mol": 627.509474, "kj/mol": 2625.499639}
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"

PATTERNS = [
    ("energy", re.compile(r"OPT\|\s+Total energy \[([^\]]+)\]\s+(" + FLOAT + ")"), "unit"),
    ("pressure", re.compile(r"OPT\|\s+Internal pressure \[([^\]]+)\]\s+(" + FLOAT + ")"), "unit"),
    ("pressure_deviation", re.compile(r"OPT\|\s+Pressure deviation(?: with constraint| without constraint)?"
                                      r" \[([^\]]+)\]\s+(" + FLOAT + ")"), "unit"),
    ("pressure_tolerance", re.compile(r"OPT\|\s+Pressure tolerance \[([^\]]+)\]\s+(" + FLOAT + ")"), "unit"),
    ("max_step_limit", re.compile(r"OPT\|\s+Convergence limit for maximum step size\s+(" + FLOAT + ")"), None),
    ("rms_step_limit", re.compile(r"OPT\|\s+Convergence limit for RMS step size\s+(" + FLOAT + ")"), None),
    ("max_grad_limit", re.compile(r"OPT\|\s+Convergence limit for maximum gradient\s+(" + FLOAT + ")"), None),
    ("rms_grad_limit", re.compile(r"OPT\|\s+Convergence limit for RMS gradient\s+(" + FLOAT + ")"), None),
    ("max_step", re.compile(r"OPT\|\s+Maximum step size\s+(" + FLOAT + r")\s*$"), None),
    ("rms_step", re.compile(r"OPT\|\s+RMS step size\s+(" + FLOAT + r")\s*$"), None),
    ("max_grad", re.compile(r"OPT\|\s+Maximum gradient\s+(" + FLOAT + r")\s*$"), None),
    ("rms_grad", re.compile(r"OPT\|\s+RMS gradient\s+(" + FLOAT + r")\s*$"), None),
]
FLAGS = [
    ("max_step_ok", re.compile(r"OPT\|\s+Maximum step size is converged\s+(YES|NO)")),
    ("rms_step_ok", re.compile(r"OPT\|\s+RMS step size is converged\s+(YES|NO)")),
    ("max_grad_ok", re.compile(r"OPT\|\s+Maximum gradient is converged\s+(YES|NO)")),
    ("rms_grad_ok", re.compile(r"OPT\|\s+RMS gradient is converged\s+(YES|NO)")),
    ("pressure_ok", re.compile(r"OPT\|\s+Pressure is converged\s+(YES|NO)")),
]
RE_STEP = re.compile(r"OPT\|\s+Step number\s+(\d+)")
RE_PROJECT = re.compile(r"GLOBAL\|\s+Project name\s+(\S+)")
RE_RUNTYPE = re.compile(r"GLOBAL\|\s+Run type\s+(\S+)")
# legacy (<= 2022) layout
RE_LEG_STEP = re.compile(r"Informations at step\s*=\s*(\d+)")
LEGACY = [
    ("energy", re.compile(r"^\s*Total Energy\s*=\s*(" + FLOAT + ")")),
    ("pressure", re.compile(r"^\s*Internal Pressure \[bar\]\s*=\s*(" + FLOAT + ")")),
    ("max_step", re.compile(r"^\s*Max\. step size\s*=\s*(" + FLOAT + ")")),
    ("rms_step", re.compile(r"^\s*RMS step size\s*=\s*(" + FLOAT + ")")),
    ("max_grad", re.compile(r"^\s*Max\. gradient\s*=\s*(" + FLOAT + ")")),
    ("rms_grad", re.compile(r"^\s*RMS gradient\s*=\s*(" + FLOAT + ")")),
    ("pressure_deviation", re.compile(r"^\s*Pressure Deviation \[bar\]\s*=\s*(" + FLOAT + ")")),
]

CSV_COLUMNS = ["step", "energy", "max_step", "rms_step", "max_grad", "rms_grad", "pressure",
               "pressure_deviation", "volume", "max_step_ok", "rms_step_ok", "max_grad_ok",
               "rms_grad_ok", "pressure_ok"]


def to_float(text):
    return float(text.replace("D", "E").replace("d", "e"))


def parse_opt(path):
    """Return (steps, info) parsed from a CP2K output file."""
    if not os.path.isfile(path):
        raise SystemExit(f"ERROR: output file not found: {path}")
    steps, info, cur = [], {"energy_unit": "hartree", "pressure_unit": "bar"}, None
    status = None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if "project" not in info:
                m = RE_PROJECT.search(line)
                if m:
                    info["project"] = m.group(1)
            if "run_type" not in info:
                m = RE_RUNTYPE.search(line)
                if m:
                    info["run_type"] = m.group(1)
            m = RE_STEP.search(line)
            if m:
                cur = {"step": int(m.group(1))}
                steps.append(cur)
                continue
            m = RE_LEG_STEP.search(line)
            if m:
                cur = {"step": int(m.group(1)), "legacy": True}
                steps.append(cur)
                continue
            if "GEOMETRY OPTIMIZATION COMPLETED" in line:
                status = "converged"
            elif "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" in line:
                status = "max_iter_reached"
            if cur is None:
                continue
            if cur.get("legacy"):
                for key, rx in LEGACY:
                    m = rx.match(line)
                    if m:
                        cur[key] = to_float(m.group(1))
                        break
                continue
            if "OPT|" not in line:
                continue
            for key, rx in FLAGS:
                m = rx.search(line)
                if m:
                    cur[key] = m.group(1) == "YES"
                    break
            else:
                for key, rx, has_unit in PATTERNS:
                    m = rx.search(line)
                    if m:
                        if has_unit:
                            cur[key] = to_float(m.group(2))
                            if key == "energy":
                                info["energy_unit"] = m.group(1).strip().lower()
                            elif key == "pressure":
                                info["pressure_unit"] = m.group(1).strip()
                        else:
                            cur[key] = to_float(m.group(1))
                        break
    info["status"] = status
    return steps, info


def read_cell_file(path):
    """Return {step: volume} from a CP2K <project>-1.cell file."""
    vols = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tok = s.split()
            if len(tok) < 12:
                raise SystemExit(f"ERROR: {path}:{ln}: expected 12 columns (Step Time Ax..Cz Volume), got {len(tok)}")
            try:
                vols[int(tok[0])] = float(tok[11])
            except ValueError:
                raise SystemExit(f"ERROR: {path}:{ln}: cannot read numbers from '{s[:60]}'")
    return vols


def find_cell_file(outpath, project):
    d = os.path.dirname(os.path.abspath(outpath))
    if project:
        cand = os.path.join(d, f"{project}-1.cell")
        if os.path.isfile(cand):
            return cand
    cands = sorted(glob.glob(os.path.join(d, "*-1.cell")))
    return cands[0] if len(cands) == 1 else None


def build_parser():
    p = argparse.ArgumentParser(
        description="Plot energy, gradients, step sizes, pressure deviation and cell volume "
                    "versus optimization step for CP2K GEO_OPT / CELL_OPT.",
        epilog="Examples:\n"
               "  python plot_opt.py h2o_geo_opt_bfgs.out\n"
               "  python plot_opt.py si_cell_opt.out --cell si_cell-1.cell -o si.png\n"
               "  python plot_opt.py opt.out --energy-unit kcal/mol --csv steps.csv --no-plot",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("output", help="CP2K output file of a GEO_OPT or CELL_OPT run")
    p.add_argument("--cell", metavar="FILE",
                   help="<PROJECT>-1.cell file for the volume panel (default: auto-detect next to the output)")
    p.add_argument("--no-cell", action="store_true", help="do not look for a .cell file")
    p.add_argument("--energy-unit", default="meV", choices=["hartree", "eV", "meV", "kcal/mol", "kJ/mol"],
                   help="unit of the relative energy axis (default meV)")
    p.add_argument("-o", "--outfile", help="image file (default: <output stem>_opt.png)")
    p.add_argument("--dpi", type=int, default=150, help="image resolution (default 150)")
    p.add_argument("--csv", metavar="FILE", help="also write the parsed steps to a CSV file")
    p.add_argument("--no-plot", action="store_true", help="only print the table / write CSV")
    p.add_argument("--show", action="store_true", help="open an interactive window instead of saving")
    return p


def print_table(steps, info):
    print(f"{'step':>5} {'energy[' + info['energy_unit'] + ']':>20} {'max_grad':>10} {'rms_grad':>10} "
          f"{'max_step':>10} {'rms_step':>10} {'P_dev':>10} {'volume':>10}")

    def f(v, spec):
        return "-" if v is None else format(v, spec)
    for s in steps:
        print(f"{s['step']:>5} {f(s.get('energy'), '.10f'):>20} {f(s.get('max_grad'), '.2e'):>10} "
              f"{f(s.get('rms_grad'), '.2e'):>10} {f(s.get('max_step'), '.2e'):>10} "
              f"{f(s.get('rms_step'), '.2e'):>10} {f(s.get('pressure_deviation'), '.1f'):>10} "
              f"{f(s.get('volume'), '.3f'):>10}")


def plot(steps, info, args, outfile):
    try:
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise SystemExit(f"ERROR: plotting needs numpy and matplotlib ({exc}); use --no-plot or "
                         "'pip install numpy matplotlib'")

    def series(key):
        xs = [s["step"] for s in steps if s.get(key) is not None]
        ys = [s[key] for s in steps if s.get(key) is not None]
        return np.array(xs), np.array(ys)

    def limit(key):
        vals = [s[key] for s in steps if s.get(key) is not None]
        return vals[-1] if vals else None

    panels = ["energy", "grad", "step"]
    if any(s.get("pressure_deviation") is not None for s in steps):
        panels.append("pressure")
    if any(s.get("volume") is not None for s in steps):
        panels.append("volume")

    fig, axes = plt.subplots(len(panels), 1, figsize=(6.4, 2.1 * len(panels) + 0.6), sharex=True)
    axes = np.atleast_1d(axes)
    blue, orange, grey = "#1f5aa6", "#d9822b", "0.45"

    def draw(ax, key, okkey, color, label):
        x, y = series(key)
        if len(x) == 0:
            return
        ax.plot(x, y, "-", color=color, lw=1.2, label=label)
        ok = np.array([bool(s.get(okkey)) for s in steps if s.get(key) is not None])
        ax.plot(x[ok], y[ok], "o", color=color, ms=4)
        ax.plot(x[~ok], y[~ok], "o", mfc="white", mec=color, ms=4)
        lim = limit(key + "_limit")
        if lim:
            ax.axhline(lim, ls="--", lw=0.9, color=color, alpha=0.8)

    for ax, name in zip(axes, panels):
        if name == "energy":
            x, y = series("energy")
            fac = ENERGY_FACTORS[args.energy_unit.lower()] / ENERGY_FACTORS.get(info["energy_unit"], 1.0)
            if len(y):
                ax.plot(x, (y - y.min()) * fac, "o-", color=blue, ms=3.5, lw=1.2)
            ax.set_ylabel(f"E - E_min ({args.energy_unit})")
            title = f"{info.get('project', '')} {info.get('run_type', '')}".strip()
            status = {"converged": "converged", "max_iter_reached": "MAX_ITER reached"}.get(
                info["status"], "not finished")
            ax.set_title(f"{title}  [{status}]", fontsize=10)
        elif name == "grad":
            draw(ax, "max_grad", "max_grad_ok", blue, "max gradient")
            draw(ax, "rms_grad", "rms_grad_ok", orange, "RMS gradient")
            ax.set_yscale("log")
            ax.set_ylabel("gradient (Ha/bohr)")
            ax.legend(fontsize=8, frameon=False)
        elif name == "step":
            draw(ax, "max_step", "max_step_ok", blue, "max step")
            draw(ax, "rms_step", "rms_step_ok", orange, "RMS step")
            ax.set_yscale("log")
            ax.set_ylabel("step size (bohr)")
            ax.legend(fontsize=8, frameon=False)
        elif name == "pressure":
            x, y = series("pressure_deviation")
            ok = np.array([bool(s.get("pressure_ok")) for s in steps if s.get("pressure_deviation") is not None])
            y = np.abs(y)
            ax.plot(x, y, "-", color=blue, lw=1.2)
            ax.plot(x[ok], y[ok], "o", color=blue, ms=4)
            ax.plot(x[~ok], y[~ok], "o", mfc="white", mec=blue, ms=4)
            tol = limit("pressure_tolerance")
            if tol:
                ax.axhline(tol, ls="--", lw=0.9, color=blue)
            if np.all(y > 0):
                ax.set_yscale("log")
            ax.set_ylabel(f"|P_int - P_ext| ({info['pressure_unit']})")
        elif name == "volume":
            x, y = series("volume")
            ax.plot(x, y, "o-", color=blue, ms=3.5, lw=1.2)
            ax.set_ylabel("volume (A^3)")
        ax.grid(alpha=0.3, which="both")
    axes[-1].set_xlabel("optimization step")
    fig.text(0.99, 0.005, "dashed: convergence limit; filled: criterion met", ha="right",
             fontsize=7, color=grey)
    fig.tight_layout()
    if args.show:
        plt.show()
    else:
        fig.savefig(outfile, dpi=args.dpi)
        print(f"saved {outfile}")


def main(argv=None):
    args = build_parser().parse_args(argv)
    steps, info = parse_opt(args.output)
    if not steps:
        raise SystemExit(f"ERROR: no OPT| blocks found in {args.output} "
                         "(is RUN_TYPE GEO_OPT/CELL_OPT and has the first step finished?)")
    if not args.no_cell:
        cell = args.cell or find_cell_file(args.output, info.get("project"))
        if args.cell and not os.path.isfile(args.cell):
            raise SystemExit(f"ERROR: cell file not found: {args.cell}")
        if cell:
            vols = read_cell_file(cell)
            for s in steps:
                if s["step"] in vols:
                    s["volume"] = vols[s["step"]]
            print(f"cell file: {cell} ({len(vols)} rows)")
    print_table(steps, info)
    last = steps[-1]
    print(f"status: {info['status'] or 'not finished / interrupted'}; last step {last['step']}")
    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(CSV_COLUMNS)
            for s in steps:
                w.writerow(["" if s.get(c) is None else s.get(c) for c in CSV_COLUMNS])
        print(f"written: {args.csv}")
    if not args.no_plot:
        out = args.outfile or os.path.splitext(args.output)[0] + "_opt.png"
        plot(steps, info, args, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
