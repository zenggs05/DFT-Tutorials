#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conv_cutoff.py - CUTOFF / REL_CUTOFF convergence driver for CP2K (two-step method).

CP2K 多重网格截断能的两步收敛测试：在模板中用 @SET 变量替换 CUTOFF 和 REL_CUTOFF，
生成扫描目录 cutoff_<值>Ry/ 或 relcutoff_<值>Ry/；可在本地依次（或并行）运行，或逐个 sbatch 提交；
随后从每个输出中提取 "ENERGY| Total FORCE_EVAL"、SCF 是否收敛、各级网格的 "count for grid"
以及（如有）"STRESS| 1/3 Trace"，写成 CSV 并作图，给出满足阈值的最小截断。

Usage examples
--------------
  # step 1: scan CUTOFF at fixed REL_CUTOFF = 60 Ry, run locally with 4 MPI ranks, then collect
  python conv_cutoff.py run -t template.inp --scan cutoff --values 150 200 250 300 400 500 600 \\
         --rel-cutoff 60 --extra si_bulk8.xyz --cp2k "mpirun -np 4 cp2k.psmp"

  # step 2: scan REL_CUTOFF at the CUTOFF chosen in step 1, submit every point with sbatch
  python conv_cutoff.py gen -t template.inp --scan relcutoff --values 40 50 60 70 80 100 --cutoff 400 \\
         --extra si_bulk8.xyz
  python conv_cutoff.py submit --scan relcutoff --job job_cp2k_slurm.sh
  python conv_cutoff.py collect --scan relcutoff --threshold 0.03        # after the jobs finished

  # only rebuild CSV + figure from existing directories
  python conv_cutoff.py collect --scan cutoff --natoms 8 -o cutoff_scan.png

Sub-commands
------------
  gen      write <scan>_<value>Ry/input.inp (+ copies of --extra files)
  run      gen + run every point locally (--cp2k, --npar) + collect
  submit   submit every existing <scan>_*Ry directory with sbatch/qsub (--job, --submit-cmd)
  collect  parse <scan>_*Ry/output.out -> <scan>_scan.csv and <scan>_scan.png

Template requirements
---------------------
The template must contain the two preprocessor lines (any values)
    @SET CUTOFF 400
    @SET REL_CUTOFF 60
and use ${CUTOFF} / ${REL_CUTOFF} in &MGRID (see CP2K/examples/04-cutoff-convergence-si/template.inp).
Other variable names can be given with --var-cutoff / --var-rel.  Keep PRINT_LEVEL MEDIUM (or higher)
so that the "count for grid" lines of the MULTIGRID INFO block are printed, use a tight EPS_SCF
(1e-8) and EPS_DEFAULT (1e-12), and switch off .wfn restart files for the scan.

Output lines parsed (CP2K 2023.x - 2026.2):
  ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]   (2024.x and older: "[a.u.]:")
  *** SCF run converged in N steps ***   /   SCF run NOT converged
  count for grid        1:            912          cutoff [a.u.]          150.00
  STRESS| 1/3 Trace   (only if &FORCE_EVAL/&PRINT/&STRESS_TENSOR is on)
  - Atoms:   (number of atoms, for meV/atom)

Dependencies: Python >= 3.9; numpy + matplotlib only for the figure (collect without them: --no-plot).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import csv
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys
import time

HARTREE_MEV = 27211.386245988
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
RE_ENERGY = re.compile(r"ENERGY\|\s+Total FORCE_EVAL\s+\(\s*\S+\s*\)\s+energy\s+\[([^\]]+)\]:?\s+(" + FLOAT + ")")
RE_GRID = re.compile(r"count for grid\s+(\d+):\s+(\d+)\s+cutoff \[a\.u\.\]\s+(" + FLOAT + ")")
RE_STRESS = re.compile(r"STRESS\|\s+1/3 Trace\s+(" + FLOAT + ")")
RE_NATOMS = re.compile(r"^\s*-\s+Atoms:\s+(\d+)\s*$")
RE_SCF_OK = re.compile(r"SCF run converged in\s+\d+\s+steps|outer SCF loop converged in")
RE_SCF_BAD = re.compile(r"SCF run NOT converged|outer SCF loop FAILED")
RE_DIR = re.compile(r"^(cutoff|relcutoff)_(" + FLOAT + r")Ry$")


def fmt_value(v):
    return f"{v:g}"


def dir_name(scan, value):
    return f"{scan}_{fmt_value(value)}Ry"


# ----------------------------------------------------------------------------- gen
def make_input(template_text, var_cut, var_rel, cutoff, rel_cutoff, template_name):
    out, found = [], {var_cut: False, var_rel: False}
    rx = re.compile(r"^(\s*)@SET\s+(\S+)\s+.*$", re.IGNORECASE)
    for line in template_text.splitlines():
        m = rx.match(line)
        if m and m.group(2) in found:
            val = cutoff if m.group(2) == var_cut else rel_cutoff
            out.append(f"{m.group(1)}@SET {m.group(2)} {fmt_value(val)}")
            found[m.group(2)] = True
        else:
            out.append(line)
    missing = [k for k, v in found.items() if not v]
    if missing:
        raise SystemExit(f"ERROR: {template_name} has no '@SET {missing[0]} <value>' line "
                         "(the scan variables must be defined with @SET, see -h)")
    return "\n".join(out) + "\n"


def cmd_gen(args):
    if not os.path.isfile(args.template):
        raise SystemExit(f"ERROR: template not found: {args.template}")
    with open(args.template, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    for f in args.extra:
        if not os.path.isfile(f):
            raise SystemExit(f"ERROR: extra file not found: {f}")
    if args.scan == "cutoff" and args.rel_cutoff is None:
        raise SystemExit("ERROR: --scan cutoff needs the fixed --rel-cutoff (e.g. 60)")
    if args.scan == "relcutoff" and args.cutoff is None:
        raise SystemExit("ERROR: --scan relcutoff needs the fixed --cutoff chosen in step 1")
    dirs = []
    for v in args.values:
        cut = v if args.scan == "cutoff" else args.cutoff
        rel = v if args.scan == "relcutoff" else args.rel_cutoff
        d = dir_name(args.scan, v)
        os.makedirs(d, exist_ok=True)
        dest = os.path.join(d, args.inp_name)
        if os.path.exists(os.path.join(d, args.out_name)) and not args.force:
            print(f"skip {d}: {args.out_name} already exists (use --force to regenerate)")
            continue
        with open(dest, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(make_input(text, args.var_cutoff, args.var_rel, cut, rel, args.template))
        for f in args.extra:
            shutil.copy2(f, d)
        dirs.append(d)
        print(f"written {dest}  (CUTOFF {fmt_value(cut)} Ry, REL_CUTOFF {fmt_value(rel)} Ry)")
    return dirs


# ----------------------------------------------------------------------------- run / submit
def cmd_run(args):
    dirs = cmd_gen(args)
    base = shlex.split(args.cp2k)
    if shutil.which(base[0]) is None:
        raise SystemExit(f"ERROR: '{base[0]}' not found in PATH (set --cp2k, load the CP2K environment first)")
    running = []
    for d in dirs:
        cmd = base + ["-i", args.inp_name, "-o", args.out_name]
        print(f"[{time.strftime('%H:%M:%S')}] running in {d}: {' '.join(cmd)}")
        running.append((d, subprocess.Popen(cmd, cwd=d)))
        while len([p for _, p in running if p.poll() is None]) >= args.npar:
            time.sleep(2)
    failed = [d for d, p in running if p.wait() != 0]
    if failed:
        print(f"WARNING: non-zero exit status in {', '.join(failed)} (see their {args.out_name})", file=sys.stderr)
    return cmd_collect(args)


def cmd_submit(args):
    if not os.path.isfile(args.job):
        raise SystemExit(f"ERROR: job script not found: {args.job}")
    dirs = find_dirs(args.scan)
    if not dirs:
        raise SystemExit(f"ERROR: no {args.scan}_*Ry directories (run 'gen' first)")
    submit = shlex.split(args.submit_cmd)
    if shutil.which(submit[0]) is None:
        raise SystemExit(f"ERROR: '{submit[0]}' not found (are you on the cluster login node?)")
    job = os.path.abspath(args.job)
    for _, d in dirs:
        if os.path.exists(os.path.join(d, args.out_name)) and not args.force:
            print(f"skip {d}: {args.out_name} exists")
            continue
        if os.path.basename(submit[0]) == "qsub":
            cmd = submit + ["-v", f"INP={args.inp_name},OUT={args.out_name}", job]
        else:
            cmd = submit + [f"--export=ALL,INP={args.inp_name},OUT={args.out_name}", job]
        print(f"{d}: {' '.join(cmd)}")
        subprocess.run(cmd, cwd=d, check=True)
    print(f"after all jobs have finished:  python {os.path.basename(sys.argv[0])} collect --scan {args.scan}")
    return 0


# ----------------------------------------------------------------------------- collect
def find_dirs(scan):
    found = []
    for d in glob.glob(f"{scan}_*Ry"):
        m = RE_DIR.match(os.path.basename(d))
        if m and m.group(1) == scan and os.path.isdir(d):
            found.append((float(m.group(2)), d))
    return sorted(found)


def parse_out(path):
    res = {"energy": None, "unit": None, "grids": {}, "pressure": None, "natoms": None,
           "scf": "unknown", "ended": False}
    ok = bad = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = RE_ENERGY.search(line)
            if m:
                res["energy"], res["unit"] = float(m.group(2)), m.group(1).strip().lower()
                continue
            m = RE_GRID.search(line)
            if m:
                res["grids"][int(m.group(1))] = (int(m.group(2)), float(m.group(3)))
                continue
            m = RE_STRESS.search(line)
            if m:
                res["pressure"] = float(m.group(1))
                continue
            if res["natoms"] is None:
                m = RE_NATOMS.match(line)
                if m:
                    res["natoms"] = int(m.group(1))
            if RE_SCF_OK.search(line):
                ok += 1
            elif RE_SCF_BAD.search(line):
                bad += 1
            if "PROGRAM ENDED AT" in line:
                res["ended"] = True
    res["scf"] = "NO" if bad else ("yes" if ok else "unknown")
    return res


def cmd_collect(args):
    dirs = find_dirs(args.scan)
    if not dirs:
        raise SystemExit(f"ERROR: no {args.scan}_*Ry directories in {os.getcwd()}")
    rows, ngrid = [], 0
    for val, d in dirs:
        f = os.path.join(d, args.out_name)
        if not os.path.isfile(f):
            print(f"skip {d}: no {args.out_name}", file=sys.stderr)
            continue
        r = parse_out(f)
        if r["energy"] is None:
            print(f"skip {d}: no 'ENERGY| Total FORCE_EVAL' line (SCF aborted or job unfinished)", file=sys.stderr)
            continue
        if r["unit"] not in ("hartree", "a.u."):
            raise SystemExit(f"ERROR: {f}: energy printed in [{r['unit']}]; remove ENERGY_UNIT or use hartree")
        if not r["grids"]:
            print(f"note {d}: no 'count for grid' lines (PRINT_LEVEL LOW?)", file=sys.stderr)
        ngrid = max(ngrid, max(r["grids"]) if r["grids"] else 0)
        r["value"], r["dir"] = val, d
        rows.append(r)
    if not rows:
        raise SystemExit("ERROR: no usable output found")
    natoms = args.natoms or rows[-1]["natoms"]
    if not natoms:
        raise SystemExit("ERROR: number of atoms unknown; give --natoms")
    ref = rows[-1]["energy"]  # highest scanned value
    csv_name = args.csv or f"{args.scan}_scan.csv"
    header = [f"{args.scan}_Ry", "energy_Ha", "dE_meV_per_atom", "scf_converged", "pressure_bar"] + \
             [f"grid{i}_count" for i in range(1, ngrid + 1)]
    with open(csv_name, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            r["dE"] = (r["energy"] - ref) / natoms * HARTREE_MEV
            w.writerow([fmt_value(r["value"]), f"{r['energy']:.10f}", f"{r['dE']:.4f}", r["scf"],
                        "" if r["pressure"] is None else f"{r['pressure']:.4f}"] +
                       [r["grids"].get(i, ("", 0))[0] for i in range(1, ngrid + 1)])

    # table on screen
    print(f"{'value':>8} {'energy_Ha':>18} {'dE_meV/atom':>12} {'scf':>4} {'P_bar':>12}  grid counts")
    for r in rows:
        counts = " ".join(str(r["grids"][i][0]) for i in sorted(r["grids"]))
        p = "" if r["pressure"] is None else f"{r['pressure']:.2f}"
        print(f"{fmt_value(r['value']):>8} {r['energy']:>18.10f} {r['dE']:>12.4f} {r['scf']:>4} {p:>12}  {counts}")
    print(f"written: {csv_name}  (natoms = {natoms}, reference = {fmt_value(rows[-1]['value'])} Ry = last row)")

    # convergence statement
    good = [abs(r["dE"]) < args.threshold for r in rows]
    for i, r in enumerate(rows[:-1]):
        if all(good[i:]):
            print(f"|dE| < {args.threshold:g} meV/atom from {args.scan.upper()} = {fmt_value(r['value'])} Ry on "
                  "(relative to the highest value scanned; extend the scan if this is the last-but-one point)")
            break
    else:
        print(f"no scanned value is converged to {args.threshold:g} meV/atom yet - extend the scan")
    for r in rows:
        if r["grids"] and r["grids"].get(1, (1, 0))[0] == 0:
            print(f"WARNING: count for grid 1 is 0 at {fmt_value(r['value'])} Ry - the finest grid is unused; "
                  "increase REL_CUTOFF instead of CUTOFF")
            break
    if any(r["scf"] == "NO" for r in rows):
        print("WARNING: some points have SCF NOT converged - do not use them for the decision")
    if not args.no_plot:
        plot(rows, args, ngrid)
    return 0


def plot(rows, args, ngrid):
    try:
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        print(f"WARNING: no plot ({exc}); install numpy + matplotlib or use --no-plot", file=sys.stderr)
        return
    has_p = any(r["pressure"] is not None for r in rows)
    npan = 2 + (1 if has_p else 0) if ngrid else 1 + (1 if has_p else 0)
    fig, axes = plt.subplots(npan, 1, figsize=(5.6, 2.4 * npan + 0.4), sharex=True)
    axes = np.atleast_1d(axes)
    x = np.array([r["value"] for r in rows])
    de = np.array([abs(r["dE"]) for r in rows])
    bad = np.array([r["scf"] == "NO" for r in rows])
    ax = axes[0]
    mask = de > 0
    ax.semilogy(x[mask], de[mask], "o-", color="#1f5aa6", label="|E - E_ref|")
    if np.any(bad & mask):
        ax.semilogy(x[bad & mask], de[bad & mask], "x", color="#c0392b", ms=9, label="SCF not converged")
    ax.axhline(args.threshold, ls="--", color="0.4", lw=1, label=f"{args.threshold:g} meV/atom")
    ax.set_ylabel("|dE| (meV/atom)")
    ax.legend(frameon=False, fontsize=8)
    k = 1
    if has_p:
        xp = np.array([r["value"] for r in rows if r["pressure"] is not None])
        yp = np.array([r["pressure"] for r in rows if r["pressure"] is not None])
        axes[k].plot(xp, yp, "o-", color="#d9822b")
        axes[k].set_ylabel("STRESS 1/3 Trace (bar)")
        k += 1
    if ngrid:
        for g in range(1, ngrid + 1):
            axes[k].plot(x, [r["grids"].get(g, (0, 0))[0] for r in rows], "o-", ms=3, label=f"grid {g}")
        axes[k].set_ylabel("count for grid")
        axes[k].legend(frameon=False, fontsize=7, ncol=min(ngrid, 4))
    for a in axes:
        a.grid(alpha=0.3, which="both")
    axes[-1].set_xlabel(("CUTOFF" if args.scan == "cutoff" else "REL_CUTOFF") + " (Ry)")
    fig.tight_layout()
    if args.show:
        plt.show()
    else:
        out = args.outfile or f"{args.scan}_scan.png"
        fig.savefig(out, dpi=args.dpi)
        print(f"saved {out}")


# ----------------------------------------------------------------------------- CLI
def build_parser():
    p = argparse.ArgumentParser(
        description="Two-step CUTOFF / REL_CUTOFF convergence test for CP2K: generate, run or submit, "
                    "collect (energy, SCF flag, count for grid, pressure) and plot.",
        epilog="Typical workflow:\n"
               "  python conv_cutoff.py run -t template.inp --scan cutoff --values 200 300 400 500 600 "
               "--rel-cutoff 60 --extra si.xyz\n"
               "  python conv_cutoff.py run -t template.inp --scan relcutoff --values 40 50 60 80 100 "
               "--cutoff 400 --extra si.xyz\n"
               "Run 'python conv_cutoff.py <command> -h' for the options of each command.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--scan", choices=["cutoff", "relcutoff"], required=True,
                        help="which parameter is scanned (directory prefix cutoff_ / relcutoff_)")
    common.add_argument("--inp-name", default="input.inp", help="input file name in each directory (default input.inp)")
    common.add_argument("--out-name", default="output.out", help="output file name in each directory (default output.out)")
    common.add_argument("--force", action="store_true", help="regenerate / resubmit points that already have output")

    gen = argparse.ArgumentParser(add_help=False)
    gen.add_argument("-t", "--template", required=True, help="template input with '@SET CUTOFF' and '@SET REL_CUTOFF'")
    gen.add_argument("--values", type=float, nargs="+", required=True, help="scanned values in Ry")
    gen.add_argument("--cutoff", type=float, help="fixed CUTOFF (Ry) for --scan relcutoff")
    gen.add_argument("--rel-cutoff", type=float, help="fixed REL_CUTOFF (Ry) for --scan cutoff")
    gen.add_argument("--extra", nargs="*", default=[], help="files copied into every directory (xyz, basis ...)")
    gen.add_argument("--var-cutoff", default="CUTOFF", help="@SET variable used for CUTOFF (default CUTOFF)")
    gen.add_argument("--var-rel", default="REL_CUTOFF", help="@SET variable used for REL_CUTOFF (default REL_CUTOFF)")

    col = argparse.ArgumentParser(add_help=False)
    col.add_argument("--natoms", type=int, help="atoms in the cell (default: read '- Atoms:' from the output)")
    col.add_argument("--threshold", type=float, default=0.03,
                     help="convergence threshold in meV/atom (default 0.03, about 1e-6 Ha/atom)")
    col.add_argument("--csv", help="CSV file name (default <scan>_scan.csv)")
    col.add_argument("-o", "--outfile", help="figure file name (default <scan>_scan.png)")
    col.add_argument("--dpi", type=int, default=200, help="figure resolution (default 200)")
    col.add_argument("--no-plot", action="store_true", help="skip the figure")
    col.add_argument("--show", action="store_true", help="show the figure interactively instead of saving")

    sub.add_parser("gen", parents=[common, gen], help="generate the scan directories only")
    r = sub.add_parser("run", parents=[common, gen, col], help="generate, run locally and collect")
    r.add_argument("--cp2k", default="cp2k.psmp",
                   help='CP2K command, may include the MPI launcher (default "cp2k.psmp"; '
                        'e.g. "mpirun -np 8 cp2k.psmp")')
    r.add_argument("--npar", type=int, default=1, help="number of points run at the same time (default 1)")
    s = sub.add_parser("submit", parents=[common], help="submit every scan directory as a batch job")
    s.add_argument("--job", required=True, help="job script, e.g. ../scripts/job_cp2k_slurm.sh (reads INP/OUT)")
    s.add_argument("--submit-cmd", default="sbatch",
                   help='submission command (default "sbatch"; "qsub" passes INP/OUT with -v)')
    sub.add_parser("collect", parents=[common, col], help="parse outputs, write CSV and plot")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "gen":
        cmd_gen(args)
        return 0
    if args.command == "run":
        return cmd_run(args)
    if args.command == "submit":
        return cmd_submit(args)
    return cmd_collect(args)


if __name__ == "__main__":
    sys.exit(main())
