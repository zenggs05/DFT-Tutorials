#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conv_kpoints.py - k-point mesh convergence test for CP2K (Monkhorst-Pack).

按晶格长度和 k 点密度自动给出 n1×n2×n3 网格序列（非周期方向固定为 1），把模板中
&KPOINTS 的 SCHEME 行改写为 "SCHEME MONKHORST-PACK n1 n2 n3"，生成 kp_<n1>x<n2>x<n3>/ 目录，
可在本地运行或用 sbatch/qsub 提交；最后汇总每原子能量、相对最密网格的能量差、费米能级、
压力（STRESS| 1/3 Trace）和不可约 k 点数，写成 CSV 并作图。

Usage examples
--------------
  # meshes from "length" values R (n_i = max(1, ceil(R * |b_i|)), |b_i| without 2*pi, R in Angstrom)
  python conv_kpoints.py run -t kpoint_template.inp --lengths 10 15 20 25 30 40 \\
         --cp2k "mpirun -np 8 cp2k.psmp"

  # meshes from a k-spacing in 1/Angstrom (2*pi included, like VASP KSPACING)
  python conv_kpoints.py gen -t slab.inp --kspacing 0.5 0.4 0.3 0.25 0.2 --extra slab.xyz
  python conv_kpoints.py submit --job job_cp2k_slurm.sh
  python conv_kpoints.py collect --threshold 1.0

  # explicit meshes; cell taken from a structure file (needs ase) instead of the template
  python conv_kpoints.py gen -t t.inp --grids 2x2x2 4x4x4 6x6x6 8x8x8 --structure POSCAR

  # just print the meshes that would be generated
  python conv_kpoints.py mesh -t si_primitive.inp --lengths 10 20 30 40 50

Sub-commands: mesh | gen | run | submit | collect   (see "<command> -h")

Template requirements
---------------------
* a &KPOINTS section with a SCHEME line inside &FORCE_EVAL/&DFT (the line is replaced; any
  ${VAR} on it is fine, e.g. the @SET NK template of CP2K/examples/06-periodic-surface-kpoints);
* the cell is read from &SUBSYS/&CELL (ABC [+ ALPHA_BETA_GAMMA] or A/B/C vectors; units
  [angstrom]/[bohr] allowed).  If the cell uses ${variables} or CELL_FILE_NAME, give
  --cell a b c [alpha beta gamma] or --structure FILE;
* PERIODIC of &CELL decides which directions get n = 1 (override with --nonperiodic).
* 2026.2 and older: k-points require diagonalization (OT + k-points only in trunk), so the
  template must use &DIAGONALIZATION (+ &MIXING, and &SMEAR/ADDED_MOS for metals).
* For the pressure column switch on STRESS_TENSOR ANALYTICAL and &FORCE_EVAL/&PRINT/&STRESS_TENSOR.

Output lines parsed (CP2K 2023.x - 2026.2): "ENERGY| Total FORCE_EVAL ... energy",
"Fermi energy:" (Hartree, printed with &SMEAR), "STRESS| 1/3 Trace", "BRILLOUIN| List of Kpoints"
(number of k-points actually used), "- Atoms:", SCF convergence lines.

Dependencies: Python >= 3.9; numpy + matplotlib for the figure; ase only for --structure.

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import csv
import glob
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time

HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903
FLOAT = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
RE_ENERGY = re.compile(r"ENERGY\|\s+Total FORCE_EVAL\s+\(\s*\S+\s*\)\s+energy\s+\[([^\]]+)\]:?\s+(" + FLOAT + ")")
RE_FERMI = re.compile(r"^\s*Fermi energy:\s+(" + FLOAT + ")")
RE_STRESS = re.compile(r"STRESS\|\s+1/3 Trace\s+(" + FLOAT + ")")
RE_NKP = re.compile(r"BRILLOUIN\|\s+List of Kpoints.*?(\d+)\s*$")
RE_NATOMS = re.compile(r"^\s*-\s+Atoms:\s+(\d+)\s*$")
RE_SCF_OK = re.compile(r"SCF run converged in\s+\d+\s+steps|outer SCF loop converged in")
RE_SCF_BAD = re.compile(r"SCF run NOT converged|outer SCF loop FAILED")
RE_DIR = re.compile(r"^kp_(\d+)x(\d+)x(\d+)$")
AXES = "XYZ"


# ----------------------------------------------------------------------------- geometry
def cross(u, v):
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def norm(u):
    return math.sqrt(dot(u, u))


def cell_from_abc(abc, angles):
    a, b, c = abc
    al, be, ga = (math.radians(x) for x in angles)
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    return [[a, 0.0, 0.0], [b * math.cos(ga), b * math.sin(ga), 0.0],
            [cx, cy, math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))]]


def recip_lengths(cell):
    """|b_i| without the 2*pi factor (1/Angstrom)."""
    a1, a2, a3 = cell
    vol = dot(a1, cross(a2, a3))
    if abs(vol) < 1e-8:
        raise SystemExit("ERROR: cell volume is zero - check the cell")
    return [norm(cross(a2, a3)) / abs(vol), norm(cross(a3, a1)) / abs(vol), norm(cross(a1, a2)) / abs(vol)]


# ----------------------------------------------------------------------------- template
def strip_comment(line):
    for i, ch in enumerate(line):
        if ch in "#!":
            return line[:i]
    return line


def scan_template(text):
    """Return (kpoint scheme line index, cell, periodic string, errors)."""
    lines = text.splitlines()
    stack, scheme_idx = [], None
    vec, abc, angles, periodic = {}, None, [90.0, 90.0, 90.0], "XYZ"
    problems = []
    for i, line in enumerate(lines):
        s = strip_comment(line).strip()
        if not s:
            continue
        if s.startswith("&"):
            tok = s[1:].split()
            if tok and tok[0].upper() == "END":
                if stack:
                    stack.pop()
            elif tok:
                stack.append(tok[0].upper())
            continue
        tok = s.split()
        key = tok[0].upper()
        if stack and stack[-1] == "KPOINTS" and key == "SCHEME" and "DFT" in stack and scheme_idx is None:
            scheme_idx = i
        if len(stack) >= 2 and stack[-1] == "CELL" and stack[-2] == "SUBSYS":
            if "${" in s or s.startswith("$"):
                problems.append(f"cell line uses a variable: {s}")
                continue
            vals = tok[1:]
            fac = 1.0
            if vals and vals[0].startswith("["):
                unit = vals[0].strip("[]").lower()
                fac = {"angstrom": 1.0, "bohr": BOHR_ANG, "nm": 10.0}.get(unit)
                if fac is None and key != "ALPHA_BETA_GAMMA":
                    problems.append(f"unsupported unit [{unit}]")
                    fac = 1.0
                vals = vals[1:]
            try:
                if key in ("A", "B", "C"):
                    vec[key] = [float(x) * fac for x in vals[:3]]
                elif key == "ABC":
                    abc = [float(x) * fac for x in vals[:3]]
                elif key == "ALPHA_BETA_GAMMA":
                    angles = [float(x) for x in vals[:3]]
                elif key == "PERIODIC":
                    periodic = vals[0].upper() if vals else "XYZ"
                elif key == "CELL_FILE_NAME":
                    problems.append("cell given by CELL_FILE_NAME")
            except ValueError:
                problems.append(f"cannot read cell line: {s}")
    cell = None
    if len(vec) == 3:
        cell = [vec["A"], vec["B"], vec["C"]]
    elif abc:
        cell = cell_from_abc(abc, angles)
    return lines, scheme_idx, cell, periodic, problems


def get_cell(args, text):
    lines, scheme_idx, cell, periodic, problems = scan_template(text)
    if args.cell:
        if len(args.cell) == 3:
            cell = cell_from_abc(args.cell, [90.0, 90.0, 90.0])
        elif len(args.cell) == 6:
            cell = cell_from_abc(args.cell[:3], args.cell[3:])
        else:
            raise SystemExit("ERROR: --cell needs 3 (a b c) or 6 (a b c alpha beta gamma) numbers")
    elif args.structure:
        try:
            from ase.io import read
        except ImportError:
            raise SystemExit("ERROR: --structure needs ASE ('pip install ase'); or give --cell a b c [al be ga]")
        if not os.path.isfile(args.structure):
            raise SystemExit(f"ERROR: structure file not found: {args.structure}")
        atoms = read(args.structure)
        if atoms.cell.rank < 3:
            raise SystemExit(f"ERROR: {args.structure} has no 3D cell")
        cell = [list(map(float, v)) for v in atoms.cell[:]]
    elif cell is None or problems:
        why = "; ".join(problems) if problems else "no &SUBSYS/&CELL found"
        raise SystemExit(f"ERROR: cannot determine the cell from the template ({why}); use --cell or --structure")
    if args.nonperiodic is not None:
        periodic = "".join(a for a in AXES if a not in args.nonperiodic.upper()) or "NONE"
    return lines, scheme_idx, cell, periodic


def meshes(args, cell, periodic):
    rl = recip_lengths(cell)
    per = [ax in periodic for ax in AXES] if periodic != "NONE" else [False] * 3
    out = []
    if args.grids:
        for g in args.grids:
            try:
                n = [int(x) for x in g.lower().split("x")]
            except ValueError:
                n = []
            if len(n) != 3 or min(n) < 1:
                raise SystemExit(f"ERROR: bad mesh '{g}' (use e.g. 4x4x4)")
            out.append((tuple(n), f"grid {g}"))
    elif args.kspacing:
        for ks in args.kspacing:
            n = [max(1, math.ceil(2 * math.pi * rl[i] / ks - 1e-6)) if per[i] else 1 for i in range(3)]
            out.append((tuple(n), f"kspacing {ks:g} 1/A"))
    else:
        for L in args.lengths:
            n = [max(1, math.ceil(L * rl[i] - 1e-6)) if per[i] else 1 for i in range(3)]
            out.append((tuple(n), f"length {L:g} A"))
    uniq, seen = [], set()
    for n, why in out:
        if n not in seen:
            seen.add(n)
            uniq.append((n, why))
    if args.even_only:
        uniq = [(tuple(x + 1 if x > 1 and x % 2 else x for x in n), w) for n, w in uniq]
        seen, tmp = set(), []
        for n, w in uniq:
            if n not in seen:
                seen.add(n)
                tmp.append((n, w))
        uniq = tmp
    return sorted(uniq, key=lambda t: (t[0][0] * t[0][1] * t[0][2], t[0])), rl, per


def mesh_name(n):
    return f"kp_{n[0]}x{n[1]}x{n[2]}"


def cmd_mesh(args, quiet=False):
    if not os.path.isfile(args.template):
        raise SystemExit(f"ERROR: template not found: {args.template}")
    with open(args.template, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    lines, scheme_idx, cell, periodic = get_cell(args, text)
    ms, rl, per = meshes(args, cell, periodic)
    if not quiet:
        lens = [norm(v) for v in cell]
        print("cell lengths |a| |b| |c| = {:.4f} {:.4f} {:.4f} A; periodic: {}".format(*lens, periodic))
        print("reciprocal |b_i| (no 2pi) = {:.4f} {:.4f} {:.4f} 1/A".format(*rl))
        print(f"{'mesh':>12} {'n1*n2*n3':>9}  {'k-spacing (2pi/A, per axis)':<30} from")
        for n, why in ms:
            ksp = " ".join(f"{2 * math.pi * rl[i] / n[i]:.3f}" if per[i] else "  -  " for i in range(3))
            print(f"{'x'.join(map(str, n)):>12} {n[0] * n[1] * n[2]:>9}  {ksp:<30} {why}")
    return lines, scheme_idx, ms


# ----------------------------------------------------------------------------- gen / run / submit
def cmd_gen(args):
    lines, scheme_idx, ms = cmd_mesh(args)
    if scheme_idx is None:
        raise SystemExit(f"ERROR: {args.template}: no SCHEME line inside &FORCE_EVAL/&DFT/&KPOINTS")
    for f in args.extra:
        if not os.path.isfile(f):
            raise SystemExit(f"ERROR: extra file not found: {f}")
    indent = re.match(r"\s*", lines[scheme_idx]).group(0)
    dirs = []
    for n, _ in ms:
        d = mesh_name(n)
        if os.path.exists(os.path.join(d, args.out_name)) and not args.force:
            print(f"skip {d}: {args.out_name} exists (use --force)")
            continue
        os.makedirs(d, exist_ok=True)
        new = list(lines)
        new[scheme_idx] = f"{indent}SCHEME MONKHORST-PACK {n[0]} {n[1]} {n[2]}"
        with open(os.path.join(d, args.inp_name), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(new) + "\n")
        for f in args.extra:
            shutil.copy2(f, d)
        dirs.append(d)
        print(f"written {d}/{args.inp_name}")
    return dirs


def cmd_run(args):
    dirs = cmd_gen(args)
    base = shlex.split(args.cp2k)
    if shutil.which(base[0]) is None:
        raise SystemExit(f"ERROR: '{base[0]}' not found in PATH (set --cp2k, load the CP2K environment first)")
    procs = []
    for d in dirs:
        cmd = base + ["-i", args.inp_name, "-o", args.out_name]
        print(f"[{time.strftime('%H:%M:%S')}] running in {d}: {' '.join(cmd)}")
        procs.append((d, subprocess.Popen(cmd, cwd=d)))
        while len([p for _, p in procs if p.poll() is None]) >= args.npar:
            time.sleep(2)
    failed = [d for d, p in procs if p.wait() != 0]
    if failed:
        print(f"WARNING: non-zero exit status in {', '.join(failed)}", file=sys.stderr)
    return cmd_collect(args)


def find_dirs():
    out = []
    for d in glob.glob("kp_*x*x*"):
        m = RE_DIR.match(os.path.basename(d))
        if m and os.path.isdir(d):
            n = tuple(int(x) for x in m.groups())
            out.append((n[0] * n[1] * n[2], n, d))
    return sorted(out)


def cmd_submit(args):
    if not os.path.isfile(args.job):
        raise SystemExit(f"ERROR: job script not found: {args.job}")
    dirs = find_dirs()
    if not dirs:
        raise SystemExit("ERROR: no kp_*x*x* directories (run 'gen' first)")
    submit = shlex.split(args.submit_cmd)
    if shutil.which(submit[0]) is None:
        raise SystemExit(f"ERROR: '{submit[0]}' not found (are you on the cluster login node?)")
    job = os.path.abspath(args.job)
    for _, _, d in dirs:
        if os.path.exists(os.path.join(d, args.out_name)) and not args.force:
            print(f"skip {d}: {args.out_name} exists")
            continue
        if os.path.basename(submit[0]) == "qsub":
            cmd = submit + ["-v", f"INP={args.inp_name},OUT={args.out_name}", job]
        else:
            cmd = submit + [f"--export=ALL,INP={args.inp_name},OUT={args.out_name}", job]
        print(f"{d}: {' '.join(cmd)}")
        subprocess.run(cmd, cwd=d, check=True)
    print(f"after all jobs have finished:  python {os.path.basename(sys.argv[0])} collect")
    return 0


# ----------------------------------------------------------------------------- collect
def parse_out(path):
    r = {"energy": None, "unit": None, "fermi": None, "pressure": None, "nkp": None, "natoms": None}
    ok = bad = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = RE_ENERGY.search(line)
            if m:
                r["energy"], r["unit"] = float(m.group(2)), m.group(1).strip().lower()
                continue
            m = RE_FERMI.match(line)
            if m:
                r["fermi"] = float(m.group(1)) * HARTREE_EV
                continue
            m = RE_STRESS.search(line)
            if m:
                r["pressure"] = float(m.group(1))
                continue
            m = RE_NKP.search(line)
            if m:
                r["nkp"] = int(m.group(1))
                continue
            if r["natoms"] is None:
                m = RE_NATOMS.match(line)
                if m:
                    r["natoms"] = int(m.group(1))
            if RE_SCF_OK.search(line):
                ok += 1
            elif RE_SCF_BAD.search(line):
                bad += 1
    r["scf"] = "NO" if bad else ("yes" if ok else "unknown")
    return r


def cmd_collect(args):
    dirs = find_dirs()
    if not dirs:
        raise SystemExit(f"ERROR: no kp_*x*x* directories in {os.getcwd()}")
    rows = []
    for prod, n, d in dirs:
        f = os.path.join(d, args.out_name)
        if not os.path.isfile(f):
            print(f"skip {d}: no {args.out_name}", file=sys.stderr)
            continue
        r = parse_out(f)
        if r["energy"] is None:
            print(f"skip {d}: no final energy (SCF aborted or job unfinished)", file=sys.stderr)
            continue
        if r["unit"] not in ("hartree", "a.u."):
            raise SystemExit(f"ERROR: {f}: energy printed in [{r['unit']}]; use the default hartree")
        r.update(mesh=n, prod=prod, dir=d)
        rows.append(r)
    if not rows:
        raise SystemExit("ERROR: no usable output found")
    natoms = args.natoms or rows[-1]["natoms"]
    if not natoms:
        raise SystemExit("ERROR: number of atoms unknown; give --natoms")
    ref = rows[-1]["energy"]  # densest mesh
    for r in rows:
        r["e_atom_eV"] = r["energy"] * HARTREE_EV / natoms
        r["dE"] = (r["energy"] - ref) * HARTREE_EV * 1000.0 / natoms
    csv_name = args.csv or "kpoint_scan.csv"
    with open(csv_name, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["mesh", "n1", "n2", "n3", "nkpts", "energy_Ha", "energy_per_atom_eV",
                    "dE_meV_per_atom", "fermi_eV", "pressure_bar", "scf_converged"])
        for r in rows:
            n = r["mesh"]
            w.writerow(["x".join(map(str, n)), n[0], n[1], n[2], "" if r["nkp"] is None else r["nkp"],
                        f"{r['energy']:.10f}", f"{r['e_atom_eV']:.6f}", f"{r['dE']:.4f}",
                        "" if r["fermi"] is None else f"{r['fermi']:.4f}",
                        "" if r["pressure"] is None else f"{r['pressure']:.4f}", r["scf"]])
    print(f"{'mesh':>10} {'nkpts':>6} {'energy_Ha':>18} {'dE_meV/atom':>12} {'E_F_eV':>9} {'P_bar':>11} {'scf':>4}")
    for r in rows:
        print(f"{'x'.join(map(str, r['mesh'])):>10} {r['nkp'] if r['nkp'] is not None else '-':>6} "
              f"{r['energy']:>18.10f} {r['dE']:>12.4f} "
              f"{'-' if r['fermi'] is None else format(r['fermi'], '.4f'):>9} "
              f"{'-' if r['pressure'] is None else format(r['pressure'], '.2f'):>11} {r['scf']:>4}")
    print(f"written: {csv_name}  (natoms = {natoms}, reference = densest mesh {'x'.join(map(str, rows[-1]['mesh']))})")
    good = [abs(r["dE"]) < args.threshold for r in rows]
    for i, r in enumerate(rows[:-1]):
        if all(good[i:]):
            print(f"|dE| < {args.threshold:g} meV/atom from mesh {'x'.join(map(str, r['mesh']))} on "
                  "(relative to the densest mesh scanned)")
            break
    else:
        print(f"no mesh converged to {args.threshold:g} meV/atom yet - add denser meshes")
    if any(r["scf"] == "NO" for r in rows):
        print("WARNING: some meshes have SCF NOT converged - do not use them")
    if not args.no_plot:
        plot(rows, args)
    return 0


def plot(rows, args):
    try:
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        print(f"WARNING: no plot ({exc}); install numpy + matplotlib or use --no-plot", file=sys.stderr)
        return
    has_f = any(r["fermi"] is not None for r in rows)
    has_p = any(r["pressure"] is not None for r in rows)
    npan = 1 + has_f + has_p
    fig, axes = plt.subplots(npan, 1, figsize=(5.8, 2.4 * npan + 0.5), sharex=True)
    axes = np.atleast_1d(axes)
    x = np.arange(len(rows))
    labels = ["x".join(map(str, r["mesh"])) for r in rows]
    de = np.abs(np.array([r["dE"] for r in rows]))
    m = de > 0
    axes[0].semilogy(x[m], de[m], "o-", color="#1f5aa6")
    axes[0].axhline(args.threshold, ls="--", color="0.4", lw=1, label=f"{args.threshold:g} meV/atom")
    axes[0].set_ylabel("|E - E_dense| (meV/atom)")
    axes[0].legend(frameon=False, fontsize=8)
    k = 1
    if has_f:
        xs = [i for i, r in enumerate(rows) if r["fermi"] is not None]
        axes[k].plot(xs, [rows[i]["fermi"] for i in xs], "o-", color="#2a9d5c")
        axes[k].set_ylabel("Fermi energy (eV)")
        k += 1
    if has_p:
        xs = [i for i, r in enumerate(rows) if r["pressure"] is not None]
        axes[k].plot(xs, [rows[i]["pressure"] for i in xs], "o-", color="#d9822b")
        axes[k].set_ylabel("STRESS 1/3 Trace (bar)")
    for a in axes:
        a.grid(alpha=0.3, which="both")
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    axes[-1].set_xlabel("Monkhorst-Pack mesh")
    fig.tight_layout()
    if args.show:
        plt.show()
    else:
        out = args.outfile or "kpoint_scan.png"
        fig.savefig(out, dpi=args.dpi)
        print(f"saved {out}")


# ----------------------------------------------------------------------------- CLI
def build_parser():
    p = argparse.ArgumentParser(
        description="k-point convergence test for CP2K: build Monkhorst-Pack meshes from the cell, "
                    "generate/run/submit kp_<n1>x<n2>x<n3>/ directories, collect energy per atom, "
                    "Fermi level and pressure, and plot.",
        epilog="Run 'python conv_kpoints.py <command> -h' for the options of each command.\n"
               "Example: python conv_kpoints.py run -t template.inp --lengths 10 20 30 40 "
               "--cp2k \"mpirun -np 8 cp2k.psmp\"",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    io = argparse.ArgumentParser(add_help=False)
    io.add_argument("--inp-name", default="input.inp", help="input file name in each directory (default input.inp)")
    io.add_argument("--out-name", default="output.out", help="output file name in each directory (default output.out)")
    io.add_argument("--force", action="store_true", help="regenerate / resubmit meshes that already have output")

    mesh = argparse.ArgumentParser(add_help=False)
    mesh.add_argument("-t", "--template", required=True, help="template input with &KPOINTS/SCHEME")
    g = mesh.add_mutually_exclusive_group()
    g.add_argument("--lengths", type=float, nargs="+", default=[10, 15, 20, 25, 30, 40, 50],
                   help="length values R in Angstrom: n_i = max(1, ceil(R*|b_i|)) (default 10 15 20 25 30 40 50)")
    g.add_argument("--kspacing", type=float, nargs="+",
                   help="k-spacings in 1/Angstrom incl. 2*pi: n_i = max(1, ceil(2*pi*|b_i|/kspacing))")
    g.add_argument("--grids", nargs="+", help="explicit meshes, e.g. 2x2x2 4x4x4 6x6x6")
    mesh.add_argument("--cell", type=float, nargs="+", metavar="X",
                      help="cell as a b c [alpha beta gamma] (Angstrom, degrees) instead of reading the template")
    mesh.add_argument("--structure", help="read the cell from a structure file (POSCAR, cif, extxyz ...; needs ase)")
    mesh.add_argument("--nonperiodic", metavar="AXES",
                      help="axes with n = 1, e.g. Z for a slab (default: from &CELL PERIODIC)")
    mesh.add_argument("--even-only", action="store_true",
                      help="round odd n > 1 up to the next even number")
    gen = argparse.ArgumentParser(add_help=False)
    gen.add_argument("--extra", nargs="*", default=[], help="files copied into every directory (xyz ...)")

    col = argparse.ArgumentParser(add_help=False)
    col.add_argument("--natoms", type=int, help="atoms in the cell (default: '- Atoms:' from the output)")
    col.add_argument("--threshold", type=float, default=1.0,
                     help="convergence threshold in meV/atom (default 1.0)")
    col.add_argument("--csv", help="CSV file name (default kpoint_scan.csv)")
    col.add_argument("-o", "--outfile", help="figure file name (default kpoint_scan.png)")
    col.add_argument("--dpi", type=int, default=200, help="figure resolution (default 200)")
    col.add_argument("--no-plot", action="store_true", help="skip the figure")
    col.add_argument("--show", action="store_true", help="show the figure interactively instead of saving")

    sub.add_parser("mesh", parents=[mesh], help="only print the meshes")
    sub.add_parser("gen", parents=[io, mesh, gen], help="generate kp_* directories")
    r = sub.add_parser("run", parents=[io, mesh, gen, col], help="generate, run locally and collect")
    r.add_argument("--cp2k", default="cp2k.psmp",
                   help='CP2K command incl. MPI launcher (default "cp2k.psmp", e.g. "mpirun -np 8 cp2k.psmp")')
    r.add_argument("--npar", type=int, default=1, help="number of meshes run at the same time (default 1)")
    s = sub.add_parser("submit", parents=[io], help="submit every kp_* directory as a batch job")
    s.add_argument("--job", required=True, help="job script (reads INP/OUT), e.g. job_cp2k_slurm.sh")
    s.add_argument("--submit-cmd", default="sbatch", help='"sbatch" (default) or "qsub"')
    sub.add_parser("collect", parents=[io, col], help="parse outputs, write CSV and plot")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "mesh":
        cmd_mesh(args)
        return 0
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
