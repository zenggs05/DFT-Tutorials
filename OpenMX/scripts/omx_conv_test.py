#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_conv_test.py -- convergence tests for OpenMX: cutoff energy, k-grid, k-spacing, basis set.
                    截断能 / k 点 / 基组收敛测试驱动：由模板批量生成输入、依次调用 mpirun 运行
                    OpenMX、解析 .out 中的 Utot，输出 CSV 与 ΔE (meV/atom) 曲线，并给出满足阈值的
                    推荐值。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Modes
-----
  ecut      values are scf.energycutoff in Ry          e.g. 150 180 220 260 300 350
            (an active scf.Ngrid line is commented out, since it overrides the cutoff)
  kgrid     values are k-grids: "8" (= 8 8 8) or "8x8x1"  e.g. 4 6 8 10 12 14
  kspacing  values are k spacings in 1/Ang (2*pi included) e.g. 0.40 0.30 0.25 0.20 0.15;
            n_i = ceil(|b_i| / spacing), duplicated grids are removed
  basis     values are PAO strings (Si7.0-s2p1d1 Si7.0-s2p2d1 Si7.0-s3p3d2) for --species,
            or the levels quick / standard / precise of the DFT_DATA19 table
            (needs omx_make_input.py in the same folder) applied to all species
            (or only to --species)

How it works
------------
1. Every value gives one input <workdir>/<prefix>_<tag>.dat created from the
   template by editing the keyword (case-insensitive) or, if the template
   contains the placeholder @VALUE@, by replacing that text.  System.Name is set
   to <prefix>_<tag>; a relative DATA.PATH is made absolute (the runs are done
   inside <workdir>); level.of.fileout is set to 0 unless --keep-fileout.
2. Each input is run as   <mpirun> <openmx> <input> -nt <nt>  > <name>.std
   (OMP_NUM_THREADS = nt).  Finished runs are skipped unless --rerun.
   --dry-run only writes the inputs; --collect only parses existing outputs.
3. Utot, SCF convergence, grid and time are read from <name>.out
   (via omx_parse_out.py of the same folder when available).
4. dE = (Utot - Utot_ref) / Natoms in meV/atom; the reference is the most
   accurate value (largest cutoff, densest grid, last basis) unless --ref.
   The recommended value is the cheapest one from which on every more
   accurate value stays within --thresh (default 1 meV/atom) of the reference.

Usage examples
--------------
  python omx_conv_test.py ecut Si_scf.dat --values 120 150 180 220 260 300 350 400
  python omx_conv_test.py kgrid Si_scf.dat --values 4 6 8 10 12 14 16 --mpirun "mpirun -np 8" --nt 2
  python omx_conv_test.py kgrid slab.dat --values 4x4x1 6x6x1 8x8x1 10x10x1 12x12x1
  python omx_conv_test.py kspacing Si_scf.dat --values 0.4 0.3 0.25 0.2 0.15 0.12
  python omx_conv_test.py basis Si_scf.dat --species Si --values Si7.0-s2p1d1 Si7.0-s2p2d1 Si7.0-s3p3d2
  python omx_conv_test.py basis Si_scf.dat --values quick standard precise
  python omx_conv_test.py ecut Si_scf.dat --values 150 220 300 --dry-run      # inputs only
  python omx_conv_test.py ecut Si_scf.dat --values 150 220 300 --collect      # parse only

Environment: MPIRUN (default "mpirun -np 4") and OPENMX (default "openmx") are
used when --mpirun / --openmx are not given.

Dependencies: python >= 3.9, numpy, matplotlib (plot), an OpenMX executable.
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
import time

HA2MEV = 27211.386245988
FINISHED_MSG = "The calculation was normally finished."
HERE = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# optional sibling modules
# ---------------------------------------------------------------------------
def _sibling(name):
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    try:
        return __import__(name)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# OpenMX input editing (keywords are case-insensitive, '#' starts a comment)
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


def kw_comment_out(lines, key):
    low = key.lower()
    hit = False
    for i, line in enumerate(lines):
        tok = line.split("#", 1)[0].split()
        if tok and tok[0].lower() == low:
            lines[i] = "# " + line if line.endswith("\n") else "# " + line + "\n"
            hit = True
    return hit


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


def read_cell(lines):
    """Unit vectors in Angstrom (rows) or None."""
    rng = block_range(lines, "Atoms.UnitVectors")
    if rng is None:
        return None
    unit = (kw_get(lines, "Atoms.UnitVectors.Unit") or ["Ang"])[0].lower()
    rows = []
    for line in lines[rng[0] + 1: rng[1]]:
        tok = line.split("#", 1)[0].split()
        if len(tok) >= 3:
            rows.append([float(x) for x in tok[:3]])
    if len(rows) != 3:
        return None
    if unit == "au":
        rows = [[x * 0.529177210903 for x in r] for r in rows]
    return rows


def species_rows(lines):
    rng = block_range(lines, "Definition.of.Atomic.Species")
    if rng is None:
        sys.exit("ERROR: <Definition.of.Atomic.Species block not found in the template")
    out = []
    for i in range(rng[0] + 1, rng[1]):
        tok = lines[i].split("#", 1)[0].split()
        if len(tok) >= 3:
            out.append((i, tok))
    return out


# ---------------------------------------------------------------------------
# value handling
# ---------------------------------------------------------------------------
def parse_kgrid(text):
    parts = re.split(r"[x, ]+", text.strip().lower())
    try:
        nums = [int(p) for p in parts if p]
    except ValueError:
        sys.exit("ERROR: bad k-grid '%s' (use 8 or 8x8x1)" % text)
    if len(nums) == 1:
        nums = nums * 3
    if len(nums) != 3 or min(nums) < 1:
        sys.exit("ERROR: bad k-grid '%s' (use 8 or 8x8x1)" % text)
    return nums


def kgrid_from_spacing(cell, spacing):
    import numpy as np
    b = 2.0 * math.pi * np.linalg.inv(np.array(cell)).T
    return [max(1, int(math.ceil(np.linalg.norm(b[i]) / spacing - 1e-8))) for i in range(3)]


def safe_tag(text):
    return re.sub(r"[^A-Za-z0-9_-]", "p", text)


def make_cases(args, lines):
    """List of dicts: label, x (numeric for plots/sorting), edits."""
    cases = []
    if args.mode == "ecut":
        vals = sorted({float(v) for v in args.values})
        for v in vals:
            lab = ("%g" % v)
            cases.append({"label": lab, "x": v, "tag": "ecut" + safe_tag(lab), "ecut": v})
    elif args.mode == "kgrid":
        for v in args.values:
            g = parse_kgrid(v)
            lab = "x".join(map(str, g))
            cases.append({"label": lab, "x": g[0] * g[1] * g[2], "tag": "k" + lab, "kgrid": g})
    elif args.mode == "kspacing":
        cell = read_cell(lines)
        if cell is None:
            sys.exit("ERROR: kspacing mode needs <Atoms.UnitVectors in the template")
        seen = {}
        for v in sorted({float(x) for x in args.values}, reverse=True):
            g = kgrid_from_spacing(cell, v)
            key = tuple(g)
            if key in seen:
                print("NOTE: spacing %g gives the same grid %s as %g; skipped" % (v, "x".join(map(str, g)),
                                                                                seen[key]), file=sys.stderr)
                continue
            seen[key] = v
            lab = "%g" % v
            cases.append({"label": lab, "x": v, "tag": "ks" + safe_tag(lab), "kgrid": g})
    elif args.mode == "basis":
        spec = species_rows(lines)
        names = [t[0] for _, t in spec]
        if args.species and args.species not in names:
            sys.exit("ERROR: species '%s' not in Definition.of.Atomic.Species (%s)" % (args.species, names))
        mk = None
        for n, v in enumerate(args.values):
            if v.lower() in ("quick", "standard", "precise"):
                mk = mk or _sibling("omx_make_input")
                if mk is None:
                    sys.exit("ERROR: level names need omx_make_input.py next to this script")
                paos = {}
                for _, t in spec:
                    if args.species and t[0] != args.species:
                        continue
                    pao = mk.recommended_basis(t[2], v.lower())
                    if pao is None:
                        sys.exit("ERROR: VPS %s of species %s is not in the DFT_DATA19 table" % (t[2], t[0]))
                    paos[t[0]] = pao
            else:
                if not args.species:
                    if len(names) != 1:
                        sys.exit("ERROR: PAO strings need --species (template has %s)" % names)
                    target = names[0]
                else:
                    target = args.species
                paos = {target: v}
            cases.append({"label": v, "x": n, "tag": "b" + safe_tag(v), "pao": paos})
    if len(cases) < 2:
        sys.exit("ERROR: need at least two different values")
    return cases


def build_input(template_lines, case, args, prefix, data_path_abs):
    lines = list(template_lines)
    name = "%s_%s" % (prefix, case["tag"])
    text = "".join(lines)
    if "@VALUE@" in text and args.mode in ("ecut", "kgrid", "kspacing"):
        val = ("%g" % case["ecut"]) if "ecut" in case else " ".join(map(str, case["kgrid"]))
        lines = text.replace("@VALUE@", val).splitlines(True)
    else:
        if "ecut" in case:
            kw_set(lines, "scf.energycutoff", "%g" % case["ecut"])
            if kw_comment_out(lines, "scf.Ngrid"):
                print("NOTE: scf.Ngrid commented out in %s (it overrides scf.energycutoff)" % name,
                      file=sys.stderr)
        if "kgrid" in case:
            kw_set(lines, "scf.Kgrid", " ".join(map(str, case["kgrid"])))
        if "pao" in case:
            for i, tok in species_rows(lines):
                if tok[0] in case["pao"]:
                    lines[i] = "  %-6s %-22s %s\n" % (tok[0], case["pao"][tok[0]], tok[2])
    kw_set(lines, "System.Name", name)
    if data_path_abs:
        kw_set(lines, "DATA.PATH", data_path_abs)
    if not args.keep_fileout:
        kw_set(lines, "level.of.fileout", "0")
    return name, lines


# ---------------------------------------------------------------------------
# running and parsing
# ---------------------------------------------------------------------------
def run_one(workdir, name, args):
    dat = name + ".dat"
    cmd = shlex.split(args.mpirun) + [args.openmx, dat, "-nt", str(args.nt)]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(args.nt)
    print(">>> (%s) %s" % (workdir, " ".join(cmd)), flush=True)
    t0 = time.time()
    try:
        with open(os.path.join(workdir, name + ".std"), "w") as fh:
            rc = subprocess.run(cmd, cwd=workdir, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode
    except FileNotFoundError as exc:
        sys.exit("ERROR: cannot start '%s': %s (check --mpirun/--openmx or load your modules)"
                 % (cmd[0], exc))
    print("    return code %d, %.1f s" % (rc, time.time() - t0), flush=True)
    return rc


def finished(workdir, name):
    std = os.path.join(workdir, name + ".std")
    if not os.path.isfile(std):
        return False
    with open(std, errors="replace") as fh:
        return FINISHED_MSG in fh.read()


def collect(workdir, name, parser_mod):
    out = os.path.join(workdir, name + ".out")
    row = {"utot": None, "conv": None, "grid": "", "time": None, "finished": finished(workdir, name)}
    if not os.path.isfile(out):
        return row
    if parser_mod is not None:
        try:
            r = parser_mod.parse_out(out)
        except Exception as exc:  # malformed / truncated output
            print("WARNING: could not parse %s: %s" % (out, exc), file=sys.stderr)
            return row
        row["utot"] = r.get("utot_Ha")
        row["conv"] = r.get("scf_converged")
        row["grid"] = "x".join(map(str, r.get("ngrid", []))) if r.get("ngrid") else ""
        row["time"] = r.get("elapsed_s")
    else:
        with open(out, errors="replace") as fh:
            text = fh.read()
        m = re.findall(r"^\s*Utot\.\s+([-+0-9.eE]+)\s*$", text, re.M)
        row["utot"] = float(m[-1]) if m else None
        m = re.search(r"Num\. of grids of a-, b-, and c-axes\s*=\s*(\d+),\s*(\d+),\s*(\d+)", text)
        row["grid"] = "x".join(m.groups()) if m else ""
        m = re.findall(r"Elapsed\.Time\.\s+([0-9.]+)", text)
        row["time"] = float(m[-1]) if m else None
    return row


def recommend(rows, thresh):
    """Cheapest index from which all following entries are within thresh (ref excluded)."""
    ok = [r["dE"] is not None and abs(r["dE"]) < thresh and r["conv"] is not False for r in rows]
    best = None
    for i in range(len(rows)):
        tail = [ok[j] for j in range(i, len(rows)) if not rows[j]["is_ref"]]
        if tail and all(tail):
            best = i
            break
    return best


def plot(rows, args, png):
    try:
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not installed; no plot", file=sys.stderr)
        return
    good = [r for r in rows if r["dE"] is not None]
    if len(good) < 2:
        print("WARNING: fewer than two finished runs; no plot", file=sys.stderr)
        return
    numeric = args.mode in ("ecut", "kspacing")
    xs = [r["x"] for r in good] if numeric else list(range(len(good)))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True)
    ax1.plot(xs, [r["dE"] for r in good], "o-", color="#1f77b4")
    ax1.axhspan(-args.thresh, args.thresh, color="#2ca02c", alpha=0.15, lw=0,
                label="+/- %g meV/atom" % args.thresh)
    ax1.set_ylabel("E - E(ref)  (meV/atom)")
    ax1.legend(frameon=False)
    absd = [max(abs(r["dE"]), 1e-4) for r in good if not r["is_ref"]]
    xa = [x for x, r in zip(xs, good) if not r["is_ref"]]
    ax2.semilogy(xa, absd, "s-", color="#d62728")
    ax2.axhline(args.thresh, color="#2ca02c", ls="--", lw=1)
    ax2.set_ylabel("|E - E(ref)|  (meV/atom)")
    xlabel = {"ecut": "scf.energycutoff (Ry)", "kgrid": "scf.Kgrid",
              "kspacing": "k spacing (1/Ang, 2pi included)", "basis": "basis"}[args.mode]
    ax2.set_xlabel(xlabel)
    if not numeric:
        ax2.set_xticks(xs)
        ax2.set_xticklabels([r["label"] for r in good], rotation=30, ha="right")
    if args.mode == "kspacing":
        ax2.invert_xaxis()
    fig.suptitle("OpenMX %s convergence (%s)" % (args.mode, os.path.basename(args.template)))
    fig.tight_layout()
    fig.savefig(png, dpi=150)
    print("wrote %s" % png)
    if args.show:
        plt.show()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Run and analyse OpenMX convergence tests (scf.energycutoff, scf.Kgrid, "
                    "k spacing, PAO basis) from one template input.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python omx_conv_test.py ecut Si_scf.dat --values 120 150 180 220 260 300 350\n"
               "  python omx_conv_test.py kgrid Si_scf.dat --values 4 6 8 10 12 14 --nt 2\n"
               "  python omx_conv_test.py kspacing Si_scf.dat --values 0.4 0.3 0.2 0.15\n"
               "  python omx_conv_test.py basis Si_scf.dat --values quick standard precise\n"
               "  python omx_conv_test.py ecut Si_scf.dat --values 150 220 300 --collect")
    p.add_argument("mode", choices=["ecut", "kgrid", "kspacing", "basis"], help="quantity to scan")
    p.add_argument("template", help="template OpenMX input (.dat), optionally containing @VALUE@")
    p.add_argument("--values", nargs="+", required=True, help="values to test (see the mode list in -h)")
    p.add_argument("--species", default=None, help="basis mode: species label whose PAO is changed")
    p.add_argument("--workdir", default=None, help="directory for inputs/outputs (default conv_<mode>)")
    p.add_argument("--prefix", default=None, help="System.Name prefix (default: System.Name of the template)")
    p.add_argument("--mpirun", default=os.environ.get("MPIRUN", "mpirun -np 4"),
                   help="MPI launcher with options (default $MPIRUN or 'mpirun -np 4'); '' for serial")
    p.add_argument("--openmx", default=os.environ.get("OPENMX", "openmx"),
                   help="OpenMX executable (default $OPENMX or 'openmx')")
    p.add_argument("--nt", type=int, default=1, help="OpenMP threads per MPI process, passed as -nt (default 1)")
    p.add_argument("--data-path", default=None, help="override DATA.PATH in all generated inputs")
    p.add_argument("--keep-fileout", action="store_true", help="keep level.of.fileout of the template")
    p.add_argument("--dry-run", action="store_true", help="only write the inputs")
    p.add_argument("--collect", action="store_true", help="only parse existing outputs (no runs)")
    p.add_argument("--rerun", action="store_true", help="run again even if a finished output exists")
    p.add_argument("--thresh", type=float, default=1.0, help="tolerance in meV/atom (default 1.0)")
    p.add_argument("--ref", default=None,
                   help="label of the reference value (default: the most accurate = last after sorting)")
    p.add_argument("--natoms", type=int, default=None, help="atoms per cell (default: Atoms.Number)")
    p.add_argument("--csv", default=None, help="CSV file (default <workdir>/<mode>_conv.csv)")
    p.add_argument("--plot", default=None, help="PNG file (default <workdir>/<mode>_conv.png)")
    p.add_argument("--no-plot", action="store_true", help="do not draw the plot")
    p.add_argument("--show", action="store_true", help="show the plot window (interactive backend)")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if not os.path.isfile(args.template):
        sys.exit("ERROR: template not found: %s" % args.template)
    with open(args.template, encoding="utf-8-sig", errors="replace") as fh:
        tlines = fh.read().replace("\r\n", "\n").splitlines(True)
    if tlines and not tlines[-1].endswith("\n"):
        tlines[-1] += "\n"

    natoms = args.natoms
    if natoms is None:
        v = kw_get(tlines, "Atoms.Number")
        try:
            natoms = int(float(v[0])) if v else 0
        except ValueError:
            natoms = 0
    if not natoms or natoms < 1:
        sys.exit("ERROR: Atoms.Number not found in the template; give --natoms")
    if args.mode in ("kgrid", "kspacing"):
        solver = (kw_get(tlines, "scf.EigenvalueSolver") or ["?"])[0].lower()
        if solver != "band":
            print("WARNING: scf.EigenvalueSolver is '%s'; scf.Kgrid is only used by Band" % solver,
                  file=sys.stderr)

    prefix = args.prefix or (kw_get(tlines, "System.Name") or ["omx"])[0]
    workdir = args.workdir or ("conv_" + args.mode)
    os.makedirs(workdir, exist_ok=True)

    dp = args.data_path
    if dp is None:
        cur = kw_get(tlines, "DATA.PATH")
        if cur:
            dp = cur[0]
        else:
            dp = "../DFT_DATA19"          # OpenMX default, relative to the run directory
            print("NOTE: DATA.PATH not set in the template; OpenMX would use ../DFT_DATA19", file=sys.stderr)
    if not os.path.isabs(dp):
        dp = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.template)), dp))
    if not os.path.isdir(os.path.join(dp, "VPS")) and not args.dry_run:
        print("WARNING: %s does not contain VPS/; check DATA.PATH" % dp, file=sys.stderr)

    cases = make_cases(args, tlines)
    names = []
    for c in cases:
        name, lines = build_input(tlines, c, args, prefix, dp)
        with open(os.path.join(workdir, name + ".dat"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write("".join(lines))
        names.append(name)
    print("wrote %d inputs to %s/" % (len(names), workdir))
    if args.dry_run:
        return 0

    if not args.collect:
        for name in names:
            if finished(workdir, name) and not args.rerun:
                print("skip %s (finished output exists; --rerun to repeat)" % name)
                continue
            run_one(workdir, name, args)

    parser_mod = _sibling("omx_parse_out")
    rows = []
    for c, name in zip(cases, names):
        r = collect(workdir, name, parser_mod)
        r.update({"label": c["label"], "x": c["x"], "name": name})
        rows.append(r)
    done = [r for r in rows if r["utot"] is not None]
    if not done:
        sys.exit("ERROR: no Utot found in any output of %s/" % workdir)

    if args.ref is not None:
        refs = [r for r in rows if r["label"] == args.ref]
        if not refs:
            sys.exit("ERROR: --ref %s is not one of %s" % (args.ref, [r["label"] for r in rows]))
        ref = refs[0]
    else:
        ref = rows[-1]
    if ref["utot"] is None:
        sys.exit("ERROR: the reference run %s has no Utot (not finished?)" % ref["name"])
    for r in rows:
        r["is_ref"] = r is ref
        r["dE"] = None if r["utot"] is None else (r["utot"] - ref["utot"]) * HA2MEV / natoms

    best = recommend(rows, args.thresh)
    csv_path = args.csv or os.path.join(workdir, "%s_conv.csv" % args.mode)
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["value", "name", "Utot_Ha", "dE_meV_per_atom", "scf_converged", "finished",
                    "grid", "elapsed_s", "reference"])
        for r in rows:
            w.writerow([r["label"], r["name"], r["utot"], r["dE"], r["conv"], r["finished"],
                        r["grid"], r["time"], "yes" if r["is_ref"] else ""])
    print("\n%-14s %20s %14s %8s %10s %12s" % (args.mode, "Utot (Ha)", "dE (meV/at)", "SCF", "grid", "time (s)"))
    for r in rows:
        print("%-14s %20s %14s %8s %10s %12s%s" % (
            r["label"], "%.10f" % r["utot"] if r["utot"] is not None else "missing",
            "%.3f" % r["dE"] if r["dE"] is not None else "-",
            {True: "ok", False: "NOCONV", None: "?"}[r["conv"]], r["grid"] or "-",
            "%.1f" % r["time"] if r["time"] else "-", "   <- reference" if r["is_ref"] else ""))
    print("wrote %s" % csv_path)
    if best is not None:
        print("\nRecommended (|dE| < %g meV/atom from here on): %s %s"
              % (args.thresh, args.mode, rows[best]["label"]))
    else:
        print("\nNo value meets |dE| < %g meV/atom for all more accurate settings; extend the range."
              % args.thresh)
    if any(r["conv"] is False for r in rows):
        print("WARNING: some runs did not reach SCF convergence; their energies are unreliable.")
    if not args.no_plot:
        plot(rows, args, args.plot or os.path.join(workdir, "%s_conv.png" % args.mode))
    return 0


if __name__ == "__main__":
    sys.exit(main())
