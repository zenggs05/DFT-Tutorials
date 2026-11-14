#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nmr_shift.py - convert Gaussian 16 GIAO isotropic shieldings into chemical shifts.

读取 Gaussian 16 NMR 计算（NMR=GIAO/CSGT/IGAIM）输出中的 "Isotropic =" 屏蔽常数，换算为化学位移：
  * 参比法：delta = sigma_ref - sigma，sigma_ref 来自同一水平计算的 TMS（--ref-log tms.log，自动对
    所有同种原子取平均）或手工给出（--sigma-ref H=31.80 --sigma-ref C=184.0）；
  * 线性回归缩放：delta = (intercept - sigma) / (-slope)（--scale H=-1.0381,31.7427，参数需与你的
    计算协议完全一致，例如 Cheshire 数据库中的条目）；
  * 等价核平均（甲基 H、对称等价核）、多构象 Boltzmann 平均（按 SCF 能量、自由能或外部能量文件）；
  * 与实验值比较：每个元素的 MAE、RMSE、最大偏差，并可用 --fit 对 sigma vs delta_exp 做线性拟合。

Usage examples
--------------
  python nmr_shift.py camphor_nmr.log --ref-log tms_nmr.log
  python nmr_shift.py camphor_nmr.log --sigma-ref H=31.80 --sigma-ref C=184.00 --elements H C
  python nmr_shift.py camphor_nmr.log --ref-log tms_nmr.log --assign camphor_assign.txt
  python nmr_shift.py conf*.log --ref-log tms.log --assign assign.txt --weight gibbs --temp 298.15
  python nmr_shift.py conf*.log --ref-log tms.log --weight file --energies conf_G.csv
  python nmr_shift.py mol.log --scale H=-1.0381,31.7427 --equiv 12-14 --equiv 15-17
  python nmr_shift.py mol.log --ref-log tms.log --assign assign.txt --fit --csv shifts.csv

Assignment file (--assign), one group per line, '#' starts a comment:
    # label   atoms (Gaussian numbering)   [experimental shift / ppm]
    C1        7                            57.7
    H8        12-14                        0.84       <- fast-rotating methyl: averaged
    H3exo     26                           2.35       <- diastereotopic H: keep separate
  Atoms are written without spaces (12-14 or 3,5,9). All atoms of a group must be the same
  element; the group value is the arithmetic mean.

Energy file (--energies), CSV or whitespace separated: <log file name> <energy in Hartree>.

Notes
-----
* The last shielding block of each log is used (Opt + --Link1-- NMR jobs are fine).
* The reference must be computed at exactly the same level (geometry level, NMR functional/basis,
  solvent model).  The script compares the NMR routes of the reference and target logs and warns
  if they differ.
* Conformers must share the same atom numbering.  Boltzmann weights: w_i ~ exp(-(E_i - E_min)/RT).
* Exchangeable protons (OH/NH/COOH) are usually excluded from statistics - simply leave them out of
  the assignment file or give no experimental value.

Dependencies: Python >= 3.9 standard library (numpy is not required).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import csv
import math
import os
import re
import sys

HARTREE_TO_KCAL = 627.509474
R_KCAL = 1.987204259e-3          # kcal mol^-1 K^-1

RE_DASH = re.compile(r"^ -{10,}\s*$")
RE_ISO = re.compile(r"^\s*(\d+)\s+([A-Za-z]{1,2})\s+Isotropic\s*=\s*(\S+)\s+Anisotropy\s*=\s*(\S+)")
RE_SCF = re.compile(r"SCF Done:\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+)")
RE_GIBBS = re.compile(r"Sum of electronic and thermal Free Energies=\s*(-?\d+\.\d+)")
RE_NMR_KW = re.compile(r"(^|[\s#/,])nmr(?=$|[\s=(,])", re.IGNORECASE)
IGNORED_ROUTE = re.compile(r"^(#\w*|geom=.*|guess=.*|nosymm|pop=.*|iop\(.*|gfinput|gfprint|test|"
                           r"nmr(=.*)?|scf=.*|int(egral)?=.*|output=.*|formcheck.*)$", re.IGNORECASE)


# ---------------------------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------------------------
def parse_log(path):
    """Return dict(sigma={atom: (element, iso)}, scf, gibbs, nmr_route, normal)."""
    out = {"file": path, "sigma": {}, "scf": None, "gibbs": None, "nmr_route": None, "normal": 0,
           "error": False}
    block = None
    capturing, buf, prev, route = False, [], "", ""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if capturing:
                if RE_DASH.match(line):
                    capturing = False
                    route = "".join(buf).strip()
                else:
                    buf.append(line[1:] if line.startswith(" ") else line)
                prev = line
                continue
            if line.startswith(" #") and RE_DASH.match(prev):
                capturing, buf, prev = True, [line[1:]], line
                continue
            prev = line
            if "Magnetic shielding tensor" in line:
                block = {}
                out["sigma"] = block
                if RE_NMR_KW.search(route):
                    out["nmr_route"] = route
                continue
            if block is not None and "Isotropic" in line:
                m = RE_ISO.match(line)
                if m:
                    try:
                        iso = float(m.group(3))
                    except ValueError:
                        iso = float("nan")       # '**********' overflow
                    block[int(m.group(1))] = (m.group(2).capitalize(), iso)
                continue
            m = RE_SCF.search(line)
            if m:
                out["scf"] = float(m.group(2))
                continue
            m = RE_GIBBS.search(line)
            if m:
                out["gibbs"] = float(m.group(1))
                continue
            if "Normal termination of Gaussian" in line:
                out["normal"] += 1
            elif "Error termination" in line:
                out["error"] = True
    return out


def route_key(route):
    """Normalised set of route tokens that define the level of theory."""
    if not route:
        return None
    toks = []
    depth, cur = 0, ""
    for ch in route:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch.isspace() and depth == 0:
            if cur:
                toks.append(cur)
            cur = ""
        else:
            cur += ch
    if cur:
        toks.append(cur)
    return frozenset(t.lower() for t in toks if not IGNORED_ROUTE.match(t))


def parse_atoms(text):
    """'12-14,20' -> [12, 13, 14, 20]."""
    atoms = []
    for part in re.split(r"[,\s]+", text.strip()):
        if not part:
            continue
        m = re.fullmatch(r"(\d+)-(\d+)", part)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            if b < a:
                raise ValueError("bad atom range '%s'" % part)
            atoms += list(range(a, b + 1))
        elif re.fullmatch(r"\d+", part):
            atoms.append(int(part))
        else:
            raise ValueError("bad atom specification '%s'" % part)
    return atoms


def read_assign(path):
    """Return list of (label, atoms, exp or None)."""
    groups = []
    with open(path, "r", encoding="utf-8") as fh:
        for n, raw in enumerate(fh, 1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            tok = line.split()
            if len(tok) not in (2, 3):
                raise ValueError("%s:%d: expected '<label> <atoms> [exp]' with atoms written without "
                                 "spaces (e.g. 12-14 or 3,5)" % (path, n))
            exp = None
            atom_txt = tok[1]
            if len(tok) == 3:
                try:
                    exp = float(tok[2])
                except ValueError:
                    raise ValueError("%s:%d: experimental shift '%s' is not a number" % (path, n, tok[2]))
            try:
                groups.append((tok[0], parse_atoms(atom_txt), exp))
            except ValueError as exc:
                raise ValueError("%s:%d: %s" % (path, n, exc))
    if not groups:
        raise ValueError("%s contains no assignment lines" % path)
    return groups


def read_energies(path):
    table = {}
    with open(path, "r", encoding="utf-8") as fh:
        for n, raw in enumerate(fh, 1):
            line = raw.split("#", 1)[0].strip()
            if not line:
                continue
            tok = [t for t in re.split(r"[,\s;]+", line) if t]
            if len(tok) < 2:
                raise ValueError("%s:%d: need '<log> <energy>'" % (path, n))
            try:
                table[os.path.basename(tok[0])] = float(tok[1])
            except ValueError:
                if n == 1:
                    continue           # header line
                raise ValueError("%s:%d: energy is not a number" % (path, n))
    return table


def parse_pairs(items, what, nvals):
    """['H=31.8', 'C=184'] -> {'H': (31.8,)}; for --scale 'H=-1.04,31.7' -> {'H': (-1.04, 31.7)}."""
    out = {}
    for it in items:
        m = re.fullmatch(r"\s*([A-Za-z]{1,2})\s*[=:]\s*(.+)", it)
        if not m:
            raise ValueError("bad %s '%s'" % (what, it))
        try:
            vals = tuple(float(v) for v in re.split(r"[,\s]+", m.group(2).strip()) if v)
        except ValueError:
            raise ValueError("bad number in %s '%s'" % (what, it))
        if len(vals) != nvals:
            raise ValueError("%s '%s' needs %d value(s)" % (what, it, nvals))
        out[m.group(1).capitalize()] = vals
    return out


# ---------------------------------------------------------------------------------------------
# calculations
# ---------------------------------------------------------------------------------------------
def boltzmann(energies, temp):
    emin = min(energies)
    x = [math.exp(-(e - emin) * HARTREE_TO_KCAL / (R_KCAL * temp)) for e in energies]
    s = sum(x)
    return [v / s for v in x]


def linfit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    slope = sxy / sxx
    icpt = my - slope * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (slope * x + icpt)) ** 2 for x, y in zip(xs, ys))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return slope, icpt, r2


def stats(errors):
    n = len(errors)
    mae = sum(abs(e) for e in errors) / n
    rmse = math.sqrt(sum(e * e for e in errors) / n)
    return mae, rmse, max(errors, key=abs)


def build_parser():
    p = argparse.ArgumentParser(
        description="GIAO isotropic shieldings (Gaussian 16) -> chemical shifts: TMS reference or linear "
                    "scaling, equivalent-nucleus and Boltzmann averaging, comparison with experiment (MAE).  "
                    "屏蔽常数换算化学位移：参比/线性缩放、等价核与构象平均、与实验比较。",
        epilog="Examples:\n"
               "  python nmr_shift.py mol_nmr.log --ref-log tms_nmr.log\n"
               "  python nmr_shift.py mol_nmr.log --sigma-ref H=31.80 --sigma-ref C=184.0 --elements H C\n"
               "  python nmr_shift.py conf*.log --ref-log tms.log --assign assign.txt --weight gibbs\n"
               "  python nmr_shift.py mol.log --scale H=-1.0381,31.7427 --equiv 12-14\n"
               "Assignment file line: <label> <atoms e.g. 12-14,20> [exp. shift]",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("logs", nargs="+", help="NMR output(s); several logs = conformers of one molecule")
    g = p.add_argument_group("reference / 参比")
    g.add_argument("--ref-log", metavar="LOG", help="NMR output of the reference (e.g. TMS) at the same level")
    g.add_argument("--sigma-ref", action="append", default=[], metavar="EL=SIGMA",
                   help="reference shielding for an element, e.g. H=31.80 (overrides --ref-log)")
    g.add_argument("--scale", action="append", default=[], metavar="EL=SLOPE,INTERCEPT",
                   help="linear scaling delta=(intercept-sigma)/(-slope) for an element (overrides references)")
    a = p.add_argument_group("averaging / 平均")
    a.add_argument("--assign", metavar="FILE", help="assignment file: label, atoms, [experimental shift]")
    a.add_argument("--equiv", action="append", default=[], metavar="ATOMS",
                   help="average these atoms, e.g. 12-14 (repeatable; without --assign)")
    a.add_argument("--elements", nargs="+", metavar="EL", help="only list these elements (default: all)")
    a.add_argument("--weight", choices=["scf", "gibbs", "file", "equal"], default="scf",
                   help="Boltzmann weights from: last SCF energy (default), Gibbs free energy in the "
                        "same log, --energies file, or equal weights")
    a.add_argument("--energies", metavar="FILE", help="'<log name> <energy/Hartree>' per line (--weight file)")
    a.add_argument("--temp", type=float, default=298.15, help="temperature in K (default 298.15)")
    o = p.add_argument_group("output / 输出")
    o.add_argument("--fit", action="store_true",
                   help="fit sigma_calc = slope*delta_exp + intercept per element (needs experimental values)")
    o.add_argument("--csv", metavar="FILE", help="write the result table to CSV")
    o.add_argument("--show-sigma", action="store_true", help="also print the shielding of every conformer")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        sigma_ref = {k: v[0] for k, v in parse_pairs(args.sigma_ref, "--sigma-ref", 1).items()}
        scale = parse_pairs(args.scale, "--scale", 2)
        groups = read_assign(args.assign) if args.assign else None
    except (ValueError, OSError) as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 2
    for el, (slope, _) in scale.items():
        if slope == 0:
            print("ERROR: --scale %s: slope must not be 0" % el, file=sys.stderr)
            return 2

    confs = []
    for log in args.logs:
        if not os.path.isfile(log):
            print("ERROR: %s not found" % log, file=sys.stderr)
            return 1
        try:
            d = parse_log(log)
        except OSError as exc:
            print("ERROR: cannot read %s: %s" % (log, exc), file=sys.stderr)
            return 1
        if not d["sigma"]:
            print("ERROR: %s has no 'Magnetic shielding tensor' block (not an NMR job, or it failed)" % log,
                  file=sys.stderr)
            return 1
        if d["error"] or d["normal"] == 0:
            print("WARNING: %s did not end normally - shieldings may be incomplete" % log, file=sys.stderr)
        confs.append(d)
    natoms = {len(c["sigma"]) for c in confs}
    if len(natoms) > 1:
        print("ERROR: conformers have different numbers of atoms %s" % sorted(natoms), file=sys.stderr)
        return 1
    elems0 = {k: v[0] for k, v in confs[0]["sigma"].items()}
    for c in confs[1:]:
        if {k: v[0] for k, v in c["sigma"].items()} != elems0:
            print("ERROR: %s: atom order differs from %s" % (c["file"], confs[0]["file"]), file=sys.stderr)
            return 1

    # reference shieldings
    ref_info = ""
    if args.ref_log:
        if not os.path.isfile(args.ref_log):
            print("ERROR: %s not found" % args.ref_log, file=sys.stderr)
            return 1
        ref = parse_log(args.ref_log)
        if not ref["sigma"]:
            print("ERROR: %s has no shielding block" % args.ref_log, file=sys.stderr)
            return 1
        by_el = {}
        for el, iso in ref["sigma"].values():
            by_el.setdefault(el, []).append(iso)
        for el, vals in by_el.items():
            sigma_ref.setdefault(el, sum(vals) / len(vals))
        ref_info = "reference %s: " % args.ref_log + ", ".join(
            "%s %.4f (n=%d, spread %.3f)" % (el, sum(v) / len(v), len(v), max(v) - min(v))
            for el, v in sorted(by_el.items()))
        kr, kt = route_key(ref["nmr_route"]), route_key(confs[0]["nmr_route"])
        if kr is not None and kt is not None and kr != kt:
            print("WARNING: NMR routes differ between reference and target - the reference must be computed "
                  "at exactly the same level:\n  ref:    %s\n  target: %s" % (ref["nmr_route"], confs[0]["nmr_route"]),
                  file=sys.stderr)

    # Boltzmann weights
    if len(confs) == 1 or args.weight == "equal":
        weights = [1.0 / len(confs)] * len(confs)
        energies = [None] * len(confs)
    else:
        if args.weight == "file":
            if not args.energies:
                print("ERROR: --weight file needs --energies FILE", file=sys.stderr)
                return 2
            try:
                table = read_energies(args.energies)
            except (ValueError, OSError) as exc:
                print("ERROR: %s" % exc, file=sys.stderr)
                return 2
            energies = [table.get(os.path.basename(c["file"])) for c in confs]
        else:
            energies = [c[args.weight] for c in confs]
        missing = [c["file"] for c, e in zip(confs, energies) if e is None]
        if missing:
            print("ERROR: no %s energy for: %s" % (args.weight, ", ".join(missing)), file=sys.stderr)
            return 1
        weights = boltzmann(energies, args.temp)
        if args.weight == "scf":
            print("NOTE: weights from SCF energies of the NMR step; Gibbs free energies (--weight gibbs/file) "
                  "are usually preferred", file=sys.stderr)

    atoms_all = sorted(elems0)
    avg = {}
    for a in atoms_all:
        avg[a] = sum(w * c["sigma"][a][1] for w, c in zip(weights, confs))

    # groups
    if groups is None:
        used = set()
        groups = []
        for spec in args.equiv:
            try:
                at = parse_atoms(spec)
            except ValueError as exc:
                print("ERROR: --equiv: %s" % exc, file=sys.stderr)
                return 2
            groups.append(("%s%s" % (elems0.get(at[0], "?"), "+".join(map(str, at))), at, None))
            used.update(at)
        for a in atoms_all:
            if a not in used:
                groups.append(("%s%d" % (elems0[a], a), [a], None))
        groups.sort(key=lambda g: g[1][0])
    for label, at, _ in groups:
        bad = [x for x in at if x not in elems0]
        if bad:
            print("ERROR: group %s: atom(s) %s not in the molecule (1..%d)" % (label, bad, len(elems0)),
                  file=sys.stderr)
            return 1
        if len({elems0[x] for x in at}) > 1:
            print("ERROR: group %s mixes elements %s" % (label, sorted({elems0[x] for x in at})), file=sys.stderr)
            return 1

    # shifts
    rows = []
    missing_ref = set()
    for label, at, exp in groups:
        el = elems0[at[0]]
        if args.elements and el not in [e.capitalize() for e in args.elements]:
            continue
        sig = sum(avg[x] for x in at) / len(at)
        spread = max(avg[x] for x in at) - min(avg[x] for x in at)
        if el in scale:
            slope, icpt = scale[el]
            delta, method = (icpt - sig) / (-slope), "scale"
        elif el in sigma_ref:
            delta, method = sigma_ref[el] - sig, "ref"
        else:
            delta, method = None, "-"
            missing_ref.add(el)
        err = delta - exp if (delta is not None and exp is not None) else None
        rows.append({"label": label, "element": el, "atoms": ",".join(map(str, at)), "n": len(at),
                     "sigma": sig, "spread": spread, "delta": delta, "exp": exp, "error": err, "method": method})

    # print
    print("=" * 86)
    for c, w, e in zip(confs, weights, energies):
        extra = "  E=%.8f" % e if e is not None else ""
        print("conformer %-40s weight %6.3f%s" % (os.path.basename(c["file"]), w, extra))
    if ref_info:
        print(ref_info)
    if sigma_ref:
        print("sigma_ref used: " + ", ".join("%s=%.4f" % kv for kv in sorted(sigma_ref.items())))
    if scale:
        print("linear scaling: " + ", ".join("%s slope=%.4f intercept=%.4f" % (k, v[0], v[1])
                                             for k, v in sorted(scale.items())))
    print("-" * 86)
    print("%-12s %-3s %-14s %3s %10s %7s %9s %9s %8s" % ("label", "el", "atoms", "n", "sigma", "spread",
                                                           "delta", "exp", "error"))
    for r in rows:
        print("%-12s %-3s %-14s %3d %10.3f %7.3f %9s %9s %8s" % (
            r["label"][:12], r["element"], r["atoms"][:14], r["n"], r["sigma"], r["spread"],
            "%.2f" % r["delta"] if r["delta"] is not None else "-",
            "%.2f" % r["exp"] if r["exp"] is not None else "-",
            "%+.2f" % r["error"] if r["error"] is not None else "-"))
    if missing_ref:
        print("NOTE: no reference for element(s) %s - give --ref-log, --sigma-ref or --scale"
              % ", ".join(sorted(missing_ref)))
    if any(r["n"] > 1 and r["spread"] > 1.0 and r["element"] == "H" for r in rows):
        print("NOTE: some averaged H groups spread by > 1 ppm (normal for methyl H; check the grouping otherwise)")

    print("-" * 86)
    for el in sorted({r["element"] for r in rows}):
        er = [r["error"] for r in rows if r["element"] == el and r["error"] is not None]
        if er:
            mae, rmse, mx = stats(er)
            print("%-3s N=%3d  MAE=%.3f  RMSE=%.3f  max|err|=%.3f ppm" % (el, len(er), mae, rmse, abs(mx)))
        if args.fit:
            pts = [(r["exp"], r["sigma"]) for r in rows if r["element"] == el and r["exp"] is not None]
            if len(pts) >= 3:
                fit = linfit([p[0] for p in pts], [p[1] for p in pts])
                if fit and fit[0] != 0:
                    slope, icpt, r2 = fit
                    sc = [((icpt - s) / (-slope)) - e for e, s in pts]
                    mae, rmse, _ = stats(sc)
                    print("    fit %s: sigma = %.4f * delta_exp + %.4f (R^2=%.4f); scaled MAE=%.3f RMSE=%.3f "
                          "(in-sample - validate on other molecules)" % (el, slope, icpt, r2, mae, rmse))
            else:
                print("    fit %s: needs >= 3 experimental values" % el)
    if not any(r["error"] is not None for r in rows):
        print("(no experimental values: give them in the --assign file to get MAE)")

    if args.show_sigma:
        print("-" * 86)
        print("%-6s %-3s " % ("atom", "el") + " ".join("%12s" % os.path.basename(c["file"])[:12] for c in confs))
        for a in atoms_all:
            print("%-6d %-3s " % (a, elems0[a]) + " ".join("%12.4f" % c["sigma"][a][1] for c in confs))

    if args.csv:
        keys = ["label", "element", "atoms", "n", "sigma", "spread", "delta", "exp", "error", "method"]
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow({k: ("" if r[k] is None else (round(r[k], 4) if isinstance(r[k], float) else r[k]))
                            for k in keys})
        print("CSV written to %s" % args.csv, file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
