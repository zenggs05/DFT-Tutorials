#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
log2xyz.py - extract geometries from Gaussian 16 opt / scan / IRC output into (multi-frame) xyz,
             and optionally write a restart .gjf.

从 Gaussian 16 的 .log/.out 中提取结构：最后一个、能量最低、全部、每个扫描点的收敛结构，
或按反应坐标排好序的 IRC 路径；能量写在 xyz 注释行。也可以为中断/失败的作业生成续算 gjf。

Usage examples
--------------
  python log2xyz.py opt.log                      # last geometry      -> opt.xyz
  python log2xyz.py opt.log --lowest             # lowest energy      -> opt_lowest.xyz
  python log2xyz.py opt.log --all -o traj.xyz    # every geometry (optimization trajectory)
  python log2xyz.py scan.log --converged         # optimized structure of every scan point
  python log2xyz.py irc.log --irc                # IRC path: reverse end -> TS -> forward end
  python log2xyz.py opt.log --list               # only print a table of all frames
  python log2xyz.py opt.log --index 12 -o -      # frame 12 to stdout
  python log2xyz.py failed.log --gjf restart.gjf --maxcycles 200 --guess-read

Notes
-----
* A "frame" is one geometry printed in an orientation block. Input and Standard orientation of
  the same step are merged; --orient chooses which one is written (auto = Input orientation if
  present, otherwise Standard; NoSymm jobs only print Input orientation).
* The energy of a frame is the most specific energy printed for that geometry
  (Counterpoise > ONIOM extrapolated > CCSD(T) > MP2/double hybrid > TD total > SCF Done).
* IRC (--irc): Gaussian prints one gradient geometry per IRC point; the frame written for point n
  is the last geometry printed before "Point Number: n". With the default HPC integrator this is
  the predictor geometry, which is very close to (not identical with) the corrected point that
  appears in the "Summary of reaction path following" table.
* The comment line uses key=value pairs (energy in Hartree), readable by ASE (extxyz).
* The restart .gjf re-uses the first user route of the log (Geom=Check/AllCheck/Connectivity are
  removed because explicit coordinates are written), the echoed Link 0 lines, charge/multiplicity,
  atom labels / freeze flags / ONIOM layer columns from the "Symbolic Z-matrix" echo, and any
  ModRedundant lines. Gen/GenECP basis blocks are NOT in the log and must be appended by hand.

Dependencies: Python >= 3.9 standard library only.

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import os
import re
import sys

HARTREE_TO_KCAL = 627.509474

SYMBOLS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
    "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl "
    "Mc Lv Ts Og"
).split()

RE_DASH = re.compile(r"^ -{10,}\s*$")
RE_SCF = re.compile(r"SCF Done:\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+)")
RE_MP2 = re.compile(r"EUMP2\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_DH = re.compile(r"E2\(([^)]+)\)\s*=\s*\S+\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_CCSDT = re.compile(r"\bCCSD\(T\)\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_ONIOM = re.compile(r"ONIOM: extrapolated energy\s*=\s*(-?\d+\.\d+)")
RE_CP = re.compile(r"Counterpoise corrected energy\s*=\s*(-?\d+\.\d+)")
RE_TD = re.compile(r"Total Energy, E\(([^)]*)\)\s*=\s*(-?\d+\.\d+)")
RE_CHARGE = re.compile(r"Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)")
RE_SCANPT = re.compile(r"Step number\s+\d+\s+out of a maximum of\s+\d+\s+on scan point\s+(\d+)\s+out of\s+(\d+)")
RE_IRCPT = re.compile(r"Point Number:\s+(\d+)\s+Path Number:\s+(\d+)")
RE_TS_REF = re.compile(r"Energies reported relative to the TS energy of\s+(-?\d+\.\d+)")
RE_OPT_KW = re.compile(r"(^|[\s#/,])opt(?=$|[\s=(,])", re.IGNORECASE)

ENERGY_PRIORITY = ("cp", "oniom", "ccsdt", "mp2", "td", "scf")
ENERGY_LABEL = {"cp": "Counterpoise", "oniom": "ONIOM", "ccsdt": "CCSD(T)", "mp2": "MP2",
                "td": "TD", "scf": "SCF"}


def dfloat(text):
    return float(text.replace("D", "E").replace("d", "e"))


def symbol(z, ghost="X"):
    if z == 0:
        return ghost
    return SYMBOLS[z] if 0 < z < len(SYMBOLS) else "X"


# --------------------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------------------
def read_orientation(lines_iter):
    """Return (znums, coords) of an orientation block; the header line was already consumed."""
    znums, coords = [], []
    dashes = 0
    for line in lines_iter:
        if RE_DASH.match(line):
            dashes += 1
            if dashes == 3:
                break
            continue
        if dashes == 2:
            tok = line.split()
            if len(tok) >= 5:
                try:
                    z = int(tok[1])
                    xyz = [float(t) for t in tok[-3:]]
                except ValueError:
                    continue
                znums.append(z)
                coords.append(xyz)
    return znums, coords


class LogData:
    """Everything log2xyz needs from one output file."""

    def __init__(self):
        self.frames = []          # list of dicts
        self.link0 = []           # echoed Link 0 lines of the first job
        self.routes = []          # (job, route)
        self.charge = None
        self.mult = None
        self.symbolic = []        # atom lines of the first "Symbolic Z-matrix" echo
        self.modred = []          # ModRedundant lines of the first job
        self.irc_dirs = {}        # path number -> 'forward' / 'reverse'
        self.irc_summary = []     # (rel_energy, rx_coord)
        self.irc_ts_energy = None
        self.normal = 0
        self.error = False


def parse_log(path):
    data = LogData()
    job = 1
    job_kind = "other"
    cur = None
    prev = ""
    capturing_route = False
    route_buf = []
    in_symbolic = False
    in_modred = False
    in_irc_summary = False
    current_path = None
    header_done = False

    def new_frame():
        fr = {"job": job, "coords": {}, "znums": None, "energies": {}, "converged": False,
              "scan_point": None, "irc": None}
        data.frames.append(fr)
        return fr

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        it = iter(fh)
        for line in it:
            line = line.rstrip("\r\n")

            if capturing_route:
                if RE_DASH.match(line):
                    capturing_route = False
                    route = "".join(route_buf).strip()
                    data.routes.append((job, route))
                    job_kind = "opt" if RE_OPT_KW.search(route) else "other"
                    header_done = True
                else:
                    route_buf.append(line[1:] if line.startswith(" ") else line)
                prev = line
                continue
            if line.startswith(" #") and RE_DASH.match(prev):
                capturing_route = True
                route_buf = [line[1:]]
                prev = line
                continue
            prev = line

            if not header_done and job == 1 and line.startswith(" %"):
                data.link0.append(line.strip())
                continue

            if in_symbolic:
                if not line.strip():
                    in_symbolic = False
                elif RE_CHARGE.search(line):
                    if data.charge is None:
                        m = RE_CHARGE.search(line)
                        data.charge, data.mult = int(m.group(1)), int(m.group(2))
                else:
                    data.symbolic.append(line.strip())
                continue
            if in_modred:
                if not line.strip():
                    in_modred = False
                else:
                    data.modred.append(line.strip())
                continue
            if in_irc_summary:
                tok = line.split()
                if len(tok) == 3:
                    try:
                        data.irc_summary.append((float(tok[1]), float(tok[2])))
                        continue
                    except ValueError:
                        pass
                if data.irc_summary and RE_DASH.match(line):
                    in_irc_summary = False
                continue

            if "orientation:" in line and ("Input" in line or "Standard" in line or "Z-Matrix" in line):
                kind = "standard" if "Standard" in line else "input"
                znums, coords = read_orientation(it)
                if not znums:
                    continue
                if (cur is None or cur["energies"] or kind in cur["coords"] or cur["job"] != job
                        or len(cur["znums"]) != len(znums)):
                    cur = new_frame()
                    cur["znums"] = znums
                cur["coords"][kind] = coords
                continue
            if "Symbolic Z-matrix:" in line and job == 1 and not data.symbolic:
                in_symbolic = True
                continue
            if "The following ModRedundant input section has been read" in line and job == 1:
                in_modred = True
                continue
            if data.charge is None:
                m = RE_CHARGE.search(line)
                if m:
                    data.charge, data.mult = int(m.group(1)), int(m.group(2))
                    continue
            if "Proceeding to internal job step number" in line:
                job += 1
                cur = None
                continue
            if "Normal termination of Gaussian" in line:
                data.normal += 1
                continue
            if "Error termination" in line:
                data.error = True
                continue

            energy = None
            m = RE_SCF.search(line)
            if m:
                energy = ("scf", float(m.group(2)))
            elif "EUMP2" in line and RE_MP2.search(line):
                energy = ("mp2", dfloat(RE_MP2.search(line).group(1)))
            elif "E2(" in line and RE_DH.search(line):
                energy = ("mp2", dfloat(RE_DH.search(line).group(3)))
            elif "CCSD(T)" in line and RE_CCSDT.search(line):
                energy = ("ccsdt", dfloat(RE_CCSDT.search(line).group(1)))
            elif "ONIOM: extrapolated energy" in line and RE_ONIOM.search(line):
                energy = ("oniom", float(RE_ONIOM.search(line).group(1)))
            elif "Counterpoise corrected energy" in line and RE_CP.search(line):
                energy = ("cp", float(RE_CP.search(line).group(1)))
            elif "Total Energy, E(" in line and RE_TD.search(line):
                energy = ("td", float(RE_TD.search(line).group(2)))
            if energy:
                if cur is not None:
                    cur["energies"][energy[0]] = energy[1]
                continue

            if "Stationary point found" in line:
                if cur is not None and job_kind == "opt":
                    cur["converged"] = True
                continue
            m = RE_SCANPT.search(line)
            if m:
                if cur is not None:
                    cur["scan_point"] = int(m.group(1))
                continue
            m = RE_IRCPT.search(line)
            if m:
                point, pathno = int(m.group(1)), int(m.group(2))
                current_path = pathno
                if cur is not None and cur["irc"] is None:
                    if point == 0 and any(f["irc"] and f["irc"][1] == 0 for f in data.frames):
                        continue
                    cur["irc"] = (pathno, point)
                continue
            if "Calculation of FORWARD path complete" in line and current_path:
                data.irc_dirs[current_path] = "forward"
                continue
            if "Calculation of REVERSE path complete" in line and current_path:
                data.irc_dirs[current_path] = "reverse"
                continue
            m = RE_TS_REF.search(line)
            if m:
                data.irc_ts_energy = float(m.group(1))
                continue
            if "Summary of reaction path following" in line:
                in_irc_summary = True
                data.irc_summary = []
                continue

    merge_duplicates(data)
    return data


def frame_coords(fr, orient):
    if orient == "standard":
        return fr["coords"].get("standard") or fr["coords"].get("input")
    if orient == "input":
        return fr["coords"].get("input") or fr["coords"].get("standard")
    return fr["coords"].get("input") or fr["coords"].get("standard")


def best_energy(fr):
    for key in ENERGY_PRIORITY:
        if key in fr["energies"]:
            return fr["energies"][key], ENERGY_LABEL[key]
    return None, None


def merge_duplicates(data, tol=1e-6):
    """Merge consecutive frames with identical geometry (e.g. final re-print after an optimization)."""
    merged = []
    for fr in data.frames:
        if merged:
            last = merged[-1]
            a, b = frame_coords(last, "auto"), frame_coords(fr, "auto")
            same_kind = set(last["coords"]) == set(fr["coords"])
            if (same_kind and last["znums"] == fr["znums"] and a and b
                    and all(abs(p - q) <= tol for ra, rb in zip(a, b) for p, q in zip(ra, rb))
                    and (not last["energies"] or not fr["energies"]
                         or best_energy(last)[0] == best_energy(fr)[0])):
                for k, v in fr["energies"].items():
                    last["energies"].setdefault(k, v)
                last["converged"] = last["converged"] or fr["converged"]
                last["scan_point"] = last["scan_point"] or fr["scan_point"]
                last["irc"] = last["irc"] or fr["irc"]
                continue
        merged.append(fr)
    data.frames = merged


# --------------------------------------------------------------------------------------
# frame selection
# --------------------------------------------------------------------------------------
def select_frames(data, args):
    frames = data.frames
    if not frames:
        raise ValueError("no orientation blocks (geometries) found")
    idx_all = list(range(len(frames)))
    if args.index is not None:
        i = args.index - 1 if args.index > 0 else len(frames) + args.index
        if not 0 <= i < len(frames):
            raise ValueError("--index %d out of range (1..%d)" % (args.index, len(frames)))
        return [i], "frame%d" % (i + 1)
    if args.all:
        return idx_all, "all"
    if args.lowest:
        cand = [i for i in idx_all if best_energy(frames[i])[0] is not None]
        if not cand:
            raise ValueError("no frame has an energy")
        return [min(cand, key=lambda i: best_energy(frames[i])[0])], "lowest"
    if args.converged:
        sel = [i for i in idx_all if frames[i]["converged"]]
        if not sel:
            raise ValueError("no converged optimization ('Stationary point found') in this file")
        return sel, "converged"
    if args.irc:
        return irc_order(data), "irc"
    return [len(frames) - 1], ""


def irc_order(data):
    tagged = [(i, fr["irc"]) for i, fr in enumerate(data.frames) if fr["irc"]]
    if not tagged:
        raise ValueError("no 'Point Number: n  Path Number: p' lines found - is this an IRC job?")
    ts = [i for i, t in tagged if t[1] == 0]
    ts_idx = ts[0] if ts else 0
    dirs = dict(data.irc_dirs)
    paths = sorted({t[0] for _, t in tagged})
    for p in paths:                      # Gaussian default: path 1 forward, path 2 reverse
        dirs.setdefault(p, "forward" if p == 1 and "forward" not in dirs.values() else "reverse")
    rev, fwd = [], []
    for i, (p, n) in tagged:
        if n == 0:
            continue
        (fwd if dirs.get(p) == "forward" else rev).append((n, i))
    order = [i for n, i in sorted(rev, reverse=True)] + [ts_idx] + [i for n, i in sorted(fwd)]
    return order


# --------------------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------------------
def comment_line(data, i, fr, source, extra=None):
    e, etype = best_energy(fr)
    parts = ["frame=%d" % (i + 1)]
    if e is not None:
        parts += ["energy=%.10f" % e, "energy_type=%s" % etype, "unit=Hartree"]
    if fr["converged"]:
        parts.append("converged=T")
    if fr["scan_point"]:
        parts.append("scan_point=%d" % fr["scan_point"])
    if fr["irc"]:
        parts.append("irc_path=%d irc_point=%d" % fr["irc"])
    if extra:
        parts.append(extra)
    parts.append("source=%s" % os.path.basename(source).replace(" ", "_"))
    return " ".join(parts)


def write_xyz(data, indices, orient, fh, source, irc_mode=False):
    extras = {}
    if irc_mode and indices and data.irc_summary and len(data.irc_summary) == len(indices) and all(
            data.frames[i]["irc"] for i in indices):
        for k, i in enumerate(indices):
            rel, rx = data.irc_summary[k]
            extras[i] = "rx_coord=%.5f" % rx
            if data.irc_ts_energy is not None:
                extras[i] += " energy_irc_summary=%.6f" % (data.irc_ts_energy + rel)
    for i in indices:
        fr = data.frames[i]
        coords = frame_coords(fr, orient)
        atoms = [(z, c) for z, c in zip(fr["znums"], coords) if z >= 0]
        fh.write("%d\n" % len(atoms))
        fh.write(comment_line(data, i, fr, source, extras.get(i)) + "\n")
        for z, (x, y, zc) in atoms:
            fh.write("%-3s %14.8f %14.8f %14.8f\n" % (symbol(z), x, y, zc))


def print_list(data, source):
    energies = [best_energy(f)[0] for f in data.frames]
    valid = [e for e in energies if e is not None]
    emin = min(valid) if valid else None
    print("# %s : %d frame(s), %d normal termination(s)%s" %
          (source, len(data.frames), data.normal, ", ERROR termination" if data.error else ""))
    print("%6s %4s %18s %10s  %s" % ("frame", "job", "energy/Eh", "dE/kcal", "flags"))
    for i, fr in enumerate(data.frames):
        e, et = best_energy(fr)
        flags = []
        if fr["converged"]:
            flags.append("converged")
        if fr["scan_point"]:
            flags.append("scan_pt=%d" % fr["scan_point"])
        if fr["irc"]:
            flags.append("irc=%d:%d" % fr["irc"])
        if e is not None:
            flags.append(et)
        print("%6d %4d %18s %10s  %s" % (
            i + 1, fr["job"], "%.8f" % e if e is not None else "-",
            "%.3f" % ((e - emin) * HARTREE_TO_KCAL) if e is not None else "-", " ".join(flags)))


# --------------------------------------------------------------------------------------
# restart .gjf
# --------------------------------------------------------------------------------------
def split_route(route):
    """Split a route line into top-level tokens (parentheses kept together)."""
    tokens, buf, depth = [], "", 0
    for ch in route:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(depth - 1, 0)
        if ch in " \t" and depth == 0:
            if buf:
                tokens.append(buf)
            buf = ""
        else:
            buf += ch
    if buf:
        tokens.append(buf)
    return tokens


def set_maxcycles(token, n):
    """Return an Opt token with MaxCycles=n set."""
    m = re.match(r"(?i)^opt\s*(=)?\s*(\((.*)\)|(.*))?$", token)
    if not m:
        return token
    inner = m.group(3) if m.group(3) is not None else (m.group(4) or "")
    opts = [o for o in re.split(r"[,\s]+", inner) if o]
    opts = [o for o in opts if not o.lower().startswith("maxcyc")]
    opts.append("MaxCycles=%d" % n)
    return "Opt=(%s)" % ",".join(opts)


def restart_route(data, args, warnings):
    if args.route:
        return args.route
    user = [r for j, r in data.routes if not r.upper().startswith("#N GEOM=ALLCHECK GUESS=TCHECK")]
    if not user:
        raise ValueError("no route section found in the log; give one with --route")
    route = user[0]
    tokens = split_route(route)
    out = []
    has_opt = False
    for tok in tokens:
        key = re.split(r"[=(]", tok, maxsplit=1)[0].lower()
        if key == "geom":
            low = tok.lower()
            if any(w in low for w in ("connect", "check", "allcheck", "modconnect")):
                continue
        if key == "guess" and "read" in tok.lower():
            continue
        if key == "opt":
            has_opt = True
            if "restart" in tok.lower():
                warnings.append("original route used Opt=Restart; the restart option was removed")
                tok = re.sub(r"(?i),?\s*restart", "", tok).replace("()", "").replace("=,", "=")
            if args.maxcycles:
                tok = set_maxcycles(tok, args.maxcycles)
        out.append(tok)
    if args.maxcycles and not has_opt:
        warnings.append("--maxcycles ignored: the route has no Opt keyword")
    if args.guess_read:
        out.append("Guess=Read")
    if args.add_route:
        out.append(args.add_route)
    low = " ".join(out).lower()
    if re.search(r"(^|[/\s])(gen|genecp)(\s|$)", low) or "pseudo=read" in low or "extrabasis" in low:
        warnings.append("route uses Gen/GenECP/ExtraBasis/Pseudo=Read: append the basis (and ECP) "
                        "blocks from the original input by hand")
    return " ".join(out)


def restart_atom_lines(data, fr, orient, warnings):
    coords = frame_coords(fr, orient)
    atoms = [(z, c) for z, c in zip(fr["znums"], coords) if z >= 0]
    lines = []
    sym = data.symbolic
    usable = len(sym) == len(atoms) and all(_cart_index(s.split()) is not None for s in sym)
    if sym and not usable:
        warnings.append("input geometry was not in plain Cartesian form (Z-matrix?); writing element "
                        "symbols only - re-apply freeze flags / fragment / ONIOM labels by hand")
    for k, (z, (x, y, zc)) in enumerate(atoms):
        xyz = "%14.8f %14.8f %14.8f" % (x, y, zc)
        if usable:
            tok = sym[k].split()
            j = _cart_index(tok)
            prefix = " ".join(tok[:j])
            suffix = " ".join(tok[j + 3:])
            lines.append(("%-16s %s %s" % (prefix, xyz, suffix)).rstrip())
        else:
            lines.append("%-4s %s" % (symbol(z, ghost="Bq"), xyz))
    return lines


def _is_coord(t):
    try:
        float(t.replace("D", "E"))
    except ValueError:
        return False
    return "." in t or "e" in t.lower()


def _cart_index(tok):
    """Index of the x coordinate in an echoed Cartesian atom line, or None."""
    for j in range(1, len(tok) - 2):
        if _is_coord(tok[j]) and _is_coord(tok[j + 1]) and _is_coord(tok[j + 2]):
            if all(_is_int(t) for t in tok[1:j]):
                return j
    return None


def _is_int(t):
    return re.fullmatch(r"[+-]?\d+", t) is not None


def write_restart(data, idx, args, path_out, source):
    warnings = []
    fr = data.frames[idx]
    route = restart_route(data, args, warnings)
    stem = os.path.splitext(os.path.basename(path_out))[0]
    old_chk = None
    link0 = []
    for l0 in data.link0:
        key = l0.split("=")[0].lower()
        if key in ("%chk", "%oldchk"):
            if key == "%chk":
                old_chk = l0.split("=", 1)[1].strip()
            continue
        if key in ("%mem", "%nprocshared", "%nproc", "%cpu", "%gpucpu", "%lindaworkers", "%nosave"):
            link0.append(l0)
    if args.guess_read:
        if old_chk:
            link0.insert(0, "%%oldchk=%s" % old_chk)
        else:
            warnings.append("--guess-read: no %chk found in the log, Guess=Read will look in the new chk")
    link0.insert(0, "%%chk=%s.chk" % stem)
    if not route.startswith("#"):
        route = "#p " + route
    title = re.sub(r"[@#!_\\\-]", " ", "Restart of %s from frame %d" % (os.path.basename(source), idx + 1))
    e, et = best_energy(fr)
    if e is not None:
        title += " E %.8f" % e
    if data.charge is None:
        raise ValueError("charge/multiplicity not found in the log")
    lines = link0 + [route, "", title, "", "%d %d" % (data.charge, data.mult)]
    lines += restart_atom_lines(data, fr, args.orient, warnings)
    lines.append("")
    if data.modred:
        lines += data.modred
        lines.append("")
        if any(re.search(r"\sS\s", " %s " % m) for m in data.modred):
            warnings.append("ModRedundant scan (S) lines copied unchanged: the scan restarts from the "
                            "written geometry - reduce the number of steps if part of it is done")
    lines += ["", ""]
    with open(path_out, "w", newline="\n", encoding="ascii", errors="replace") as fh:
        fh.write("\n".join(lines))
    print("restart input written to %s (frame %d)" % (path_out, idx + 1), file=sys.stderr)
    for w in warnings:
        print("WARNING: " + w, file=sys.stderr)


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description="Extract geometries from Gaussian 16 output (opt/scan/IRC) into xyz; optionally write "
                    "a restart .gjf.  从 Gaussian 输出中提取结构（最后/最低能量/全部/扫描点/IRC 路径），"
                    "可生成续算输入。",
        epilog="Examples:\n"
               "  python log2xyz.py opt.log                    # last geometry -> opt.xyz\n"
               "  python log2xyz.py opt.log --lowest\n"
               "  python log2xyz.py opt.log --all -o traj.xyz\n"
               "  python log2xyz.py scan.log --converged       # one frame per scan point\n"
               "  python log2xyz.py irc.log --irc              # ordered IRC path\n"
               "  python log2xyz.py opt.log --list\n"
               "  python log2xyz.py failed.log --gjf restart.gjf --maxcycles 200 --guess-read",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("logs", nargs="+", help="Gaussian .log/.out file(s) / 输出文件")
    sel = p.add_mutually_exclusive_group()
    sel.add_argument("--last", action="store_true", help="last geometry (default) / 最后一个结构（默认）")
    sel.add_argument("--lowest", action="store_true", help="lowest-energy geometry / 能量最低的结构")
    sel.add_argument("--all", action="store_true", help="every geometry, multi-frame xyz / 全部结构")
    sel.add_argument("--converged", "--scan", dest="converged", action="store_true",
                     help="converged geometries ('Stationary point found'), e.g. every relaxed-scan point "
                          "/ 每个收敛点（柔性扫描每个点）")
    sel.add_argument("--irc", action="store_true",
                     help="IRC path ordered reverse end -> TS -> forward end / 按反应坐标排序的 IRC 路径")
    sel.add_argument("--index", type=int, metavar="N",
                     help="frame N (1-based, negative counts from the end) / 第 N 个结构")
    p.add_argument("--orient", choices=["auto", "input", "standard"], default="auto",
                   help="orientation to write (default auto: Input if present) / 取向")
    p.add_argument("-o", "--output", metavar="FILE",
                   help="output xyz ('-' = stdout; default <log stem>[_mode].xyz) / 输出文件")
    p.add_argument("--list", action="store_true", help="only print a table of frames / 只列出所有结构")
    g = p.add_argument_group("restart input / 续算输入")
    g.add_argument("--gjf", metavar="FILE",
                   help="also write a restart .gjf from the selected (single) frame / 生成续算 gjf")
    g.add_argument("--route", help="route line for the restart input (default: original route)")
    g.add_argument("--maxcycles", type=int, help="set Opt MaxCycles in the restart route")
    g.add_argument("--guess-read", action="store_true",
                   help="add Guess=Read and %%oldchk=<original chk> / 读取原 chk 的波函数作为初猜")
    g.add_argument("--add-route", metavar="TEXT", help="extra route keywords, e.g. 'SCF=XQC'")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if len(args.logs) > 1 and args.output and args.output != "-":
        print("ERROR: -o/--output can only be used with a single log file", file=sys.stderr)
        return 2
    if args.gjf and len(args.logs) > 1:
        print("ERROR: --gjf can only be used with a single log file", file=sys.stderr)
        return 2
    status = 0
    for log in args.logs:
        if not os.path.isfile(log):
            print("ERROR: %s not found" % log, file=sys.stderr)
            status = 1
            continue
        try:
            data = parse_log(log)
        except OSError as exc:
            print("ERROR: cannot read %s: %s" % (log, exc), file=sys.stderr)
            status = 1
            continue
        if not data.frames:
            print("ERROR: %s contains no geometry (orientation block); is it a Gaussian output?" % log,
                  file=sys.stderr)
            status = 1
            continue
        if args.list:
            print_list(data, log)
            continue
        try:
            indices, tag = select_frames(data, args)
        except ValueError as exc:
            print("ERROR: %s: %s" % (log, exc), file=sys.stderr)
            status = 1
            continue
        stem = os.path.splitext(log)[0]
        out = args.output or (stem + ("_" + tag if tag else "") + ".xyz")
        if out == "-":
            write_xyz(data, indices, args.orient, sys.stdout, log, irc_mode=(tag == "irc"))
        else:
            if os.path.abspath(out) == os.path.abspath(log):
                print("ERROR: output would overwrite the log file", file=sys.stderr)
                return 2
            with open(out, "w", newline="\n", encoding="ascii", errors="replace") as fh:
                write_xyz(data, indices, args.orient, fh, log, irc_mode=(tag == "irc"))
            print("%s: %d frame(s) -> %s" % (log, len(indices), out), file=sys.stderr)
        if data.error:
            print("NOTE: %s ended with an error termination" % log, file=sys.stderr)
        elif data.normal == 0:
            print("NOTE: %s has no 'Normal termination' (running, killed or truncated)" % log, file=sys.stderr)
        if args.gjf:
            if len(indices) != 1:
                print("NOTE: several frames selected; the restart input uses the last of them", file=sys.stderr)
            try:
                write_restart(data, indices[-1], args, args.gjf, log)
            except ValueError as exc:
                print("ERROR: %s" % exc, file=sys.stderr)
                status = 1
    return status


if __name__ == "__main__":
    sys.exit(main())
