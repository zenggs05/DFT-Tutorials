#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scan_irc_profile.py - energy profiles of Gaussian 16 PES scans and IRC paths.

从 Gaussian 16 输出中提取：
  * 柔性扫描 Opt=ModRedundant（一维、二维）：每个扫描点收敛结构（Stationary point found）的能量，
    扫描坐标值由收敛结构按 ModRedundant 中 S 行的原子编号直接计算（键长/键角/二面角）；
  * 刚性扫描（Scan 关键词）：读取 "Summary of the potential surface scan" 汇总表，没有时逐点读 SCF Done；
  * IRC：读取 "Summary of reaction path following"（相对 TS 的能量与反应坐标），作业未完成时逐点读取
    "Point Number / Path Number / NET REACTION COORDINATE"；正反两个方向（同一个 log 或 Forward/Reverse
    两个 log）自动合并成一条曲线。
以 kcal/mol（或 kJ/mol、eV）相对最低点作图，写 CSV，并报告扫描的最高点（TS 初猜）与局部极大值。

Usage examples
--------------
  python scan_irc_profile.py butane_scan.log                       # 1D relaxed scan -> png + table
  python scan_irc_profile.py sn2_scan.log --write-max ts_guess.xyz # geometry of the highest point
  python scan_irc_profile.py scan2d.log --csv scan2d.csv           # 2D scan -> contour plot
  python scan_irc_profile.py irc.log                               # both IRC directions in one log
  python scan_irc_profile.py irc_fwd.log irc_rev.log -o irc.png    # forward + reverse logs merged
  python scan_irc_profile.py irc.log --ref ts --unit kj
  python scan_irc_profile.py gic_scan.log --coord "B 1 5" --coord "B 1 6"   # coordinates by hand

Notes
-----
* Energy of a point = the most specific energy of its geometry (Counterpoise > ONIOM > MP2/double
  hybrid > TD total > SCF Done).  Unconverged relaxed-scan points (no "Stationary point found") are
  skipped and reported - their energies are not trustworthy.
* Coordinates are computed from the geometry (Angstrom / degree), so they do not depend on the
  layout of the summary table; dihedrals are unwrapped to a continuous curve (--no-unwrap to keep
  (-180, 180]).  Scans that are not described by ModRedundant lines (GIC, Z-matrix variables) fall
  back to the point index unless --coord is given.
* IRC: the reaction coordinate is in Bohr amu^1/2 as printed by Gaussian; reverse-path points are
  negative.  Endpoints are only approximate minima - optimise them before reporting energies.

Dependencies: Python >= 3.9 standard library; matplotlib for figures (--no-plot works without it).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import csv
import math
import os
import re
import sys

HARTREE = {"kcal": 627.509474, "kj": 2625.49964, "ev": 27.211386}
UNIT_LABEL = {"kcal": "kcal/mol", "kj": "kJ/mol", "ev": "eV"}

SYMBOLS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce "
    "Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn"
).split()

RE_DASH = re.compile(r"^ -{10,}\s*$")
RE_SCF = re.compile(r"SCF Done:\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+)")
RE_MP2 = re.compile(r"EUMP2\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_DH = re.compile(r"E2\(([^)]+)\)\s*=\s*\S+\s+E\(([^)]+)\)\s*=\s*(-?\d+\.\d+[DdEe][+-]\d+)")
RE_ONIOM = re.compile(r"ONIOM: extrapolated energy\s*=\s*(-?\d+\.\d+)")
RE_CP = re.compile(r"Counterpoise corrected energy\s*=\s*(-?\d+\.\d+)")
RE_TD = re.compile(r"Total Energy, E\(([^)]*)\)\s*=\s*(-?\d+\.\d+)")
RE_SCANPT = re.compile(r"on scan point\s+(\d+)\s+out of\s+(\d+)")
RE_IRCPT = re.compile(r"Point Number:\s+(\d+)\s+Path Number:\s+(\d+)")
RE_NETRC = re.compile(r"NET REACTION COORDINATE UP TO THIS POINT\s*=\s*(-?\d+\.\d+)")
RE_TSREF = re.compile(r"Energies reported relative to the TS energy of\s+(-?\d+\.\d+)")
RE_OPT_KW = re.compile(r"(^|[\s#/,])opt(?=$|[\s=(,])", re.IGNORECASE)
RE_IRC_KW = re.compile(r"(^|[\s#/,])irc(?=$|[\s=(,])", re.IGNORECASE)
RE_SCAN_KW = re.compile(r"(^|[\s#/,])scan(?=$|[\s=(,])", re.IGNORECASE)
ENERGY_PRIORITY = ("cp", "oniom", "mp2", "td", "scf")

COLOR_SERIES = "#2a78d6"
COLOR_ACCENT = "#d0582a"
COLOR_GRID = "#e4e3df"
COLOR_TEXT = "#0b0b0b"
COLOR_MUTED = "#8a8984"
COLOR_SURFACE = "#fcfcfb"


def dfloat(text):
    return float(text.replace("D", "E").replace("d", "e"))


def best(energies):
    for k in ENERGY_PRIORITY:
        if k in energies:
            return energies[k], k.upper()
    return None, None


# ---------------------------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------------------------
def _sub(a, b):
    return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def internal(coords, atoms):
    """Distance (2 atoms), angle (3) or dihedral (4) from 1-based atom indices."""
    p = [coords[i - 1] for i in atoms]
    if len(p) == 2:
        return math.sqrt(_dot(_sub(p[0], p[1]), _sub(p[0], p[1])))
    if len(p) == 3:
        u, v = _sub(p[0], p[1]), _sub(p[2], p[1])
        c = _dot(u, v) / math.sqrt(_dot(u, u) * _dot(v, v))
        return math.degrees(math.acos(max(-1.0, min(1.0, c))))
    b1, b2, b3 = _sub(p[1], p[0]), _sub(p[2], p[1]), _sub(p[3], p[2])
    n1, n2 = _cross(b1, b2), _cross(b2, b3)
    m1 = _cross(n1, [x / math.sqrt(_dot(b2, b2)) for x in b2])
    return math.degrees(math.atan2(_dot(m1, n2), _dot(n1, n2)))


COORD_NAME = {2: "R", 3: "A", 4: "D"}
COORD_UNIT = {2: "Angstrom", 3: "deg", 4: "deg"}


def parse_coord_spec(text):
    """'B 1 2', 'D 1 2 3 4', '1,2,3' -> tuple of atoms."""
    tok = [t for t in re.split(r"[\s,]+", text.strip()) if t]
    if tok and tok[0].upper() in ("B", "R", "A", "D"):
        tok = tok[1:]
    try:
        atoms = tuple(int(t) for t in tok)
    except ValueError:
        raise ValueError("bad --coord '%s' (use e.g. 'B 1 2', 'A 1 2 3', 'D 1 2 3 4')" % text)
    if len(atoms) not in (2, 3, 4) or min(atoms) < 1:
        raise ValueError("--coord needs 2, 3 or 4 positive atom numbers: '%s'" % text)
    return atoms


def modred_scan_coords(lines):
    """Atom tuples of the ModRedundant lines that contain an 'S' (scan) action."""
    out = []
    for ln in lines:
        tok = ln.split()
        if "S" not in [t.upper() for t in tok]:
            continue
        if tok and tok[0].upper() in ("B", "A", "D", "L", "R"):
            tok = tok[1:]
        atoms = []
        for t in tok:
            if re.fullmatch(r"\d+", t):
                atoms.append(int(t))
            else:
                break
        if len(atoms) in (2, 3, 4):
            out.append(tuple(atoms))
    return out


# ---------------------------------------------------------------------------------------------
# parsing
# ---------------------------------------------------------------------------------------------
def read_orientation(it):
    znums, coords, dashes = [], [], 0
    for line in it:
        if RE_DASH.match(line):
            dashes += 1
            if dashes == 3:
                break
            continue
        if dashes == 2:
            tok = line.split()
            if len(tok) >= 5:
                try:
                    znums.append(int(tok[1]))
                    coords.append([float(t) for t in tok[-3:]])
                except ValueError:
                    continue
    return znums, coords


def parse_log(path):
    d = {"file": path, "kind": None, "routes": [], "modred": [], "points": [], "unconverged": [],
         "scan_total": None, "rigid_table": [], "rigid_header": [], "irc_summary": [], "irc_ts": None,
         "irc_points": [], "irc_dirs": {}, "frames": [], "normal": 0, "error": False}
    capturing, route_buf, prev = False, [], ""
    in_modred = in_rigid = in_irc = False
    geom, znums, energies = None, None, {}
    scan_point = None
    netrc = None
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        it = (ln.rstrip("\r\n") for ln in fh)
        for line in it:
            if capturing:
                if RE_DASH.match(line):
                    capturing = False
                    d["routes"].append("".join(route_buf).strip())
                else:
                    route_buf.append(line[1:] if line.startswith(" ") else line)
                prev = line
                continue
            if line.startswith(" #") and RE_DASH.match(prev):
                capturing, route_buf, prev = True, [line[1:]], line
                continue
            prev = line
            if in_modred:
                if not line.strip():
                    in_modred = False
                else:
                    d["modred"].append(line.strip())
                continue
            if in_rigid:
                tok = line.split()
                if not tok:
                    if d["rigid_table"]:
                        in_rigid = False
                    continue
                if tok[0] == "N" and not d["rigid_header"]:
                    d["rigid_header"] = tok
                    continue
                if set(line.strip()) <= set("- "):          # '----  ---------' separator
                    continue
                if re.fullmatch(r"\d+", tok[0]) and len(tok) >= 3:
                    try:
                        d["rigid_table"].append([dfloat(t) for t in tok])
                        continue
                    except ValueError:
                        pass
                in_rigid = False
                continue
            if in_irc:
                tok = line.split()
                if len(tok) >= 3 and re.fullmatch(r"\d+", tok[0]):
                    try:
                        d["irc_summary"].append((float(tok[1]), float(tok[2])))
                        continue
                    except ValueError:
                        pass
                if d["irc_summary"] and RE_DASH.match(line):
                    in_irc = False
                continue

            if "orientation:" in line and ("Input" in line or "Standard" in line or "Z-Matrix" in line):
                z, c = read_orientation(it)
                if z:
                    if energies:                 # a new geometry starts
                        energies = {}
                    geom, znums = c, z
                continue
            m = RE_SCF.search(line)
            if m:
                energies["scf"] = float(m.group(2))
                d["frames"].append((znums, geom, float(m.group(2))))
                continue
            if "EUMP2" in line and RE_MP2.search(line):
                energies["mp2"] = dfloat(RE_MP2.search(line).group(1))
                continue
            if "E2(" in line and RE_DH.search(line):
                energies["mp2"] = dfloat(RE_DH.search(line).group(3))
                continue
            if "ONIOM: extrapolated energy" in line and RE_ONIOM.search(line):
                energies["oniom"] = float(RE_ONIOM.search(line).group(1))
                continue
            if "Counterpoise corrected energy" in line and RE_CP.search(line):
                energies["cp"] = float(RE_CP.search(line).group(1))
                continue
            if "Total Energy, E(" in line and RE_TD.search(line):
                energies["td"] = float(RE_TD.search(line).group(2))
                continue
            if "The following ModRedundant input section has been read" in line and not d["modred"]:
                in_modred = True
                continue
            m = RE_SCANPT.search(line)
            if m:
                if scan_point is not None and int(m.group(1)) != scan_point and \
                        not any(p["scan_point"] == scan_point for p in d["points"]):
                    d["unconverged"].append(scan_point)
                scan_point = int(m.group(1))
                d["scan_total"] = int(m.group(2))
                continue
            if "Stationary point found" in line:
                e, et = best(energies)
                if geom is not None and e is not None and not any(r.upper().startswith("#N GEOM=ALLCHECK")
                                                                  for r in d["routes"][-1:]):
                    d["points"].append({"scan_point": scan_point, "energy": e, "etype": et,
                                        "coords": geom, "znums": znums})
                continue
            if "Summary of the potential surface scan" in line:
                in_rigid = True
                d["rigid_table"], d["rigid_header"] = [], []
                continue
            m = RE_NETRC.search(line)
            if m:
                netrc = float(m.group(1))
                if d["irc_points"]:
                    d["irc_points"][-1]["rc"] = netrc
                continue
            m = RE_IRCPT.search(line)
            if m:
                e, et = best(energies)
                d["irc_points"].append({"point": int(m.group(1)), "path": int(m.group(2)), "energy": e,
                                        "rc": None, "coords": geom, "znums": znums})
                continue
            if "Calculation of FORWARD path complete" in line and d["irc_points"]:
                d["irc_dirs"][d["irc_points"][-1]["path"]] = 1.0
                continue
            if "Calculation of REVERSE path complete" in line and d["irc_points"]:
                d["irc_dirs"][d["irc_points"][-1]["path"]] = -1.0
                continue
            m = RE_TSREF.search(line)
            if m:
                d["irc_ts"] = float(m.group(1))
                continue
            if "Summary of reaction path following" in line:
                in_irc = True
                d["irc_summary"] = []
                continue
            if "Normal termination of Gaussian" in line:
                d["normal"] += 1
                continue
            if "Error termination" in line:
                d["error"] = True
    route = " ".join(d["routes"])
    if RE_IRC_KW.search(route):
        d["kind"] = "irc"
    elif RE_OPT_KW.search(route) and (d["scan_total"] or modred_scan_coords(d["modred"])):
        d["kind"] = "relaxed"
    elif RE_SCAN_KW.search(route) or d["rigid_table"]:
        d["kind"] = "rigid"
    elif d["points"] and d["scan_total"]:
        d["kind"] = "relaxed"
    return d


# ---------------------------------------------------------------------------------------------
# profile building
# ---------------------------------------------------------------------------------------------
RE_ENERGY_HDR = re.compile(r"^(SCF|MP\d|CC|QCI|CI|BD|E\()", re.IGNORECASE)


def rigid_nvar(header, ncol):
    """Number of scan-variable columns of the rigid-scan summary (N, variables..., energies...)."""
    if header:
        for i, h in enumerate(header[1:]):
            if RE_ENERGY_HDR.match(h):
                return max(i, 1)
    return max(ncol - 2, 1)


def unwrap(values):
    out = []
    for v in values:
        if out:
            while v - out[-1] > 180.0:
                v -= 360.0
            while v - out[-1] < -180.0:
                v += 360.0
        out.append(v)
    return out


def scan_profile(d, coord_specs, do_unwrap):
    """Return (rows, names); row = dict(x=[...], energy, point, coords, znums)."""
    rows = []
    if d["kind"] == "relaxed":
        pts = sorted(d["points"], key=lambda p: (p["scan_point"] or 0))
        uniq = {}
        for p in pts:                         # keep the last converged structure of each point
            uniq[p["scan_point"]] = p
        pts = [uniq[k] for k in sorted(uniq, key=lambda k: k or 0)]
        for p in pts:
            rows.append({"point": p["scan_point"], "energy": p["energy"], "coords": p["coords"],
                         "znums": p["znums"]})
    else:
        if d["rigid_table"]:
            nvar = rigid_nvar(d["rigid_header"], len(d["rigid_table"][0]))
            for r in d["rigid_table"]:
                rows.append({"point": int(r[0]), "energy": r[-1], "table_x": r[1:1 + nvar], "coords": None,
                             "znums": None})
            frames = d["frames"]
            if len(frames) == len(rows):
                for r, (z, c, _) in zip(rows, frames):
                    r["coords"], r["znums"] = c, z
        else:
            for i, (z, c, e) in enumerate(d["frames"], 1):
                rows.append({"point": i, "energy": e, "coords": c, "znums": z})
    specs = list(coord_specs) or modred_scan_coords(d["modred"])
    names = []
    if rows and specs and all(r["coords"] for r in rows):
        natom = min(len(r["coords"]) for r in rows)
        if max(max(s) for s in specs) > natom:
            raise ValueError("coordinate atom number larger than the number of atoms (%d)" % natom)
        for s in specs:
            vals = [internal(r["coords"], s) for r in rows]
            if len(s) == 4 and do_unwrap:
                vals = unwrap(vals)
            for r, v in zip(rows, vals):
                r.setdefault("x", []).append(v)
            names.append("%s(%s) / %s" % (COORD_NAME[len(s)], ",".join(map(str, s)), COORD_UNIT[len(s)]))
    elif rows and rows[0].get("table_x"):
        hdr = d["rigid_header"]
        for r in rows:
            r["x"] = list(r["table_x"])
        names = (hdr[1:1 + len(rows[0]["x"])] if hdr else ["var%d" % (i + 1) for i in range(len(rows[0]["x"]))])
    else:
        for r in rows:
            r["x"] = [float(r["point"])]
        names = ["scan point"]
    return rows, names[:2] if len(names) > 2 else names


def irc_profile(logs_data, flip):
    """Merge the IRC data of one or more logs; return rows sorted by reaction coordinate."""
    rows, ts_energies = [], []
    for d in logs_data:
        route = " ".join(d["routes"]).lower()
        if d["irc_summary"] and d["irc_ts"] is not None:
            ts_energies.append(d["irc_ts"])
            for rel, rc in d["irc_summary"]:
                rows.append({"x": [rc], "energy": d["irc_ts"] + rel, "source": "summary", "file": d["file"],
                             "coords": None, "znums": None})
            continue
        pts = d["irc_points"]
        if not pts:
            continue
        dirs = dict(d["irc_dirs"])
        only_rev = re.search(r"irc\s*=\s*\(?[^\s]*\breverse\b", route) is not None
        for p in pts:
            sign = dirs.get(p["path"], -1.0 if (p["path"] == 2 or only_rev) else 1.0)
            if p["point"] == 0 or p["rc"] is None:
                if p["point"] == 0 and p["energy"] is not None:
                    ts_energies.append(p["energy"])
                    rows.append({"x": [0.0], "energy": p["energy"], "source": "points", "file": d["file"],
                                 "coords": p["coords"], "znums": p["znums"]})
                continue
            if p["energy"] is None:
                continue
            rows.append({"x": [sign * abs(p["rc"])], "energy": p["energy"], "source": "points",
                         "file": d["file"], "coords": p["coords"], "znums": p["znums"]})
    if flip:
        for r in rows:
            r["x"][0] = -r["x"][0]
    # drop duplicate TS points (one per log / path)
    rows.sort(key=lambda r: r["x"][0])
    merged = []
    for r in rows:
        if merged and abs(r["x"][0] - merged[-1]["x"][0]) < 1e-6 and abs(r["energy"] - merged[-1]["energy"]) < 1e-5:
            continue
        merged.append(r)
    for i, r in enumerate(merged, 1):
        r["point"] = i
    return merged, ts_energies


def local_maxima(rows):
    e = [r["energy"] for r in rows]
    return [i for i in range(1, len(e) - 1) if e[i] > e[i - 1] and e[i] >= e[i + 1]]


# ---------------------------------------------------------------------------------------------
# output
# ---------------------------------------------------------------------------------------------
def write_xyz(row, path, comment):
    if not row.get("coords"):
        raise ValueError("no geometry stored for this point")
    with open(path, "w", newline="\n", encoding="ascii") as fh:
        fh.write("%d\n%s\n" % (len(row["coords"]), comment))
        for z, (x, y, c) in zip(row["znums"], row["coords"]):
            sym = SYMBOLS[z] if 0 < z < len(SYMBOLS) else "X"
            fh.write("%-3s %14.8f %14.8f %14.8f\n" % (sym, x, y, c))


def plot_1d(rows, rel, xname, unit, title, out, show, dpi, mark_idx, is_irc):
    try:
        import matplotlib
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("ERROR: matplotlib is not installed (pip install matplotlib) - use --no-plot", file=sys.stderr)
        return False
    x = [r["x"][0] for r in rows]
    fig, ax = plt.subplots(figsize=(6.4, 4.2), facecolor=COLOR_SURFACE)
    ax.set_facecolor(COLOR_SURFACE)
    ax.grid(True, color=COLOR_GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.plot(x, rel, color=COLOR_SERIES, linewidth=2, marker="o", markersize=3.5, zorder=3)
    if mark_idx is not None:
        ax.plot([x[mark_idx]], [rel[mark_idx]], linestyle="none", marker="o", markersize=9,
                markerfacecolor="none", markeredgecolor=COLOR_ACCENT, markeredgewidth=2, zorder=4)
        ax.annotate("%s %.2f" % ("TS" if is_irc else "max", rel[mark_idx]), (x[mark_idx], rel[mark_idx]),
                    xytext=(6, 6), textcoords="offset points", color=COLOR_TEXT, fontsize=9)
    ax.set_xlabel(xname, color=COLOR_TEXT)
    ax.set_ylabel("relative energy (%s)" % UNIT_LABEL[unit], color=COLOR_TEXT)
    ax.set_title(title, loc="left", fontsize=10, color=COLOR_TEXT)
    ax.tick_params(colors=COLOR_MUTED)
    fig.tight_layout()
    return _finish(fig, plt, out, show, dpi)


def plot_2d(rows, rel, names, unit, title, out, show, dpi):
    try:
        import matplotlib
        if not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except ImportError:
        print("ERROR: matplotlib is not installed (pip install matplotlib) - use --no-plot", file=sys.stderr)
        return False
    x = [r["x"][0] for r in rows]
    y = [r["x"][1] for r in rows]
    fig, ax = plt.subplots(figsize=(6.0, 5.0), facecolor=COLOR_SURFACE)
    try:
        tri = mtri.Triangulation(x, y)
        cs = ax.tricontourf(tri, rel, levels=20, cmap="viridis")
        ax.tricontour(tri, rel, levels=10, colors="white", linewidths=0.5, alpha=0.6)
    except (ValueError, RuntimeError):
        cs = ax.scatter(x, y, c=rel, cmap="viridis", s=40)
    ax.plot(x, y, linestyle="none", marker=".", color="black", markersize=2)
    fig.colorbar(cs, ax=ax, label="relative energy (%s)" % UNIT_LABEL[unit])
    ax.set_xlabel(names[0])
    ax.set_ylabel(names[1])
    ax.set_title(title, loc="left", fontsize=10)
    fig.tight_layout()
    return _finish(fig, plt, out, show, dpi)


def _finish(fig, plt, out, show, dpi):
    if show:
        plt.show()
    else:
        fig.savefig(out, dpi=dpi, facecolor=COLOR_SURFACE, bbox_inches="tight")
        print("figure written to %s" % out, file=sys.stderr)
    plt.close(fig)
    return True


def reference(rows, ref, is_irc):
    e = [r["energy"] for r in rows]
    if ref == "min":
        return min(e)
    if ref == "first":
        return e[0]
    if ref == "ts":
        if not is_irc:
            raise ValueError("--ref ts only makes sense for IRC")
        zero = [r["energy"] for r in rows if abs(r["x"][0]) < 1e-8]
        return zero[0] if zero else max(e)
    try:
        return float(ref)
    except ValueError:
        raise ValueError("--ref must be min, first, ts or an energy in Hartree")


def report(rows, names, rel, unit, kind, d_list, args):
    print("%6s " % "#" + " ".join("%16s" % n.split(" /")[0] for n in names) +
          " %18s %12s" % ("E / Hartree", "rel/" + UNIT_LABEL[unit]))
    for r, v in zip(rows, rel):
        print("%6s " % r["point"] + " ".join("%16.5f" % x for x in r["x"]) + " %18.8f %12.3f" % (r["energy"], v))
    print("-" * 70)
    imin = min(range(len(rows)), key=lambda i: rows[i]["energy"])
    imax = max(range(len(rows)), key=lambda i: rows[i]["energy"])

    def xs(i):
        return ", ".join("%.4f" % x for x in rows[i]["x"])

    print("lowest  point #%s at (%s): %.3f %s" % (rows[imin]["point"], xs(imin), rel[imin], UNIT_LABEL[unit]))
    print("highest point #%s at (%s): %.3f %s" % (rows[imax]["point"], xs(imax), rel[imax], UNIT_LABEL[unit]))
    if kind == "irc":
        f = HARTREE[unit]
        ts = [r for r in rows if abs(r["x"][0]) < 1e-8]
        if ts:
            ets = ts[0]["energy"]
            print("reverse end (x=%.3f): %.3f %s below TS; forward end (x=%.3f): %.3f %s below TS"
                  % (rows[0]["x"][0], (ets - rows[0]["energy"]) * f, UNIT_LABEL[unit],
                     rows[-1]["x"][0], (ets - rows[-1]["energy"]) * f, UNIT_LABEL[unit]))
        mono = all(rows[i]["energy"] <= rows[i + 1]["energy"] + 1e-7 for i in range(len(rows) - 1)
                   if rows[i + 1]["x"][0] <= 0) and \
            all(rows[i]["energy"] >= rows[i + 1]["energy"] - 1e-7 for i in range(len(rows) - 1)
                if rows[i]["x"][0] >= 0)
        print("energy decreases monotonically from the TS on both sides: %s" % ("yes" if mono else
              "NO - check the path (wrong TS, mode not connected, or a second barrier)"))
        print("NOTE: optimise both endpoints (Opt + Freq) before using their energies.")
    else:
        if len(names) == 1:
            lm = local_maxima(rows)
            if lm:
                print("interior local maxima (TS guesses): " + "; ".join(
                    "#%s at %s = %.3f" % (rows[i]["point"], xs(i), rel[i]) for i in lm))
            else:
                print("no interior local maximum: the barrier lies outside the scanned range or at its end")
        print("TS initial guess: geometry of point #%s (use --write-max FILE.xyz, or "
              "log2xyz.py <log> --converged and take that frame), then Opt=(TS,CalcFC,NoEigenTest) Freq."
              % rows[imax]["point"])
        for d in d_list:
            if d["unconverged"]:
                print("WARNING: %s: scan point(s) %s never reached 'Stationary point found' and were skipped"
                      % (d["file"], ",".join(map(str, sorted(set(d["unconverged"]))))))
            if d["kind"] == "relaxed" and d["scan_total"] and len(d["points"]) < d["scan_total"]:
                print("NOTE: %s: %d of %d scan points converged (job unfinished or failed?)"
                      % (d["file"], len({p["scan_point"] for p in d["points"]}), d["scan_total"]))
    return imax


def build_parser():
    p = argparse.ArgumentParser(
        description="Energy profiles of Gaussian 16 relaxed/rigid scans (1D, 2D) and IRC paths, in "
                    "kcal/mol relative to the minimum; reports the highest scan point as TS guess.  "
                    "提取柔性/刚性扫描与 IRC 能量曲线，作图并报告最高点。",
        epilog="Examples:\n"
               "  python scan_irc_profile.py scan.log\n"
               "  python scan_irc_profile.py scan.log --write-max ts_guess.xyz --csv scan.csv\n"
               "  python scan_irc_profile.py irc.log\n"
               "  python scan_irc_profile.py irc_fwd.log irc_rev.log -o irc.png\n"
               "  python scan_irc_profile.py gic_scan.log --coord \"B 1 5\"",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("logs", nargs="+", help="Gaussian .log/.out file(s); several IRC logs are merged")
    p.add_argument("--mode", choices=["auto", "scan", "irc"], default="auto",
                   help="force scan or IRC parsing (default: from the route)")
    p.add_argument("--coord", action="append", default=[], metavar="SPEC",
                   help="scan coordinate from atoms, e.g. 'B 1 2', 'D 1 2 3 4' (repeat for 2D); "
                        "default: the S lines of the ModRedundant input")
    p.add_argument("--no-unwrap", action="store_true", help="keep dihedrals in (-180,180]")
    p.add_argument("--ref", default="min", help="zero of energy: min (default), first, ts (IRC) or Hartree value")
    p.add_argument("--unit", choices=["kcal", "kj", "ev"], default="kcal", help="energy unit (default kcal)")
    p.add_argument("--flip", action="store_true", help="IRC: invert the sign of the reaction coordinate")
    p.add_argument("--separate", action="store_true", help="IRC: do not merge several logs")
    p.add_argument("-o", "--output", help="figure file (default <first log stem>_profile.png)")
    p.add_argument("--csv", metavar="FILE", help="write the profile to CSV")
    p.add_argument("--write-max", metavar="XYZ", help="write the geometry of the highest point (scans)")
    p.add_argument("--no-plot", action="store_true", help="no figure (matplotlib not needed)")
    p.add_argument("--show", action="store_true", help="interactive window instead of a file")
    p.add_argument("--dpi", type=int, default=150, help="figure resolution (default 150)")
    return p


def process(group, args, kind):
    stem = os.path.splitext(group[0]["file"])[0]
    if kind == "irc":
        rows, ts_e = irc_profile(group, args.flip)
        names = ["reaction coordinate / Bohr amu^1/2"]
        if len(set(round(e, 5) for e in ts_e)) > 1:
            print("WARNING: TS energies differ between the merged logs: %s" % ", ".join("%.6f" % e for e in ts_e))
        if any(r["source"] == "points" for r in rows):
            print("NOTE: no 'Summary of reaction path following' (unfinished job?) - points read one by one")
    else:
        specs = [parse_coord_spec(c) for c in args.coord]
        rows, names = [], []
        for d in group:
            r, names = scan_profile(d, specs, not args.no_unwrap)
            rows += r
        if len(names) == 2 and len(rows) < 3:
            names = names[:1]
    if not rows:
        print("ERROR: %s: no %s points found" % (group[0]["file"], "IRC" if kind == "irc" else "scan"),
              file=sys.stderr)
        return 1
    e0 = reference(rows, args.ref, kind == "irc")
    f = HARTREE[args.unit]
    rel = [(r["energy"] - e0) * f for r in rows]
    print("=" * 70)
    print("%s  [%s, %d points]" % (", ".join(d["file"] for d in group), kind, len(rows)))
    imax = report(rows, names, rel, args.unit, kind, group, args)
    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["point"] + [n for n in names] + ["energy_hartree", "relative_" + args.unit])
            for r, v in zip(rows, rel):
                w.writerow([r["point"]] + ["%.6f" % x for x in r["x"]] + ["%.10f" % r["energy"], "%.4f" % v])
        print("CSV written to %s" % args.csv, file=sys.stderr)
    if args.write_max:
        if kind == "irc":
            print("NOTE: --write-max is for scans; for IRC use log2xyz.py <log> --irc", file=sys.stderr)
        else:
            try:
                write_xyz(rows[imax], args.write_max, "highest scan point %s E=%.10f Hartree source=%s"
                          % (rows[imax]["point"], rows[imax]["energy"], os.path.basename(group[0]["file"])))
                print("geometry of the highest point written to %s" % args.write_max, file=sys.stderr)
            except ValueError as exc:
                print("ERROR: --write-max: %s" % exc, file=sys.stderr)
    if not args.no_plot:
        out = args.output or stem + "_profile.png"
        title = os.path.basename(group[0]["file"]) + (" (+%d)" % (len(group) - 1) if len(group) > 1 else "")
        if len(names) == 2:
            ok = plot_2d(rows, rel, names, args.unit, title, out, args.show, args.dpi)
        else:
            mark = None
            if kind == "irc":
                zero = [i for i, r in enumerate(rows) if abs(r["x"][0]) < 1e-8]
                mark = zero[0] if zero else None
            else:
                mark = imax
            ok = plot_1d(rows, rel, names[0], args.unit, title, out, args.show, args.dpi, mark, kind == "irc")
        if not ok:
            return 1
    return 0


def main(argv=None):
    args = build_parser().parse_args(argv)
    data = []
    for log in args.logs:
        if not os.path.isfile(log):
            print("ERROR: %s not found" % log, file=sys.stderr)
            return 1
        try:
            d = parse_log(log)
        except OSError as exc:
            print("ERROR: cannot read %s: %s" % (log, exc), file=sys.stderr)
            return 1
        if not d["routes"]:
            print("ERROR: %s has no route section - is it a Gaussian output?" % log, file=sys.stderr)
            return 1
        if args.mode == "irc":
            d["kind"] = "irc"
        elif args.mode == "scan" and d["kind"] in (None, "irc"):
            d["kind"] = "relaxed" if d["points"] else "rigid"
        if d["kind"] is None:
            print("ERROR: %s is neither a scan nor an IRC job (use --mode to force)" % log, file=sys.stderr)
            return 1
        data.append(d)
    irc = [d for d in data if d["kind"] == "irc"]
    scans = [d for d in data if d["kind"] != "irc"]
    groups = []
    if irc:
        groups += [([d], "irc") for d in irc] if args.separate else [(irc, "irc")]
    groups += [([d], "scan") for d in scans]
    if len(groups) > 1 and (args.output or args.csv or args.write_max):
        print("ERROR: -o/--csv/--write-max need a single profile (one scan log, or IRC logs merged)",
              file=sys.stderr)
        return 2
    rc = 0
    for group, kind in groups:
        try:
            rc |= process(group, args, kind)
        except ValueError as exc:
            print("ERROR: %s" % exc, file=sys.stderr)
            rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
