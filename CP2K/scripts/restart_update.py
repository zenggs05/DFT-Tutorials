#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
restart_update.py - put the latest CP2K coordinates/cell back into an input file or an xyz file.

从 CP2K 的 <PROJECT>-1.restart（或 -pos-1.xyz 轨迹的最后一帧 / 指定帧，可配合 -1.cell）
读取最新的坐标和晶胞，然后：
  (1) 写成 xyz / 扩展 xyz（带 Lattice），供可视化或后续计算使用；
  (2) 写回一个已有的 CP2K 输入（替换 &SUBSYS 中的 &COORD 与 &CELL 的 A/B/C，
      或在使用 &TOPOLOGY/COORD_FILE_NAME 时写出新的 xyz 并改文件名），
      其余设置原样保留，便于在优化中断后续算、或把优化结构交给能带/振动计算。

Usage examples
--------------
  # last geometry of an optimization -> xyz
  python restart_update.py h2o_opt-1.restart --xyz h2o_opt_final.xyz

  # write the optimized coordinates + cell into a new input based on an old one
  python restart_update.py si_cell-1.restart -i si_band.inp -o si_band_opt.inp

  # last frame of a trajectory, cell taken from the matching line of the .cell file
  python restart_update.py si_cell-pos-1.xyz --cell-file si_cell-1.cell -i si.inp -o si_new.inp

  # frame 10 (0-based) of a trajectory, extended xyz with the cell
  python restart_update.py md-pos-1.xyz --frame 10 --cell-file md-1.cell --xyz frame10.xyz --extxyz

  # continue the same run: copy coordinates AND add an &EXT_RESTART block
  python restart_update.py h2o_opt-1.restart -i h2o_opt.inp -o h2o_opt_cont.inp --ext-restart

Notes
-----
* Coordinates are written in Angstrom (Cartesian).  Atom labels are kept exactly as in the
  restart / trajectory (they are CP2K kind names, e.g. Fe1/Fe2 for an AFM setup); extra
  columns of &COORD lines (molecule name, residue index) are kept as well.
* The cell is written as A/B/C vectors in Angstrom; ABC, ALPHA_BETA_GAMMA, A, B, C,
  CELL_FILE_NAME and CELL_FILE_FORMAT at the top level of &CELL are replaced, while PERIODIC,
  SYMMETRY, MULTIPLE_UNIT_CELL and &CELL_REF are kept.
* A restart file contains much more than the structure (velocities, thermostat state, the
  optimizer history, the wavefunction file name ...).  To *continue the same run* the cleanest
  way is &EXT_RESTART (option --ext-restart adds it).  This script is meant for moving a
  structure into a new or modified input.
* Ordinary .xyz trajectories carry no cell; for CELL_OPT/NPT use --cell-file <PROJECT>-1.cell
  (columns: Step Time Ax Ay Az Bx By Bz Cx Cy Cz Volume) or an extended xyz with Lattice="...".
* With SPLIT_RESTART_FILE the coordinates are not stored in the .restart text; use the
  -pos-1.xyz trajectory instead.

Dependencies: Python >= 3.9, standard library only.

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import collections
import math
import os
import re
import sys

BOHR_ANG = 0.529177210903
LENGTH_UNITS = {"angstrom": 1.0, "bohr": BOHR_ANG, "nm": 10.0, "pm": 0.01, "a.u.": BOHR_ANG, "au": BOHR_ANG}
CELL_KEYS_REPLACED = {"ABC", "ALPHA_BETA_GAMMA", "A", "B", "C", "CELL_FILE_NAME", "CELL_FILE_FORMAT"}


class Structure:
    def __init__(self):
        self.labels = []      # kind labels
        self.positions = []   # Cartesian, Angstrom
        self.extra = []       # extra tokens of &COORD lines
        self.cell = None      # [[ax,ay,az],[bx,by,bz],[cx,cy,cz]] in Angstrom
        self.source = ""
        self.step = None


# ----------------------------------------------------------------------------- helpers
def strip_comment(line):
    """Remove CP2K comments (# or ! outside quotes)."""
    out, quote = [], None
    for ch in line:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "\"'":
            quote = ch
        elif ch in "#!":
            break
        out.append(ch)
    return "".join(out)


def parse_unit_token(tokens, default="angstrom"):
    """Pop a leading [unit] token. Returns (factor_to_angstrom, remaining tokens)."""
    if tokens and tokens[0].startswith("["):
        unit = tokens[0].strip("[]").lower()
        if unit not in LENGTH_UNITS:
            raise SystemExit(f"ERROR: unsupported length unit [{unit}] (use angstrom or bohr)")
        return LENGTH_UNITS[unit], tokens[1:]
    return LENGTH_UNITS[default], tokens


def cell_from_abc(abc, angles):
    a, b, c = abc
    al, be, ga = (math.radians(x) for x in angles)
    va = [a, 0.0, 0.0]
    vb = [b * math.cos(ga), b * math.sin(ga), 0.0]
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return [va, vb, [cx, cy, cz]]


def frac_to_cart(frac, cell):
    return [sum(frac[k] * cell[k][i] for k in range(3)) for i in range(3)]


def section_name(stripped):
    """Return ('open', NAME) / ('close', NAME or '') / None for a line."""
    if not stripped.startswith("&"):
        return None
    tok = stripped[1:].split()
    if not tok:
        return None
    name = tok[0].upper()
    if name == "END":
        return ("close", tok[1].upper() if len(tok) > 1 else "")
    return ("open", name)


# ----------------------------------------------------------------------------- readers
def read_restart(path):
    """Read cell and coordinates of the first FORCE_EVAL/SUBSYS of a .restart (input-format) file."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    st = Structure()
    st.source = path
    stack, nsubsys = [], 0
    coord_scaled, coord_fac = False, 1.0
    raw_coords, cell_vec, abc, angles = [], {}, None, [90.0, 90.0, 90.0]
    in_first = False
    for line in lines:
        s = strip_comment(line).strip()
        if not s:
            continue
        sec = section_name(s)
        if sec:
            if sec[0] == "open":
                stack.append(sec[1])
                if sec[1] == "SUBSYS" and "FORCE_EVAL" in stack:
                    nsubsys += 1
                    in_first = nsubsys == 1
            else:
                if stack:
                    stack.pop()
            continue
        if not in_first:
            continue
        if stack and stack[-1] == "COORD" and len(stack) >= 2 and stack[-2] == "SUBSYS":
            tok = s.split()
            key = tok[0].upper()
            if key == "UNIT":
                unit = tok[1].lower() if len(tok) > 1 else "angstrom"
                if unit not in LENGTH_UNITS:
                    raise SystemExit(f"ERROR: {path}: unsupported &COORD UNIT {unit}")
                coord_fac = LENGTH_UNITS[unit]
                continue
            if key == "SCALED":
                coord_scaled = len(tok) == 1 or tok[1].upper() in ("T", "TRUE", ".TRUE.", "YES", "ON")
                continue
            if len(tok) < 4:
                raise SystemExit(f"ERROR: {path}: malformed &COORD line: {line.strip()}")
            try:
                xyz = [float(t) for t in tok[1:4]]
            except ValueError:
                raise SystemExit(f"ERROR: {path}: cannot read coordinates in line: {line.strip()}")
            raw_coords.append((tok[0], xyz, tok[4:]))
        elif stack and stack[-1] == "CELL" and len(stack) >= 2 and stack[-2] == "SUBSYS":
            tok = s.split()
            key = tok[0].upper()
            if key in ("A", "B", "C"):
                fac, rest = parse_unit_token(tok[1:])
                cell_vec[key] = [float(x) * fac for x in rest[:3]]
            elif key == "ABC":
                fac, rest = parse_unit_token(tok[1:])
                abc = [float(x) * fac for x in rest[:3]]
            elif key == "ALPHA_BETA_GAMMA":
                rest = [t for t in tok[1:] if not t.startswith("[")]  # angles are in degrees
                angles = [float(x) for x in rest[:3]]
            elif key == "CELL_FILE_NAME":
                print(f"WARNING: {path}: &CELL uses CELL_FILE_NAME; cell not read", file=sys.stderr)
    if nsubsys > 1:
        print(f"WARNING: {path}: {nsubsys} &SUBSYS sections found; using the first one", file=sys.stderr)
    if len(cell_vec) == 3:
        st.cell = [cell_vec["A"], cell_vec["B"], cell_vec["C"]]
    elif abc:
        st.cell = cell_from_abc(abc, angles)
    if not raw_coords:
        raise SystemExit(f"ERROR: {path}: no &COORD block found in &FORCE_EVAL/&SUBSYS "
                         "(SPLIT_RESTART_FILE used? then read the -pos-1.xyz trajectory instead)")
    for label, xyz, extra in raw_coords:
        if coord_scaled:
            if st.cell is None:
                raise SystemExit(f"ERROR: {path}: SCALED coordinates but no cell found")
            pos = frac_to_cart(xyz, st.cell)
        else:
            pos = [x * coord_fac for x in xyz]
        st.labels.append(label)
        st.positions.append(pos)
        st.extra.append(extra)
    return st


RE_LATTICE = re.compile(r'Lattice\s*=\s*"([^"]+)"', re.IGNORECASE)
RE_STEP = re.compile(r"\bi\s*=\s*(\d+)")


def read_xyz_frame(path, frame=-1):
    """Read one frame of a (multi-frame) xyz / extended xyz file without loading the whole file."""
    keep = collections.deque(maxlen=1) if frame < 0 else None
    want = None if frame < 0 else frame
    nframes = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            if not head.strip():
                continue
            try:
                nat = int(head.split()[0])
            except ValueError:
                raise SystemExit(f"ERROR: {path}: expected an atom count, got '{head.strip()[:40]}' "
                                 f"(frame {nframes})")
            comment = fh.readline()
            body = [fh.readline() for _ in range(nat)]
            if len(body) < nat or (nat and not body[-1]):
                print(f"WARNING: {path}: last frame is incomplete (job still running?); ignored", file=sys.stderr)
                break
            if keep is not None:
                keep.append((nframes, comment, body))
            elif nframes == want:
                keep = [(nframes, comment, body)]
                break
            nframes += 1
    if not keep:
        raise SystemExit(f"ERROR: {path}: frame {frame} not found ({nframes} complete frames)")
    idx, comment, body = keep[-1] if frame < 0 else keep[0]
    st = Structure()
    st.source = f"{path} (frame {idx})"
    m = RE_STEP.search(comment)
    st.step = int(m.group(1)) if m else None
    m = RE_LATTICE.search(comment)
    if m:
        vals = [float(x) for x in m.group(1).split()]
        if len(vals) == 9:
            st.cell = [vals[0:3], vals[3:6], vals[6:9]]
    for ln in body:
        tok = ln.split()
        if len(tok) < 4:
            raise SystemExit(f"ERROR: {path}: malformed atom line '{ln.strip()}' in frame {idx}")
        st.labels.append(tok[0])
        st.positions.append([float(t) for t in tok[1:4]])
        st.extra.append([])
    return st


def read_cell_file(path, step=None):
    """Return (cell, step) from a <PROJECT>-1.cell file (matching step if given, else the last line)."""
    rows = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for ln, line in enumerate(fh, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            tok = s.split()
            if len(tok) < 11:
                raise SystemExit(f"ERROR: {path}:{ln}: expected Step Time Ax Ay Az Bx By Bz Cx Cy Cz [Volume]")
            rows.append((int(tok[0]), [float(x) for x in tok[2:11]]))
    if not rows:
        raise SystemExit(f"ERROR: {path}: no data lines")
    chosen = rows[-1]
    if step is not None:
        match = [r for r in rows if r[0] == step]
        if match:
            chosen = match[-1]
        else:
            print(f"WARNING: step {step} not in {path}; using its last line (step {rows[-1][0]})", file=sys.stderr)
    v = chosen[1]
    return [v[0:3], v[3:6], v[6:9]], chosen[0]


# ----------------------------------------------------------------------------- writers
def coord_lines(st, indent):
    out = []
    for lab, pos, extra in zip(st.labels, st.positions, st.extra):
        tail = ("  " + " ".join(extra)) if extra else ""
        out.append(f"{indent}{lab:<4s} {pos[0]:18.10f} {pos[1]:18.10f} {pos[2]:18.10f}{tail}")
    return out


def cell_lines(cell, indent):
    return [f"{indent}{n} [angstrom] {v[0]:18.10f} {v[1]:18.10f} {v[2]:18.10f}" for n, v in zip("ABC", cell)]


def write_xyz(st, path, extxyz=False):
    with open(path, "w", encoding="ascii", errors="replace", newline="\n") as fh:
        fh.write(f"{len(st.labels)}\n")
        if extxyz and st.cell:
            lat = " ".join(f"{x:.10f}" for v in st.cell for x in v)
            fh.write(f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T"\n')
        else:
            fh.write(f"from {os.path.basename(st.source.split(' (')[0])}"
                     + (f", step {st.step}" if st.step is not None else "") + "\n")
        for lab, pos in zip(st.labels, st.positions):
            fh.write(f"{lab:<4s} {pos[0]:18.10f} {pos[1]:18.10f} {pos[2]:18.10f}\n")


def update_input(template, st, out_inp, xyz_name, ext_restart, force):
    """Return the text of an updated input file."""
    with open(template, "r", encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    # locate sections in the first FORCE_EVAL/SUBSYS
    stack, spans = [], {}
    subsys_open = None
    for i, line in enumerate(lines):
        s = strip_comment(line).strip()
        sec = section_name(s)
        if not sec:
            continue
        if sec[0] == "open":
            stack.append((sec[1], i))
            path = [n for n, _ in stack]
            if sec[1] == "SUBSYS" and "FORCE_EVAL" in path and subsys_open is None:
                subsys_open = i
        else:
            if not stack:
                continue
            name, start = stack.pop()
            path = [n for n, _ in stack]
            if subsys_open is not None and "SUBSYS" not in spans:
                if name == "SUBSYS" and start == subsys_open:
                    spans["SUBSYS"] = (start, i)
                elif path and path[-1] == "SUBSYS" and name in ("COORD", "CELL", "TOPOLOGY") and name not in spans:
                    spans[name] = (start, i)
    if "SUBSYS" not in spans:
        raise SystemExit(f"ERROR: {template}: no &FORCE_EVAL/&SUBSYS section found "
                         "(an @INCLUDE'd subsys is not supported)")
    indent_sub = re.match(r"\s*", lines[spans["SUBSYS"][0]]).group(0) + "  "
    edits = []  # (start, end_exclusive, new_lines)
    natoms_template = None

    # ---- coordinates
    if "COORD" in spans:
        a, b = spans["COORD"]
        old = [strip_comment(x).strip() for x in lines[a + 1:b]]
        old_atoms = [x for x in old if x and x.split()[0].upper() not in ("UNIT", "SCALED")]
        natoms_template = len(old_atoms)
        ind = re.match(r"\s*", lines[a]).group(0) + "  "
        edits.append((a + 1, b, coord_lines(st, ind)))
        coord_mode = "&COORD replaced"
    elif "TOPOLOGY" in spans:
        a, b = spans["TOPOLOGY"]
        name_idx, fmt_idx, fmt_val = None, None, None
        for k in range(a + 1, b):
            tok = strip_comment(lines[k]).split()
            if tok and tok[0].upper() == "COORD_FILE_NAME":
                name_idx = k
            elif tok and tok[0].upper() == "COORD_FILE_FORMAT":
                fmt_idx, fmt_val = k, (tok[1].upper() if len(tok) > 1 else "")
        if name_idx is None:
            ind = re.match(r"\s*", lines[spans["SUBSYS"][0]]).group(0) + "  "
            edits.append((spans["SUBSYS"][0] + 1, spans["SUBSYS"][0] + 1,
                          [f"{ind}&COORD"] + coord_lines(st, ind + "  ") + [f"{ind}&END COORD"]))
            coord_mode = "&COORD inserted"
        else:
            if fmt_val not in (None, "XYZ") and not force:
                raise SystemExit(f"ERROR: {template}: COORD_FILE_FORMAT {fmt_val} would be replaced by XYZ, "
                                 "which can lose connectivity/molecule information; use --force to do it anyway")
            ind = re.match(r"\s*", lines[name_idx]).group(0)
            edits.append((name_idx, name_idx + 1, [f"{ind}COORD_FILE_NAME {xyz_name}"]))
            if fmt_idx is None:
                edits.append((name_idx + 1, name_idx + 1, [f"{ind}COORD_FILE_FORMAT XYZ"]))
            elif fmt_val != "XYZ":
                edits.append((fmt_idx, fmt_idx + 1, [f"{ind}COORD_FILE_FORMAT XYZ"]))
            coord_mode = f"COORD_FILE_NAME -> {xyz_name}"
    else:
        s0 = spans["SUBSYS"][0]
        edits.append((s0 + 1, s0 + 1, [f"{indent_sub}&COORD"] + coord_lines(st, indent_sub + "  ")
                      + [f"{indent_sub}&END COORD"]))
        coord_mode = "&COORD inserted"
    if natoms_template is not None and natoms_template != len(st.labels) and not force:
        raise SystemExit(f"ERROR: template has {natoms_template} atoms in &COORD but the structure has "
                         f"{len(st.labels)}; use --force if this is intended")

    # ---- cell
    cell_mode = "cell unchanged (no cell information)"
    if st.cell is not None:
        if "CELL" in spans:
            a, b = spans["CELL"]
            depth, remove = 0, []
            for k in range(a + 1, b):
                s = strip_comment(lines[k]).strip()
                sec = section_name(s)
                if sec:
                    depth += 1 if sec[0] == "open" else -1
                    continue
                tok = s.split()
                if depth == 0 and tok:
                    key = tok[0].upper()
                    if key in CELL_KEYS_REPLACED:
                        remove.append(k)
                    elif key == "MULTIPLE_UNIT_CELL" and tok[1:4] != ["1", "1", "1"]:
                        print("WARNING: &CELL has MULTIPLE_UNIT_CELL != 1 1 1; the A/B/C written now are the "
                              "full cell of the restart - check the result", file=sys.stderr)
            for k in remove:
                edits.append((k, k + 1, []))
            ind = re.match(r"\s*", lines[a]).group(0) + "  "
            edits.append((a + 1, a + 1, cell_lines(st.cell, ind)))
            cell_mode = "&CELL A/B/C replaced"
        else:
            s0 = spans["SUBSYS"][0]
            edits.append((s0 + 1, s0 + 1, [f"{indent_sub}&CELL"] + cell_lines(st.cell, indent_sub + "  ")
                          + [f"{indent_sub}&END CELL"]))
            cell_mode = "&CELL inserted (PERIODIC default XYZ)"

    # apply edits from the bottom up
    new = list(lines)
    for start, end, repl in sorted(edits, key=lambda e: (e[0], e[1]), reverse=True):
        new[start:end] = repl
    text = "\n".join(new) + "\n"
    if ext_restart:
        if re.search(r"^\s*&EXT_RESTART", text, re.IGNORECASE | re.MULTILINE):
            print("WARNING: the template already has &EXT_RESTART; not added again", file=sys.stderr)
        elif not st.source.endswith(".restart"):
            print("WARNING: --ext-restart needs a .restart file as source; ignored", file=sys.stderr)
        else:
            text += ("\n&EXT_RESTART\n"
                     f"  RESTART_FILE_NAME {os.path.basename(st.source)}\n"
                     "&END EXT_RESTART\n")
    return text, coord_mode, cell_mode


# ----------------------------------------------------------------------------- CLI
def build_parser():
    p = argparse.ArgumentParser(
        description="Take coordinates/cell from a CP2K *-1.restart file (or a frame of *-pos-1.xyz) and "
                    "write them to an xyz file and/or into a copy of an existing input file.",
        epilog="Examples:\n"
               "  python restart_update.py opt-1.restart --xyz final.xyz\n"
               "  python restart_update.py opt-1.restart -i old.inp -o new.inp\n"
               "  python restart_update.py cell-pos-1.xyz --cell-file cell-1.cell -i old.inp -o new.inp\n"
               "  python restart_update.py opt-1.restart -i opt.inp -o opt_cont.inp --ext-restart",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("source", help="<PROJECT>-1.restart (input format) or an xyz / extended-xyz trajectory")
    p.add_argument("-i", "--input", metavar="INP", help="existing CP2K input to update (it is not modified)")
    p.add_argument("-o", "--output", metavar="NEW_INP",
                   help="file name of the updated input (default: <input stem>_updated.inp)")
    p.add_argument("--xyz", metavar="FILE",
                   help="write the structure as xyz; with -i and COORD_FILE_NAME this is the file the new "
                        "input points to (default then: <new input stem>.xyz)")
    p.add_argument("--extxyz", action="store_true", help="write an extended xyz with Lattice=\"...\"")
    p.add_argument("--frame", type=int, default=-1,
                   help="trajectory frame, 0-based; -1 = last complete frame (default)")
    p.add_argument("--cell-file", metavar="CELL",
                   help="<PROJECT>-1.cell file providing the cell for a trajectory frame "
                        "(the line with the same step i = ... is used, else the last line)")
    p.add_argument("--no-cell", action="store_true", help="do not change the cell of the input")
    p.add_argument("--ext-restart", action="store_true",
                   help="also append &EXT_RESTART RESTART_FILE_NAME <source> (source must be a .restart)")
    p.add_argument("--force", action="store_true",
                   help="allow a different number of atoms / overwriting files / non-XYZ COORD_FILE_FORMAT")
    return p


def guard_overwrite(path, force, protected=()):
    if os.path.abspath(path) in [os.path.abspath(x) for x in protected if x]:
        raise SystemExit(f"ERROR: refusing to overwrite the input/source file {path}")
    if os.path.exists(path) and not force:
        raise SystemExit(f"ERROR: {path} exists; choose another name or use --force")


def main(argv=None):
    args = build_parser().parse_args(argv)
    if not os.path.isfile(args.source):
        raise SystemExit(f"ERROR: source file not found: {args.source}")
    if args.input and not os.path.isfile(args.input):
        raise SystemExit(f"ERROR: input file not found: {args.input}")

    if args.frame < -1:
        raise SystemExit("ERROR: only --frame -1 (last) or a non-negative index is supported")
    low = args.source.lower()
    if low.endswith(".xyz") or low.endswith(".extxyz"):
        st = read_xyz_frame(args.source, args.frame)
        if args.cell_file:
            st.cell, cstep = read_cell_file(args.cell_file, st.step)
            print(f"cell from {args.cell_file} (step {cstep})")
    else:
        if args.cell_file:
            print("WARNING: --cell-file ignored for a .restart source (the restart contains the cell)",
                  file=sys.stderr)
        st = read_restart(args.source)
    if args.no_cell:
        st.cell = None
    print(f"read {len(st.labels)} atoms from {st.source}" + (f", step {st.step}" if st.step is not None else ""))
    if st.cell:
        vol = abs(st.cell[0][0] * (st.cell[1][1] * st.cell[2][2] - st.cell[1][2] * st.cell[2][1])
                  - st.cell[0][1] * (st.cell[1][0] * st.cell[2][2] - st.cell[1][2] * st.cell[2][0])
                  + st.cell[0][2] * (st.cell[1][0] * st.cell[2][1] - st.cell[1][1] * st.cell[2][0]))
        lens = [math.sqrt(sum(x * x for x in v)) for v in st.cell]
        print("cell |a| |b| |c| = {:.6f} {:.6f} {:.6f} Angstrom, volume = {:.4f} Angstrom^3".format(*lens, vol))
        if any(abs(x) > 1e-6 for x in (st.cell[0][1], st.cell[0][2], st.cell[1][2])):
            print("NOTE: the cell is not in CP2K's lower-triangular orientation (a along x, b in xy); "
                  "fine for ENERGY/GEO_OPT, but CELL_OPT expects it (see &CELL/CANONICALIZE in 2026.2)")
    elif not args.no_cell:
        print("no cell information found (plain xyz without --cell-file)")

    protected = [args.source, args.input]
    if args.input:
        out_inp = args.output or os.path.splitext(args.input)[0] + "_updated.inp"
        guard_overwrite(out_inp, args.force, protected)
        xyz_name = args.xyz or os.path.splitext(out_inp)[0] + ".xyz"
        text, cmode, cellmode = update_input(args.input, st, out_inp, os.path.basename(xyz_name),
                                             args.ext_restart, args.force)
        if "COORD_FILE_NAME" in cmode:
            guard_overwrite(xyz_name, args.force, protected)
            write_xyz(st, xyz_name, args.extxyz)
            print(f"written: {xyz_name}")
        with open(out_inp, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        print(f"written: {out_inp}  ({cmode}; {cellmode})")
        if args.xyz and "COORD_FILE_NAME" not in cmode:
            guard_overwrite(args.xyz, args.force, protected)
            write_xyz(st, args.xyz, args.extxyz)
            print(f"written: {args.xyz}")
        print("Check the new input with:  cp2k.psmp -c " + os.path.basename(out_inp))
    else:
        xyz = args.xyz or os.path.splitext(os.path.basename(args.source))[0] + "_last.xyz"
        guard_overwrite(xyz, args.force, protected)
        write_xyz(st, xyz, args.extxyz)
        print(f"written: {xyz}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
