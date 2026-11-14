#!/usr/bin/env python3
"""Static pre-check of CP2K inputs written by generators (Multiwfn, ASE, cp2kgen,
old templates) before they are run with a given CP2K version.

Chapter 17 of the CP2K tutorial. This is NOT a replacement for `cp2k.psmp -c`
(which is the only authoritative syntax check of your CP2K build); it points at
the version-dependent constructs that the generators most often get wrong and
explains how to fix them.

Checks
  ERROR  2026.2+: &MOTION/&CELL_OPT/TYPE (keyword removed in 2026.2)
  ERROR  2026.2+: &DFT/&PRINT/&PDOS (moved to &DFT/&PRINT/&DOS/&PDOS)
  ERROR  2026.2+: &REFTRAJ EVAL_ENERGY_FORCES / EVAL_FORCES (replaced by EVAL)
  ERROR  <2026.2: 2026.2-only syntax (&DOS with &PDOS/&LDOS/&R_LDOS/&CURVE/&BROADEN)
  ERROR  <2024.1: REFTRAJ EVAL (the keyword exists since CP2K 2024.1)
  WARN   OT together with a k-point mesh (unsupported up to and including 2026.2)
  WARN   CELL PERIODIC and POISSON PERIODIC differ (both default to XYZ)
  WARN   basis set and potential of a KIND with different -qN
  WARN   HFX CUTOFF_RADIUS larger than half of the smallest cell width
  NOTE   @INCLUDE files are not followed; PRINT_LEVEL LOW (cp2kdata needs MEDIUM)

Limitations: lines inside all @IF/@ELSE branches are scanned; ${VAR} values are not
substituted; cell and periodicity of several &FORCE_EVAL sections are merged.

Usage:
    python inp_migration_check.py generated.inp
    python inp_migration_check.py --target 2025.2 old_*.inp
Exit status 1 if at least one ERROR was found. Python >= 3.8, standard library only.
"""
import argparse
import math
import re
import sys

NEW_DOS_SUBSECTIONS = {"PDOS", "LDOS", "R_LDOS", "CURVE", "BROADEN"}
ON_VALUES = {"", "ON", "T", "TRUE", ".TRUE.", "YES", "Y"}


def parse_version(text):
    m = re.match(r"^(\d{4})\.(\d+)$", text.strip())
    if not m:
        raise argparse.ArgumentTypeError("version must look like 2026.2")
    return int(m.group(1)), int(m.group(2))


def logical_lines(path):
    """Yield (line_number, text) with comments removed and continuations joined."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        buf, start = "", None
        for n, raw in enumerate(fh, 1):
            line = re.split(r"[#!]", raw.rstrip("\n"), maxsplit=1)[0].rstrip()
            if start is None:
                start = n
            if line.endswith("\\"):
                buf += line[:-1] + " "
                continue
            buf += line
            if buf.strip():
                yield start, buf.strip()
            buf, start = "", None
        if buf.strip():
            yield start, buf.strip()


def to_float(token):
    try:
        return float(token.replace("d", "e").replace("D", "e"))
    except ValueError:
        return None


def strip_unit(tokens):
    """Return (unit or None, remaining tokens) for a value list like [angstrom] 10 10 10."""
    if tokens and tokens[0].startswith("[") and tokens[0].endswith("]"):
        return tokens[0][1:-1].lower(), tokens[1:]
    return None, tokens


def cell_widths(cell):
    """Perpendicular widths (A) of the cell, or None if they cannot be determined."""
    if "A" in cell and "B" in cell and "C" in cell:
        a, b, c = cell["A"], cell["B"], cell["C"]

        def cross(u, v):
            return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]

        def norm(u):
            return math.sqrt(sum(x * x for x in u))

        vol = abs(sum(x * y for x, y in zip(a, cross(b, c))))
        if vol < 1e-8:
            return None
        return [vol / norm(cross(b, c)), vol / norm(cross(c, a)), vol / norm(cross(a, b))]
    if "ABC" in cell and "ALPHA_BETA_GAMMA" not in cell:
        return list(cell["ABC"])
    return None


def check_file(path, target):
    msgs = []
    new = target >= (2026, 2)

    def add(level, n, text):
        msgs.append((level, n, text))

    stack = []                  # list of (SECTION_NAME, param)
    cell, cell_line, cell_periodic, poisson_periodic = {}, None, None, None
    kpoints_line, ot_line = None, None
    hfx_radius = None           # (value_in_A, line)
    kind = None                 # dict for the current &KIND
    mult_cell = None

    for n, line in logical_lines(path):
        tokens = line.split()
        head = tokens[0].upper()
        path_names = [s[0] for s in stack]

        if head.startswith("@"):
            if head == "@INCLUDE":
                add("NOTE", n, f"@INCLUDE {' '.join(tokens[1:])}: included file is not checked")
            continue

        if head.startswith("&END"):
            if stack:
                name = stack.pop()[0]
                if name == "KIND" and kind is not None:
                    qb = re.search(r"-q(\d+)$", kind.get("BASIS_SET", ""))
                    qp = re.search(r"-q(\d+)$", kind.get("POTENTIAL", ""))
                    if qb and qp and qb.group(1) != qp.group(1):
                        add("WARN", kind["line"],
                            f"KIND {kind['name']}: basis -q{qb.group(1)} but potential "
                            f"-q{qp.group(1)}; valence electrons must match")
                    kind = None
            continue

        if head.startswith("&"):
            name = head[1:]
            param = " ".join(tokens[1:]).upper()
            parent = path_names
            # 1) old PDOS syntax directly under &DFT/&PRINT
            if name == "PDOS" and parent[-2:] == ["DFT", "PRINT"] and new:
                add("ERROR", n, "&DFT/&PRINT/&PDOS was replaced in 2026.2: write "
                    "&PRINT / &DOS / &PDOS ... &END PDOS / &END DOS "
                    "(CP2K stops with 'unknown subsection PDOS of section PRINT')")
            # 2) new DOS syntax used for an old target
            if (name in NEW_DOS_SUBSECTIONS and not new
                    and any(parent[i:i + 3] == ["DFT", "PRINT", "DOS"] for i in range(len(parent)))):
                add("ERROR", n, f"&{name} inside &DFT/&PRINT/&DOS is 2026.2 syntax; for "
                    f"{target[0]}.{target[1]} use &DFT/&PRINT/&PDOS")
            if name == "OT" and param in ON_VALUES:
                ot_line = n
            if name == "KIND":
                kind = {"name": param, "line": n}
            stack.append((name, param))
            continue

        if not stack:
            continue
        sec = stack[-1][0]
        values = tokens[1:]

        # 3) CELL_OPT TYPE
        if sec == "CELL_OPT" and head == "TYPE" and new:
            extra = ""
            if values and values[0].upper() != "DIRECT_CELL_OPT":
                extra = (f" The old value {values[0]} is no longer available: 2026.2 always "
                         "uses DIRECT_CELL_OPT, check the optimizer settings.")
            add("ERROR", n, "&CELL_OPT TYPE was removed in 2026.2 (CP2K stops with "
                "\"The specified keyword 'TYPE' is not available anymore\"); delete the line."
                + extra)
        # 4) REFTRAJ evaluation keywords
        if sec == "REFTRAJ":
            if head in ("EVAL_ENERGY_FORCES", "EVAL_FORCES") and new:
                add("ERROR", n, f"{head} was replaced in 2026.2 by the keyword EVAL "
                    "(e.g. EVAL ENERGY_FORCES); check the allowed values in the manual")
            if head == "EVAL" and target < (2024, 1):
                add("ERROR", n, "REFTRAJ EVAL exists only since CP2K 2024.1; older versions use "
                    "EVAL_ENERGY_FORCES / EVAL_FORCES")
        # k-points
        if sec == "KPOINTS" and head == "SCHEME" and values and values[0].upper() not in ("NONE", "GAMMA"):
            kpoints_line = n
        # periodicity
        if sec == "CELL" and head == "PERIODIC" and values:
            cell_periodic = (values[0].upper(), n)
        if sec == "POISSON" and head == "PERIODIC" and values:
            poisson_periodic = (values[0].upper(), n)
        # cell geometry (only plain Angstrom values are evaluated)
        if sec == "CELL" and len(stack) >= 2 and stack[-2][0] == "SUBSYS":
            unit, vals = strip_unit(values)
            if head in ("ABC", "A", "B", "C") and unit in (None, "angstrom"):
                nums = [to_float(v) for v in vals[:3]]
                if len(nums) == 3 and None not in nums:
                    cell[head] = nums
                    cell_line = n
            if head == "ALPHA_BETA_GAMMA":
                cell[head] = True
            if head == "MULTIPLE_UNIT_CELL":
                nums = [to_float(v) for v in vals[:3]]
                if len(nums) == 3 and None not in nums:
                    mult_cell = nums
        # HFX truncation radius
        if sec == "INTERACTION_POTENTIAL" and head == "CUTOFF_RADIUS":
            unit, vals = strip_unit(values)
            r = to_float(vals[0]) if vals else None
            if r is not None and unit in (None, "angstrom"):
                hfx_radius = (r, n)
        # basis / potential per KIND
        if sec == "KIND" and kind is not None and head in ("BASIS_SET", "POTENTIAL") and values:
            # BASIS_SET may carry a type first (ORB, AUX_FIT, RI_AUX ...): keep the ORB / single value
            if head == "BASIS_SET" and len(values) > 1:
                if values[0].upper() != "ORB":
                    continue
                kind[head] = values[1]
            else:
                kind[head] = values[-1]
        if sec == "GLOBAL" and head == "PRINT_LEVEL" and values and values[0].upper() == "LOW":
            add("NOTE", n, "PRINT_LEVEL LOW: cp2kdata refuses LOW outputs, use MEDIUM if you "
                "want to parse the output with cp2kdata/dpdata")

    # file-level checks
    if kpoints_line and ot_line and target <= (2026, 2):
        add("WARN", ot_line, f"OT together with k-points (line {kpoints_line}) is not supported "
            "up to CP2K 2026.2; use &DIAGONALIZATION + &MIXING (+ &SMEAR, ADDED_MOS)")
    # both keywords default to XYZ when they are not given
    cp = cell_periodic[0] if cell_periodic else "XYZ"
    pp = poisson_periodic[0] if poisson_periodic else "XYZ"
    if cp != pp:
        if poisson_periodic is None:
            hint = ("add &POISSON PERIODIC NONE + POISSON_SOLVER MT/WAVELET" if cp == "NONE"
                    else f"add &POISSON PERIODIC {cp} with a matching POISSON_SOLVER")
            add("WARN", cell_periodic[1], f"CELL PERIODIC {cp} but no &POISSON PERIODIC "
                f"(default XYZ); {hint}")
        else:
            where = f"line {cell_periodic[1]}" if cell_periodic else "default"
            add("WARN", poisson_periodic[1], f"POISSON PERIODIC {pp} differs from "
                f"CELL PERIODIC {cp} ({where})")
    if hfx_radius:
        widths = cell_widths(cell)
        if widths:
            if mult_cell:
                widths = [w * m for w, m in zip(widths, mult_cell)]
            half = min(widths) / 2.0
            if hfx_radius[0] > half + 1e-6:
                add("WARN", hfx_radius[1], f"CUTOFF_RADIUS {hfx_radius[0]:.2f} A > half of the "
                    f"smallest cell width {half:.2f} A (cell line {cell_line}); enlarge the cell "
                    "or reduce the radius (k-point runs: judge with the Born-von Karman supercell)")
    return msgs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="CP2K input files")
    ap.add_argument("--target", type=parse_version, default=(2026, 2),
                    help="CP2K version the input will be run with (default 2026.2)")
    args = ap.parse_args()

    n_err = 0
    for f in args.inputs:
        try:
            msgs = check_file(f, args.target)
        except OSError as exc:
            print(f"{f}: cannot read ({exc})")
            n_err += 1
            continue
        for level, n, text in sorted(msgs, key=lambda m: m[1]):
            print(f"{f}:{n}: {level}: {text}")
        n_err += sum(1 for m in msgs if m[0] == "ERROR")
        if not msgs:
            print(f"{f}: no version issues found for {args.target[0]}.{args.target[1]}")
    print(f"-- {n_err} error(s). Always confirm with: cp2k.psmp -c <input>")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())
