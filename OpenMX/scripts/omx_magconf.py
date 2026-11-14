#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_magconf.py -- write initial spin configurations into an OpenMX input (.dat).
                  根据价电子数与目标初始磁矩，自动改写 Atoms.SpeciesAndCoordinates 坐标行的
                  上/下自旋初始电子数（第 6、7 列），支持 FM / AFM / 逐原子自定义；非共线计算
                  额外写第 8-13 列（自旋 Euler 角、轨道 Euler 角、约束开关、+U 轨道极化开关），
                  并检查每个原子的电子数守恒。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Rules implemented (OpenMX manual Sec. 8.2 "Atoms.SpeciesAndCoordinates", Sec. 37)
---------------------------------------------------------------------------------
* columns 6/7 = initial up/down electrons; up + down must equal the valence
  charge of the species' VPS (e.g. Fe_PBE19H 16, Fe_PBE19S 14, Ni_PBE19H 18);
  up - down = initial spin moment (muB).  Charged cells are handled with
  scf.system.charge, never by changing these numbers.
* non-collinear (scf.SpinPolarization nc): columns 8,9 = Euler angles
  theta, phi (deg) of the spin moment, 10,11 = angles of the orbital moment,
  12 = 1/0 switch for the constraint / Zeeman terms, 13 = on/off orbital
  polarization for DFT+U.  As in the official Cr2_CNC.dat example a negative
  moment (down > up) is allowed together with the angles.
* OpenMX reports moments in its own convention (sign opposite to the physical
  magnetic moment of electrons); only relative signs matter for FM/AFM.

Valence charges come from DATA.PATH/VPS/<vps>.vps ("valence.electron") when the
directory exists (DATA.PATH of the input or --data-path), otherwise from the
DFT_DATA19 table in omx_make_input.py of the same folder, or from --valence.

Usage examples
--------------
  # ferromagnetic bcc Fe, 3 muB per Fe (species names or element symbols are accepted)
  python omx_magconf.py fe.dat --mode fm --moment Fe=3.0 -o fe_fm.dat
  # AFM: alternate signs over the Ni atoms in file order, O stays non-magnetic
  python omx_magconf.py nio.dat --mode afm --moment Ni=1.8 -o nio_afm.dat
  # AFM with an explicit sign pattern for the selected atoms 1-4
  python omx_magconf.py cell.dat --mode afm --moment Cr=3 --atoms 1-4 --signs "+ - - +" --inplace
  # per-atom list (VASP MAGMOM style)
  python omx_magconf.py mno.dat --mode custom --magmom "2*4.5 2*-4.5 4*0" -o mno_afm.dat
  # non-collinear: spins in the xz plane, constraint switch on for magnetic atoms
  python omx_magconf.py cr2.dat --mode afm --moment Cr=2 --nc --theta 90 --phi 0 --constrain
  # only check the electron count of an existing file
  python omx_magconf.py nio.dat --check

Dependencies: python >= 3.9, standard library only.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys

TOL = 1.0e-6


# ---------------------------------------------------------------------------
# valence charges
# ---------------------------------------------------------------------------
def _table_valence():
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    try:
        from omx_make_input import valence_of  # type: ignore
        return valence_of
    except Exception:
        return None


def _vps_file_valence(vps: str, data_path: str | None):
    if not data_path:
        return None
    path = os.path.join(data_path, "VPS", vps + ".vps")
    if not os.path.isfile(path):
        return None
    with open(path, encoding="ascii", errors="replace") as fh:
        for line in fh:
            tok = line.split()
            if len(tok) >= 2 and tok[0].lower() == "valence.electron":
                try:
                    return float(tok[1])
                except ValueError:
                    return None
    return None


def element_of(species: str, vps: str | None) -> str:
    """Element symbol from the VPS name (Fe_PBE19H -> Fe) or from the species label."""
    src = vps if vps else species
    m = re.match(r"([A-Z][a-z]?)", src)
    return m.group(1) if m else species


# ---------------------------------------------------------------------------
# input file handling
# ---------------------------------------------------------------------------
def find_block(lines, name):
    """(start, end) line indices of <name ... name> (case-insensitive), or None."""
    start = end = None
    low = name.lower()
    for i, line in enumerate(lines):
        s = line.split("#", 1)[0].strip().lower()
        if start is None and s.startswith("<" + low) and s[1:].split()[0] == low:
            start = i
        elif start is not None and s.endswith(low + ">"):
            end = i
            break
    if start is None or end is None:
        return None
    return start, end


def keyword_value(lines, key):
    low = key.lower()
    for line in lines:
        tok = line.split("#", 1)[0].split()
        if len(tok) >= 2 and tok[0].lower() == low:
            return tok[1]
    return None


def set_keyword(lines, key, value):
    low = key.lower()
    for i, line in enumerate(lines):
        tok = line.split("#", 1)[0].split()
        if tok and tok[0].lower() == low:
            lines[i] = "%-27s %s\n" % (key, value)
            return False
    lines.append("%-27s %s\n" % (key, value))
    return True


def parse_index_list(text: str, natoms: int):
    out = []
    for part in text.replace(",", " ").split():
        if "-" in part:
            a, b = part.split("-", 1)
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    for i in out:
        if i < 1 or i > natoms:
            sys.exit("ERROR: atom index %d out of range 1..%d" % (i, natoms))
    return out


def expand_list(text: str):
    vals = []
    for tok in text.replace(",", " ").split():
        if "*" in tok:
            n, v = tok.split("*", 1)
            vals.extend([float(v)] * int(n))
        else:
            vals.append(float(tok))
    return vals


def parse_signs(text: str):
    out = []
    for tok in text.replace(",", " ").split():
        if tok in ("+", "u", "up"):
            out.append(1.0)
        elif tok in ("-", "d", "dn", "down"):
            out.append(-1.0)
        else:
            out.append(1.0 if float(tok) >= 0 else -1.0)
    return out


def parse_pairs(items, what):
    res = {}
    for it in items or []:
        if "=" not in it:
            sys.exit("ERROR: %s expects KEY=VALUE, got '%s'" % (what, it))
        k, v = it.split("=", 1)
        try:
            res[k.strip()] = float(v)
        except ValueError:
            sys.exit("ERROR: %s: '%s' is not a number" % (what, v))
    return res


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Set the initial up/down electrons (and non-collinear Euler angles) of every atom "
                    "in an OpenMX input, with an electron-count check.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python omx_magconf.py fe.dat --mode fm --moment Fe=3.0 -o fe_fm.dat\n"
               "  python omx_magconf.py nio.dat --mode afm --moment Ni=1.8 -o nio_afm.dat\n"
               "  python omx_magconf.py mno.dat --mode custom --magmom '2*4.5 2*-4.5 4*0' --inplace\n"
               "  python omx_magconf.py cr2.dat --mode afm --moment Cr=2 --nc --theta 90 --constrain\n"
               "  python omx_magconf.py nio.dat --check")
    p.add_argument("input", help="OpenMX input file (.dat)")
    out = p.add_mutually_exclusive_group()
    out.add_argument("-o", "--output", help="write the modified input here (default: stdout)")
    out.add_argument("--inplace", action="store_true", help="overwrite the input (a .bak copy is kept)")
    out.add_argument("--block-only", action="store_true", help="print only the coordinate block")
    p.add_argument("--mode", choices=["fm", "afm", "custom", "none"], default="fm",
                   help="fm: same sign; afm: alternating/explicit signs; custom: --magmom list; "
                        "none: non-magnetic (up = down) (default fm)")
    p.add_argument("--moment", action="append", metavar="KEY=MU",
                   help="initial moment in muB for a species label or element (repeatable), e.g. Fe=3.0")
    p.add_argument("--atoms", default=None,
                   help="restrict fm/afm to these atom indices, e.g. '1-4,9' (default: all atoms "
                        "whose species/element appears in --moment)")
    p.add_argument("--signs", default=None,
                   help="afm sign pattern for the selected atoms in file order, e.g. '+ - - +' or "
                        "'1 -1 -1 1' (default: alternate per element)")
    p.add_argument("--magmom", default=None, help="custom: per-atom moments, e.g. '2*3.0 2*-3.0 4*0'")
    p.add_argument("--valence", action="append", metavar="KEY=N",
                   help="override the valence charge of a species label or element (repeatable)")
    p.add_argument("--data-path", default=None,
                   help="DFT_DATA19 directory to read valence.electron from the VPS files "
                        "(default: DATA.PATH of the input, if it exists)")
    p.add_argument("--nc", action="store_true",
                   help="non-collinear: write columns 8-13 and set scf.SpinPolarization nc")
    p.add_argument("--theta", type=float, default=0.0, help="nc: spin Euler angle theta in deg (default 0)")
    p.add_argument("--phi", type=float, default=0.0, help="nc: spin Euler angle phi in deg (default 0)")
    p.add_argument("--angles", action="append", metavar="I:THETA:PHI",
                   help="nc: per-atom spin angles, e.g. 2:180:0 (repeatable)")
    p.add_argument("--constrain", action="store_true",
                   help="nc: column 12 = 1 for magnetic atoms (constraint/Zeeman), 0 otherwise")
    p.add_argument("--orbpol", action="store_true",
                   help="nc: column 13 = on for magnetic atoms (DFT+U orbital polarization)")
    p.add_argument("--no-set-spin", action="store_true",
                   help="do not change scf.SpinPolarization (default: set on/nc when moments are nonzero)")
    p.add_argument("--check", action="store_true",
                   help="only report the current occupations and the electron-count check")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if not os.path.isfile(args.input):
        sys.exit("ERROR: input file not found: %s" % args.input)
    with open(args.input, encoding="utf-8-sig", errors="replace") as fh:   # strips a UTF-8 BOM
        lines = fh.read().replace("\r\n", "\n").splitlines(True)

    # ---- species and valence
    sp_blk = find_block(lines, "Definition.of.Atomic.Species")
    if sp_blk is None:
        sys.exit("ERROR: <Definition.of.Atomic.Species block not found in %s" % args.input)
    species_vps = {}
    for line in lines[sp_blk[0] + 1: sp_blk[1]]:
        tok = line.split("#", 1)[0].split()
        if len(tok) >= 3:
            species_vps[tok[0]] = tok[2]
    data_path = args.data_path or keyword_value(lines, "DATA.PATH")
    if data_path and not os.path.isabs(data_path):
        data_path = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.input)), data_path))
    if data_path and not os.path.isdir(data_path):
        data_path = None
    val_over = parse_pairs(args.valence, "--valence")
    table = _table_valence()
    valence = {}
    for sp, vps in species_vps.items():
        el = element_of(sp, vps)
        v = val_over.get(sp, val_over.get(el))
        src = "--valence"
        if v is None:
            v = _vps_file_valence(vps, data_path)
            src = "VPS file"
        if v is None and table is not None:
            v = table(vps, None)
            src = "DFT_DATA19 table"
        if v is None:
            sys.exit("ERROR: unknown valence charge of species %s (%s); give --valence %s=N or --data-path"
                     % (sp, vps, sp))
        valence[sp] = (v, src)

    # ---- coordinates
    blk = find_block(lines, "Atoms.SpeciesAndCoordinates")
    if blk is None:
        sys.exit("ERROR: <Atoms.SpeciesAndCoordinates block not found in %s" % args.input)
    atoms = []
    for k in range(blk[0] + 1, blk[1]):
        tok = lines[k].split("#", 1)[0].split()
        if not tok:
            continue
        if len(tok) < 7:
            sys.exit("ERROR: line %d has fewer than 7 columns: %s" % (k + 1, lines[k].strip()))
        sp = tok[1]
        if sp not in species_vps:
            sys.exit("ERROR: line %d uses species '%s' that is not defined in "
                     "Definition.of.Atomic.Species" % (k + 1, sp))
        try:
            up, dn = float(tok[5]), float(tok[6])
        except ValueError:
            sys.exit("ERROR: line %d: columns 6-7 are not numbers: %s" % (k + 1, lines[k].strip()))
        atoms.append({"line": k, "tok": tok, "species": sp, "el": element_of(sp, species_vps[sp]),
                      "up": up, "down": dn})
    natoms = len(atoms)
    if natoms == 0:
        sys.exit("ERROR: no atoms in the coordinate block")
    nat_kw = keyword_value(lines, "Atoms.Number")
    if nat_kw is not None and int(float(nat_kw)) != natoms:
        print("WARNING: Atoms.Number = %s but the block has %d atoms" % (nat_kw, natoms), file=sys.stderr)

    if args.check:
        bad = 0
        print(" idx species  valence     up     down    sum   moment  status")
        for i, a in enumerate(atoms, 1):
            v = valence[a["species"]][0]
            ok = abs(a["up"] + a["down"] - v) < TOL
            bad += 0 if ok else 1
            print("%4d %-8s %7.2f %7.3f %7.3f %7.3f %7.3f  %s"
                  % (i, a["species"], v, a["up"], a["down"], a["up"] + a["down"], a["up"] - a["down"],
                     "ok" if ok else "MISMATCH"))
        tot_v = sum(valence[a["species"]][0] for a in atoms)
        tot_m = sum(a["up"] - a["down"] for a in atoms)
        print("total valence electrons %.3f, initial moment %.3f muB, mismatching atoms: %d"
              % (tot_v, tot_m, bad))
        return 1 if bad else 0

    # ---- target moments
    moments = [0.0] * natoms
    if args.mode == "custom":
        if not args.magmom:
            sys.exit("ERROR: --mode custom needs --magmom")
        moments = expand_list(args.magmom)
        if len(moments) != natoms:
            sys.exit("ERROR: --magmom has %d values, the input has %d atoms" % (len(moments), natoms))
    elif args.mode in ("fm", "afm"):
        mom = parse_pairs(args.moment, "--moment")
        if not mom:
            sys.exit("ERROR: --mode %s needs at least one --moment KEY=MU" % args.mode)
        for key in mom:
            if not any(key in (a["species"], a["el"]) for a in atoms):
                sys.exit("ERROR: --moment %s matches no species label or element in the input" % key)
        selected = set(parse_index_list(args.atoms, natoms)) if args.atoms else None
        chosen = []
        for i, a in enumerate(atoms, 1):
            m = mom.get(a["species"], mom.get(a["el"]))
            if m is None or (selected is not None and i not in selected):
                continue
            chosen.append((i, a, m))
        if not chosen:
            sys.exit("ERROR: no atom selected for the moments")
        if args.mode == "fm":
            for i, a, m in chosen:
                moments[i - 1] = m
        else:
            if args.signs:
                signs = parse_signs(args.signs)
                if len(signs) != len(chosen):
                    sys.exit("ERROR: --signs has %d entries but %d atoms are selected"
                             % (len(signs), len(chosen)))
                for (i, a, m), s in zip(chosen, signs):
                    moments[i - 1] = s * abs(m)
            else:
                count = {}
                for i, a, m in chosen:
                    n = count.get(a["el"], 0)
                    moments[i - 1] = abs(m) if n % 2 == 0 else -abs(m)
                    count[a["el"]] = n + 1
                if any(c % 2 for c in count.values()):
                    print("WARNING: an odd number of atoms of some element was selected; the AFM "
                          "configuration is not compensated", file=sys.stderr)

    # ---- nc angles
    angles = [(args.theta, args.phi)] * natoms
    for item in args.angles or []:
        try:
            i, t, f = item.split(":")
            angles[int(i) - 1] = (float(t), float(f))
        except (ValueError, IndexError):
            sys.exit("ERROR: --angles expects I:THETA:PHI with 1 <= I <= %d, got '%s'" % (natoms, item))
    if (args.angles or args.constrain or args.orbpol) and not args.nc:
        print("WARNING: --angles/--constrain/--orbpol are only written with --nc", file=sys.stderr)

    # ---- rewrite lines
    report = []
    for i, a in enumerate(atoms):
        v = valence[a["species"]][0]
        m = moments[i]
        if abs(m) > v + TOL:
            sys.exit("ERROR: atom %d (%s): |moment| %.3f exceeds the valence charge %.3f"
                     % (i + 1, a["species"], abs(m), v))
        up, dn = 0.5 * (v + m), 0.5 * (v - m)
        if abs(up + dn - v) > TOL or up < -TOL or dn < -TOL:
            sys.exit("ERROR: electron count check failed for atom %d" % (i + 1))
        tok = a["tok"]
        cols = "  %4s  %-6s %14s %14s %14s   %8.4f %8.4f" % (tok[0], tok[1], tok[2], tok[3], tok[4], up, dn)
        if args.nc:
            th, ph = angles[i]
            mag = abs(m) > TOL
            cols += "   %7.2f %7.2f   %7.2f %7.2f   %d  %s" % (
                th, ph, th, ph, 1 if (args.constrain and mag) else 0, "on" if (args.orbpol and mag) else "off")
        lines[a["line"]] = cols + "\n"
        report.append((i + 1, a["species"], v, up, dn, m))

    total_m = sum(r[5] for r in report)
    if not args.no_set_spin and not args.block_only:
        target = "nc" if args.nc else ("on" if any(abs(r[5]) > TOL for r in report) else None)
        cur = (keyword_value(lines, "scf.SpinPolarization") or "off").lower()
        if target and cur != target:
            set_keyword(lines, "scf.SpinPolarization", target)
            print("NOTE: scf.SpinPolarization %s -> %s" % (cur, target), file=sys.stderr)
        if args.nc and cur != "nc":
            print("NOTE: nc calculations cost about 4x a collinear run; SOC additionally needs "
                  "scf.SpinOrbit.Coupling on", file=sys.stderr)

    print(" idx species  valence      up    down   moment", file=sys.stderr)
    for r in report:
        print("%4d %-8s %7.2f %7.3f %7.3f %7.3f" % r, file=sys.stderr)
    print("electrons per atom conserved; total valence %.3f, total initial moment %.3f muB"
          % (sum(r[2] for r in report), total_m), file=sys.stderr)
    srcs = sorted({src for _, src in valence.values()})
    print("valence source: %s" % ", ".join(srcs), file=sys.stderr)

    if args.block_only:
        sys.stdout.write("".join(lines[blk[0]: blk[1] + 1]))
    elif args.inplace:
        shutil.copy2(args.input, args.input + ".bak")
        with open(args.input, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("".join(lines))
        print("updated %s (backup %s.bak)" % (args.input, args.input), file=sys.stderr)
    elif args.output:
        with open(args.output, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("".join(lines))
        print("wrote %s" % args.output, file=sys.stderr)
    else:
        sys.stdout.write("".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
