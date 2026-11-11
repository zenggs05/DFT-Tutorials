#!/usr/bin/env python3
"""Inspect the OpenMX VPS/PAO database (DFT_DATA19) and check species in .dat files.

DFT-Tutorials, OpenMX chapter 1 (01-overview.md).  Standard library only, Python >= 3.9.

Sub-commands
------------
  vps        list pseudopotentials: element, AtomSpecies (Z), valence.electron, xc.type
  pao        list pseudo-atomic orbital files: radial.cutoff.pao, PAO.Lmax, PAO.Mul
  check-dat  check an OpenMX input file: PAO/VPS files exist, basis specification is
             valid, number of basis functions, initial up+down electrons vs. valence
             electrons, total electron count, a few common pitfalls

Database directory: for "vps" and "pao" the first match of --data, the
environment variable OMX_DATA and ../DFT_DATA19 is used. For "check-dat" the
directory is --data if given, otherwise the DATA.PATH keyword of the file
(default ../DFT_DATA19, resolved from the current directory, as OpenMX does);
if that is invalid an error is reported and OMX_DATA is used for the rest.

Examples
--------
  python list_dft_data.py vps --data ~/software/openmx4.0/DFT_DATA19 -e H C O Fe
  python list_dft_data.py pao --data ~/software/openmx4.0/DFT_DATA19 -e Fe --csv
  python list_dft_data.py check-dat Methane.dat
"""
import argparse
import os
import re
import sys
from pathlib import Path

L_OF = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}


# --------------------------------------------------------------------------- helpers
def element_of(name):
    m = re.match(r"([A-Z][a-z]?)", name)
    return m.group(1) if m else name


def matches_element(name, elements):
    if not elements:
        return True
    return any(re.match(r"%s(?![a-z])" % re.escape(e), name) for e in elements)


def scan_keys(path, keys):
    """Return {key: first value} for 'key value' lines; stops when all are found."""
    found = {}
    pats = {k: re.compile(r"^\s*%s\s+(\S+)" % re.escape(k), re.I) for k in keys}
    try:
        with open(path, "r", encoding="ascii", errors="replace") as fh:
            for line in fh:
                for k, p in pats.items():
                    if k not in found:
                        m = p.match(line)
                        if m:
                            found[k] = m.group(1)
                if len(found) == len(keys):
                    break
    except OSError as exc:
        print(f"warning: cannot read {path}: {exc}", file=sys.stderr)
    return found


def vps_info(path):
    d = scan_keys(path, ["AtomSpecies", "valence.electron", "xc.type"])
    val = d.get("valence.electron")
    return {
        "Z": d.get("AtomSpecies", "?"),
        "valence": float(val) if val is not None else None,
        "xc": d.get("xc.type", "?"),
    }


def pao_info(path):
    d = scan_keys(path, ["radial.cutoff.pao", "PAO.Lmax", "PAO.Mul"])
    return {
        "rc": d.get("radial.cutoff.pao", "?"),
        "Lmax": int(d["PAO.Lmax"]) if "PAO.Lmax" in d else None,
        "Mul": int(d["PAO.Mul"]) if "PAO.Mul" in d else None,
    }


def resolve_one(cand):
    """Return Path(cand) if it contains PAO/ and VPS/, else None."""
    if not cand:
        return None
    p = Path(os.path.expanduser(cand))
    if (p / "PAO").is_dir() and (p / "VPS").is_dir():
        return p
    return None


def resolve_data(arg):
    """Database for the vps/pao listings: --data, then $OMX_DATA, then ../DFT_DATA19."""
    for cand in (arg, os.environ.get("OMX_DATA"), "../DFT_DATA19"):
        p = resolve_one(cand)
        if p is not None:
            return p
    return None


def parse_basis(spec):
    """'s2p2d1' or 's2>1p2>1' -> list of (l, primitive, contracted)."""
    out = []
    pos = 0
    for m in re.finditer(r"([spdfg])(\d+)(?:>(\d+))?", spec):
        if m.start() != pos:
            return None
        prim = int(m.group(2))
        cont = int(m.group(3)) if m.group(3) else prim
        out.append((L_OF[m.group(1)], prim, cont))
        pos = m.end()
    return out if out and pos == len(spec) else None


# --------------------------------------------------------------------------- vps / pao
def cmd_vps(args):
    data = resolve_data(args.data)
    if data is None:
        sys.exit("error: DFT_DATA19 not found; use --data or set OMX_DATA")
    rows = []
    for f in sorted((data / "VPS").glob("*.vps")):
        if not matches_element(f.stem, args.element):
            continue
        info = vps_info(f)
        rows.append((element_of(f.stem), f.stem, info["Z"],
                     "?" if info["valence"] is None else f"{info['valence']:g}", info["xc"]))
    header = ("element", "vps", "Z", "valence", "xc.type")
    print_table(header, rows, args.csv)


def cmd_pao(args):
    data = resolve_data(args.data)
    if data is None:
        sys.exit("error: DFT_DATA19 not found; use --data or set OMX_DATA")
    rows = []
    for f in sorted((data / "PAO").glob("*.pao")):
        if not matches_element(f.stem, args.element):
            continue
        info = pao_info(f)
        lmax = "?" if info["Lmax"] is None else "spdfg"[: info["Lmax"] + 1]
        rows.append((element_of(f.stem), f.stem, info["rc"], lmax,
                     "?" if info["Mul"] is None else str(info["Mul"])))
    header = ("element", "pao", "rc(Bohr)", "l available", "radial fns per l")
    print_table(header, rows, args.csv)


def print_table(header, rows, csv):
    if csv:
        print(",".join(header))
        for r in rows:
            print(",".join(str(x) for x in r))
        return
    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)] if rows else [len(h) for h in header]
    fmt = "  ".join("%-" + str(w) + "s" for w in widths)
    print(fmt % header)
    print(fmt % tuple("-" * w for w in widths))
    for r in rows:
        print(fmt % tuple(str(x) for x in r))
    print(f"\n{len(rows)} entries")


# --------------------------------------------------------------------------- check-dat
def read_dat(path):
    raw = Path(path).read_bytes()
    problems = []
    if raw.startswith(b"\xef\xbb\xbf"):
        problems.append("file starts with a UTF-8 BOM")
    if raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        problems.append("file is UTF-16 encoded (convert with iconv -f UTF-16 -t UTF-8)")
    if b"\r\n" in raw:
        problems.append("CRLF line endings (convert with dos2unix)")
    bad = [i + 1 for i, line in enumerate(raw.split(b"\n")) if any(b > 127 for b in line)]
    if bad:
        problems.append(f"non-ASCII bytes on line(s) {bad[:10]}{' ...' if len(bad) > 10 else ''}")
    text = raw.decode("ascii", errors="replace").replace("\r", "")
    lines = [ln.split("#", 1)[0].strip() for ln in text.split("\n")]
    return lines, problems


def keyword(lines, key):
    """First occurrence wins (as in the OpenMX input parser)."""
    for ln in lines:
        tok = ln.split()
        if len(tok) >= 2 and tok[0].lower() == key.lower():
            return tok[1:]
    return None


def block(lines, name):
    start = "<" + name.lower()
    end = name.lower() + ">"
    out, inside = [], False
    for ln in lines:
        tok = ln.split()
        if not tok:
            continue
        if not inside and tok[0].lower() == start:
            inside = True
            continue
        if inside:
            if tok[0].lower() == end:
                return out
            out.append(tok)
    return out if inside else None


def cmd_check(args):
    lines, problems = read_dat(args.datfile)
    errors, warnings = [], list(problems)

    dp = keyword(lines, "DATA.PATH")
    dp_val = dp[0] if dp else None
    if dp_val is None:
        warnings.append("DATA.PATH not set: OpenMX will use ../DFT_DATA19 relative to the run directory")
    elif not (dp_val.startswith(("/", "~")) or os.path.isabs(dp_val)):
        warnings.append(f"DATA.PATH '{dp_val}' is relative (it is resolved from the directory where openmx runs)")
    if args.data:
        data = resolve_one(args.data)
        if data is None:
            errors.append(f"--data '{args.data}' does not contain PAO/ and VPS/")
    else:
        effective = dp_val or "../DFT_DATA19"
        data = resolve_one(effective)
        if data is None:
            errors.append(f"DATA.PATH '{effective}' does not contain PAO/ and VPS/: openmx would stop "
                          "when reading the basis/pseudopotential files (edit DATA.PATH)")
            data = resolve_one(os.environ.get("OMX_DATA"))
            if data is not None:
                warnings.append(f"remaining checks use OMX_DATA={data}")

    for ln in lines:
        tok = ln.split()
        if tok and tok[0].lower() == "system.currentdirectory":
            warnings.append("'System.CurrentDirectory' (two r) is not the official keyword; "
                            "use System.CurrrentDirectory (three r)")

    xc = (keyword(lines, "scf.XcType") or ["?"])[0]
    spin = (keyword(lines, "scf.SpinPolarization") or ["off"])[0].lower()
    charge = float((keyword(lines, "scf.system.charge") or ["0"])[0])
    nspec_kw = keyword(lines, "Species.Number")
    natom_kw = keyword(lines, "Atoms.Number")

    species = block(lines, "Definition.of.Atomic.Species")
    if not species:
        errors.append("block <Definition.of.Atomic.Species not found or empty")
        species = []
    if nspec_kw is None or int(nspec_kw[0]) != len(species):
        errors.append(f"Species.Number = {nspec_kw[0] if nspec_kw else 'missing'} but the block has {len(species)} line(s)")

    print(f"file          : {args.datfile}")
    print(f"DFT_DATA19    : {data if data else 'NOT FOUND'}")
    print(f"scf.XcType    : {xc}    scf.SpinPolarization: {spin}    scf.system.charge: {charge:g}\n")

    spec_info = {}
    print("species  PAO spec            basis  valence  VPS            checks")
    for tok in species:
        if len(tok) < 3:
            errors.append(f"malformed species line: {' '.join(tok)}")
            continue
        name, paospec, vps = tok[0], tok[1], tok[2]
        notes = []
        pao_file, _, spec = paospec.rpartition("-")
        shells = parse_basis(spec) if pao_file else None
        nbasis = sum(c * (2 * l + 1) for l, p, c in shells) if shells else None
        if shells is None:
            errors.append(f"{name}: cannot parse basis specification '{paospec}' (expected e.g. C6.0-s2p2d1)")
        valence = None
        if data is not None:
            pfile = data / "PAO" / (pao_file + ".pao")
            vfile = data / "VPS" / (vps + ".vps")
            if not pfile.is_file():
                errors.append(f"{name}: PAO file not found: {pfile}")
            elif shells:
                info = pao_info(pfile)
                for l, p, c in shells:
                    if info["Lmax"] is not None and l > info["Lmax"]:
                        errors.append(f"{name}: l={'spdfg'[l]} requested but {pfile.name} has PAO.Lmax={info['Lmax']}")
                    if info["Mul"] is not None and p > info["Mul"]:
                        errors.append(f"{name}: {p} radial functions for l={'spdfg'[l]} but PAO.Mul={info['Mul']}")
                    if c > p:
                        errors.append(f"{name}: contraction {p}>{c} is invalid")
                notes.append(f"rc={info['rc']}")
            if not vfile.is_file():
                errors.append(f"{name}: VPS file not found: {vfile}")
            else:
                vi = vps_info(vfile)
                valence = vi["valence"]
                vxc = vi["xc"].upper()
                if xc.upper().startswith("GGA") and vxc != "GGA":
                    warnings.append(f"{name}: scf.XcType={xc} but {vps}.vps has xc.type={vi['xc']}")
                if xc.upper() in ("LDA", "LSDA-CA", "LSDA-PW") and vxc != "LDA":
                    warnings.append(f"{name}: scf.XcType={xc} but {vps}.vps has xc.type={vi['xc']}")
        spec_info[name.lower()] = (valence, nbasis)
        print(f"{name:8s} {paospec:18s} {str(nbasis):>5s}  {'?' if valence is None else f'{valence:g}':>7s}  {vps:14s} {' '.join(notes)}")

    atoms = block(lines, "Atoms.SpeciesAndCoordinates")
    if not atoms:
        errors.append("block <Atoms.SpeciesAndCoordinates not found or empty")
        atoms = []
    if natom_kw is None or int(natom_kw[0]) != len(atoms):
        errors.append(f"Atoms.Number = {natom_kw[0] if natom_kw else 'missing'} but the block has {len(atoms)} line(s)")

    nel, nbas, mom = 0.0, 0, 0.0
    for i, tok in enumerate(atoms, start=1):
        if len(tok) < 7:
            errors.append(f"atom line {i}: expected at least 7 columns (index species x y z up down)")
            continue
        if tok[0] != str(i):
            errors.append(f"atom line {i}: sequential number is '{tok[0]}' (must be {i})")
        sp = tok[1].lower()
        if sp not in spec_info:
            errors.append(f"atom {i}: undefined species name '{tok[1]}'")
            continue
        up, dn = float(tok[5]), float(tok[6])
        val, nb = spec_info[sp]
        if val is not None:
            nel += val
            if abs(up + dn - val) > 1e-6:
                warnings.append(f"atom {i} ({tok[1]}): up+down = {up + dn:g} but valence electrons = {val:g}")
        if nb:
            nbas += nb
        mom += up - dn

    print(f"\natoms: {len(atoms)}   basis functions: {nbas}   valence electrons: {nel:g}"
          f"   electrons after scf.system.charge: {nel - charge:g}")
    print(f"initial spin moment (sum of up-down): {mom:g} muB")
    if spin == "off" and abs(mom) > 1e-6:
        warnings.append("scf.SpinPolarization is off but the initial up/down electrons are not balanced")

    print()
    for w in warnings:
        print("WARNING:", w)
    for e in errors:
        print("ERROR  :", e)
    print("RESULT : " + ("OK" if not errors else f"{len(errors)} error(s)"))
    return 1 if errors else 0


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, helptext in (("vps", "list VPS files"), ("pao", "list PAO files")):
        p = sub.add_parser(name, help=helptext)
        p.add_argument("--data", help="DFT_DATA19 directory (contains PAO/ and VPS/)")
        p.add_argument("-e", "--element", nargs="*", default=[], help="element symbols to show, e.g. -e H C Fe")
        p.add_argument("--csv", action="store_true", help="CSV output")
    p = sub.add_parser("check-dat", help="check species/atoms of an OpenMX .dat file")
    p.add_argument("datfile")
    p.add_argument("--data", help="DFT_DATA19 directory (overrides DATA.PATH of the file)")
    args = ap.parse_args()
    if args.cmd == "vps":
        cmd_vps(args)
    elif args.cmd == "pao":
        cmd_pao(args)
    else:
        sys.exit(cmd_check(args))


if __name__ == "__main__":
    main()
