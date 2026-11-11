#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vasp_incar_writer.py -- write INCAR.static for static (DOS) runs of
                        non-magnetic / magnetic / DFT+U systems.
                        为非磁、磁性、DFT+U 三类静态（DOS）计算写 INCAR.static。

Author: Gaoshan Zeng
Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

This file creates VASP INCAR.static files for some common cases found with
transition-metal oxides.  Only electronic control parameters are specified
here; ionic-relaxation tags are handled by the job-submission script writers.

Usage examples
--------------
    python vasp_incar_writer.py                       # non-magnetic, PBEsol
    python vasp_incar_writer.py -m 2*0_2_-2_6*0       # spin-polarised, MAGMOM given
    python vasp_incar_writer.py -m 2_-2_2*0 -hu       # + DFT+U, default U table
    python vasp_incar_writer.py -m 2_-2_2*0 -hu Ni=6.2,Mn=3.9
    python vasp_incar_writer.py -xc PBE --encut 520 --ncore 8 -o INCAR
    python vasp_incar_writer.py -h

Behaviour kept from the original script: reads POSCAR for SYSTEM, writes
INCAR.static, ISTART=1 / ICHARG=1 (continue from a relaxation), EDIFF=1e-7,
NEDOS=3001, ISMEAR=-5, LORBIT=11, ENCUT=600, PBEsol (GGA = PS) as default.

Fixed in this version (2026-11):
* the script crashed when -m was omitted (None.replace) and had a syntax
  error in the -name argument;
* -hu was compared as a string with the bool False, so the DFT+U block was
  never written; -hu now accepts no value / True / False / 'El=U,...';
* MAGMOM was written commented out; it is now active and its length is
  checked against the number of atoms (N*value syntax understood);
* LDAUL/LDAUU/LDAUJ are built per species from the POSCAR element line
  (were hard-coded to three species);
* outdated parallel tags NPAR/LPLANE/LSCALU replaced by NCORE (+ optional KPAR);
* GGA = PS was hard-coded; the functional is now an option (-xc).

Dependencies: python >= 3.9 (standard library only).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

__version__ = "2.0.0"

# Dudarev U_eff (eV) used when -hu is given without explicit values:
# Materials-Project-style values for oxides; check literature for your system.
DEFAULT_U = {"Co": 3.32, "Cr": 3.7, "Fe": 5.3, "Mn": 3.9, "Mo": 4.38,
             "Ni": 6.2, "V": 3.25, "W": 6.2}

F_BLOCK = set("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
              "Ac Th Pa U Np Pu Am Cm".split())
D_BLOCK = set("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
              "Hf Ta W Re Os Ir Pt Au Hg".split())

# functional name -> INCAR lines
FUNCTIONALS = {
    "PBESOL": ["GGA = PS # PBEsol exchange-correlation functional"],
    "PBE": ["GGA = PE # PBE exchange-correlation functional"],
    "RPBE": ["GGA = RP # RPBE exchange-correlation functional"],
    "REVPBE": ["GGA = RE # revPBE exchange-correlation functional"],
    "PW91": ["GGA = 91 # PW91 exchange-correlation functional"],
    "AM05": ["GGA = AM # AM05 exchange-correlation functional"],
    "R2SCAN": ["METAGGA = R2SCAN # r2SCAN meta-GGA (VASP >= 6.2)",
               "LMIXTAU = .TRUE. # mix kinetic-energy density (helps meta-GGA SCF)"],
    "SCAN": ["METAGGA = SCAN # SCAN meta-GGA",
             "LMIXTAU = .TRUE. # mix kinetic-energy density (helps meta-GGA SCF)"],
}

TRUE_WORDS = {"true", "t", "yes", "y", "1", "on", ".true."}
FALSE_WORDS = {"false", "f", "no", "n", "0", "off", ".false.", "none", ""}


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def read_poscar_header(path: Path):
    """Return (first-line name, species list or None, counts list)."""
    if not path.is_file():
        die(f"'{path}' not found - run the script in a directory with a POSCAR or use -p")
    lines = path.read_text(errors="replace").splitlines()
    if len(lines) < 7:
        die(f"'{path}' is too short to be a POSCAR")
    name = lines[0].strip().replace(" ", "")
    tok = lines[5].split()
    species = None
    i = 5
    if tok and not all(t.isdigit() for t in tok):
        species = [re.match(r"[A-Z][a-z]?", t).group(0) if re.match(r"[A-Z][a-z]?", t) else t for t in tok]
        i = 6
    try:
        counts = [int(x) for x in lines[i].split()]
    except ValueError:
        die(f"cannot read the atom counts on line {i + 1} of '{path}'")
    if species is not None and len(species) != len(counts):
        die(f"'{path}': {len(species)} element symbols but {len(counts)} atom counts")
    return name, species, counts


def expand_magmom(text: str) -> list:
    """'2*0_2_-2_6*0' (or spaces) -> list of floats."""
    vals = []
    for tok in re.split(r"[_\s,]+", text.strip()):
        if not tok:
            continue
        try:
            if "*" in tok:
                n, v = tok.split("*", 1)
                vals += [float(v)] * int(n)
            else:
                vals.append(float(tok))
        except ValueError:
            die(f"cannot understand MAGMOM item '{tok}' (use e.g. 2*0_2_-2_6*0)")
    return vals


def parse_hubbard(value: str | None, species):
    """Return dict element -> U, or {} for no DFT+U."""
    if value is None:
        return {}
    v = value.strip()
    if v.lower() in FALSE_WORDS:
        return {}
    if species is None:
        die("DFT+U needs the element line (VASP 5 format) in the POSCAR")
    if v.lower() in TRUE_WORDS:
        u = {e: DEFAULT_U[e] for e in species if e in DEFAULT_U}
        if not u:
            die(f"none of {species} has a default U; give values, e.g. -hu Cu=4.0")
        print("note: using default U values " + ", ".join(f"{k}={x}" for k, x in u.items()))
        return u
    u = {}
    for item in re.split(r"[,_;\s]+", v):
        if not item:
            continue
        m = re.fullmatch(r"([A-Z][a-z]?)\s*[=:]\s*([-+\d.eE]+)", item)
        if not m:
            die(f"cannot parse -hu item '{item}' (expected Element=U, e.g. Ni=6.2)")
        u[m.group(1)] = float(m.group(2))
    missing = [k for k in u if k not in species]
    if missing:
        die(f"-hu elements {missing} are not in the POSCAR {species}")
    return u


def build_text(a, name, species, counts, magmom, hubu) -> str:
    istart, icharg = (0, 2) if a.scratch else (1, 1)
    lines = [
        f"System = {name}",
        "#Start parameters",
        f"ISTART = {istart} # job : 0-new 1-cont 2-samecut",
        f"ICHARG = {icharg} # charge: 1-file 2-atom 11-const",
        "",
        "# Convergence tolerance",
        f"EDIFF = {a.ediff}",
        "SYMPREC = 1e-3",
        f"NELM = {a.nelm}",
        "",
        "# Density of States Info",
        f"NEDOS = {a.nedos}",
    ]
    if a.emin is not None:
        lines.append(f"EMIN = {a.emin:g} # absolute energy (eV), not relative to E-fermi")
    if a.emax is not None:
        lines.append(f"EMAX = {a.emax:g} # absolute energy (eV), not relative to E-fermi")
    lines += [
        "LORBIT = 11",
        "",
        "# BZ Integrations",
        "ISMEAR = -5 # tetrahedron method with Bloechl corrections (needs >= 4 k-points)",
        "#SIGMA = 0.05 # use ISMEAR = 0 and this SIGMA for Gamma-only or very few k-points",
        "",
        "# Calculation Details",
        "LWAVE = .FALSE. # Do we want to write the wavefunctions",
        "LCHARG = .TRUE. # Do we want to write the charge density",
        f"PREC = {a.prec}",
        f"ENCUT = {a.encut:g}",
        "LREAL = .FALSE.",
        "LASPH = .TRUE.",
        "",
        "# Exchange-correlation",
    ] + FUNCTIONALS[a.functional] + [
        "",
        "# Relaxation parameters",
        "NSW = 0 # number of ionic steps",
        "#ISIF = 3 # calculate force and relax ions",
        "#IBRION = 1 # Newtonian method for mini_opt",
        "#POTIM = 0.3 # quasi_newton scalar",
        "#EDIFFG = -0.0005",
        "",
        "# Parallelization (NPAR/LPLANE are outdated; NCORE ~ sqrt(cores per node))",
    ]
    if a.ncore:
        lines.append(f"NCORE = {a.ncore}")
    if a.kpar:
        lines.append(f"KPAR = {a.kpar}")
    text = "\n".join(lines) + "\n"

    if magmom is not None:
        text += "\n#Magnetism\nISPIN = 2\n"
        text += f"MAGMOM = {compress(magmom)}\n"

    if hubu:
        ldaul = [(3 if e in F_BLOCK else 2) if hubu.get(e, 0.0) != 0.0 else -1 for e in species]
        lmaxmix = 6 if 3 in ldaul else 4
        text += ("\n#LSDA + U\n"
                 "LDAU = .TRUE.\n"
                 "LDAUTYPE = 2 # Dudarev: only U - J matters\n"
                 f"LDAUL = {' '.join(str(x) for x in ldaul)}\n"
                 f"LDAUU = {' '.join(f'{hubu.get(e, 0.0):g}' for e in species)}\n"
                 f"LDAUJ = {' '.join('0' for _ in species)}\n"
                 f"LMAXMIX = {lmaxmix}\n"
                 "LDAUPRINT = 1\n")
    elif species and set(species) & (D_BLOCK | F_BLOCK):
        text += f"\nLMAXMIX = {6 if set(species) & F_BLOCK else 4} # d/f elements: keep l-components when reading CHGCAR\n"
    return text


def compress(vals) -> str:
    out, i = [], 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[j + 1] == vals[i]:
            j += 1
        v = f"{vals[i]:g}"
        out.append(f"{j - i + 1}*{v}" if j > i else v)
        i = j + 1
    return " ".join(out)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="vasp_incar_writer.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Write INCAR.static for non-magnetic / magnetic / DFT+U static runs (author: Gaoshan Zeng).\n"
                     "为非磁/磁性/DFT+U 静态计算写 INCAR.static（作者：Gaoshan Zeng）。"),
        epilog=("Examples / 示例:\n"
                "  python vasp_incar_writer.py -m 2*0_2_-2_6*0\n"
                "  python vasp_incar_writer.py -m 2_-2_2*0 -hu Ni=6.2\n"
                "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License"))
    # original options (kept)
    p.add_argument("-name", dest="system_name", default=None,
                   help="SYSTEM name (default: first line of POSCAR without spaces)")
    p.add_argument("-m", dest="magmom", default=None,
                   help="MAGMOM separated by '_' (or spaces), e.g. 2*0_2_-2_6*0; "
                        "if omitted ISPIN = 1 (non-magnetic). 磁矩，下划线分隔")
    p.add_argument("-hu", dest="hubU", nargs="?", const="True", default="False",
                   help="DFT+U (Dudarev). '-hu' alone or '-hu True' uses default U values "
                        f"({', '.join(f'{k}={v}' for k, v in DEFAULT_U.items())}); "
                        "'-hu Ni=6.2,Mn=3.9' sets them explicitly; default False")
    # new options
    p.add_argument("-p", "--poscar", default="POSCAR", help="structure file (default POSCAR)")
    p.add_argument("-o", "--output", default="INCAR.static", help="output file (default INCAR.static)")
    p.add_argument("-xc", "--functional", default="PBEsol", type=str.upper,
                   choices=[k for k in FUNCTIONALS],
                   help="exchange-correlation functional (default PBEsol as in the original script)")
    p.add_argument("--encut", type=float, default=600, help="ENCUT in eV (default 600)")
    p.add_argument("--prec", default="Accurate", choices=("Accurate", "High", "Normal"),
                   help="PREC (default Accurate; the original used High)")
    p.add_argument("--ediff", default="1e-7", help="EDIFF (default 1e-7)")
    p.add_argument("--nelm", type=int, default=40, help="NELM (default 40)")
    p.add_argument("--nedos", type=int, default=3001, help="NEDOS (default 3001)")
    p.add_argument("--emin", type=float, default=-12.0, help="EMIN in eV (absolute; default -12)")
    p.add_argument("--emax", type=float, default=12.0, help="EMAX in eV (absolute; default 12)")
    p.add_argument("--no-window", action="store_true", help="do not write EMIN/EMAX (let VASP choose)")
    p.add_argument("--ncore", type=int, default=4, help="NCORE (default 4; 0 = omit)")
    p.add_argument("--kpar", type=int, help="KPAR (optional)")
    p.add_argument("--scratch", action="store_true", help="ISTART=0, ICHARG=2 instead of continuing")
    p.add_argument("-f", "--force", action="store_true", help="overwrite an existing output file")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    poscar = Path(a.poscar)
    name, species, counts = read_poscar_header(poscar)
    if a.system_name:
        name = a.system_name.strip()
    if a.no_window:
        a.emin = a.emax = None
    natoms = sum(counts)

    magmom = None
    if a.magmom is not None and a.magmom.strip():
        magmom = expand_magmom(a.magmom)
        if len(magmom) != natoms:
            die(f"MAGMOM has {len(magmom)} values but the POSCAR has {natoms} atoms")
    hubu = parse_hubbard(a.hubU, species)
    if hubu and magmom is None:
        print("note: DFT+U without -m: ISPIN stays 1 (give -m for magnetic oxides)")

    out = Path(a.output)
    if out.exists() and not a.force and out.name != "INCAR.static":
        die(f"{out} exists (use -f to overwrite)")
    text = build_text(a, name, species, counts, magmom, hubu)
    with open(out, "w", newline="\n") as fh:
        fh.write(text)
    kind = "DFT+U " if hubu else ""
    kind += "magnetic" if magmom is not None else "non-magnetic"
    print(f"wrote {out}: {kind}, {a.functional}, {natoms} atoms")
    return 0


if __name__ == "__main__":
    sys.exit(main())
