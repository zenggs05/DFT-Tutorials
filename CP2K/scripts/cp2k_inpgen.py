#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cp2k_inpgen.py - generate a complete CP2K Quickstep (GPW) input from an xyz / cif / POSCAR file.

从结构文件（xyz、cif、POSCAR 等 ASE 可读格式）生成完整的 CP2K 输入：
自动写 &CELL / PERIODIC（孤立分子自动建盒并居中，配 &POISSON PERIODIC NONE + MT）；
按元素在 BASIS_MOLOPT_UZH 与 POTENTIAL_UZH 中查找基组和赝势，并核对两者的 -qN 是否一致、
名字是否真的存在；可选 OT 或对角化 + Fermi-Dirac 展宽、Monkhorst-Pack k 点、DFT-D3(BJ)、
UKS / 多重度 / 电荷 / 初始磁矩、表面偶极校正，以及 RUN_TYPE ENERGY / ENERGY_FORCE /
GEO_OPT / CELL_OPT / MD 的 &MOTION 段；--version 决定 DOS/PDOS（2026.2 起 &DOS + &PDOS +
&CURVE）和 CELL_OPT（2026.2 起删除 TYPE）的写法。

Usage examples
--------------
  python cp2k_inpgen.py h2o.xyz                                   # isolated molecule, OT, ENERGY
  python cp2k_inpgen.py h2o.xyz --run-type GEO_OPT --vacuum 6
  python cp2k_inpgen.py Si.cif --kpoints 4 4 4 --smear 300 --run-type CELL_OPT
  python cp2k_inpgen.py POSCAR --kpoints 6 6 1 --d3bj --surface-dipole --run-type GEO_OPT
  python cp2k_inpgen.py fe.cif --uks --multiplicity 5 --magmom Fe=4 --smear 500 --kpoints 8 8 8
  python cp2k_inpgen.py water64.xyz --box 12.42 12.42 12.42 --periodic XYZ --run-type MD \\
         --temperature 330 --steps 20000 --d3bj
  python cp2k_inpgen.py TiO2.cif --pdos --version 2025.2         # legacy &PRINT/&PDOS syntax
  python cp2k_inpgen.py Cu.cif --q Cu=11 --basis TZVP --data-dir /opt/cp2k/data

What is generated
-----------------
* &GLOBAL (PROJECT, RUN_TYPE, PRINT_LEVEL MEDIUM)
* &FORCE_EVAL/&DFT: BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH, POTENTIAL_FILE_NAME POTENTIAL_UZH,
  &MGRID CUTOFF/REL_CUTOFF, &QS EPS_DEFAULT, &SCF (OT + OUTER_SCF, or DIAGONALIZATION + Broyden
  MIXING [+ SMEAR, ADDED_MOS]), &KPOINTS, &XC PBE [+ DFT-D3(BJ)], &POISSON for molecules,
  optional DOS/PDOS print, &SUBSYS (&CELL, &COORD or &TOPOLOGY, one &KIND per element)
* &MOTION for GEO_OPT / CELL_OPT / MD (CSVR thermostat for NVT)
Only PBE is generated (basis/potential names *-PBE-*).  For hybrids/ADMM, DFT+U, SCAN ... edit
the result (see the CP2K chapters of DFT-Tutorials).  Always check the result with
"cp2k.psmp -c <input>" before submitting.

Basis / pseudopotential lookup
------------------------------
Data files are searched in --data-dir, $CP2K_DATA_DIR and ../data next to this script.
For each element the valence q is the intersection of the q values offered by
<BASIS>-MOLOPT-PBE-GTH-q* in BASIS_MOLOPT_UZH and GTH-PBE-q* in POTENTIAL_UZH; when several are
possible the largest (semicore) one is taken and the alternatives are printed (--q El=N to choose).
Without data files only H C N O F Si P S Cl Br I (unambiguous valences) work unless --q is given.

Dependencies: Python >= 3.9, ase (structure reading), numpy (pulled in by ase).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

import argparse
import math
import os
import re
import sys

# valence of the common GTH-PBE pseudopotentials that exist with a single q (used only when no
# POTENTIAL_UZH file is available); everything else needs the data files or --q
FALLBACK_Q = {"H": 1, "C": 4, "N": 5, "O": 6, "F": 7, "Si": 4, "P": 5, "S": 6, "Cl": 7, "Br": 7, "I": 7}
RE_Q = re.compile(r"-q(\d+)$")
ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$")


# ----------------------------------------------------------------------------- data files
def find_data_file(name, data_dir):
    dirs = []
    if data_dir:
        dirs.append(data_dir)
    if os.environ.get("CP2K_DATA_DIR"):
        dirs.extend(os.environ["CP2K_DATA_DIR"].split(os.pathsep))
    dirs.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "data"))
    for d in dirs:
        p = os.path.join(d, name)
        if os.path.isfile(p):
            return os.path.normpath(p)
    return None


def read_names(path):
    """Return {element: set(names incl. aliases)} from a CP2K basis or potential file."""
    table = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.split("#", 1)[0].strip()
            if not s:
                continue
            tok = s.split()
            if len(tok) >= 2 and ELEMENT_RE.match(tok[0]) and not re.match(r"^[-+]?\d", tok[1]):
                table.setdefault(tok[0], set()).update(tok[1:])
    return table


def q_of(name):
    m = RE_Q.search(name)
    return int(m.group(1)) if m else None


# ----------------------------------------------------------------------------- geometry
def standard_cell(cellpar):
    """CP2K convention: a along x, b in the xy plane (lower-triangular cell matrix)."""
    a, b, c, al, be, ga = cellpar
    al, be, ga = (math.radians(x) for x in (al, be, ga))
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return [[a, 0.0, 0.0], [b * math.cos(ga), b * math.sin(ga), 0.0], [cx, cy, cz]]


def parse_version(text):
    if text.lower() in ("trunk", "master", "dev"):
        return (9999, 9)
    m = re.match(r"^(\d{4})\.(\d+)$", text.strip())
    if not m:
        raise SystemExit(f"ERROR: --version must look like 2026.2, 2025.1 or 'trunk' (got {text})")
    return (int(m.group(1)), int(m.group(2)))


def parse_pairs(items, cast, what):
    out = {}
    for it in items or []:
        if "=" not in it:
            raise SystemExit(f"ERROR: {what} entries must look like El=value (got {it})")
        k, v = it.split("=", 1)
        try:
            out[k.strip()] = cast(v.strip())
        except ValueError:
            raise SystemExit(f"ERROR: bad value in {what}: {it}")
    return out


# ----------------------------------------------------------------------------- kinds
def choose_kinds(elements, args):
    qset = parse_pairs(args.q, int, "--q")
    user_basis = parse_pairs(args.kind_basis, str, "--kind-basis")
    user_pot = parse_pairs(args.kind_potential, str, "--kind-potential")
    bfile = find_data_file(args.basis_file, args.data_dir)
    pfile = find_data_file(args.potential_file, args.data_dir)
    btab = read_names(bfile) if bfile else None
    ptab = read_names(pfile) if pfile else None
    notes = []
    if btab is None or ptab is None:
        notes.append(f"WARNING: {args.basis_file}/{args.potential_file} not found (use --data-dir or set "
                     "CP2K_DATA_DIR); names were NOT checked against the data files")
    else:
        notes.append(f"basis file: {bfile}")
        notes.append(f"potential file: {pfile}")
    kinds = {}
    for el in elements:
        if el in user_basis and el in user_pot:
            bas, pot = user_basis[el], user_pot[el]
        else:
            if el in qset:
                q = qset[el]
            elif btab is not None and ptab is not None:
                bq = {q_of(n) for n in btab.get(el, ()) if n.startswith(f"{args.basis}-MOLOPT-PBE-GTH-q")}
                pq = {q_of(n) for n in ptab.get(el, ()) if n.startswith("GTH-PBE-q")}
                common = sorted(x for x in bq & pq if x)
                if not common:
                    raise SystemExit(
                        f"ERROR: no common q for {el}: {args.basis}-MOLOPT-PBE-GTH-q{sorted(bq) or '?'} in "
                        f"{args.basis_file}, GTH-PBE-q{sorted(pq) or '?'} in {args.potential_file}. "
                        "Try another --basis (DZVP/TZVP/TZV2P), --q, or --kind-basis/--kind-potential.")
                q = common[-1]
                if len(common) > 1:
                    notes.append(f"{el}: q options {common}, using q{q} (choose with --q {el}=N)")
            elif el in FALLBACK_Q:
                q = FALLBACK_Q[el]
            else:
                raise SystemExit(f"ERROR: valence of {el} unknown without the data files; give --q {el}=N "
                                 "(see POTENTIAL_UZH) or --data-dir")
            bas = user_basis.get(el, f"{args.basis}-MOLOPT-PBE-GTH-q{q}")
            pot = user_pot.get(el, f"GTH-PBE-q{q}")
        qb, qp = q_of(bas), q_of(pot)
        if qb is not None and qp is not None and qb != qp:
            msg = f"{el}: basis {bas} (q{qb}) and potential {pot} (q{qp}) have different -qN"
            if not args.no_check:
                raise SystemExit("ERROR: " + msg + " (use --no-check to write it anyway)")
            notes.append("WARNING: " + msg)
        if btab is not None and bas not in btab.get(el, ()):
            msg = f"{el}: basis {bas} not found in {args.basis_file}"
            if not args.no_check:
                raise SystemExit("ERROR: " + msg)
            notes.append("WARNING: " + msg)
        if ptab is not None and pot not in ptab.get(el, ()):
            msg = f"{el}: potential {pot} not found in {args.potential_file}"
            if not args.no_check:
                raise SystemExit("ERROR: " + msg)
            notes.append("WARNING: " + msg)
        kinds[el] = (bas, pot, qp if qp is not None else qb)
    return kinds, notes


# ----------------------------------------------------------------------------- writer
def build_input(atoms, args, kinds, version, notes):
    new_syntax = version >= (2026, 2)
    symbols = atoms.get_chemical_symbols()
    periodic = args.periodic
    L = []
    w = L.append

    # electrons / spin
    nval = [kinds[s][2] for s in symbols]
    nelec = sum(nval) - args.charge if all(v is not None for v in nval) else None
    uks = args.uks or (args.multiplicity is not None and args.multiplicity != 1)
    mult = args.multiplicity
    if nelec is not None:
        if nelec % 2 == 1 and not uks:
            raise SystemExit(f"ERROR: {nelec} valence electrons (odd) - add --uks (and --multiplicity), "
                             "otherwise CP2K stops with 'Use the LSD option for an odd number of electrons'. "
                             "If the system should be closed-shell, check for missing/extra H atoms.")
        if uks and mult is None:
            mult = 2 if nelec % 2 else 1
        if mult is not None and (nelec + mult - 1) % 2 != 0:
            raise SystemExit(f"ERROR: multiplicity {mult} is impossible with {nelec} electrons")
    elif uks and mult is None:
        mult = 1
    smear = args.smear is not None
    kp = args.kpoints is not None and args.kpoints != [1, 1, 1]
    scf = args.scf or ("diag" if (smear or kp) else "ot")
    if scf == "ot" and (smear or kp):
        if version < (9999, 0):
            raise SystemExit("ERROR: OT cannot be combined with k-points or smearing in CP2K <= 2026.2; "
                             "use --scf diag (default when --kpoints/--smear are given)")
    added = None
    if scf == "diag":
        if args.added_mos is not None:
            added = args.added_mos
        elif nelec is not None:
            added = max(10, math.ceil(0.2 * math.ceil(nelec / 2)))
        else:
            added = 20

    w(f"# CP2K input generated by cp2k_inpgen.py for CP2K {args.version}")
    w(f"# structure: {os.path.basename(args.structure)}, {len(symbols)} atoms, "
      f"{nelec if nelec is not None else '?'} valence electrons")
    w("# Check before running:  cp2k.psmp -c " + os.path.basename(args.output or (args.project + ".inp")))
    w("")
    w("&GLOBAL")
    w(f"  PROJECT     {args.project}")
    w(f"  RUN_TYPE    {args.run_type}")
    w("  PRINT_LEVEL MEDIUM")
    w("&END GLOBAL")
    w("")
    w("&FORCE_EVAL")
    w("  METHOD QUICKSTEP")
    if args.run_type == "CELL_OPT" or args.stress:
        w("  STRESS_TENSOR ANALYTICAL")
    w("  &DFT")
    w(f"    BASIS_SET_FILE_NAME {args.basis_file}")
    w(f"    POTENTIAL_FILE_NAME {args.potential_file}")
    if args.charge:
        w(f"    CHARGE {args.charge}")
    if uks:
        w("    UKS")
        w(f"    MULTIPLICITY {mult}")
    if args.surface_dipole:
        w("    SURFACE_DIPOLE_CORRECTION T   ! asymmetric slab: compensate the dipole")
        w("    SURF_DIP_DIR              Z   ! surface normal must be along z")
    w("    &MGRID")
    w(f"      CUTOFF     {args.cutoff:g}            ! Ry, converge it (conv_cutoff.py)")
    w(f"      REL_CUTOFF {args.rel_cutoff:g}")
    w("    &END MGRID")
    w("    &QS")
    w("      EPS_DEFAULT 1.0E-12")
    if args.run_type == "MD":
        w("      EXTRAPOLATION       ASPC")
        w("      EXTRAPOLATION_ORDER 3")
    w("    &END QS")
    if periodic == "NONE":
        w("    &POISSON")
        w("      PERIODIC       NONE          ! must match &CELL PERIODIC")
        w("      POISSON_SOLVER MT            ! box >= 2x the extent of the electron density")
        w("    &END POISSON")
    if kp:
        n = args.kpoints
        w("    &KPOINTS")
        w(f"      SCHEME MONKHORST-PACK {n[0]} {n[1]} {n[2]}")
        w("    &END KPOINTS")
    w("    &SCF")
    w("      SCF_GUESS ATOMIC")
    w(f"      EPS_SCF   {args.eps_scf:.1E}")
    if scf == "ot":
        w("      MAX_SCF   30")
        w("      &OT ON")
        w("        MINIMIZER      DIIS")
        w("        PRECONDITIONER FULL_SINGLE_INVERSE")
        w("      &END OT")
        w("      &OUTER_SCF ON")
        w(f"        EPS_SCF {args.eps_scf:.1E}")
        w("        MAX_SCF 10")
        w("      &END OUTER_SCF")
    else:
        w("      MAX_SCF   200")
        w(f"      ADDED_MOS {added}             ! empty states (needed for smearing / DOS)")
        w("      &DIAGONALIZATION ON")
        w("        ALGORITHM STANDARD")
        w("      &END DIAGONALIZATION")
        w("      &MIXING ON")
        w("        METHOD  BROYDEN_MIXING")
        w("        ALPHA   0.3")
        w("        NBUFFER 8")
        w("      &END MIXING")
        if smear:
            w("      &SMEAR ON")
            w("        METHOD FERMI_DIRAC")
            w(f"        ELECTRONIC_TEMPERATURE [K] {args.smear:g}")
            w("      &END SMEAR")
    w("      &PRINT")
    w("        &RESTART")
    w("          BACKUP_COPIES 1")
    w("        &END RESTART")
    w("      &END PRINT")
    w("    &END SCF")
    w("    &XC")
    w("      &XC_FUNCTIONAL PBE")
    w("      &END XC_FUNCTIONAL")
    if args.d3bj:
        w("      &VDW_POTENTIAL")
        w("        DISPERSION_FUNCTIONAL PAIR_POTENTIAL")
        w("        &PAIR_POTENTIAL")
        w("          TYPE DFTD3(BJ)              ! always write TYPE explicitly")
        w("          PARAMETER_FILE_NAME dftd3.dat")
        w("          REFERENCE_FUNCTIONAL PBE")
        w("        &END PAIR_POTENTIAL")
        w("      &END VDW_POTENTIAL")
    w("    &END XC")
    if args.pdos:
        if new_syntax:
            w("    &PRINT")
            w("      &DOS                          ! CP2K >= 2026.2 syntax")
            if scf == "ot":
                w("        NLUMO   20                  ! virtual states computed after the OT SCF")
            w("        DELTA_E 0.0005              ! grid spacing in hartree")
            w("        &PDOS")
            w("          COMPONENTS F")
            w("        &END PDOS")
            w("        &CURVE")
            w("          ENERGY_UNIT EV")
            w("          ENERGY_ZERO AUTO            ! FERMI with smearing, otherwise HOCO")
            w("          &BROADEN")
            w("            TYPE  GAUSSIAN")
            w("            WIDTH [eV] 0.1            ! FWHM")
            w("          &END BROADEN")
            w("        &END CURVE")
            w("      &END DOS")
            w("    &END PRINT")
        elif kp:
            notes.append("WARNING: legacy &PRINT/&PDOS (< 2026.2) is not available with k-points; "
                         "PDOS block omitted (use a supercell at Gamma or CP2K 2026.2)")
        else:
            w("    &PRINT")
            w("      &PDOS                         ! legacy syntax, CP2K <= 2026.1")
            w("        NLUMO      20")
            w("        COMPONENTS F")
            w("      &END PDOS")
            w("    &END PRINT")
    w("  &END DFT")

    # ---------------- SUBSYS
    w("  &SUBSYS")
    w("    &CELL")
    cp = args._cellpar
    w(f"      ABC              {cp[0]:.8f} {cp[1]:.8f} {cp[2]:.8f}")
    w(f"      ALPHA_BETA_GAMMA {cp[3]:.6f} {cp[4]:.6f} {cp[5]:.6f}")
    w(f"      PERIODIC {periodic}")
    if args.run_type == "CELL_OPT":
        w("      ! &CELL_REF with a slightly larger cell keeps the number of grid points fixed")
    w("    &END CELL")
    pos = args._positions
    if args.coord_file:
        w("    &TOPOLOGY")
        w(f"      COORD_FILE_NAME   {args.coord_file}")
        w("      COORD_FILE_FORMAT XYZ")
        w("    &END TOPOLOGY")
    else:
        w("    &COORD")
        for s, p in zip(symbols, pos):
            w(f"      {s:<3s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}")
        w("    &END COORD")
    magmom = parse_pairs(args.magmom, float, "--magmom")
    for el in sorted(set(symbols), key=symbols.index):
        bas, pot, _ = kinds[el]
        w(f"    &KIND {el}")
        w(f"      BASIS_SET {bas}")
        w(f"      POTENTIAL {pot}")
        if el in magmom:
            if not uks:
                notes.append(f"WARNING: --magmom {el} ignored without --uks")
            else:
                w(f"      MAGNETIZATION {magmom[el]:g}")
        w("    &END KIND")
    w("  &END SUBSYS")
    if args.run_type == "CELL_OPT" or args.stress:
        w("  &PRINT")
        w("    &STRESS_TENSOR ON")
        w("    &END STRESS_TENSOR")
        w("  &END PRINT")
    w("&END FORCE_EVAL")

    # ---------------- MOTION
    if args.run_type in ("GEO_OPT", "CELL_OPT", "MD"):
        w("")
        w("&MOTION")
        if args.run_type == "GEO_OPT":
            w("  &GEO_OPT")
            w(f"    OPTIMIZER {args.optimizer}")
            w(f"    MAX_ITER  {args.max_iter}")
            w(f"    MAX_FORCE {args.max_force:.1E}         ! hartree/bohr")
            w("  &END GEO_OPT")
        elif args.run_type == "CELL_OPT":
            w("  &CELL_OPT")
            if not new_syntax:
                w("    TYPE DIRECT_CELL_OPT                ! removed in CP2K 2026.2 (always direct)")
            w(f"    OPTIMIZER {args.optimizer}")
            w(f"    MAX_ITER  {args.max_iter}")
            w(f"    MAX_FORCE {args.max_force:.1E}")
            w(f"    EXTERNAL_PRESSURE [bar] {args.pressure:g}")
            w("    PRESSURE_TOLERANCE [bar] 100")
            w("  &END CELL_OPT")
        else:
            w("  &MD")
            w(f"    ENSEMBLE    {args.ensemble}")
            w(f"    STEPS       {args.steps}")
            w(f"    TIMESTEP    {args.timestep:g}             ! fs")
            w(f"    TEMPERATURE {args.temperature:g}")
            if args.ensemble == "NVT":
                w("    &THERMOSTAT")
                w("      TYPE   CSVR")
                w("      REGION GLOBAL")
                w("      &CSVR")
                w(f"        TIMECON {args.timecon:g}              ! fs")
                w("      &END CSVR")
                w("    &END THERMOSTAT")
            w("  &END MD")
        w("  &PRINT")
        w("    &TRAJECTORY")
        w("      FORMAT XYZ")
        w("    &END TRAJECTORY")
        if args.run_type == "CELL_OPT":
            w("    &CELL ON")
            w("    &END CELL")
        if args.run_type == "MD":
            w("    &VELOCITIES ON")
            w("    &END VELOCITIES")
        w("    &RESTART")
        w("      BACKUP_COPIES 1")
        w("    &END RESTART")
        w("  &END PRINT")
        w("&END MOTION")
    info = dict(nelec=nelec, uks=uks, mult=mult, scf=scf, added=added)
    return "\n".join(L) + "\n", info


# ----------------------------------------------------------------------------- structure
def prepare_structure(args):
    try:
        from ase.io import read
    except ImportError:
        raise SystemExit("ERROR: cp2k_inpgen.py needs ASE to read structures: pip install ase")
    if not os.path.isfile(args.structure):
        raise SystemExit(f"ERROR: structure file not found: {args.structure}")
    try:
        atoms = read(args.structure, format=args.format) if args.format else read(args.structure)
    except Exception as exc:  # ase raises many different exception types
        raise SystemExit(f"ERROR: ASE could not read {args.structure}: {exc}")
    if len(atoms) == 0:
        raise SystemExit(f"ERROR: {args.structure} contains no atoms")
    has_cell = atoms.cell.rank == 3
    if args.periodic is None:
        args.periodic = "XYZ" if (has_cell and any(atoms.pbc)) or args.box else "NONE"
    import numpy as np
    pos = atoms.get_positions()
    if args.box:
        a, b, c = args.box
        cellpar = [a, b, c, 90.0, 90.0, 90.0]
        if args.periodic == "NONE" or not has_cell:
            pos = pos - pos.mean(axis=0) + np.array([a, b, c]) / 2.0
        else:
            raise SystemExit("ERROR: --box is meant for structures without a cell (xyz); the file has a cell")
    elif args.periodic == "NONE":
        ext = pos.max(axis=0) - pos.min(axis=0)
        box = [math.ceil(2.0 * max(max(2.0 * e, e + 2.0 * args.vacuum), 1.0)) / 2.0 for e in ext]
        if args.cubic_box:
            box = [max(box)] * 3
        cellpar = [box[0], box[1], box[2], 90.0, 90.0, 90.0]
        center = (pos.max(axis=0) + pos.min(axis=0)) / 2.0
        pos = pos - center + np.array(box) / 2.0
    else:
        if not has_cell:
            raise SystemExit(f"ERROR: {args.structure} has no cell; give --box a b c or --periodic NONE")
        cellpar = list(atoms.cell.cellpar())
        frac = atoms.get_scaled_positions(wrap=False)
        pos = np.array(frac) @ np.array(standard_cell(cellpar))
    args._cellpar = cellpar
    args._positions = pos.tolist()
    return atoms


def build_parser():
    p = argparse.ArgumentParser(
        description="Generate a CP2K (GPW, PBE) input from a structure file, with basis/pseudopotential "
                    "lookup in BASIS_MOLOPT_UZH / POTENTIAL_UZH and -qN consistency check.",
        epilog="Examples:\n"
               "  python cp2k_inpgen.py h2o.xyz --run-type GEO_OPT\n"
               "  python cp2k_inpgen.py Si.cif --kpoints 4 4 4 --run-type CELL_OPT\n"
               "  python cp2k_inpgen.py slab.vasp --kpoints 4 4 1 --smear 300 --d3bj --surface-dipole\n"
               "  python cp2k_inpgen.py box.xyz --box 12 12 12 --periodic XYZ --run-type MD --d3bj\n"
               "Always run 'cp2k.psmp -c <input>' on the result before submitting.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("structure", help="xyz / cif / POSCAR / extxyz ... (anything ASE can read)")
    p.add_argument("--format", help="ASE format name if it cannot be guessed (e.g. vasp)")
    p.add_argument("-o", "--output", help="input file to write (default <project>.inp)")
    p.add_argument("--project", help="PROJECT name (default: structure file stem)")
    p.add_argument("--force", action="store_true", help="overwrite existing files")
    p.add_argument("--version", default="2026.2",
                   help="target CP2K version, e.g. 2026.2 (default), 2025.2, 2024.1 or trunk; "
                        "decides DOS/PDOS and CELL_OPT syntax")
    g = p.add_argument_group("system")
    g.add_argument("--run-type", default="ENERGY", choices=["ENERGY", "ENERGY_FORCE", "GEO_OPT", "CELL_OPT", "MD"],
                   help="RUN_TYPE (default ENERGY)")
    g.add_argument("--periodic", choices=["XYZ", "NONE"],
                   help="XYZ for crystals/slabs/liquids, NONE for isolated molecules (default: from the file)")
    g.add_argument("--box", type=float, nargs=3, metavar=("A", "B", "C"),
                   help="orthorhombic box in Angstrom for structures without a cell (molecule is centred)")
    g.add_argument("--vacuum", type=float, default=5.0,
                   help="isolated molecule: box = max(2*extent, extent + 2*vacuum) per axis (default 5 A)")
    g.add_argument("--cubic-box", action="store_true", help="isolated molecule: use a cubic box")
    g.add_argument("--charge", type=int, default=0, help="total charge (default 0)")
    g.add_argument("--uks", action="store_true", help="spin-polarized (UKS) calculation")
    g.add_argument("--multiplicity", type=int, help="2S+1 (implies --uks if != 1)")
    g.add_argument("--magmom", nargs="+", metavar="El=M",
                   help="initial MAGNETIZATION per element kind with --uks, e.g. Fe=4")
    g.add_argument("--surface-dipole", action="store_true",
                   help="SURFACE_DIPOLE_CORRECTION T, SURF_DIP_DIR Z (slab with vacuum along z)")
    g.add_argument("--coord-file", metavar="XYZ",
                   help="write coordinates to this xyz file and use &TOPOLOGY instead of an inline &COORD")
    b = p.add_argument_group("basis / pseudopotential")
    b.add_argument("--basis", default="DZVP", choices=["DZVP", "TZVP", "TZV2P"],
                   help="MOLOPT-UZH family (default DZVP -> DZVP-MOLOPT-PBE-GTH-qN)")
    b.add_argument("--q", nargs="+", metavar="El=N", help="choose the valence q per element, e.g. Cu=11 Na=9")
    b.add_argument("--kind-basis", nargs="+", metavar="El=NAME", help="explicit basis name per element")
    b.add_argument("--kind-potential", nargs="+", metavar="El=NAME", help="explicit potential name per element")
    b.add_argument("--basis-file", default="BASIS_MOLOPT_UZH", help="BASIS_SET_FILE_NAME (default BASIS_MOLOPT_UZH)")
    b.add_argument("--potential-file", default="POTENTIAL_UZH", help="POTENTIAL_FILE_NAME (default POTENTIAL_UZH)")
    b.add_argument("--data-dir", help="CP2K data directory (default: $CP2K_DATA_DIR, then ../data)")
    b.add_argument("--no-check", action="store_true", help="write names even if the checks fail (warn only)")
    s = p.add_argument_group("SCF / numerics")
    s.add_argument("--cutoff", type=float, default=400, help="CUTOFF in Ry (default 400)")
    s.add_argument("--rel-cutoff", type=float, default=60, help="REL_CUTOFF in Ry (default 60)")
    s.add_argument("--eps-scf", type=float, default=1e-6, help="EPS_SCF (default 1e-6)")
    s.add_argument("--scf", choices=["ot", "diag"],
                   help="SCF solver (default: ot; diag automatically with --kpoints or --smear)")
    s.add_argument("--smear", type=float, metavar="T_K",
                   help="Fermi-Dirac smearing with this electronic temperature in K (implies diag)")
    s.add_argument("--added-mos", type=int, help="ADDED_MOS for diag (default max(10, 20%% of occupied states))")
    s.add_argument("--kpoints", type=int, nargs=3, metavar=("N1", "N2", "N3"),
                   help="Monkhorst-Pack mesh (implies diag in CP2K <= 2026.2)")
    s.add_argument("--d3bj", action="store_true", help="add DFT-D3(BJ) dispersion (PBE reference)")
    s.add_argument("--stress", action="store_true", help="analytical stress + STRESS_TENSOR print (auto for CELL_OPT)")
    s.add_argument("--pdos", action="store_true", help="add a DOS/PDOS print block (syntax per --version)")
    m = p.add_argument_group("GEO_OPT / CELL_OPT / MD")
    m.add_argument("--optimizer", default="BFGS", choices=["BFGS", "LBFGS", "CG"], help="optimizer (default BFGS)")
    m.add_argument("--max-iter", type=int, default=200, help="MAX_ITER (default 200)")
    m.add_argument("--max-force", type=float, default=4.5e-4, help="MAX_FORCE in hartree/bohr (default 4.5e-4)")
    m.add_argument("--pressure", type=float, default=1.01325, help="CELL_OPT EXTERNAL_PRESSURE in bar (default 1.01325)")
    m.add_argument("--ensemble", default="NVT", choices=["NVT", "NVE"], help="MD ensemble (default NVT, CSVR)")
    m.add_argument("--temperature", type=float, default=300.0, help="MD temperature in K (default 300)")
    m.add_argument("--timestep", type=float, default=0.5, help="MD time step in fs (default 0.5)")
    m.add_argument("--steps", type=int, default=1000, help="MD steps (default 1000)")
    m.add_argument("--timecon", type=float, default=100.0, help="CSVR time constant in fs (default 100)")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    version = parse_version(args.version)
    if args.project is None:
        stem = os.path.splitext(os.path.basename(args.structure))[0] or "cp2k"
        args.project = re.sub(r"[^A-Za-z0-9_.-]", "_", stem)
    out = args.output or args.project + ".inp"
    for f in [out] + ([args.coord_file] if args.coord_file else []):
        if os.path.exists(f) and not args.force:
            raise SystemExit(f"ERROR: {f} exists (use -o another name or --force)")
    if args.kpoints and min(args.kpoints) < 1:
        raise SystemExit("ERROR: --kpoints values must be >= 1")
    atoms = prepare_structure(args)
    if args.periodic == "NONE" and args.kpoints:
        raise SystemExit("ERROR: k-points make no sense for PERIODIC NONE")
    if args.periodic == "NONE" and args.run_type == "CELL_OPT":
        raise SystemExit("ERROR: CELL_OPT needs a periodic cell")
    elements = list(dict.fromkeys(atoms.get_chemical_symbols()))
    kinds, notes = choose_kinds(elements, args)
    text, info = build_input(atoms, args, kinds, version, notes)
    with open(out, "w", encoding="ascii", errors="replace", newline="\n") as fh:
        fh.write(text)
    if args.coord_file:
        with open(args.coord_file, "w", encoding="ascii", newline="\n") as fh:
            fh.write(f"{len(atoms)}\ncoordinates for {args.project} (cp2k_inpgen.py)\n")
            for s, p in zip(atoms.get_chemical_symbols(), args._positions):
                fh.write(f"{s:<3s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n")
        print(f"written: {args.coord_file}")
    print(f"written: {out}")
    cp = args._cellpar
    print(f"  atoms {len(atoms)}, cell {cp[0]:.4f} {cp[1]:.4f} {cp[2]:.4f} / {cp[3]:.2f} {cp[4]:.2f} {cp[5]:.2f}, "
          f"PERIODIC {args.periodic}")
    for el, (bas, pot, q) in kinds.items():
        print(f"  {el:<3s} {bas:<28s} {pot}")
    print(f"  valence electrons {info['nelec'] if info['nelec'] is not None else '?'}"
          + (f", UKS multiplicity {info['mult']}" if info["uks"] else "")
          + f", SCF {info['scf'].upper()}" + (f", ADDED_MOS {info['added']}" if info["added"] else ""))
    for n in notes:
        print("  " + n)
    print(f"next: cp2k.psmp -c {out}   (syntax check), then converge CUTOFF and k-points")
    return 0


if __name__ == "__main__":
    sys.exit(main())
