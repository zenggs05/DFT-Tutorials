#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_any.py -- convert crystal / molecular structures between CIF, POSCAR,
                  xyz/extxyz, Gaussian .gjf, CP2K &CELL/&COORD and OpenMX .dat,
                  optionally reducing to the primitive / standard conventional
                  cell, building a supercell, adding vacuum and fixing layers.
                  万能结构转换：CIF/POSCAR/xyz/gjf/CP2K/OpenMX 互转，可选原胞或标准胞、
                  超胞、加真空、按层固定原子（Selective dynamics / FIXED_ATOMS /
                  MD.Fixed.XYZ / Gaussian 冻结码）。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Input formats
-------------
* everything ASE can read (CIF, POSCAR/CONTCAR, vasprun.xml, OUTCAR, xyz,
  extxyz, .gjf/.com, Gaussian .log, PDB, ...); --in-format forces the format;
* CP2K input/restart files (*.inp, *.restart, *.inc): &CELL (A/B/C or
  ABC + ALPHA_BETA_GAMMA, optional [angstrom]/[bohr] units, PERIODIC) and
  &COORD (Cartesian or SCALED, UNIT); kind names such as Fe1/O_w are mapped
  to elements through &KIND/ELEMENT or their leading letters;
* OpenMX input files (*.dat, *.dat#): Atoms.UnitVectors(.Unit),
  Atoms.SpeciesAndCoordinates(.Unit Ang|AU|FRAC) and the species table
  (element taken from the VPS name, e.g. Fe_PBE19 -> Fe).
Selective-dynamics flags read by ASE (FixAtoms/FixScaled/FixCartesian) are
kept through all operations except --primitive/--conventional and general
(9-integer) supercells.

Output formats (from the file name or --to)
-------------------------------------------
  vasp    POSCAR, CONTCAR, *.vasp, *POSCAR*   VASP 5 format, Direct coordinates,
                                              "Selective dynamics" if atoms are fixed
  cif     *.cif                               via ASE
  xyz     *.xyz                               plain xyz (no cell!)
  extxyz  *.extxyz                            extended xyz (cell + pbc)
  gjf     *.gjf, *.com                        Gaussian input, Tv lines for periodic
                                              cells, freeze code -1 for fixed atoms
  cp2k    *.inc, *.inp                        &CELL + &COORD for @INCLUDE in &SUBSYS
                                              (+ <name>_fixed.inc with &FIXED_ATOMS)
  openmx  *.dat                               Atoms.* blocks (+ <MD.Fixed.XYZ)
  any other ASE format with --to ase:<format>

Operation order: read -> --primitive/--conventional -> --supercell ->
--vacuum/--center -> --wrap -> --sort -> fixing -> write.

Usage examples
--------------
    python convert_any.py Si.cif -o POSCAR --primitive
    python convert_any.py Si.cif -t vasp cp2k openmx gjf         # Si.vasp, Si.inc, Si.dat, Si.gjf
    python convert_any.py POSCAR -o slab.vasp --supercell 3 3 1 --vacuum 15 --fix-layers 2
    python convert_any.py slab.vasp -o slab.inc --fix-below 8.0  # CP2K + slab_fixed.inc
    python convert_any.py h2o.xyz -o h2o.gjf --route "B3LYP/6-31G(d) Opt Freq" --nproc 8 --mem 16GB
    python convert_any.py h2o.xyz -o POSCAR --box 12              # molecule in a 12 A padded box
    python convert_any.py cp2k_opt-1.restart -o CONTCAR.vasp --sort
    python convert_any.py Si.dat# -o relaxed.cif                  # OpenMX restart input -> CIF
    python convert_any.py OUTCAR --index -1 -o final.xyz

Dependencies: Python >= 3.9, numpy, ase; spglib for --primitive/--conventional.
"""

import argparse
import os
import re
import sys

BOHR = 0.529177210903          # Angstrom (CODATA 2018)

FORMATS = ("vasp", "cif", "xyz", "extxyz", "gjf", "cp2k", "openmx")
EXT = {"vasp": ".vasp", "cif": ".cif", "xyz": ".xyz", "extxyz": ".extxyz",
       "gjf": ".gjf", "cp2k": ".inc", "openmx": ".dat"}

# Valence electrons of a few DFT_DATA19 *_PBE19 pseudopotentials (OpenMX
# initial spin densities). Other elements: --valence El=N (see the VPS table
# of the OpenMX DFT_DATA19 web page; many metals include semicore states).
OPENMX_VALENCE = {"H": 1.0, "B": 3.0, "C": 4.0, "N": 5.0, "O": 6.0, "F": 7.0,
                  "Al": 3.0, "Si": 4.0, "P": 5.0, "S": 6.0, "Cl": 7.0}
OPENMX_PAO_HINT = {"H": "H6.0-s2p1", "C": "C6.0-s2p2d1", "N": "N6.0-s2p2d1",
                   "O": "O6.0-s2p2d1", "Si": "Si7.0-s2p2d1"}


def die(msg, code=1):
    print("convert_any.py: error: " + msg, file=sys.stderr)
    sys.exit(code)


def need(module, hint):
    try:
        return __import__(module)
    except ImportError:
        die("python package '%s' is required: %s" % (module, hint))


# ==========================================================================
# readers
# ==========================================================================
def element_from_label(label):
    m = re.match(r"([A-Za-z]{1,2})", label)
    if not m:
        return None
    from ase.data import chemical_symbols
    s = m.group(1)
    cand = [s[0].upper() + s[1:].lower(), s[0].upper()]
    for c in cand:
        if c in chemical_symbols[1:]:
            return c
    return None


def strip_comment(line, chars="#!"):
    for c in chars:
        k = line.find(c)
        if k >= 0:
            line = line[:k]
    return line.strip()


def unit_factor(token, default=1.0):
    """'[angstrom]' / '[bohr]' / '[nm]' / '[pm]' -> factor to Angstrom."""
    t = token.strip("[]").lower()
    table = {"angstrom": 1.0, "ang": 1.0, "a": 1.0, "bohr": BOHR, "au": BOHR,
             "a.u.": BOHR, "nm": 10.0, "pm": 0.01}
    if not t:
        return default
    if t not in table:
        die("unsupported CP2K length unit [%s]" % t)
    return table[t]


def read_cp2k(path):
    """Minimal parser for &SUBSYS/&CELL and &COORD of a CP2K input/restart."""
    import numpy as np
    from ase import Atoms
    with open(path, encoding="utf-8", errors="replace") as fh:
        raw = [strip_comment(ln) for ln in fh]
    stack, vec, abc, angles, periodic = [], {}, None, (90.0, 90.0, 90.0), "XYZ"
    coords, scaled, cunit, kind_elem, cur_kind = [], False, 1.0, {}, None
    for ln in raw:
        if not ln:
            continue
        w = ln.split()
        key = w[0].upper()
        if key.startswith("&END"):
            if stack:
                stack.pop()
            continue
        if key.startswith("&"):
            sec = key[1:]
            stack.append(sec)
            if sec == "KIND" and len(w) > 1:
                cur_kind = w[1]
            continue
        if not stack:
            continue
        sec = stack[-1]
        if sec == "CELL" and "SUBSYS" in stack:
            fac, vals = 1.0, w[1:]
            if vals and vals[0].startswith("["):
                fac, vals = (unit_factor(vals[0]) if key != "ALPHA_BETA_GAMMA" else 1.0), vals[1:]
            if key in ("A", "B", "C"):
                vec[key] = [float(x) * fac for x in vals[:3]]
            elif key == "ABC":
                abc = [float(x) * fac for x in vals[:3]]
            elif key == "ALPHA_BETA_GAMMA":
                angles = tuple(float(x) for x in vals[:3])
            elif key == "PERIODIC":
                periodic = vals[0].upper() if vals else "XYZ"
        elif sec == "COORD" and "SUBSYS" in stack:
            if key == "SCALED":
                scaled = (len(w) == 1) or w[1].upper() in ("T", "TRUE", ".TRUE.", "YES", "ON")
            elif key == "UNIT":
                cunit = unit_factor(w[1])
            elif len(w) >= 4:
                try:
                    coords.append((w[0], [float(x) for x in w[1:4]]))
                except ValueError:
                    die("%s: cannot parse &COORD line: %s" % (path, ln))
        elif sec == "KIND" and key == "ELEMENT" and cur_kind and len(w) > 1:
            kind_elem[cur_kind] = w[1]
    if not coords:
        die("%s: no &SUBSYS/&COORD block found (coordinates in an external "
            "COORD_FILE_NAME? convert that file instead)" % path)
    if len(vec) == 3:
        cell = np.array([vec["A"], vec["B"], vec["C"]])
    elif abc is not None:
        from ase.geometry import cellpar_to_cell
        cell = cellpar_to_cell(list(abc) + list(angles))
    else:
        cell = np.zeros((3, 3))
    symbols = []
    for name, _ in coords:
        el = kind_elem.get(name) or element_from_label(name)
        if el is None:
            die("%s: cannot map kind '%s' to an element" % (path, name))
        symbols.append(el)
    xyz = np.array([c for _, c in coords])
    pbc = [ax in periodic for ax in "XYZ"] if periodic != "NONE" else [False] * 3
    if scaled:
        return Atoms(symbols, scaled_positions=xyz, cell=cell, pbc=pbc)
    return Atoms(symbols, positions=xyz * cunit, cell=cell, pbc=pbc)


def omx_block(lines, name):
    start, out = "<" + name.lower(), []
    inside = False
    for ln in lines:
        low = ln.lower()
        if not inside and low.startswith(start):
            inside = True
            continue
        if inside:
            if low.startswith(name.lower() + ">"):
                return out
            out.append(ln)
    return out if inside else None


def omx_keyword(lines, name, default=None):
    for ln in lines:
        w = ln.split()
        if len(w) >= 2 and w[0].lower() == name.lower():
            return w[1]
    return default


def read_openmx(path):
    import numpy as np
    from ase import Atoms
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = [strip_comment(ln, "#") for ln in fh]
    lines = [ln for ln in lines if ln]
    species = {}
    for ln in omx_block(lines, "Definition.of.Atomic.Species") or []:
        w = ln.split()
        if len(w) >= 3:
            species[w[0]] = element_from_label(w[2].split("_")[0]) or element_from_label(w[0])
    coord = omx_block(lines, "Atoms.SpeciesAndCoordinates")
    if not coord:
        die("%s: no <Atoms.SpeciesAndCoordinates block" % path)
    cunit = (omx_keyword(lines, "Atoms.SpeciesAndCoordinates.Unit", "Ang") or "Ang").upper()
    symbols, xyz = [], []
    for ln in coord:
        w = ln.split()
        if len(w) < 5:
            continue
        try:
            xyz.append([float(x) for x in w[2:5]])
        except ValueError:
            die("%s: cannot parse coordinate line: %s" % (path, ln))
        el = species.get(w[1]) or element_from_label(w[1])
        if el is None:
            die("%s: cannot map species '%s' to an element" % (path, w[1]))
        symbols.append(el)
    vec = omx_block(lines, "Atoms.UnitVectors")
    cell = None
    if vec:
        vunit = (omx_keyword(lines, "Atoms.UnitVectors.Unit", "Ang") or "Ang").upper()
        cell = np.array([[float(x) for x in ln.split()[:3]] for ln in vec[:3]])
        if vunit == "AU":
            cell *= BOHR
    xyz = np.array(xyz)
    if cunit == "FRAC":
        if cell is None:
            die("%s: FRAC coordinates need <Atoms.UnitVectors" % path)
        return Atoms(symbols, scaled_positions=xyz, cell=cell, pbc=True)
    if cunit == "AU":
        xyz = xyz * BOHR
    if cell is None:
        return Atoms(symbols, positions=xyz, pbc=False)
    return Atoms(symbols, positions=xyz, cell=cell, pbc=True)


def detect_reader(path, in_format):
    if in_format in ("cp2k", "openmx"):
        return in_format
    if in_format:
        return "ase"
    name = os.path.basename(path).lower()
    if name.endswith((".inp", ".restart", ".inc")) or ".restart.bak" in name:
        return "cp2k"
    if name.endswith((".dat", ".dat#")):
        return "openmx"
    return "ase"


def read_any(path, in_format, index):
    if not os.path.isfile(path):
        die("input file not found: %s" % path)
    need("numpy", "pip install numpy")
    need("ase", "pip install ase  (or conda install -c conda-forge ase)")
    kind = detect_reader(path, in_format)
    if kind == "cp2k":
        return read_cp2k(path)
    if kind == "openmx":
        return read_openmx(path)
    from ase.io import read
    try:
        if in_format:
            atoms = read(path, index=index, format=in_format)
        else:
            atoms = read(path, index=index)
    except Exception as exc:  # ASE raises many exception types
        die("ASE cannot read %s (%s); try --in-format (e.g. vasp, cif, xyz, gaussian-in)"
            % (path, exc))
    if isinstance(atoms, list):
        if not atoms:
            die("no structure found in %s" % path)
        atoms = atoms[-1]
    return atoms


# ==========================================================================
# fixed-atom mask <-> ASE constraints
# ==========================================================================
def mask_from_constraints(atoms):
    """(N,3) bool array, True = coordinate fixed (VASP 'F')."""
    import numpy as np
    mask = np.zeros((len(atoms), 3), bool)
    for c in atoms.constraints:
        name = type(c).__name__
        try:
            if name == "FixAtoms":
                mask[c.get_indices()] = True
            elif name in ("FixScaled", "FixCartesian"):
                idx = np.atleast_1d(getattr(c, "index", getattr(c, "a", None)))
                m = np.array(c.mask, bool)
                if m.ndim == 1:
                    m = np.tile(m, (len(idx), 1))
                mask[idx] |= m
            else:
                print("note: constraint %s ignored" % name, file=sys.stderr)
        except Exception as exc:  # unknown ASE internals
            print("note: could not interpret constraint %s (%s)" % (name, exc), file=sys.stderr)
    return mask


def parse_indices(text, n):
    """'1-8,12,15-16' (1-based) -> 0-based list."""
    out = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        m = re.fullmatch(r"(\d+)(?:-(\d+))?", part)
        if not m:
            die("cannot parse --fix-indices item %r (use e.g. 1-8,12)" % part)
        a = int(m.group(1))
        b = int(m.group(2) or a)
        if a < 1 or b > n or a > b:
            die("--fix-indices %r outside 1..%d" % (part, n))
        out.extend(range(a - 1, b))
    return out


def heights(atoms, axis):
    """Height (Angstrom) of every atom along the normal of the plane of the other two vectors."""
    import numpy as np
    cell = atoms.cell.array
    if atoms.cell.rank == 3:
        j, k = [i for i in range(3) if i != axis]
        area = np.linalg.norm(np.cross(cell[j], cell[k]))
        d = abs(np.linalg.det(cell)) / area
        return atoms.get_scaled_positions(wrap=False)[:, axis] * d
    return atoms.positions[:, axis]


def layers(h, tol):
    """Group sorted heights into layers; return list of index arrays (bottom first)."""
    import numpy as np
    order = np.argsort(h)
    groups, cur = [], [order[0]]
    for i in order[1:]:
        if h[i] - h[cur[-1]] > tol:
            groups.append(cur)
            cur = [i]
        else:
            cur.append(i)
    groups.append(cur)
    return groups


# ==========================================================================
# writers
# ==========================================================================
def runs(symbols):
    out = []
    for s in symbols:
        if out and out[-1][0] == s:
            out[-1][1] += 1
        else:
            out.append([s, 1])
    return out


def write_poscar(path, atoms, mask, title):
    sel = bool(mask.any())
    frac = atoms.get_scaled_positions(wrap=False)
    groups = runs(atoms.get_chemical_symbols())
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        f.write(title + "\n   1.0\n")
        for v in atoms.cell.array:
            f.write("  %20.14f %20.14f %20.14f\n" % tuple(v))
        f.write("  " + " ".join("%4s" % s for s, _ in groups) + "\n")
        f.write("  " + " ".join("%4d" % n for _, n in groups) + "\n")
        if sel:
            f.write("Selective dynamics\n")
        f.write("Direct\n")
        for i, r in enumerate(frac):
            line = "  %18.14f %18.14f %18.14f" % tuple(r)
            if sel:
                line += "  " + " ".join("F" if x else "T" for x in mask[i])
            f.write(line + "\n")
    names = [s for s, _ in groups]
    if len(names) != len(set(names)):
        print("note: %s lists repeated species %s; the POTCAR must follow the same order "
              "(or use --sort)" % (path, " ".join(names)), file=sys.stderr)


def write_gjf(path, atoms, mask, a):
    stem = os.path.splitext(os.path.basename(path))[0]
    periodic = atoms.cell.rank > 0 and any(atoms.pbc) and not a.no_tv
    route = a.route or ("PBEPBE/6-31G(d) SCF=Tight" if periodic else "B3LYP/6-31G(d) Opt")
    if not route.lstrip().startswith("#"):
        route = "# " + route
    fixed_any = bool(mask.any())
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        f.write("%%NProcShared=%d\n%%Mem=%s\n%%Chk=%s.chk\n" % (a.nproc, a.gmem, stem))
        f.write(route + "\n\n")
        f.write("%s converted by convert_any.py\n\n" % atoms.get_chemical_formula())
        f.write("%d %d\n" % (a.charge, a.mult))
        for i, (s, p) in enumerate(zip(atoms.get_chemical_symbols(), atoms.positions)):
            if fixed_any:
                code = -1 if mask[i].any() else 0
                f.write(" %-3s %3d %15.8f %15.8f %15.8f\n" % (s, code, p[0], p[1], p[2]))
            else:
                f.write(" %-3s %15.8f %15.8f %15.8f\n" % (s, p[0], p[1], p[2]))
        if periodic:
            for i in range(3):
                if atoms.pbc[i]:
                    f.write(" Tv  %15.8f %15.8f %15.8f\n" % tuple(atoms.cell.array[i]))
        f.write("\n\n")        # blank line ends the molecule section + trailing blank line
    if fixed_any:
        print("note: %s: fixed atoms got freeze code -1 (honoured by Opt)" % path, file=sys.stderr)


def cp2k_periodic(atoms):
    s = "".join(ax for ax, p in zip("XYZ", atoms.pbc) if p)
    return s or "NONE"


def write_cp2k(path, atoms, mask):
    elems = sorted(set(atoms.get_chemical_symbols()), key=atoms.get_chemical_symbols().index)
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        f.write("! CP2K subsystem generated by convert_any.py (DFT-Tutorials)\n")
        f.write("! Use inside &FORCE_EVAL/&SUBSYS with:  @INCLUDE '%s'\n" % os.path.basename(path))
        f.write("! Add one &KIND per element, e.g. (check basis/potential files you use):\n")
        for e in elems:
            f.write("!   &KIND %s\n!     BASIS_SET DZVP-MOLOPT-SR-GTH\n!     POTENTIAL GTH-PBE\n"
                    "!   &END KIND\n" % e)
        f.write("&CELL\n")
        for name, v in zip("ABC", atoms.cell.array):
            f.write("  %s [angstrom] %16.10f %16.10f %16.10f\n" % (name, v[0], v[1], v[2]))
        f.write("  PERIODIC %s\n" % cp2k_periodic(atoms))
        f.write("&END CELL\n&COORD\n")
        for s, p in zip(atoms.get_chemical_symbols(), atoms.positions):
            f.write("  %-3s %16.10f %16.10f %16.10f\n" % (s, p[0], p[1], p[2]))
        f.write("&END COORD\n")
    if not mask.any():
        return None
    fpath = os.path.splitext(path)[0] + "_fixed.inc"
    comp = {}
    for i, m in enumerate(mask):
        if m.any():
            key = "".join(ax for ax, x in zip("XYZ", m) if x)
            comp.setdefault(key, []).append(i + 1)
    with open(fpath, "w", newline="\n", encoding="utf-8") as f:
        f.write("! Fixed atoms generated by convert_any.py: @INCLUDE inside &MOTION/&CONSTRAINT\n")
        f.write("! (only used by RUN_TYPE GEO_OPT, CELL_OPT and MD). Components refer to the\n")
        f.write("! Cartesian X/Y/Z axes (identical to the lattice flags for orthogonal cells).\n")
        for key, idx in comp.items():
            f.write("&FIXED_ATOMS\n")
            for k in range(0, len(idx), 15):
                f.write("  LIST %s\n" % " ".join(str(x) for x in idx[k:k + 15]))
            f.write("  COMPONENTS_TO_FIX %s\n&END FIXED_ATOMS\n" % key)
    return fpath


def write_openmx(path, atoms, mask, valence, unit):
    symbols = atoms.get_chemical_symbols()
    species = sorted(set(symbols), key=symbols.index)
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        f.write("# OpenMX structure blocks generated by convert_any.py (DFT-Tutorials)\n")
        f.write("# Species definition template for DFT_DATA19 (check the recommended\n")
        f.write("# PAO/VPS table of OpenMX; names must match column 2 below):\n")
        f.write("#   Species.Number  %d\n#   <Definition.of.Atomic.Species\n" % len(species))
        for s in species:
            f.write("#    %-3s  %-14s %s_PBE19\n" % (s, OPENMX_PAO_HINT.get(s, "<PAO-of-%s>" % s), s))
        f.write("#   Definition.of.Atomic.Species>\n\n")
        f.write("Atoms.Number  %d\n" % len(atoms))
        f.write("Atoms.SpeciesAndCoordinates.Unit  %s\n" % unit)
        f.write("<Atoms.SpeciesAndCoordinates\n")
        coords = atoms.get_scaled_positions(wrap=False) if unit == "FRAC" else atoms.positions
        for i, (s, r) in enumerate(zip(symbols, coords), start=1):
            nval = valence.get(s)
            ud = "  ?.?    ?.?" if nval is None else "%5.2f  %5.2f" % (nval / 2.0, nval / 2.0)
            f.write("  %4d  %-3s %15.10f %15.10f %15.10f  %s\n" % (i, s, r[0], r[1], r[2], ud))
        f.write("Atoms.SpeciesAndCoordinates>\n")
        if atoms.cell.rank == 3:
            f.write("Atoms.UnitVectors.Unit  Ang\n<Atoms.UnitVectors\n")
            for v in atoms.cell.array:
                f.write("  %16.10f %16.10f %16.10f\n" % tuple(v))
            f.write("Atoms.UnitVectors>\n")
        if mask.any():
            f.write("\n# 1 = fixed, 0 = relaxed (x y z, Cartesian); used by MD.Type optimizers and MD\n")
            f.write("<MD.Fixed.XYZ\n")
            for i, m in enumerate(mask, start=1):
                f.write("  %4d  %d %d %d\n" % (i, int(m[0]), int(m[1]), int(m[2])))
            f.write("MD.Fixed.XYZ>\n")
    missing = sorted(set(symbols) - set(valence))
    if missing:
        print("note: OpenMX valence unknown for %s -> '?.?' written; use --valence El=N"
              % ", ".join(missing), file=sys.stderr)


def write_ase(path, atoms, fmt):
    from ase.io import write
    with open(path, "w", newline="\n", encoding="utf-8") as f:
        write(f, atoms, format=fmt)


def fmt_from_name(path):
    n = os.path.basename(path).lower()
    if n in ("poscar", "contcar") or n.endswith(".vasp") or "poscar" in n or "contcar" in n:
        return "vasp"
    ext = os.path.splitext(n)[1]
    return {".cif": "cif", ".xyz": "xyz", ".extxyz": "extxyz", ".gjf": "gjf", ".com": "gjf",
            ".inc": "cp2k", ".inp": "cp2k", ".dat": "openmx"}.get(ext)


# ==========================================================================
# operations
# ==========================================================================
def standardize(atoms, primitive, symprec):
    from ase import Atoms
    spglib = need("spglib", "pip install spglib")
    cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
    res = spglib.standardize_cell(cell, to_primitive=primitive, no_idealize=False, symprec=symprec)
    if res is None:
        die("spglib could not find the symmetry (try a larger --symprec)")
    lat, pos, num = res
    try:
        sg = spglib.get_spacegroup(cell, symprec=symprec)
    except Exception:
        sg = "?"
    print("spglib (symprec=%g): space group %s -> %s cell with %d atoms"
          % (symprec, sg, "primitive" if primitive else "conventional", len(num)))
    return Atoms(numbers=num, scaled_positions=pos, cell=lat, pbc=True)


def parse_supercell(vals):
    joined = " ".join(vals).lower().replace("x", " ").replace(",", " ")
    try:
        nums = [int(x) for x in joined.split()]
    except ValueError:
        die("--supercell needs 3 or 9 integers (e.g. 2 2 1, 2x2x1)")
    if len(nums) not in (3, 9):
        die("--supercell needs 3 or 9 integers, got %d" % len(nums))
    return nums


def axis_index(text):
    t = text.lower()
    m = {"a": 0, "b": 1, "c": 2, "x": 0, "y": 1, "z": 2}
    if t not in m:
        die("axis must be a/b/c (or x/y/z), got %r" % text)
    return m[t]


def build_parser():
    p = argparse.ArgumentParser(
        prog="convert_any.py",
        description=("Convert structures between CIF/POSCAR/xyz/extxyz/gjf/CP2K/OpenMX with optional "
                     "primitive/conventional cell, supercell, vacuum and fixed layers."),
        epilog=("Examples:\n"
                "  python convert_any.py Si.cif -o POSCAR --primitive\n"
                "  python convert_any.py Si.cif -t vasp cp2k openmx gjf\n"
                "  python convert_any.py POSCAR -o slab.vasp --supercell 3 3 1 --vacuum 15 --fix-layers 2\n"
                "  python convert_any.py slab.vasp -o slab.inc --fix-below 8.0\n"
                "  python convert_any.py h2o.xyz -o h2o.gjf --route 'B3LYP/6-31G(d) Opt Freq'\n"
                "Formats: " + ", ".join(FORMATS) + ", ase:<any ASE format>"),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", help="input structure (ASE formats, CP2K .inp/.restart, OpenMX .dat)")
    p.add_argument("-o", "--output", nargs="+", default=[], metavar="FILE",
                   help="output file(s); format from the name (POSCAR, *.vasp, *.cif, *.gjf, *.inc, *.dat, ...)")
    p.add_argument("-t", "--to", nargs="+", default=[], metavar="FMT",
                   help="output format(s) written as <input-stem>_<fmt-ext>; with one -o: format of that file")
    p.add_argument("--in-format", help="input format: cp2k, openmx or any ASE format name")
    p.add_argument("--index", default="-1", help="frame of a trajectory/OUTCAR (ASE index, default -1)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--primitive", action="store_true", help="spglib standardized primitive cell")
    g.add_argument("--conventional", action="store_true", help="spglib standardized conventional cell")
    p.add_argument("--symprec", type=float, default=1e-3, help="spglib tolerance in Angstrom (default 1e-3)")
    p.add_argument("--supercell", nargs="+", metavar="N",
                   help="3 integers (2 2 1 or 2x2x1) or a 3x3 matrix given as 9 integers (row-wise)")
    p.add_argument("--vacuum", type=float, metavar="ANG",
                   help="TOTAL vacuum thickness along --vacuum-axis (half on each side)")
    p.add_argument("--vacuum-axis", default="c", help="axis for --vacuum/--center (a|b|c, default c)")
    p.add_argument("--center", action="store_true", help="center atoms along --vacuum-axis without changing the cell")
    p.add_argument("--box", type=float, default=10.0, metavar="ANG",
                   help="molecules without cell: padding added to the extent for POSCAR/CP2K/OpenMX (default 10)")
    p.add_argument("--wrap", action="store_true", help="wrap atoms into the cell")
    p.add_argument("--sort", action="store_true", help="sort atoms by element (recommended for POSCAR)")
    p.add_argument("--fix-layers", type=int, metavar="N", help="fix the N lowest layers along --fix-axis")
    p.add_argument("--fix-below", type=float, metavar="ANG",
                   help="fix atoms whose height along --fix-axis is below ANG (height from the cell origin)")
    p.add_argument("--fix-indices", metavar="LIST", help="fix atoms by 1-based index, e.g. 1-8,12")
    p.add_argument("--fix-axis", default="c", help="slab normal for fixing (a|b|c, default c)")
    p.add_argument("--layer-tol", type=float, default=0.5,
                   help="atoms closer than this (Angstrom) along the normal form one layer (default 0.5)")
    p.add_argument("--unfix", action="store_true", help="drop all fixed flags read from the input")
    gg = p.add_argument_group("Gaussian .gjf options")
    gg.add_argument("--route", help="route section (default: B3LYP/6-31G(d) Opt; periodic: PBEPBE/6-31G(d) SCF=Tight)")
    gg.add_argument("--charge", type=int, default=0, help="total charge (default 0)")
    gg.add_argument("--mult", type=int, default=1, help="spin multiplicity (default 1)")
    gg.add_argument("--nproc", type=int, default=8, help="%%NProcShared (default 8)")
    gg.add_argument("--mem", dest="gmem", default="8GB", help="%%Mem (default 8GB)")
    gg.add_argument("--no-tv", action="store_true", help="do not write Tv (translation vector) lines")
    go = p.add_argument_group("OpenMX options")
    go.add_argument("--valence", nargs="*", default=[], metavar="El=N",
                    help="valence electrons of the VPS, e.g. Ti=12 O=6 (initial up/down = N/2)")
    go.add_argument("--omx-unit", choices=("Ang", "FRAC"), default="Ang",
                    help="Atoms.SpeciesAndCoordinates.Unit (default Ang)")
    p.add_argument("-f", "--force", action="store_true", help="overwrite existing output files")
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    if not a.output and not a.to:
        die("give at least one output file (-o) or format (-t); see -h")
    if len(a.to) == 1 and len(a.output) == 1:
        targets = [(a.output[0], a.to[0])]
    elif a.to and a.output:
        die("use either several -o files (format from the names) or -t formats, or one -o with one -t")
    else:
        targets = [(o, None) for o in a.output]
        stem = os.path.splitext(os.path.basename(a.input.rstrip("#")))[0] or "structure"
        for fmt in a.to:
            ext = EXT.get(fmt, "." + fmt.split(":", 1)[-1])
            targets.append((stem + ext, fmt))
    resolved = []
    for path, fmt in targets:
        fmt = fmt or fmt_from_name(path)
        if fmt is None:
            die("cannot guess the format of %s: add -t FMT (one of %s or ase:<format>)"
                % (path, ", ".join(FORMATS)))
        if fmt not in FORMATS and not fmt.startswith("ase:"):
            die("unknown format %r (use %s or ase:<format>)" % (fmt, ", ".join(FORMATS)))
        if os.path.exists(path) and not a.force:
            die("%s exists (use -f to overwrite)" % path)
        if os.path.abspath(path) == os.path.abspath(a.input):
            die("output %s would overwrite the input" % path)
        resolved.append((path, fmt))

    import numpy as np
    index = a.index
    if re.fullmatch(r"-?\d+", index):
        index = int(index)
    atoms = read_any(a.input, a.in_format, index)
    print("read %s: %s, %d atoms, pbc=%s, cell lengths %s"
          % (a.input, atoms.get_chemical_formula(), len(atoms), atoms.pbc.tolist(),
             np.round(atoms.cell.lengths(), 4).tolist()))

    mask = np.zeros((len(atoms), 3), bool) if a.unfix else mask_from_constraints(atoms)
    from ase import Atoms
    atoms = Atoms(atoms.get_chemical_symbols(), positions=atoms.positions,
                  cell=atoms.cell, pbc=atoms.pbc)       # drop CIF/ASE extra arrays

    if a.primitive or a.conventional:
        if atoms.cell.rank != 3:
            die("--primitive/--conventional need a 3D periodic cell")
        if mask.any():
            print("note: fixed flags dropped by the cell standardization", file=sys.stderr)
        atoms = standardize(atoms, a.primitive, a.symprec)
        mask = np.zeros((len(atoms), 3), bool)

    if a.supercell:
        nums = parse_supercell(a.supercell)
        if atoms.cell.rank != 3:
            die("--supercell needs a 3D cell")
        if len(nums) == 3:
            if min(nums) < 1:
                die("--supercell repetitions must be >= 1")
            atoms = atoms.repeat(nums)
            mask = np.tile(mask, (nums[0] * nums[1] * nums[2], 1))
        else:
            from ase.build import make_supercell
            P = np.array(nums).reshape(3, 3)
            if abs(np.linalg.det(P)) < 0.5:
                die("--supercell matrix is singular")
            if mask.any():
                print("note: fixed flags dropped by the general supercell", file=sys.stderr)
            atoms = make_supercell(atoms, P)
            mask = np.zeros((len(atoms), 3), bool)
        print("supercell: %d atoms, cell lengths %s"
              % (len(atoms), np.round(atoms.cell.lengths(), 4).tolist()))

    vax = axis_index(a.vacuum_axis)
    if a.vacuum is not None:
        if a.vacuum < 0:
            die("--vacuum must be >= 0")
        if atoms.cell.rank != 3:
            atoms.center(vacuum=a.box / 2.0)
        atoms.center(vacuum=a.vacuum / 2.0, axis=vax)
        print("vacuum: %.2f A in total along axis %s -> length %.4f A"
              % (a.vacuum, "abc"[vax], atoms.cell.lengths()[vax]))
    elif a.center:
        atoms.center(axis=vax)

    needs_cell = any(f in ("vasp", "cp2k", "openmx", "cif") for _, f in resolved)
    if atoms.cell.rank != 3 and needs_cell:
        atoms.center(vacuum=a.box / 2.0)
        print("molecule without cell: placed in a box with %.1f A padding, lengths %s"
              % (a.box, np.round(atoms.cell.lengths(), 3).tolist()))

    if a.wrap and atoms.cell.rank == 3:
        atoms.wrap()

    if a.sort:
        syms = atoms.get_chemical_symbols()
        order = sorted(range(len(atoms)), key=lambda i: (syms.index(syms[i]), i))
        atoms = atoms[order]
        mask = mask[order]

    fax = axis_index(a.fix_axis)
    if a.fix_layers or a.fix_below is not None or a.fix_indices:
        h = heights(atoms, fax)
        if a.fix_layers:
            groups = layers(h, a.layer_tol)
            if a.fix_layers >= len(groups):
                die("--fix-layers %d but only %d layers found (tolerance %.2f A)"
                    % (a.fix_layers, len(groups), a.layer_tol))
            for g_ in groups[:a.fix_layers]:
                mask[g_] = True
            print("layers along %s: %d found; heights of the fixed ones: %s"
                  % ("abc"[fax], len(groups),
                     ", ".join("%.3f" % np.mean(h[g_]) for g_ in groups[:a.fix_layers])))
        if a.fix_below is not None:
            mask[h < a.fix_below] = True
        if a.fix_indices:
            mask[parse_indices(a.fix_indices, len(atoms))] = True
    nfix = int(mask.any(axis=1).sum())
    if nfix:
        print("fixed atoms: %d of %d" % (nfix, len(atoms)))

    valence = dict(OPENMX_VALENCE)
    for item in a.valence:
        try:
            el, n = item.split("=")
            valence[el] = float(n)
        except ValueError:
            die("--valence items look like Ti=12, got %r" % item)

    title = "%s converted by convert_any.py from %s" % (atoms.get_chemical_formula(),
                                                         os.path.basename(a.input))
    for path, fmt in resolved:
        d = os.path.dirname(path)
        if d:
            os.makedirs(d, exist_ok=True)
        extra = ""
        if fmt == "vasp":
            write_poscar(path, atoms, mask, title)
        elif fmt == "gjf":
            write_gjf(path, atoms, mask, a)
        elif fmt == "cp2k":
            fp = write_cp2k(path, atoms, mask)
            if fp:
                extra = " + %s (&FIXED_ATOMS for &MOTION/&CONSTRAINT)" % fp
        elif fmt == "openmx":
            write_openmx(path, atoms, mask, valence, a.omx_unit)
        elif fmt == "xyz":
            write_ase(path, atoms, "xyz")
        elif fmt == "extxyz":
            write_ase(path, atoms, "extxyz")
        elif fmt == "cif":
            write_ase(path, atoms, "cif")
        else:
            write_ase(path, atoms, fmt.split(":", 1)[1])
        if mask.any() and (fmt in ("xyz", "extxyz", "cif") or fmt.startswith("ase:")):
            extra += " (fixed flags are not stored in this format)"
        print("written %s [%s]%s" % (path, fmt, extra))
    if atoms.cell.rank == 3:
        print("volume %.3f A^3, density %.4f g/cm^3"
              % (atoms.get_volume(), atoms.get_masses().sum() / atoms.get_volume() * 1.66053907))
    return 0


if __name__ == "__main__":
    sys.exit(main())
