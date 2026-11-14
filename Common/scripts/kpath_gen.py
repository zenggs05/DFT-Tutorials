#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
kpath_gen.py -- standard high-symmetry k-path for band structures, written for
               VASP, CP2K and OpenMX at the same time.
               用 SeeK-path (HPKOT 约定) 或 ASE 生成标准高对称路径，同时输出
               VASP line-mode KPOINTS / KPOINTS_OPT、杂化泛函零权重 KPOINTS、
               CP2K &KPOINT_SET 片段和 OpenMX Band.kpath 块。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

IMPORTANT: k-point coordinates only make sense together with the cell they
refer to. With the default method (seekpath) the coordinates refer to the
*standardized primitive cell* of SeeK-path, which may be rotated / re-ordered
with respect to your input. That cell is written to POSCAR_prim: run the SCF
and band calculations with exactly this cell (convert it for CP2K/OpenMX with
convert_any.py), or use --orig-cell to get the same path expressed in the
reciprocal basis of YOUR cell (seekpath >= 2.1).

Files written to the output directory (default: kpath/)
-------------------------------------------------------
  POSCAR_prim             cell the coordinates refer to (seekpath, not --orig-cell)
  KPOINTS_band            VASP line-mode KPOINTS (ICHARG = 11 non-SCF run)
  KPOINTS_OPT             same path for VASP >= 6.3 KPOINTS_OPT (e.g. HSE06 bands)
  KPOINTS_hse             only with --ibzkpt: IBZKPT of the SCF mesh + the path
                          with zero weight (classic hybrid-functional band run)
  cp2k_kpoint_set.inc     &KPOINT_SET blocks for &DFT/&PRINT/&BAND_STRUCTURE
  openmx_band_kpath.txt   Band.dispersion / Band.Nkpath / <Band.kpath block
                          (+ <Band.KPath.UnitCell with the matching cell)
  kpath_info.txt          space group, Bravais lattice, special points, path

Usage examples
--------------
    python kpath_gen.py POSCAR                      # seekpath, 40 points per segment
    python kpath_gen.py Si.cif -n 60 -o si_kpath
    python kpath_gen.py CONTCAR --orig-cell         # path in the basis of the input cell
    python kpath_gen.py POSCAR --ibzkpt scf/IBZKPT --hse-npoints 15
    python kpath_gen.py POSCAR --density 30         # points per 1/Angstrom (CP2K/OpenMX/HSE list)
    python kpath_gen.py POSCAR --method ase         # ASE (Setyawan-Curtarolo) without seekpath
    python kpath_gen.py POSCAR --print              # also print all blocks to the terminal

Dependencies: Python >= 3.9, numpy, ase (reading the structure);
seekpath + spglib for --method seekpath (pip install seekpath).
References: Hinuma et al., Comput. Mater. Sci. 128, 140 (2017),
doi:10.1016/j.commatsci.2016.10.015 (seekpath/HPKOT);
Setyawan & Curtarolo, Comput. Mater. Sci. 49, 299 (2010),
doi:10.1016/j.commatsci.2010.05.010 (ASE bandpath).
"""

import argparse
import math
import os
import sys
import warnings


def die(msg, code=1):
    print("kpath_gen.py: error: " + msg, file=sys.stderr)
    sys.exit(code)


def need(module, hint):
    """Import an optional dependency lazily with a friendly message."""
    try:
        return __import__(module)
    except ImportError:
        die("python package '%s' is required: %s" % (module, hint))


# --------------------------------------------------------------------------
# structure input
# --------------------------------------------------------------------------
def read_structure(path, fmt=None):
    if not os.path.isfile(path):
        die("structure file not found: %s" % path)
    need("ase", "pip install ase  (or conda install -c conda-forge ase)")
    from ase.io import read
    try:
        atoms = read(path, format=fmt) if fmt else read(path)
    except Exception as exc:  # ase raises many different exception types
        die("cannot read %s with ASE (%s); give --format, e.g. --format vasp" % (path, exc))
    if isinstance(atoms, list):
        atoms = atoms[-1]
    if atoms.cell.rank != 3:
        die("%s has no 3D cell: a band path needs a periodic structure" % path)
    return atoms


# --------------------------------------------------------------------------
# path from seekpath or ase
# --------------------------------------------------------------------------
def path_seekpath(atoms, symprec, time_reversal, orig_cell):
    np = need("numpy", "pip install numpy")
    need("spglib", "pip install spglib")
    seekpath = need("seekpath", "pip install seekpath  (or use --method ase)")
    structure = (atoms.cell.array.tolist(), atoms.get_scaled_positions().tolist(),
                 atoms.numbers.tolist())
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        if orig_cell:
            if not hasattr(seekpath, "get_path_orig_cell"):
                die("--orig-cell needs seekpath >= 2.1 (pip install -U seekpath)")
            res = seekpath.get_path_orig_cell(structure, with_time_reversal=time_reversal,
                                              symprec=symprec)
        else:
            res = seekpath.get_path(structure, with_time_reversal=time_reversal,
                                    symprec=symprec)
    for w in caught:
        print("seekpath warning: %s" % w.message, file=sys.stderr)
    info = {
        "method": "seekpath (HPKOT)",
        "spacegroup": "%s (%d)" % (res.get("spacegroup_international", "?"),
                                   res.get("spacegroup_number", 0)),
        "bravais": "%s / %s" % (res.get("bravais_lattice", "?"),
                                res.get("bravais_lattice_extended", "?")),
        "notes": [],
    }
    if orig_cell:
        cell = np.array(atoms.cell.array)
        prim = None
        if res.get("is_supercell"):
            info["notes"].append("input cell is a SUPERCELL: labels are not high-symmetry "
                                 "points of its Brillouin zone (band folding)")
        info["notes"].append("coordinates refer to the reciprocal vectors of the INPUT cell")
    else:
        cell = np.array(res["primitive_lattice"])
        prim = {"cell": cell, "frac": np.array(res["primitive_positions"]),
                "numbers": list(res["primitive_types"])}
        ratio = abs(float(res.get("volume_original_wrt_prim", 1.0)))
        rot = np.array(res.get("rotation_matrix", np.eye(3)))
        changed = abs(ratio - 1.0) > 1e-3 or not np.allclose(rot, np.eye(3), atol=1e-6) \
            or not np.allclose(cell, atoms.cell.array, atol=1e-4)
        info["notes"].append("coordinates refer to the reciprocal vectors of POSCAR_prim "
                             "(seekpath standardized primitive cell)")
        if changed:
            info["notes"].append("POSCAR_prim DIFFERS from the input cell (volume ratio "
                                 "input/prim = %.3f, rotated or re-ordered axes): use "
                                 "POSCAR_prim for SCF and bands, or rerun with --orig-cell" % ratio)
    points = {k: [float(x) for x in v] for k, v in res["point_coords"].items()}
    segments = [(a, b) for a, b in res["path"]]
    return cell, points, segments, prim, info


def path_ase(atoms):
    np = need("numpy", "pip install numpy")
    from ase.dft.kpoints import parse_path_string
    try:
        bp = atoms.cell.bandpath()
    except Exception as exc:
        die("ASE could not determine the Bravais lattice/band path: %s" % exc)
    points = {}
    for lab, k in bp.special_points.items():
        points["GAMMA" if lab == "G" else lab] = [float(x) for x in k]
    segments = []
    for sub in parse_path_string(bp.path):
        sub = ["GAMMA" if s == "G" else s for s in sub]
        segments += list(zip(sub[:-1], sub[1:]))
    try:
        lat = bp.cell.get_bravais_lattice().name
    except Exception:
        lat = "?"
    info = {"method": "ASE bandpath (Setyawan-Curtarolo)", "spacegroup": "n/a (lattice only)",
            "bravais": lat,
            "notes": ["coordinates refer to the reciprocal vectors of the INPUT cell",
                      "ASE uses the lattice symmetry only: atoms that lower the symmetry "
                      "are ignored (use seekpath for crystal symmetry)"]}
    return np.array(atoms.cell.array), points, segments, None, info


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def subpaths(segments):
    """Group (start, end) pairs into continuous sub-paths."""
    out = []
    for a, b in segments:
        if out and out[-1][-1] == a:
            out[-1].append(b)
        else:
            out.append([a, b])
    return out


def path_string(segments):
    return " | ".join("-".join(s) for s in subpaths(segments))


def recip(cell):
    import numpy as np
    return 2.0 * math.pi * np.linalg.inv(np.array(cell)).T


def seg_length(points, a, b, rec):
    import numpy as np
    d = (np.array(points[b]) - np.array(points[a])) @ rec
    return float(np.linalg.norm(d))


def seg_npts(points, a, b, rec, npts, density):
    if density:
        return max(2, int(math.ceil(seg_length(points, a, b, rec) * density)))
    return npts


def label_for(code, lab):
    if lab.upper() == "GAMMA":
        return {"vasp": "GAMMA", "cp2k": "GAMMA", "openmx": "G"}[code]
    return lab


def fmt3(v):
    return "%12.8f %12.8f %12.8f" % tuple(v)


# --------------------------------------------------------------------------
# writers
# --------------------------------------------------------------------------
def kpoints_line_mode(points, segments, npts, title):
    lines = [title, "%d" % npts, "Line-mode", "Reciprocal"]
    for i, (a, b) in enumerate(segments):
        if i:
            lines.append("")
        lines.append("%s   ! %s" % (fmt3(points[a]), label_for("vasp", a)))
        lines.append("%s   ! %s" % (fmt3(points[b]), label_for("vasp", b)))
    return "\n".join(lines) + "\n"


def explicit_path(points, segments, rec, npts, density):
    """Explicit k-point list along the path (no duplicates at joints)."""
    import numpy as np
    kpts = []
    for sub in subpaths(segments):
        for j, (a, b) in enumerate(zip(sub[:-1], sub[1:])):
            n = seg_npts(points, a, b, rec, npts, density)
            ka, kb = np.array(points[a]), np.array(points[b])
            for i in range(n):
                if j > 0 and i == 0:
                    continue        # same as the last point of the previous segment
                lab = a if i == 0 else (b if i == n - 1 else "")
                kpts.append((ka + (kb - ka) * i / (n - 1), lab))
    return kpts


def read_ibzkpt(path):
    if not os.path.isfile(path):
        die("IBZKPT file not found: %s" % path)
    with open(path, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    try:
        nk = int(lines[1].split()[0])
        mode = lines[2].strip().lower()
        if not mode.startswith(("r", "k")):     # Reciprocal (K is accepted by VASP too)
            die("%s: only IBZKPT files in reciprocal coordinates are supported" % path)
        kpts = []
        for ln in lines[3:3 + nk]:
            w = ln.split()
            kpts.append(([float(x) for x in w[:3]], float(w[3])))
    except (IndexError, ValueError):
        die("%s is not a valid IBZKPT/KPOINTS explicit list (line 2 = number of "
            "k-points, line 3 = Reciprocal, then 'kx ky kz weight')" % path)
    if len(kpts) != nk:
        die("%s: header says %d k-points but %d were read" % (path, nk, len(kpts)))
    return kpts


def kpoints_hse(ibz, path_pts, title):
    lines = [title, "%d" % (len(ibz) + len(path_pts)), "Reciprocal lattice"]
    for k, w in ibz:
        lines.append("%s %10.6f" % (fmt3(k), w))
    for k, lab in path_pts:
        s = "%s %10.6f" % (fmt3(k), 0.0)
        if lab:
            s += "  %s" % label_for("vasp", lab)
        lines.append(s)
    return "\n".join(lines) + "\n"


def cp2k_block(points, segments, npts):
    out = ["! k-path for CP2K, generated by kpath_gen.py (DFT-Tutorials)",
           "! Paste inside &FORCE_EVAL/&DFT/&PRINT/&BAND_STRUCTURE (next to ADDED_MOS, FILE_NAME).",
           "! One &KPOINT_SET per continuous part of the path; coordinates in units of",
           "! the reciprocal vectors of the cell given in &SUBSYS/&CELL (must be the same cell).",
           "! NPOINTS = number of points on each line between two special points."]
    for sub in subpaths(segments):
        out.append("&KPOINT_SET")
        out.append("  UNITS B_VECTOR")
        out.append("  NPOINTS %d" % npts)
        for lab in sub:
            out.append("  SPECIAL_POINT %-8s %s" % (label_for("cp2k", lab), fmt3(points[lab])))
        out.append("&END KPOINT_SET")
    return "\n".join(out) + "\n"


def openmx_block(points, segments, rec, npts, density, cell):
    out = ["# k-path for OpenMX, generated by kpath_gen.py (DFT-Tutorials)",
           "# Coordinates are fractions of the reciprocal vectors of the cell below",
           "# (Band.KPath.UnitCell, same unit as Atoms.UnitVectors.Unit = Ang).",
           "# If Atoms.UnitVectors is exactly this cell, the UnitCell block may be removed.",
           "Band.dispersion              on",
           "<Band.KPath.UnitCell"]
    for v in cell:
        out.append("  %s" % fmt3(v))
    out.append("Band.KPath.UnitCell>")
    out.append("Band.Nkpath                  %d" % len(segments))
    out.append("<Band.kpath")
    for a, b in segments:
        n = seg_npts(points, a, b, rec, npts, density)
        out.append("  %4d  %s   %s   %s %s" % (n, fmt3(points[a]), fmt3(points[b]),
                                              label_for("openmx", a), label_for("openmx", b)))
    out.append("Band.kpath>")
    return "\n".join(out) + "\n"


def write_poscar(path, prim):
    from ase import Atoms
    from ase.io import write
    atoms = Atoms(numbers=prim["numbers"], scaled_positions=prim["frac"],
                  cell=prim["cell"], pbc=True)
    with open(path, "w", newline="\n") as fh:
        write(fh, atoms, format="vasp", direct=True)
    return atoms


def save(path, text):
    with open(path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write(text)


# --------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        prog="kpath_gen.py",
        description=("Standard high-symmetry k-path (seekpath/HPKOT or ASE) written as VASP "
                     "line-mode KPOINTS + KPOINTS_OPT, zero-weight HSE KPOINTS, CP2K "
                     "&KPOINT_SET and OpenMX Band.kpath."),
        epilog=("Examples:\n"
                "  python kpath_gen.py POSCAR\n"
                "  python kpath_gen.py Si.cif -n 60 -o si_kpath\n"
                "  python kpath_gen.py CONTCAR --orig-cell\n"
                "  python kpath_gen.py POSCAR --ibzkpt scf/IBZKPT --hse-npoints 15\n"
                "  python kpath_gen.py POSCAR --method ase --print\n"
                "Use the SAME cell for SCF and band runs as the one the path refers to\n"
                "(POSCAR_prim by default)."),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("structure", help="any structure file ASE can read (POSCAR, CIF, extxyz, ...)")
    p.add_argument("--format", help="ASE format of the structure file if not detected (e.g. vasp, cif)")
    p.add_argument("--method", choices=("seekpath", "ase"), default="seekpath",
                   help="path convention: seekpath/HPKOT (default) or ASE Setyawan-Curtarolo")
    p.add_argument("--orig-cell", action="store_true",
                   help="seekpath: express the path in the reciprocal basis of the input cell")
    p.add_argument("--symprec", type=float, default=1e-3,
                   help="spglib symmetry tolerance in Angstrom (default 1e-3)")
    p.add_argument("--no-time-reversal", action="store_true",
                   help="seekpath: no time-reversal symmetry (adds -k segments if no inversion)")
    p.add_argument("-n", "--npoints", type=int, default=40,
                   help="points per segment for line-mode KPOINTS, CP2K and OpenMX (default 40)")
    p.add_argument("--density", type=float,
                   help="points per 1/Angstrom along each segment (OpenMX and explicit HSE list)")
    p.add_argument("--hse-npoints", type=int, default=20,
                   help="points per segment in the zero-weight list (default 20)")
    p.add_argument("--ibzkpt", metavar="FILE",
                   help="IBZKPT of the SCF run -> write KPOINTS_hse (weights + zero-weight path)")
    p.add_argument("-o", "--outdir", default="kpath", help="output directory (default kpath)")
    p.add_argument("--print", dest="show", action="store_true",
                   help="also print every block to the terminal")
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    if a.npoints < 2 or a.hse_npoints < 2:
        die("--npoints and --hse-npoints must be >= 2")
    if a.density is not None and a.density <= 0:
        die("--density must be positive")
    if a.method == "ase" and (a.orig_cell or a.no_time_reversal):
        print("note: --orig-cell/--no-time-reversal are seekpath options; ignored with --method ase",
              file=sys.stderr)

    atoms = read_structure(a.structure, a.format)
    if a.method == "seekpath":
        cell, points, segments, prim, info = path_seekpath(
            atoms, a.symprec, not a.no_time_reversal, a.orig_cell)
    else:
        cell, points, segments, prim, info = path_ase(atoms)
    if not segments:
        die("no path segments returned")
    rec = recip(cell)

    os.makedirs(a.outdir, exist_ok=True)
    written = []
    if prim is not None:
        pa = write_poscar(os.path.join(a.outdir, "POSCAR_prim"), prim)
        written.append("POSCAR_prim")
        nat = len(pa)
    else:
        nat = len(atoms)

    title = "%s k-path %s (kpath_gen.py, %s)" % (atoms.get_chemical_formula(),
                                                 path_string(segments), info["method"])
    texts = {}
    texts["KPOINTS_band"] = kpoints_line_mode(points, segments, a.npoints, title)
    texts["KPOINTS_OPT"] = texts["KPOINTS_band"]
    if a.ibzkpt:
        ibz = read_ibzkpt(a.ibzkpt)
        pts = explicit_path(points, segments, rec, a.hse_npoints, a.density)
        texts["KPOINTS_hse"] = kpoints_hse(ibz, pts, "SCF mesh from %s + %d zero-weight path points"
                                           % (os.path.basename(a.ibzkpt), len(pts)))
    texts["cp2k_kpoint_set.inc"] = cp2k_block(points, segments, a.npoints)
    texts["openmx_band_kpath.txt"] = openmx_block(points, segments, rec, a.npoints, a.density, cell)

    info_lines = ["structure      : %s (%s, %d atoms in the reference cell)"
                  % (a.structure, atoms.get_chemical_formula(), nat),
                  "method         : %s" % info["method"],
                  "space group    : %s" % info["spacegroup"],
                  "Bravais lattice: %s" % info["bravais"],
                  "path           : %s" % path_string(segments),
                  "", "special points (fractions of reciprocal vectors)       |k| (1/A)"]
    for lab in sorted(points, key=lambda s: (s != "GAMMA", s)):
        kc = [sum(points[lab][i] * rec[i][j] for i in range(3)) for j in range(3)]
        info_lines.append("  %-8s %s   %8.4f" % (lab, fmt3(points[lab]),
                                                  math.sqrt(sum(x * x for x in kc))))
    info_lines.append("")
    info_lines.append("segment lengths (1/A, 2*pi included):")
    for s, e in segments:
        info_lines.append("  %-8s -> %-8s %8.4f" % (s, e, seg_length(points, s, e, rec)))
    info_lines.append("")
    info_lines += ["NOTE: " + n for n in info["notes"]]
    if a.ibzkpt:
        info_lines.append("NOTE: KPOINTS_hse is only valid if IBZKPT was produced with the SAME "
                          "cell (and k-mesh/ISYM) used for the band run")
    texts["kpath_info.txt"] = "\n".join(info_lines) + "\n"

    for name, text in texts.items():
        save(os.path.join(a.outdir, name), text)
        written.append(name)

    print("\n".join(info_lines))
    if a.show:
        for name in ("KPOINTS_band", "KPOINTS_hse", "cp2k_kpoint_set.inc", "openmx_band_kpath.txt"):
            if name in texts:
                print("\n===== %s =====" % name)
                print(texts[name], end="")
    print("\nwritten to %s/: %s" % (a.outdir, " ".join(written)))
    print("VASP: ICHARG = 11 + KPOINTS_band (non-SCF), or VASP >= 6.3 KPOINTS_OPT in the SCF run;"
          " plot with VASP/scripts/plot_band.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
