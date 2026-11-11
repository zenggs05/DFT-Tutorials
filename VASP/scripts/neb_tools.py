#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
neb_tools.py -- helper for VASP (VTST) nudged-elastic-band calculations.
                VASP/VTST NEB 全流程辅助：插值建目录、检查原子顺序与图像间距、收集能量与切向力、
                画样条拟合能垒曲线、导出 XYZ 动画。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Sub-commands
------------
interp   IS FS -n N    linear or IDPP (ASE) interpolation -> directories 00 .. N+1
                       (same layout as VTST nebmake.pl); Selective dynamics of the
                       initial state are kept; minimum-image convention by default.
check    IS FS         same species/order/cell?  per-atom displacement, atoms that
                       move suspiciously far (wrong atom order), optimal re-ordering
                       by the Hungarian algorithm (--write-reordered).
check    DIR           image spacing (like dist.pl), shortest interatomic distance
                       per image, total path length.
collect  DIR           energy (energy(sigma->0) or TOTEN) and VTST tangent force of
                       every image, reaction coordinate from the geometries, force-
                       based cubic Hermite spline (natural cubic spline without VTST
                       forces), forward/backward barrier, reaction energy; writes
                       neb.dat, spline.dat and a PNG.
xyz      DIR           extended-XYZ movie of the band (CONTCAR, fallback POSCAR),
                       unwrapped so that atoms do not jump across the cell.

Usage examples
--------------
    python neb_tools.py interp POSCAR_IS POSCAR_FS -n 5                 # linear
    python neb_tools.py interp POSCAR_IS POSCAR_FS -n 5 --method idpp \\
           --is-outcar is/OUTCAR --fs-outcar fs/OUTCAR
    python neb_tools.py check POSCAR_IS POSCAR_FS --write-reordered POSCAR_FS_new
    python neb_tools.py check .                 # spacing of 00..NN
    python neb_tools.py collect . -o neb        # neb.dat, spline.dat, neb.png
    python neb_tools.py xyz . -o movie.xyz

Notes
-----
* The tangent force is read from the VTST line
  'NEB: projections on to tangent (spring, REAL)' (last number); the max
  NEB force from 'FORCES: max atom, RMS'.  Without VTST these lines are
  missing and a plain cubic spline through the energies is used.
* The end points 00 and NN need an OUTCAR (copy the OUTCAR of the relaxed
  initial/final states there, as for nebresults.pl).
* MPI ranks must be divisible by IMAGES; run a plain NEB first, then LCLIMB.

Dependencies: python >= 3.9, numpy, matplotlib, scipy; ase only for --method idpp.
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

import numpy as np

__version__ = "1.0.0"

FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"
C1, C2 = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def fnum(s):
    try:
        return float(str(s).replace("D", "E"))
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# POSCAR I/O
# ---------------------------------------------------------------------------
class Structure:
    def __init__(self, comment, lattice, species, counts, frac, selective=None):
        self.comment = comment
        self.lattice = np.asarray(lattice, float)
        self.species = list(species)
        self.counts = list(counts)
        self.frac = np.asarray(frac, float).reshape(-1, 3)
        self.selective = None if selective is None else np.asarray(selective, bool)

    @property
    def natoms(self):
        return int(sum(self.counts))

    @property
    def symbols(self):
        out = []
        for s, n in zip(self.species, self.counts):
            out += [s] * n
        return out

    @property
    def cart(self):
        return self.frac @ self.lattice

    def copy(self, frac=None):
        return Structure(self.comment, self.lattice.copy(), self.species, self.counts,
                         self.frac.copy() if frac is None else frac,
                         None if self.selective is None else self.selective.copy())


def read_poscar(path) -> Structure:
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"'{path}' not found")
    lines = p.read_text(errors="replace").splitlines()
    if len(lines) < 8:
        raise ValueError(f"'{path}' is too short to be a POSCAR")
    try:
        scale = [float(x) for x in lines[1].split()[:3]]
        lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
    except (ValueError, IndexError):
        raise ValueError(f"malformed lattice in '{path}'") from None
    i = 5
    tok = lines[i].split()
    if all(t.isdigit() for t in tok):
        cand = lines[0].split()
        species = cand[:len(tok)] if len(cand) >= len(tok) else [f"X{k + 1}" for k in range(len(tok))]
    else:
        species = [re.match(r"[A-Za-z]+", t).group(0)[:2] if re.match(r"[A-Za-z]+", t) else t for t in tok]
        species = [s[0].upper() + (s[1:].lower() if len(s) > 1 and s[1].islower() else "") for s in species]
        i += 1
    try:
        counts = [int(x) for x in lines[i].split()]
    except ValueError:
        raise ValueError(f"cannot read atom counts in '{path}'") from None
    if len(counts) != len(species):
        raise ValueError(f"'{path}': {len(species)} symbols but {len(counts)} counts")
    i += 1
    sel = lines[i].strip()[:1] in ("s", "S")
    if sel:
        i += 1
    cart = lines[i].strip()[:1] in ("c", "C", "k", "K")
    i += 1
    nat = sum(counts)
    pos, flags = [], []
    for j in range(nat):
        try:
            t = lines[i + j].split()
            pos.append([float(x) for x in t[:3]])
        except (IndexError, ValueError):
            raise ValueError(f"cannot read atom {j + 1} in '{path}'") from None
        if sel:
            flags.append([x.upper().startswith("T") for x in t[3:6]] if len(t) >= 6 else [True] * 3)
    pos = np.array(pos)
    if len(scale) == 1:
        s = (abs(scale[0]) / abs(np.linalg.det(lat))) ** (1 / 3) if scale[0] < 0 else scale[0]
        lat_s = lat * s
        if cart:
            pos = (pos * s) @ np.linalg.inv(lat_s)
    else:
        sv = np.array(scale[:3])
        lat_s = lat * sv
        if cart:
            pos = (pos * sv) @ np.linalg.inv(lat_s)
    return Structure(lines[0].strip(), lat_s, species, counts, pos, flags if sel else None)


def write_poscar(path, s: Structure, comment=None) -> None:
    out = [comment or s.comment or "neb_tools.py", "   1.0"]
    out += ["  " + "  ".join(f"{x:20.12f}" for x in v) for v in s.lattice]
    out.append("  " + "  ".join(f"{e:>4s}" for e in s.species))
    out.append("  " + "  ".join(f"{n:>4d}" for n in s.counts))
    if s.selective is not None:
        out.append("Selective dynamics")
    out.append("Direct")
    for k, r in enumerate(s.frac):
        row = "  " + "  ".join(f"{x:18.12f}" for x in r)
        if s.selective is not None:
            row += "  " + " ".join("T" if f else "F" for f in s.selective[k])
        out.append(row)
    Path(path).write_text("\n".join(out) + "\n")


def load_geometry(d: Path, prefer_contcar=True):
    """CONTCAR of an image if usable, else POSCAR."""
    names = ("CONTCAR", "POSCAR") if prefer_contcar else ("POSCAR",)
    for nm in names:
        f = d / nm
        if f.is_file() and f.stat().st_size > 0:
            try:
                return read_poscar(f), nm
            except (ValueError, FileNotFoundError):
                continue
    return None, None


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------
def mic_frac(df):
    return df - np.round(df)


def displacement(a: Structure, b: Structure, mic=True):
    """Cartesian displacement vectors b - a (minimum image)."""
    df = b.frac - a.frac
    if mic:
        df = mic_frac(df)
    return df @ a.lattice


def image_distance(a: Structure, b: Structure, mic=True) -> float:
    """Total distance between two images, sqrt(sum_i |dr_i|^2) (as VTST dist.pl)."""
    return float(np.sqrt(np.sum(displacement(a, b, mic) ** 2)))


def min_distance(s: Structure):
    """Shortest interatomic distance (periodic, 27 neighbouring cells)."""
    n = s.natoms
    if n < 2:
        return None, None
    shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    best, pair = np.inf, None
    cart = s.cart
    for a in range(n):
        df = s.frac - s.frac[a]
        df = mic_frac(df)
        for sh in shifts:
            v = (df + sh) @ s.lattice
            dist = np.linalg.norm(v, axis=1)
            dist[a] = np.inf if not sh.any() else dist[a]
            k = int(np.argmin(dist))
            if dist[k] < best:
                best, pair = float(dist[k]), (a + 1, k + 1)
    del cart
    return best, pair


def image_dirs(root: Path):
    ds = sorted(x for x in root.iterdir() if x.is_dir() and re.fullmatch(r"\d{2,3}", x.name))
    if len(ds) < 3:
        die(f"'{root}' does not contain image directories 00, 01, ... (found {len(ds)})")
    return ds


# ---------------------------------------------------------------------------
# interp
# ---------------------------------------------------------------------------
def check_compatible(a: Structure, b: Structure, allow_cell=False):
    if a.species != b.species or a.counts != b.counts:
        die(f"species/counts differ: IS {list(zip(a.species, a.counts))} vs FS {list(zip(b.species, b.counts))}")
    dl = np.abs(a.lattice - b.lattice).max()
    if dl > 1e-3 and not allow_cell:
        die(f"lattices of IS and FS differ by {dl:.4f} A; NEB needs the same cell "
            "(use --allow-cell-change only for variable-cell NEB)")
    return dl


def ase_idpp(images, steps, fmax):
    """Run IDPP with ASE on a list of Structure objects (modified in place)."""
    try:
        from ase import Atoms
        from ase.constraints import FixAtoms, FixScaled
    except ImportError:
        die("--method idpp needs ASE (pip install ase)")
    try:
        from ase.mep import NEB, idpp_interpolate  # ASE >= 3.23
    except ImportError:
        try:
            from ase.neb import NEB, idpp_interpolate  # older ASE
        except ImportError:
            die("this ASE version has no NEB/idpp_interpolate")
    atoms_list = []
    for s in images:
        at = Atoms(s.symbols, scaled_positions=s.frac, cell=s.lattice, pbc=True)
        if s.selective is not None:
            cons = []
            full = [k for k, f in enumerate(s.selective) if not f.any()]
            if full:
                cons.append(FixAtoms(indices=full))
            for k, f in enumerate(s.selective):
                if f.any() and not f.all():
                    cons.append(FixScaled(k, mask=~f))
            if cons:
                at.set_constraint(cons)
        atoms_list.append(at)
    neb = NEB(atoms_list)
    idpp_interpolate(neb, traj=None, log=None, fmax=fmax, steps=steps, mic=True)
    for s, at in zip(images, atoms_list):
        # keep the continuous (unwrapped) representation
        new = np.linalg.solve(s.lattice.T, at.positions.T).T
        s.frac = new


def cmd_interp(a):
    try:
        s_is, s_fs = read_poscar(a.initial), read_poscar(a.final)
    except (FileNotFoundError, ValueError) as exc:
        die(str(exc))
    check_compatible(s_is, s_fs, a.allow_cell_change)
    n = a.nimages
    if n < 1:
        die("-n must be >= 1")
    out = Path(a.outdir)
    names = [f"{k:02d}" for k in range(n + 2)]
    if any((out / nm).exists() for nm in names) and not a.force:
        die(f"image directories already exist in '{out}' (use -f to overwrite)")
    df = s_fs.frac - s_is.frac
    if not a.no_mic:
        df = mic_frac(df)
    disp = np.linalg.norm(df @ s_is.lattice, axis=1)
    images = []
    for k in range(n + 2):
        t = k / (n + 1)
        img = s_is.copy(frac=s_is.frac + t * df)
        if a.allow_cell_change:
            img.lattice = (1 - t) * s_is.lattice + t * s_fs.lattice
        img.comment = f"image {k:02d} ({a.method}) from {Path(a.initial).name} -> {Path(a.final).name}"
        images.append(img)
    if a.method == "idpp":
        ase_idpp(images, a.idpp_steps, a.idpp_fmax)
        # end points stay exactly the input structures
        images[0].frac = s_is.frac.copy()
        images[-1].frac = s_is.frac + df
    for k, img in enumerate(images):
        d = out / names[k]
        d.mkdir(parents=True, exist_ok=True)
        if not a.keep_unwrapped:
            img.frac = np.mod(img.frac, 1.0)
        write_poscar(d / "POSCAR", img)
    if a.is_outcar:
        shutil.copy(a.is_outcar, out / names[0] / "OUTCAR")
    if a.fs_outcar:
        shutil.copy(a.fs_outcar, out / names[-1] / "OUTCAR")
    print(f"wrote {names[0]} .. {names[-1]} in '{out}' ({n} intermediate images, method {a.method})")
    print(f"largest atomic displacement IS->FS: {disp.max():.3f} A (atom {int(disp.argmax()) + 1}); "
          f"total distance {np.sqrt((disp ** 2).sum()):.3f} A")
    report_band(images, a.dmin)
    print(f"INCAR: IMAGES = {n}; MPI ranks must be a multiple of {n}")
    if not a.is_outcar:
        print(f"remember: copy the OUTCAR of the relaxed IS/FS into {names[0]}/ and {names[-1]}/ for 'collect'")


def report_band(images, dmin):
    print(f"{'image':>5} {'d(prev) A':>10} {'d_min A':>8}  closest pair")
    for k, img in enumerate(images):
        dprev = image_distance(images[k - 1], img) if k else 0.0
        md, pair = min_distance(img)
        flag = "  <-- atoms too close!" if md is not None and md < dmin else ""
        print(f"{k:>5d} {dprev:>10.3f} {md if md is not None else float('nan'):>8.3f}  {pair}{flag}")


# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------
def best_order(s_is: Structure, s_fs: Structure):
    """Hungarian assignment (per species) of FS atoms to IS atoms."""
    try:
        from scipy.optimize import linear_sum_assignment
    except ImportError:
        die("scipy is needed for the re-ordering analysis")
    order = np.arange(s_is.natoms)
    start = 0
    for n in s_is.counts:
        idx = np.arange(start, start + n)
        dfr = s_fs.frac[idx][None, :, :] - s_is.frac[idx][:, None, :]
        dfr = mic_frac(dfr)
        cost = np.linalg.norm(dfr @ s_is.lattice, axis=2)
        r, c = linear_sum_assignment(cost)
        order[idx[r]] = idx[c]
        start += n
    return order


def cmd_check(a):
    targets = a.paths
    if len(targets) == 2 and all(Path(t).is_file() for t in targets):
        try:
            s_is, s_fs = read_poscar(targets[0]), read_poscar(targets[1])
        except (FileNotFoundError, ValueError) as exc:
            die(str(exc))
        dl = check_compatible(s_is, s_fs, allow_cell=True)
        print(f"species / counts : {list(zip(s_is.species, s_is.counts))}  (identical)")
        print(f"max lattice diff : {dl:.5f} A" + ("  <-- NEB needs identical cells" if dl > 1e-3 else ""))
        disp = displacement(s_is, s_fs)
        dn = np.linalg.norm(disp, axis=1)
        order = np.argsort(-dn)
        sym = s_is.symbols
        print(f"total IS->FS distance: {np.sqrt((dn ** 2).sum()):.3f} A; atoms moving > {a.move:g} A:")
        for k in order[: a.top]:
            if dn[k] < a.move:
                break
            print(f"   atom {k + 1:4d} {sym[k]:>2s}  {dn[k]:7.3f} A")
        suspicious = [k for k in range(len(dn)) if dn[k] > a.max_disp]
        new = best_order(s_is, s_fs)
        changed = [(k, new[k]) for k in range(len(new)) if new[k] != k]
        if changed:
            dn2 = np.linalg.norm(mic_frac(s_fs.frac[new] - s_is.frac) @ s_is.lattice, axis=1)
            print(f"\nA different atom order of FS reduces the total distance "
                  f"{np.sqrt((dn ** 2).sum()):.3f} -> {np.sqrt((dn2 ** 2).sum()):.3f} A:")
            for k, j in changed[: a.top]:
                print(f"   IS atom {k + 1:4d} ({sym[k]}) <- FS atom {j + 1:4d}")
            if a.write_reordered:
                s_new = s_fs.copy(frac=s_fs.frac[new])
                if s_fs.selective is not None:
                    s_new.selective = s_fs.selective[new]
                write_poscar(a.write_reordered, s_new, s_fs.comment + " (re-ordered to match IS)")
                print(f"re-ordered final state written to {a.write_reordered}")
            else:
                print("   (use --write-reordered FILE to write the re-ordered final state)")
        elif suspicious:
            print(f"\natoms {[k + 1 for k in suspicious]} move more than {a.max_disp:g} A - check the structures")
        else:
            print("\natom order looks consistent (identity is the optimal assignment)")
        md_is, p_is = min_distance(s_is)
        md_fs, p_fs = min_distance(s_fs)
        print(f"shortest distance IS {md_is:.3f} A {p_is}, FS {md_fs:.3f} A {p_fs}")
        return
    root = Path(targets[0])
    if not root.is_dir():
        die("check needs either two POSCAR files (IS FS) or one directory with 00..NN")
    imgs = []
    for d in image_dirs(root):
        s, src = load_geometry(d, prefer_contcar=not a.poscar)
        if s is None:
            die(f"no POSCAR/CONTCAR in {d}")
        imgs.append(s)
    for s in imgs[1:]:
        if s.species != imgs[0].species or s.counts != imgs[0].counts:
            die("images contain different species/atom counts")
    d = [image_distance(imgs[k - 1], imgs[k]) for k in range(1, len(imgs))]
    report_band(imgs, a.dmin)
    d = np.array(d)
    print(f"path length {d.sum():.3f} A; spacing mean {d.mean():.3f} A, min {d.min():.3f}, max {d.max():.3f}"
          + ("  <-- uneven spacing (>50 % deviation)" if d.max() > 1.5 * d.mean() or d.min() < 0.5 * d.mean() else ""))


# ---------------------------------------------------------------------------
# collect
# ---------------------------------------------------------------------------
def parse_image(d: Path, nions_hint=None):
    """Energy, tangent force, max force and EDIFFG of one image directory."""
    res = {"e0": None, "toten": None, "ftan": None, "fmax_vtst": None, "fmax": None, "ediffg": None,
           "src": None}
    oc = d / "OUTCAR"
    if oc.is_file():
        res["src"] = "OUTCAR"
        nions = nions_hint
        with open(oc, errors="replace") as fh:
            it = iter(fh)
            for line in it:
                if "NIONS" in line and nions is None:
                    m = re.search(r"NIONS\s*=\s*(\d+)", line)
                    if m:
                        nions = int(m.group(1))
                elif "EDIFFG" in line and res["ediffg"] is None:
                    m = re.search(r"EDIFFG\s*=\s*(" + FLOAT + ")", line)
                    if m:
                        res["ediffg"] = fnum(m.group(1))
                elif "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in line:
                    for _ in range(6):
                        nxt = next(it, "")
                        m = re.search(r"TOTEN\s*=\s*(" + FLOAT + ")", nxt)
                        if m:
                            res["toten"] = fnum(m.group(1))
                        m = re.search(r"energy\(sigma->0\)\s*=\s*(" + FLOAT + ")", nxt)
                        if m:
                            res["e0"] = fnum(m.group(1))
                            break
                elif "NEB: projections on to tangent" in line:
                    nums = re.findall(FLOAT, line.split(":")[-1])
                    if nums:
                        res["ftan"] = fnum(nums[-1])
                elif "FORCES: max atom, RMS" in line:
                    nums = re.findall(FLOAT, line.split("RMS")[-1])
                    if nums:
                        res["fmax_vtst"] = fnum(nums[0])
                elif "TOTAL-FORCE" in line and "POSITION" in line and nions:
                    next(it, "")
                    rows = []
                    for _ in range(nions):
                        t = next(it, "").split()
                        if len(t) >= 6:
                            rows.append([fnum(x) for x in t[3:6]])
                    if len(rows) == nions and all(None not in r for r in rows):
                        res["fmax"] = float(np.max(np.linalg.norm(np.array(rows), axis=1)))
    if res["e0"] is None and (d / "OSZICAR").is_file():
        for line in (d / "OSZICAR").read_text(errors="replace").splitlines():
            m = re.search(r"\bE0=\s*(" + FLOAT + ")", line)
            if m:
                res["e0"] = fnum(m.group(1))
                res["src"] = res["src"] or "OSZICAR"
            m = re.search(r"\bF=\s*(" + FLOAT + ")", line)
            if m and res["toten"] is None:
                pass
    return res


def fit_curve(s, e, ftan, npts):
    """Return (x, y, method) of the energy profile fit."""
    from scipy.interpolate import CubicHermiteSpline, CubicSpline
    xs = np.linspace(s[0], s[-1], npts)
    if all(f is not None for f in ftan) and len(s) >= 2:
        dydx = -np.array(ftan, float)
        return xs, CubicHermiteSpline(s, e, dydx)(xs), "force-based cubic Hermite spline (VTST tangent forces)"
    if len(s) >= 3:
        return xs, CubicSpline(s, e, bc_type="natural")(xs), "natural cubic spline through energies (no tangent forces)"
    return xs, np.interp(xs, s, e), "linear interpolation"


def cmd_collect(a):
    root = Path(a.dir)
    if not root.is_dir():
        die(f"'{root}' is not a directory")
    dirs = image_dirs(root)
    geoms = []
    for d in dirs:
        g, _ = load_geometry(d, prefer_contcar=not a.poscar)
        if g is None:
            die(f"no POSCAR/CONTCAR in {d}")
        geoms.append(g)
    nions = geoms[0].natoms
    data = [parse_image(d, nions) for d in dirs]
    key = "toten" if a.energy == "toten" else "e0"
    missing = [d.name for d, r in zip(dirs, data) if r[key] is None]
    if missing:
        hint = ""
        if dirs[0].name in missing or dirs[-1].name in missing:
            hint = " - copy the OUTCAR of the relaxed initial/final state into the end-point directories"
        die(f"no energy found for image(s) {missing}{hint}")
    E = np.array([r[key] for r in data])
    s = np.concatenate([[0.0], np.cumsum([image_distance(geoms[k - 1], geoms[k]) for k in range(1, len(geoms))])])
    ftan = [0.0 if k in (0, len(dirs) - 1) else data[k]["ftan"] for k in range(len(dirs))]
    have_vtst = all(f is not None for f in ftan[1:-1])
    if not have_vtst:
        ftan = [None] * len(dirs)
    rel = E - E[0]
    xs, ys, method = fit_curve(s, rel, ftan, a.points)
    kmax = int(np.argmax(ys))
    ts_img = int(np.argmax(rel[1:-1])) + 1 if len(rel) > 2 else 0
    ea_f = float(ys[kmax])
    ea_b = float(ys[kmax] - rel[-1])
    de = float(rel[-1])
    fm = [r["fmax_vtst"] if r["fmax_vtst"] is not None else r["fmax"] for r in data]
    ediffg = next((r["ediffg"] for r in data[1:-1] if r["ediffg"] is not None), None)

    print(f"{'image':>5} {'s (A)':>8} {'E (eV)':>16} {'E-E0 (eV)':>10} {'F_tan':>9} {'F_max':>8}  source")
    for k, d in enumerate(dirs):
        ft = ftan[k]
        print(f"{d.name:>5} {s[k]:8.3f} {E[k]:16.6f} {rel[k]:10.4f} "
              f"{(f'{ft:9.4f}' if ft is not None else '        -')} "
              f"{(f'{fm[k]:8.4f}' if fm[k] is not None else '       -')}  {data[k]['src']}")
    print()
    print(f"fit: {method}")
    print(f"highest image      : {dirs[ts_img].name} (E - E0 = {rel[ts_img]:.4f} eV)")
    print(f"spline maximum     : s = {xs[kmax]:.3f} A, E - E0 = {ea_f:.4f} eV")
    print(f"forward barrier    : {ea_f:.4f} eV")
    print(f"backward barrier   : {ea_b:.4f} eV")
    print(f"reaction energy dE : {de:.4f} eV  ({'exothermic' if de < 0 else 'endothermic'})")
    inner = [f for f in fm[1:-1] if f is not None]
    if inner:
        crit = abs(ediffg) if (ediffg is not None and ediffg < 0) else None
        status = ""
        if crit is not None:
            status = "  -> converged" if max(inner) < crit else f"  -> NOT converged (|EDIFFG| = {crit:g})"
        print(f"max force (images) : {max(inner):.4f} eV/A{status}")
    if not have_vtst:
        print("note: no VTST tangent forces found ('NEB: projections on to tangent'); "
              "barrier from the energies only")
    if ts_img in (1, len(dirs) - 2) and len(dirs) > 4:
        print("note: the highest image is next to an end point - consider more images or check the end states")

    prefix = a.output
    with open(prefix + ".dat", "w") as fh:
        fh.write("# image  s(A)  E-E0(eV)  F_tan(eV/A)  F_max(eV/A)  E(eV)\n")
        for k in range(len(dirs)):
            fh.write(f"{k:3d} {s[k]:10.5f} {rel[k]:12.6f} "
                     f"{ftan[k] if ftan[k] is not None else float('nan'):12.6f} "
                     f"{fm[k] if fm[k] is not None else float('nan'):10.5f} {E[k]:16.8f}\n")
    np.savetxt(prefix + "_spline.dat", np.column_stack([xs, ys]), fmt="%12.6f",
               header="s(A)  E-E0(eV)   " + method)
    print(f"written: {prefix}.dat, {prefix}_spline.dat", end="")
    if not a.no_plot:
        plot_profile(prefix + ".png", s, rel, xs, ys, kmax, ea_f, ea_b, de, a)
        print(f", {prefix}.png")
    else:
        print()


def plot_profile(path, s, rel, xs, ys, kmax, ea_f, ea_b, de, a):
    import matplotlib
    if not a.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    ax.plot(xs, ys, color=C1, lw=1.8)
    ax.plot(s, rel, "o", color=C1, ms=7, mec="white", mew=1.2)
    for k, (x, y) in enumerate(zip(s, rel)):
        ax.annotate(f"{k:02d}", (x, y), xytext=(0, -14), textcoords="offset points",
                    ha="center", fontsize=7, color=INK2)
    ax.axhline(0.0, color=INK2, lw=0.8, ls="--")
    ax.annotate(f"Ea = {ea_f:.2f} eV", (xs[kmax], ys[kmax]), xytext=(0, 10), textcoords="offset points",
                ha="center", fontsize=9, color=INK)
    ax.annotate(f"dE = {de:.2f} eV", (s[-1], rel[-1]), xytext=(-4, 10), textcoords="offset points",
                ha="right", fontsize=9, color=INK)
    ax.set_xlabel("reaction coordinate (A)", color=INK)
    ax.set_ylabel("E - E(IS) (eV)", color=INK)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(INK2)
    ax.tick_params(colors=INK2)
    ymin, ymax = min(ys.min(), rel.min()), max(ys.max(), rel.max())
    pad = 0.15 * (ymax - ymin if ymax > ymin else 1.0)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_title(a.title or "NEB energy profile", color=INK, fontsize=11, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=a.dpi)
    if a.show:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# xyz
# ---------------------------------------------------------------------------
def cmd_xyz(a):
    root = Path(a.dir)
    if not root.is_dir():
        die(f"'{root}' is not a directory")
    dirs = image_dirs(root)
    geoms = []
    for d in dirs:
        g, _ = load_geometry(d, prefer_contcar=not a.poscar)
        if g is None:
            die(f"no POSCAR/CONTCAR in {d}")
        geoms.append(g)
    energies = []
    for d in dirs:
        r = parse_image(d, geoms[0].natoms)
        energies.append(r["e0"])
    ref = geoms[0]
    frac_prev = ref.frac.copy()
    lines = []
    for k, (g, e) in enumerate(zip(geoms, energies)):
        # unwrap relative to the previous image so atoms move continuously
        f = frac_prev + mic_frac(g.frac - frac_prev)
        frac_prev = f
        cart = f @ g.lattice
        lat = " ".join(f"{x:.8f}" for x in g.lattice.reshape(-1))
        comment = f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T" image={k}'
        if e is not None:
            comment += f" energy={e:.8f}"
        lines.append(str(g.natoms))
        lines.append(comment)
        for sym, r in zip(g.symbols, cart):
            lines.append(f"{sym:<3s} {r[0]:14.8f} {r[1]:14.8f} {r[2]:14.8f}")
    Path(a.output).write_text("\n".join(lines) + "\n")
    print(f"wrote {len(geoms)} frames to {a.output} (extended XYZ; open with OVITO, VMD, ASE gui, ...)")


# ---------------------------------------------------------------------------
def build_parser():
    ap = argparse.ArgumentParser(
        prog="neb_tools.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("Helper for VASP/VTST NEB runs: interpolate, check, collect energies & barrier, export XYZ.\n"
                     "NEB 辅助工具：插值、检查、收集能量与能垒、导出动画。"),
        epilog=("Examples / 示例:\n"
                "  python neb_tools.py interp POSCAR_IS POSCAR_FS -n 5 --method idpp\n"
                "  python neb_tools.py check POSCAR_IS POSCAR_FS\n"
                "  python neb_tools.py check .\n"
                "  python neb_tools.py collect . -o neb\n"
                "  python neb_tools.py xyz . -o movie.xyz\n"
                "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License"))
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="cmd", metavar="COMMAND")
    sub.required = True

    p = sub.add_parser("interp", help="interpolate images 00..N+1 (linear / IDPP) / 插值生成图像")
    p.add_argument("initial", help="initial-state POSCAR/CONTCAR")
    p.add_argument("final", help="final-state POSCAR/CONTCAR")
    p.add_argument("-n", "--nimages", type=int, required=True, help="number of intermediate images (= IMAGES)")
    p.add_argument("--method", choices=("linear", "idpp"), default="linear", help="interpolation (default linear)")
    p.add_argument("-o", "--outdir", default=".", help="where to create 00..NN (default .)")
    p.add_argument("--no-mic", action="store_true", help="do not apply the minimum-image convention")
    p.add_argument("--keep-unwrapped", action="store_true", help="do not wrap coordinates into [0,1)")
    p.add_argument("--allow-cell-change", action="store_true", help="interpolate the cell as well (SS-NEB)")
    p.add_argument("--is-outcar", help="OUTCAR of the relaxed IS, copied to 00/")
    p.add_argument("--fs-outcar", help="OUTCAR of the relaxed FS, copied to NN/")
    p.add_argument("--idpp-steps", type=int, default=200, help="max IDPP optimisation steps (default 200)")
    p.add_argument("--idpp-fmax", type=float, default=0.05, help="IDPP force tolerance (default 0.05)")
    p.add_argument("--dmin", type=float, default=0.7, help="warn if atoms closer than this (A, default 0.7)")
    p.add_argument("-f", "--force", action="store_true", help="overwrite existing image directories")
    p.set_defaults(func=cmd_interp)

    p = sub.add_parser("check", help="check IS/FS pair or an image directory / 检查原子顺序与图像间距")
    p.add_argument("paths", nargs="+", help="IS FS  (two files)  or  DIR (with 00..NN)")
    p.add_argument("--move", type=float, default=0.3, help="list atoms moving more than this (A, default 0.3)")
    p.add_argument("--max-disp", type=float, default=2.0, help="displacement regarded as suspicious (A, default 2)")
    p.add_argument("--top", type=int, default=15, help="max lines in lists (default 15)")
    p.add_argument("--write-reordered", help="write the FS with the optimal atom order to this file")
    p.add_argument("--poscar", action="store_true", help="DIR mode: use POSCAR even if CONTCAR exists")
    p.add_argument("--dmin", type=float, default=0.7, help="warn if atoms closer than this (A)")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("collect", help="energies, tangent forces, spline barrier / 收集能量并拟合能垒")
    p.add_argument("dir", nargs="?", default=".", help="NEB directory with 00..NN (default .)")
    p.add_argument("--energy", choices=("e0", "toten"), default="e0", help="energy(sigma->0) (default) or TOTEN")
    p.add_argument("--poscar", action="store_true", help="use POSCAR geometries instead of CONTCAR")
    p.add_argument("-o", "--output", default="neb", help="output prefix (default neb)")
    p.add_argument("--points", type=int, default=300, help="points of the fitted curve (default 300)")
    p.add_argument("--title", help="plot title")
    p.add_argument("--no-plot", action="store_true", help="no figure")
    p.add_argument("--show", action="store_true", help="show the figure interactively")
    p.add_argument("--dpi", type=int, default=200, help="figure resolution (default 200)")
    p.set_defaults(func=cmd_collect)

    p = sub.add_parser("xyz", help="extended-XYZ movie of the band / 导出 XYZ 动画")
    p.add_argument("dir", nargs="?", default=".", help="NEB directory (default .)")
    p.add_argument("-o", "--output", default="neb_movie.xyz", help="output file (default neb_movie.xyz)")
    p.add_argument("--poscar", action="store_true", help="use POSCAR instead of CONTCAR")
    p.set_defaults(func=cmd_xyz)
    return ap


def main(argv=None) -> int:
    try:
        sys.stdout.reconfigure(errors="replace")
    except AttributeError:
        pass
    a = build_parser().parse_args(argv)
    a.func(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
