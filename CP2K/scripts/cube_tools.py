#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
cube_tools.py -- small toolbox for CP2K (Gaussian-format) cube files.
                 cube 文件工具：格式检查与信息、积分、平面平均（V_HARTREE -> 真空能级与功函数）、
                 同网格 cube 的差分（差分电荷密度）。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Sub-commands
------------
  info   CUBE ...           header, grid, cell, value range, integral
  check  CUBE ...           format check; with several files also checks that
                            they share the same grid (needed for differences)
  integrate CUBE [--range Z1 Z2 --axis c] [--abs]
                            integral = sum(values) * voxel volume; for an
                            ELECTRON_DENSITY cube this is the number of
                            (valence) electrons; positive/negative parts too
  planar CUBE [--axis c] [--fermi EV | --out CP2K.out] [-o png] [--csv]
                            plane average along a cell axis.  For a
                            V_HARTREE cube (hartree -> eV) with --fermi/--out
                            it reports the vacuum level(s) and the work
                            function  Phi = V_vac - E_F.  For a density or
                            density-difference cube it adds the charge per
                            unit length and its running integral Delta Q(z).
  diff   OUT.cube A.cube B.cube [C.cube ...] [--coeffs 1 -1 -1]
                            OUT = A - B - C ... (or any linear combination);
                            all inputs must have the same grid.

Cube conventions (as written by CP2K): lengths in bohr, two comment lines,
"natoms origin", three "n_i  voxel vector" lines, one line per atom
"Z  Z_eff  x y z", then the data with the third index running fastest.
CP2K's V_HARTREE cube is in hartree and holds the electrostatic potential
energy of an electron (same sign convention as the eigenvalues), so the
vacuum plateau lies above E_F and Phi = V_vac - E_F > 0.
Vacuum windows for the work function: by default the largest gap between
atoms along the axis is taken as vacuum, and a window of --vac-width A placed
--vac-offset A outside the outermost atom on each side is averaged.

Usage examples
--------------
    python cube_tools.py info h2o_cube-ELECTRON_DENSITY-1_0.cube
    python cube_tools.py integrate h2o_cube-ELECTRON_DENSITY-1_0.cube
    python cube_tools.py planar pt111_wf-v_hartree-1_0.cube --out pt111_workfunction.out -o vh.png --csv vh.csv
    python cube_tools.py planar pt111_wf-v_hartree-1_0.cube --fermi -5.21 --window 20 22
    python cube_tools.py check co_cu-ELECTRON_DENSITY-1_0.cube cu-ELECTRON_DENSITY-1_0.cube co-ELECTRON_DENSITY-1_0.cube
    python cube_tools.py diff drho.cube co_cu-ELECTRON_DENSITY-1_0.cube cu-ELECTRON_DENSITY-1_0.cube co-ELECTRON_DENSITY-1_0.cube
    python cube_tools.py planar drho.cube --kind density -o drho_z.png

Dependencies: python >= 3.9, numpy, matplotlib (only for figures).
"""
import argparse
import os
import re
import sys

BOHR = 0.529177210903          # Angstrom
HA2EV = 27.211386245988
NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"


def die(msg):
    sys.exit("cube_tools.py: error: " + msg)


def need_numpy():
    try:
        import numpy as np
    except ImportError:
        die("numpy is required (pip install numpy)")
    return np


# --------------------------------------------------------------------------- cube I/O
class Cube:
    def __init__(self):
        self.path = ""
        self.comments = ["", ""]
        self.origin = None      # bohr
        self.n = None           # (n1, n2, n3)
        self.vox = None         # 3x3, rows = voxel vectors (bohr)
        self.atoms = []         # (Z, Zeff, x, y, z) bohr
        self.extra = None       # MO line of Gaussian cubes with natoms < 0
        self.data = None        # (n1, n2, n3)


def _floats(line, nexp, path, what):
    toks = line.split()
    try:
        vals = [float(t.replace("D", "E")) for t in toks]
    except ValueError:
        vals = []
    if len(vals) < nexp:
        # Fortran fixed format without blanks: I5 followed by F12.6 fields
        try:
            vals = [float(line[0:5])] + [float(line[5 + 12 * i:17 + 12 * i]) for i in range(nexp - 1)]
        except ValueError:
            die(f"{path}: malformed {what} line: {line.rstrip()[:70]}")
    return vals


def _parse_values(np, text, path):
    try:
        return np.array(text.split(), dtype=float)
    except ValueError:
        pass
    vals = []
    for tok in text.split():
        try:
            vals.append(float(tok))
        except ValueError:
            m = re.fullmatch(r"([-+]?\d*\.\d+)[EeDd]?([-+]\d{2,3})", tok)
            if not m:
                die(f"{path}: cannot parse value '{tok}'")
            vals.append(float(m.group(1) + "E" + m.group(2)))
    return np.array(vals)


def read_cube(path, header_only=False):
    np = need_numpy()
    c = Cube()
    c.path = path
    try:
        fh = open(path, "r", errors="replace")
    except OSError as exc:
        die(f"cannot read {path}: {exc}")
    with fh:
        c.comments = [fh.readline().rstrip("\n"), fh.readline().rstrip("\n")]
        l3 = fh.readline()
        if not l3:
            die(f"{path}: file too short for a cube header")
        v = _floats(l3, 4, path, "atom-count/origin")
        natoms = int(v[0])
        c.origin = np.array(v[1:4])
        n, vox = [], []
        for i in range(3):
            v = _floats(fh.readline(), 4, path, f"grid axis {i + 1}")
            n.append(int(v[0]))
            vox.append(v[1:4])
        if any(k <= 0 for k in n):
            die(f"{path}: negative or zero grid counts {n} (Angstrom units are not supported)")
        c.n = tuple(n)
        c.vox = np.array(vox)
        for _ in range(abs(natoms)):
            v = _floats(fh.readline(), 5, path, "atom")
            c.atoms.append((int(v[0]), v[1], v[2], v[3], v[4]))
        if natoms < 0:
            c.extra = fh.readline().rstrip("\n")
        if header_only:
            return c
        text = fh.read()
    vals = _parse_values(np, text, path)
    ntot = n[0] * n[1] * n[2]
    if vals.size < ntot:
        die(f"{path}: {vals.size} values found, {ntot} expected ({n[0]}x{n[1]}x{n[2]}) - truncated file?")
    if vals.size > ntot:
        if natoms < 0:
            print(f"note: {path}: multi-orbital cube; only the first data set is used")
            vals = vals[:ntot] if vals.size % ntot else vals.reshape(ntot, -1)[:, 0]
        else:
            die(f"{path}: {vals.size} values found, {ntot} expected - wrong header?")
    c.data = vals.reshape(n)
    return c


def write_cube(path, c, data, comment2):
    np = need_numpy()
    with open(path, "w", newline="\n") as fh:
        fh.write("cube_tools.py (DFT-Tutorials)\n")
        fh.write(comment2.strip()[:70] + "\n")
        fh.write("%5d%12.6f%12.6f%12.6f\n" % (len(c.atoms), *c.origin))
        for i in range(3):
            fh.write("%5d%12.6f%12.6f%12.6f\n" % (c.n[i], *c.vox[i]))
        for a in c.atoms:
            fh.write("%5d%12.6f%12.6f%12.6f%12.6f\n" % a)
        flat = data.reshape(-1, c.n[2])
        for row in flat:
            for k in range(0, len(row), 6):
                fh.write("".join(" %12.5E" % x for x in row[k:k + 6]) + "\n")


# --------------------------------------------------------------------------- helpers
def voxel_volume(np, c):
    return abs(np.linalg.det(c.vox))          # bohr^3


def axis_index(axis):
    return {"a": 0, "x": 0, "b": 1, "y": 1, "c": 2, "z": 2}[axis.lower()]


def axis_geometry(np, c, ax):
    """Spacing (A) of the grid planes perpendicular to axis ax and the plane area (A^2)."""
    cellv = c.vox * np.array(c.n)[:, None] * BOHR        # rows = full cell vectors (A)
    o1, o2 = [i for i in range(3) if i != ax]
    normal = np.cross(cellv[o1], cellv[o2])
    area = np.linalg.norm(normal)
    normal /= area
    dz = abs(np.dot(c.vox[ax] * BOHR, normal))
    return dz, area, normal


def guess_kind(c):
    txt = (os.path.basename(c.path) + " " + " ".join(c.comments)).lower()
    if "hartree" in txt or "potential" in txt:
        return "potential"
    if "density" in txt or "drho" in txt or "diff" in txt:
        return "density"
    return "other"


def read_fermi_ev(path):
    ef = None
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                m = re.search(r"Fermi energy:\s*(" + NUM + r")", line)
                if m:
                    ef = float(m.group(1)) * HA2EV
    except OSError as exc:
        die(f"cannot read {path}: {exc}")
    if ef is None:
        die(f"{path}: no 'Fermi energy:' line (printed only with smearing); give --fermi in eV")
    return ef


def same_grid(np, a, b, tol=1e-5):
    probs = []
    if a.n != b.n:
        probs.append(f"grid {a.n} vs {b.n}")
    if np.abs(a.vox - b.vox).max() > tol:
        probs.append("voxel vectors differ")
    if np.abs(a.origin - b.origin).max() > tol:
        probs.append("origins differ")
    return probs


# --------------------------------------------------------------------------- commands
def cmd_info(args):
    np = need_numpy()
    for f in args.cubes:
        c = read_cube(f)
        dv = voxel_volume(np, c)
        cellv = c.vox * np.array(c.n)[:, None] * BOHR
        print(f"== {f}")
        print(f"   comments : {c.comments[0].strip()[:60]} | {c.comments[1].strip()[:60]}")
        print(f"   atoms    : {len(c.atoms)}   elements(Z): "
              + ", ".join(f"{z}x{n}" for z, n in sorted(count_z(c).items())))
        print(f"   grid     : {c.n[0]} x {c.n[1]} x {c.n[2]} = {np.prod(c.n)} points")
        for i, lab in enumerate("abc"):
            print(f"   {lab}        : {cellv[i][0]:10.4f} {cellv[i][1]:10.4f} {cellv[i][2]:10.4f} A   "
                  f"|{lab}| = {np.linalg.norm(cellv[i]):.4f} A, spacing {np.linalg.norm(c.vox[i]) * BOHR:.4f} A")
        print(f"   volume   : {dv * np.prod(c.n) * BOHR ** 3:.4f} A^3, voxel {dv:.6e} bohr^3")
        d = c.data
        print(f"   values   : min {d.min():.6e}  max {d.max():.6e}  mean {d.mean():.6e}")
        print(f"   integral : {d.sum() * dv:.6f} (value x bohr^3; electrons for a density cube)")


def count_z(c):
    out = {}
    for a in c.atoms:
        out[a[0]] = out.get(a[0], 0) + 1
    return out


def cmd_check(args):
    np = need_numpy()
    cubes, ok = [], True
    for f in args.cubes:
        c = read_cube(f)
        cubes.append(c)
        issues = []
        if not np.isfinite(c.data).all():
            issues.append("non-finite values")
        if np.abs(c.vox).max() == 0:
            issues.append("zero voxel vectors")
        if abs(np.linalg.det(c.vox)) < 1e-12:
            issues.append("voxel vectors are linearly dependent")
        print(f"{f}: {'OK' if not issues else 'PROBLEM: ' + '; '.join(issues)} "
              f"(grid {c.n[0]}x{c.n[1]}x{c.n[2]}, {len(c.atoms)} atoms)")
        ok &= not issues
    for c in cubes[1:]:
        probs = same_grid(np, cubes[0], c)
        if probs:
            ok = False
            print(f"grid mismatch {cubes[0].path} vs {c.path}: {'; '.join(probs)}")
        elif len(cubes) > 1:
            print(f"same grid as {cubes[0].path}: {c.path}")
    if len(cubes) > 1 and ok:
        print("all cubes share the same grid - differences are meaningful (atom positions must also match)")
    sys.exit(0 if ok else 1)


def cmd_integrate(args):
    np = need_numpy()
    c = read_cube(args.cube)
    dv = voxel_volume(np, c)
    d = c.data
    if args.range:
        ax = axis_index(args.axis)
        dz, _, _ = axis_geometry(np, c, ax)
        z = np.arange(c.n[ax]) * dz
        z1, z2 = args.range
        sel = (z >= z1) & (z <= z2)
        if not sel.any():
            die(f"--range {z1} {z2} contains no grid planes (axis length {c.n[ax] * dz:.3f} A)")
        d = np.take(d, np.where(sel)[0], axis=ax)
        print(f"planes {int(sel.sum())} of {c.n[ax]} along {args.axis} ({z1}-{z2} A)")
    tot = d.sum() * dv
    pos = d[d > 0].sum() * dv
    neg = d[d < 0].sum() * dv
    print(f"integral          = {tot:.6f}")
    print(f"positive part     = {pos:.6f}")
    print(f"negative part     = {neg:.6f}")
    if args.abs:
        print(f"integral of |f|   = {np.abs(d).sum() * dv:.6f}")
    print(f"(voxel volume {dv:.6e} bohr^3; for densities in e/bohr^3 the result is a number of electrons)")
    zeff = sum(a[1] for a in c.atoms)
    if zeff > 0 and not args.range:
        print(f"sum of Z_eff in the header = {zeff:.3f} (valence electrons of a neutral system)")


def cmd_planar(args):
    np = need_numpy()
    c = read_cube(args.cube)
    ax = axis_index(args.axis)
    kind = args.kind if args.kind != "auto" else guess_kind(c)
    dz, area, normal = axis_geometry(np, c, ax)
    other = tuple(i for i in range(3) if i != ax)
    nz = c.n[ax]
    z = np.arange(nz) * dz
    length = nz * dz
    avg = c.data.mean(axis=other)
    print(f"{args.cube}: axis {args.axis}, {nz} planes, spacing {dz:.4f} A, length {length:.4f} A, kind = {kind}")
    if ax == 2 and abs(normal[2]) < 0.999:
        print("note: the a-b plane normal is not along z; positions are measured along the normal")

    cols, heads = [z, avg], ["z_A", "planar_avg"]
    ef = None
    if kind == "potential":
        v = avg * HA2EV
        cols, heads = [z, v], ["z_A", "V_eV"]
        if args.fermi is not None:
            ef = args.fermi
        elif args.out:
            ef = read_fermi_ev(args.out)
        if ef is not None:
            cols.append(v - ef)
            heads.append("V_minus_EF_eV")
            print(f"E_F = {ef:.4f} eV")
        # vacuum windows
        atoms_z = np.array([np.dot(np.array(a[2:5]) * BOHR - c.origin * BOHR, normal) for a in c.atoms]) % length \
            if c.atoms else np.array([])
        windows = []
        if args.window:
            windows.append(("window", args.window[0], args.window[1]))
        elif atoms_z.size:
            za = np.sort(atoms_z)
            gaps = np.diff(np.append(za, za[0] + length))
            i = int(np.argmax(gaps))
            lo, gap = za[i], gaps[i]
            print(f"largest atom-free gap: z = {lo:.2f} -> {(lo + gap) % length:.2f} A ({gap:.2f} A thick)")
            need = 2 * (args.vac_offset + args.vac_width)
            if gap < need:
                print(f"warning: vacuum gap {gap:.2f} A < 2*(offset+width) = {need:.2f} A; windows overlap or "
                      "reach the slab - use --window or a thicker vacuum")
            windows.append(("above slab", lo + args.vac_offset, lo + args.vac_offset + args.vac_width))
            windows.append(("below slab", lo + gap - args.vac_offset - args.vac_width, lo + gap - args.vac_offset))
        for name, a, b in windows:
            sel = ((z - a) % length) <= (b - a)
            if not sel.any():
                print(f"{name}: window {a:.2f}-{b:.2f} A holds no grid plane")
                continue
            vv = v[sel].mean()
            spread = v[sel].max() - v[sel].min()
            txt = f"{name:11s}: z = {a % length:6.2f}-{b % length:6.2f} A  V_vac = {vv:9.4f} eV  (spread {spread:.3f} eV)"
            if ef is not None:
                txt += f"  Phi = V_vac - E_F = {vv - ef:.3f} eV"
            print(txt)
        if ef is None:
            print("give --fermi EV or --out CP2K.out to get the work function")
        else:
            print("check: plateaus must be flat (spread ~0.01 eV); asymmetric slabs need SURFACE_DIPOLE_CORRECTION")
    elif kind == "density":
        dv = voxel_volume(np, c)
        q_plane = c.data.sum(axis=other) * dv            # electrons per grid plane
        lam = q_plane / dz                               # e / A
        cum = np.cumsum(q_plane)
        cols += [lam, cum]
        heads += ["charge_per_length_e_per_A", "cumulative_charge_e"]
        print(f"total integral = {cum[-1]:.6f} e;  max |Delta Q(z)| = {np.abs(cum).max():.6f} e at "
              f"z = {z[int(np.argmax(np.abs(cum)))]:.2f} A (meaningful for difference densities)")
    if args.csv:
        np.savetxt(args.csv, np.column_stack(cols), delimiter=",", header=",".join(heads), comments="", fmt="%.6e")
        print(f"wrote {args.csv}")
    if args.output.lower() == "none" and not args.show:
        return
    try:
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        die("matplotlib is required for plotting (pip install matplotlib); use -o none --csv FILE")
    fig, ax_ = plt.subplots(figsize=(6.0, 3.6))
    if kind == "potential":
        y = cols[2] if ef is not None else cols[1]
        ax_.plot(z, y, color="#2a78d6", lw=1.5)
        ax_.set_ylabel(r"$\bar V(z) - E_\mathrm{F}$ (eV)" if ef is not None else r"$\bar V(z)$ (eV)")
        if ef is not None:
            ax_.axhline(0.0, color="#8a8a86", lw=0.8, ls="--")
    elif kind == "density":
        ax_.plot(z, cols[2], color="#2a78d6", lw=1.5, label=r"$\Delta\lambda(z)$ (e/$\mathrm{\AA}$)")
        ax_.plot(z, cols[3], color="#eb6834", lw=1.5, label=r"$\Delta Q(z)$ (e)")
        ax_.axhline(0.0, color="#8a8a86", lw=0.8)
        ax_.legend(frameon=False, fontsize=8)
        ax_.set_ylabel("charge")
    else:
        ax_.plot(z, avg, color="#2a78d6", lw=1.5)
        ax_.set_ylabel("planar average (cube units)")
    for a in c.atoms:
        za = np.dot(np.array(a[2:5]) * BOHR - c.origin * BOHR, normal) % length
        ax_.axvline(za, color="#d9d8d2", lw=0.5, zorder=0)
    ax_.set_xlabel(f"position along {args.axis} " + r"($\mathrm{\AA}$)")
    ax_.set_xlim(0, length)
    ax_.grid(alpha=0.25, lw=0.5)
    for sp in ("top", "right"):
        ax_.spines[sp].set_visible(False)
    fig.tight_layout()
    if args.output.lower() != "none":
        fig.savefig(args.output, dpi=args.dpi)
        print(f"wrote {args.output}")
    if args.show:
        plt.show()


def cmd_diff(args):
    np = need_numpy()
    if len(args.inputs) < 2:
        die("diff needs at least two input cubes")
    coeffs = args.coeffs or [1.0] + [-1.0] * (len(args.inputs) - 1)
    if len(coeffs) != len(args.inputs):
        die(f"--coeffs: {len(coeffs)} values for {len(args.inputs)} cubes")
    if os.path.abspath(args.out) in [os.path.abspath(p) for p in args.inputs]:
        die("output file must differ from the inputs")
    ref = read_cube(args.inputs[0])
    res = coeffs[0] * ref.data
    for f, w in zip(args.inputs[1:], coeffs[1:]):
        c = read_cube(f)
        probs = same_grid(np, ref, c)
        if probs:
            die(f"{f} vs {ref.path}: {'; '.join(probs)} - all cubes need the same cell, CUTOFF and STRIDE")
        if len(c.atoms) and len(ref.atoms):
            pos_ref = {tuple(np.round(a[2:5], 3)) for a in ref.atoms}
            missing = sum(1 for a in c.atoms if tuple(np.round(a[2:5], 3)) not in pos_ref)
            if missing:
                print(f"warning: {missing} atom(s) of {f} are not at a position of {ref.path} "
                      "(fragments should be cut from the combined structure without relaxation)")
        res = res + w * c.data
    expr = " ".join(f"{'+' if w >= 0 else '-'} {abs(w):g}*{os.path.basename(p)}" for w, p in zip(coeffs, args.inputs))
    write_cube(args.out, ref, res, "linear combination: " + expr)
    dv = voxel_volume(np, ref)
    print(f"wrote {args.out} = {expr}")
    print(f"integral of the result = {res.sum() * dv:.6f}  (should be ~0 for a charge-density difference)")
    print(f"min {res.min():.4e}  max {res.max():.4e}; typical isosurface for Delta rho: +-0.001 to 0.005 e/bohr^3")


# --------------------------------------------------------------------------- CLI
def build_parser():
    ap = argparse.ArgumentParser(
        prog="cube_tools.py",
        description="Tools for CP2K cube files: info, check, integrate, planar (work function), diff. "
                    "cube 工具：信息、格式检查、积分、平面平均与功函数、差分。",
        epilog="Run 'python cube_tools.py <command> -h' for details.\nExamples:\n"
               "  python cube_tools.py integrate h2o_cube-ELECTRON_DENSITY-1_0.cube\n"
               "  python cube_tools.py planar pt111_wf-v_hartree-1_0.cube --out pt111_workfunction.out -o vh.png\n"
               "  python cube_tools.py diff drho.cube AB.cube A.cube B.cube\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("info", help="print header information and value statistics")
    p.add_argument("cubes", nargs="+", help="cube files")
    p.set_defaults(func=cmd_info)

    p = sub.add_parser("check", help="format check; several files: same-grid check (exit 1 on problems)")
    p.add_argument("cubes", nargs="+", help="cube files")
    p.set_defaults(func=cmd_check)

    p = sub.add_parser("integrate", help="integral of the data (electrons for a density cube)")
    p.add_argument("cube", help="cube file")
    p.add_argument("--axis", default="c", choices=["a", "b", "c", "x", "y", "z"],
                   help="axis for --range (default c)")
    p.add_argument("--range", nargs=2, type=float, metavar=("Z1", "Z2"),
                   help="only integrate the planes between Z1 and Z2 (A, measured from the cube origin)")
    p.add_argument("--abs", action="store_true", help="also print the integral of |f|")
    p.set_defaults(func=cmd_integrate)

    p = sub.add_parser("planar", help="plane average along an axis; work function from V_HARTREE")
    p.add_argument("cube", help="cube file (v_hartree, ELECTRON_DENSITY, difference density ...)")
    p.add_argument("--axis", default="c", choices=["a", "b", "c", "x", "y", "z"],
                   help="average over the planes spanned by the other two axes (default c)")
    p.add_argument("--kind", choices=["auto", "potential", "density", "other"], default="auto",
                   help="data type (auto: from the file name/comments)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--fermi", type=float, metavar="EV", help="Fermi energy in eV")
    g.add_argument("--out", metavar="CP2K.out", help="read the last 'Fermi energy:' (hartree) from this output")
    p.add_argument("--window", nargs=2, type=float, metavar=("Z1", "Z2"),
                   help="explicit vacuum window in A (instead of the automatic two windows)")
    p.add_argument("--vac-offset", type=float, default=4.0,
                   help="distance of the automatic windows from the outermost atoms, A (default 4)")
    p.add_argument("--vac-width", type=float, default=2.0, help="width of the automatic windows, A (default 2)")
    p.add_argument("--csv", metavar="FILE", help="write the profile to CSV")
    p.add_argument("-o", "--output", default="planar.png", help="figure file (default planar.png; 'none' = no figure)")
    p.add_argument("--dpi", type=int, default=200, help="figure resolution (default 200)")
    p.add_argument("--show", action="store_true", help="open an interactive window")
    p.set_defaults(func=cmd_planar)

    p = sub.add_parser("diff", help="OUT = A - B - C ... (same grid required)")
    p.add_argument("out", help="output cube file")
    p.add_argument("inputs", nargs="+", help="input cubes A B [C ...]")
    p.add_argument("--coeffs", nargs="+", type=float, metavar="W",
                   help="coefficients of the inputs (default 1 -1 -1 ...)")
    p.set_defaults(func=cmd_diff)
    return ap


def main():
    args = build_parser().parse_args()
    for attr in ("cubes", "inputs"):
        for f in getattr(args, attr, []) or []:
            if not os.path.isfile(f):
                die(f"file not found: {f}")
    if getattr(args, "cube", None) and not os.path.isfile(args.cube):
        die(f"file not found: {args.cube}")
    args.func(args)


if __name__ == "__main__":
    main()
