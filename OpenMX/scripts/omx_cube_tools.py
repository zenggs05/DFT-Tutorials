#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_cube_tools.py -- Gaussian-cube utilities for OpenMX output (tden/sden/dden/v0/vhart, MO cubes).
                     平面平均势与功函数（替代 gcube2oned）、电荷密度差 rho(AB)-rho(A)-rho(B)
                     （替代 diff_gcube/add_gcube）、积分电荷、|psi|^2、Bader 前处理。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Sub-commands
------------
  info       header, grid, cell, min/max and integral of a cube
  integrate  integral of the data (electrons for tden, muB for sden), optionally
             only for a slab of fractional coordinates along one axis
  avg        planar average along a / b / c (optionally a macroscopic average)
  workfunc   work function Phi = V_vac - mu from System.Name.v0.cube (or v1.cube)
             with mu read from System.Name.out ("Chemical potential (Hartree)")
  diff       rho(AB) - rho(A) - rho(B) - ...  -> new cube + planar-averaged
             difference and the charge-displacement curve dQ(z)
  combine    any linear combination  c1*cube1 + c2*cube2 + ...  (e.g. up + down)
  psi2       |psi|^2 = re^2 + im^2 from the _r / _i MO cubes of a periodic run
  bader-prep write a cube for Henkelman's "bader" program: optional sum of
             several cubes, empty atoms (Z = 0) removed from the atom list,
             small negative densities optionally clipped

OpenMX cube conventions (checked against OpenMX 4.0 OutData.c)
---------------------------------------------------------------
* lengths in Bohr; densities in electrons/Bohr^3; potentials (v0, v1, vhart,
  vxc) in Hartree; data loop a (outer) -> b -> c (inner), 6 values per line.
* values are written with "%13.3E", i.e. only 4 significant digits; planar
  averages of v0 are therefore accurate to roughly 1e-4 relative - fine for
  work functions at the 0.01 eV level, but do not expect more.
* atom lines contain the atomic number (0 for empty atoms) and the net
  Mulliken-type charge instead of the nuclear charge.
* for work functions use v0.cube / v1.cube (Kohn-Sham potential without the
  non-local part); vhart.cube is the Hartree potential of the DIFFERENCE
  density and is not the total electrostatic potential (OpenMX manual Ch. 50).

The averaging coordinate is the distance from the cube origin measured
perpendicular to the plane of the other two grid vectors (identical to the
axis length for orthogonal cells).

Usage examples
--------------
  python omx_cube_tools.py info h2o.tden.cube
  python omx_cube_tools.py integrate h2o.tden.cube
  python omx_cube_tools.py integrate fe.sden.cube --axis c --range 0.0 0.5
  python omx_cube_tools.py avg al111.v0.cube --axis a --unit eV -o al111_v0_avg.dat --plot avg.png
  python omx_cube_tools.py workfunc al111_esm_on1_2e.v0.cube --axis a --outfile al111_esm_on1_2e.out --plot wf.png
  python omx_cube_tools.py workfunc slab.v0.cube --axis c --mu -0.1523 --left 2 6 --right 40 44
  python omx_cube_tools.py diff co_cu_AB.tden.cube co_cu_A.tden.cube co_cu_B.tden.cube -o cdd.cube --axis c --plot cdd.png
  python omx_cube_tools.py combine x.up.cube x.dn.cube --coef 1 1 -o x.tot.cube
  python omx_cube_tools.py psi2 si.lumo1_0_0_r.cube si.lumo1_0_0_i.cube -o lumo_k1.cube
  python omx_cube_tools.py bader-prep h2o.tden.cube -o h2o_bader.cube

Dependencies: Python >= 3.9, numpy; matplotlib only for --plot.
"""

import argparse
import os
import re
import sys
import warnings

BOHR2ANG = 0.529177210903
HA2EV = 27.211386245988


def die(msg):
    sys.exit("omx_cube_tools.py: error: " + msg)


# ----------------------------------------------------------------------------
# cube I/O
# ----------------------------------------------------------------------------
class Cube:
    def __init__(self, comments, origin, n, vox, atoms, data, extra=""):
        self.comments = comments      # two strings
        self.origin = origin          # (3,) Bohr
        self.n = n                    # (3,) ints
        self.vox = vox                # (3,3) Bohr, rows = step vectors
        self.atoms = atoms            # list of (Z, charge, x, y, z)
        self.data = data              # ndarray (n1, n2, n3)
        self.extra = extra            # MO index line for negative natoms (Gaussian style)

    @property
    def dv(self):
        import numpy as np
        return abs(np.linalg.det(self.vox))

    def heights(self, axis):
        """Perpendicular distance (Bohr) of each grid plane along `axis` from the origin."""
        import numpy as np
        o = [i for i in range(3) if i != axis]
        nrm = np.cross(self.vox[o[0]], self.vox[o[1]])
        nrm /= np.linalg.norm(nrm)
        step = abs(np.dot(self.vox[axis], nrm))
        return np.arange(self.n[axis]) * step, step


def read_cube(fname):
    import numpy as np

    try:
        with open(fname) as f:
            c1, c2 = f.readline().rstrip("\n"), f.readline().rstrip("\n")
            head = [f.readline().split() for _ in range(4)]
            try:
                nat = int(head[0][0])
                origin = np.array(head[0][1:4], float)
                n = [int(h[0]) for h in head[1:]]
                vox = np.array([h[1:4] for h in head[1:]], float)
            except (IndexError, ValueError):
                die("%s: malformed cube header (lines 3-6)" % fname)
            if any(v <= 0 for v in n):
                die("%s: negative grid counts (Angstrom-unit cubes) are not supported" % fname)
            atoms = []
            for _ in range(abs(nat)):
                p = f.readline().split()
                try:
                    atoms.append((int(float(p[0])), float(p[1]), float(p[2]), float(p[3]), float(p[4])))
                except (IndexError, ValueError):
                    die("%s: malformed atom line" % fname)
            extra = f.readline() if nat < 0 else ""
            text = f.read()
    except OSError as exc:
        die("cannot read %s (%s)" % (fname, exc))
    npts = n[0] * n[1] * n[2]
    data = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            data = np.fromstring(text, dtype=float, sep=" ")   # fast path (text mode)
        except (AttributeError, ValueError):
            data = None
    if data is None or data.size != npts:
        try:
            data = np.array(text.split(), dtype=float)
        except ValueError as exc:
            die("%s: non-numeric volumetric data (%s)" % (fname, exc))
    if data.size != npts:
        die("%s: expected %d values (%d x %d x %d), found %d (truncated file? multi-MO cube?)"
            % (fname, npts, n[0], n[1], n[2], data.size))
    return Cube([c1, c2], origin, n, vox, atoms, data.reshape(n), extra)


def write_cube(fname, cube, comment1, comment2="", atoms=None):
    atoms = cube.atoms if atoms is None else atoms
    with open(fname, "w") as f:
        f.write(" %s\n %s\n" % (comment1[:200], comment2[:200]))
        f.write("%5d%12.6f%12.6f%12.6f\n" % ((len(atoms),) + tuple(cube.origin)))
        for i in range(3):
            f.write("%5d%12.6f%12.6f%12.6f\n" % ((cube.n[i],) + tuple(cube.vox[i])))
        for z, q, x, y, zz in atoms:
            f.write("%5d%12.6f%12.6f%12.6f%12.6f\n" % (z, q, x, y, zz))
        n3 = cube.n[2]
        rows = cube.data.reshape(-1, n3)
        for row in rows:
            for s in range(0, n3, 6):
                f.write("".join("%13.5E" % v for v in row[s:s + 6]) + "\n")
    print("cube written to %s" % fname)


def same_grid(a, b, name_a, name_b):
    import numpy as np
    if list(a.n) != list(b.n):
        die("grid of %s %s differs from %s %s (use the same cell, cutoff/scf.Ngrid)"
            % (name_b, b.n, name_a, a.n))
    if not np.allclose(a.vox, b.vox, atol=1e-5):
        die("voxel vectors of %s differ from %s" % (name_b, name_a))
    if not np.allclose(a.origin, b.origin, atol=1e-4):
        die("grid origin of %s differs from %s: put scf.fixed.grid (Bohr, copied from the "
            "AB run) into all inputs and rerun" % (name_b, name_a))


def axis_index(s):
    m = {"a": 0, "b": 1, "c": 2, "1": 0, "2": 1, "3": 2, "x": 0, "y": 1, "z": 2}
    if s.lower() not in m:
        die("axis must be a/b/c (or 1/2/3)")
    return m[s.lower()]


def planar_average(cube, axis):
    o = tuple(i for i in range(3) if i != axis)
    return cube.data.mean(axis=o)


def macro_average(z, v, width):
    """Running average over a window `width` (same unit as z), periodic."""
    import numpy as np
    if width <= 0:
        return v
    dz = z[1] - z[0]
    nw = max(1, int(round(width / dz)))
    ext = np.concatenate([v, v, v])
    ker = np.ones(nw) / nw
    sm = np.convolve(ext, ker, mode="same")
    return sm[len(v):2 * len(v)]


def read_mu(outfile):
    pat = re.compile(r"Chemical\s+potential\s+\(Hartree\)\s+([-+]?\d+\.\d+(?:[eE][-+]?\d+)?)")
    mu = None
    try:
        with open(outfile, errors="replace") as f:
            for line in f:
                m = pat.search(line)
                if m and "Leads" not in line:
                    mu = float(m.group(1))
    except OSError as exc:
        die("cannot read %s (%s)" % (outfile, exc))
    if mu is None:
        die("'Chemical potential (Hartree)' not found in %s" % outfile)
    return mu


def make_plot(args_plot, show, xs, ys, xlabel, ylabel, hlines=(), spans=(), labels=None):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    labels = labels or [None] * len(ys)
    for y, lab in zip(ys, labels):
        ax.plot(xs, y, lw=1.1, label=lab)
    for h, style in hlines:
        ax.axhline(h, **style)
    for lo, hi in spans:
        ax.axvspan(lo, hi, color="0.9", lw=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xs[0], xs[-1])
    if any(labels):
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    if show:
        plt.show()
    else:
        fig.savefig(args_plot, dpi=300)
        print("figure written to %s" % args_plot)


# ----------------------------------------------------------------------------
# sub-commands
# ----------------------------------------------------------------------------
def cmd_info(a):
    import numpy as np
    c = read_cube(a.cube)
    cell = c.vox * np.array(c.n)[:, None] * BOHR2ANG
    print("file      : %s" % a.cube)
    print("comments  : %s | %s" % (c.comments[0].strip(), c.comments[1].strip()))
    print("atoms     : %d (%d empty atoms with Z = 0)" % (len(c.atoms), sum(1 for t in c.atoms if t[0] == 0)))
    print("grid      : %d x %d x %d   voxel volume %.6e Bohr^3" % (c.n[0], c.n[1], c.n[2], c.dv))
    print("origin    : %.6f %.6f %.6f Bohr" % tuple(c.origin))
    for i, lab in enumerate("abc"):
        print("cell %s    : %10.5f %10.5f %10.5f Ang  (|%s| = %.5f Ang)"
              % ((lab,) + tuple(cell[i]) + (lab, np.linalg.norm(cell[i]))))
    print("data      : min %.6e  max %.6e  mean %.6e" % (c.data.min(), c.data.max(), c.data.mean()))
    print("integral  : %.6f  (sum x dV; electrons for tden, muB for sden)" % (c.data.sum() * c.dv))


def cmd_integrate(a):
    import numpy as np
    c = read_cube(a.cube)
    d = np.abs(c.data) if a.abs else c.data
    if a.range:
        ax = axis_index(a.axis)
        lo, hi = a.range
        frac = np.arange(c.n[ax]) / c.n[ax]
        sel = (frac >= lo) & (frac < hi)
        if not sel.any():
            die("--range %.3f %.3f selects no grid plane" % (lo, hi))
        idx = [slice(None)] * 3
        idx[ax] = sel
        val = d[tuple(idx)].sum() * c.dv
        print("integral over %s in [%.3f, %.3f) (%d planes): %.6f" % (a.axis, lo, hi, sel.sum(), val))
    else:
        print("integral: %.6f" % (d.sum() * c.dv))
    pos, neg = c.data[c.data > 0].sum() * c.dv, c.data[c.data < 0].sum() * c.dv
    print("positive part %.6f   negative part %.6f" % (pos, neg))


def cmd_avg(a):
    import numpy as np
    c = read_cube(a.cube)
    ax = axis_index(a.axis)
    z, _ = c.heights(ax)
    v = planar_average(c, ax)
    scale = {"none": 1.0, "eV": HA2EV}[a.unit]
    v = v * scale
    z_ang = z * BOHR2ANG
    vm = macro_average(z_ang, v, a.macro) if a.macro else None
    cols, head = [z_ang, v], "z_Ang avg"
    if vm is not None:
        cols.append(vm)
        head += " macro_avg"
    out = a.out or os.path.splitext(os.path.basename(a.cube))[0] + "_avg_%s.dat" % "abc"[ax]
    np.savetxt(out, np.column_stack(cols), fmt="%14.6f", header=head + "  (unit: %s)" % a.unit)
    print("planar average along %s (%d points, length %.4f Ang) written to %s"
          % ("abc"[ax], len(z), len(z) * (z[1] - z[0]) * BOHR2ANG if len(z) > 1 else 0, out))
    print("min %.6f  max %.6f  (%s)" % (v.min(), v.max(), a.unit))
    if a.plot or a.show:
        ys, labs = [v], ["planar average"]
        if vm is not None:
            ys.append(vm)
            labs.append("macroscopic (%.2f Ang)" % a.macro)
        make_plot(a.plot or "avg.png", a.show, z_ang, ys, "distance along %s (Ang)" % "abc"[ax],
                  "average (%s)" % a.unit, labels=labs)


def cmd_workfunc(a):
    import numpy as np
    c = read_cube(a.cube)
    if "vhart" in os.path.basename(a.cube).lower():
        print("WARNING: vhart.cube is the Hartree potential of the difference density; "
              "use System.Name.v0.cube (or v1.cube) for work functions")
    if a.mu is None and not a.outfile:
        guess = re.sub(r"\.v[01]\.cube$", ".out", a.cube)
        if guess != a.cube and os.path.exists(guess):
            a.outfile = guess
        else:
            die("give --mu (Hartree) or --outfile System.Name.out")
    mu = a.mu if a.mu is not None else read_mu(a.outfile)
    ax = axis_index(a.axis)
    z, _ = c.heights(ax)
    zz = z * BOHR2ANG
    v = planar_average(c, ax) * HA2EV
    length = len(zz) * (zz[1] - zz[0])
    left = a.left or (0.05 * length, 0.15 * length)
    right = a.right or (0.85 * length, 0.95 * length)
    print("# cell length along %s : %10.4f Ang" % ("abc"[ax], length))
    print("# chemical potential mu : %.8f Hartree = %.4f eV%s"
          % (mu, mu * HA2EV, " (from %s)" % a.outfile if a.outfile and a.mu is None else ""))
    res = []
    for name, (lo, hi) in (("left ", left), ("right", right)):
        sel = (zz >= lo) & (zz <= hi)
        if sel.sum() < 2:
            die("vacuum window %s [%.2f, %.2f] Ang contains < 2 grid points" % (name, lo, hi))
        vv = v[sel]
        phi = vv.mean() - mu * HA2EV
        res.append((name, lo, hi, vv.mean(), vv.max() - vv.min(), phi))
        print("%s vacuum [%7.2f, %7.2f] Ang: V = %9.4f eV (flatness %.4f eV) -> Phi = %7.3f eV"
              % (name, lo, hi, vv.mean(), vv.max() - vv.min(), phi))
        if vv.max() - vv.min() > 0.05:
            print("  WARNING: potential not flat in this window; choose a window in the vacuum "
                  "plateau (--left/--right) or check the vacuum thickness")
    print("vacuum-level difference V_right - V_left = %.3f eV (non-zero for asymmetric/polar slabs)"
          % (res[1][3] - res[0][3]))
    if a.dat:
        np.savetxt(a.dat, np.column_stack([zz, v, v - mu * HA2EV]), fmt="%14.6f",
                   header="z_Ang V_eV V-mu_eV")
        print("profile written to %s" % a.dat)
    if a.plot or a.show:
        make_plot(a.plot or "workfunc.png", a.show, zz, [v - mu * HA2EV],
                  "distance along %s (Ang)" % "abc"[ax], r"$\bar V - \mu$ (eV)",
                  hlines=[(0.0, dict(color="0.4", ls="--", lw=0.8))], spans=[left, right])


def cmd_diff(a):
    import numpy as np
    ab = read_cube(a.cubes[0])
    d = ab.data.copy()
    for fn in a.cubes[1:]:
        c = read_cube(fn)
        same_grid(ab, c, a.cubes[0], fn)
        d -= c.data
    res = Cube(ab.comments, ab.origin, ab.n, ab.vox, ab.atoms, d)
    print("integral of difference: %.6f e  (positive %.6f, negative %.6f)"
          % (d.sum() * ab.dv, d[d > 0].sum() * ab.dv, d[d < 0].sum() * ab.dv))
    if a.out:
        write_cube(a.out, res, "charge density difference %s" % " - ".join(os.path.basename(x) for x in a.cubes),
                   "electrons/Bohr^3; atoms of the first cube")
    if not a.out and not a.axis:
        print("hint: add -o cdd.cube to save the difference and/or --axis c for the planar profile")
    if a.axis:
        ax = axis_index(a.axis)
        z, step = ab.heights(ax)
        o = tuple(i for i in range(3) if i != ax)
        lam = d.sum(axis=o) * ab.dv / step          # electrons per Bohr along the axis
        dq = np.cumsum(lam) * step                  # charge displacement curve
        lam_a = lam / BOHR2ANG                      # electrons per Angstrom
        out = a.dat or os.path.splitext(os.path.basename(a.out or "cdd.cube"))[0] + "_planar_%s.dat" % "abc"[ax]
        np.savetxt(out, np.column_stack([z * BOHR2ANG, lam_a, dq]), fmt="%14.6e",
                   header="z_Ang drho_planar_e_per_Ang dQ_cumulative_e")
        i = int(np.argmax(np.abs(dq)))
        print("planar-integrated difference along %s written to %s; max |dQ(z)| = %.4f e at z = %.3f Ang"
              % ("abc"[ax], out, dq[i], z[i] * BOHR2ANG))
        if a.plot or a.show:
            make_plot(a.plot or "cdd_planar.png", a.show, z * BOHR2ANG, [lam_a, dq],
                      "distance along %s (Ang)" % "abc"[ax], r"$\Delta\rho$ (e/Ang), $\Delta Q$ (e)",
                      hlines=[(0.0, dict(color="0.4", lw=0.6))],
                      labels=[r"$\Delta\rho(z)$ plane-integrated", r"$\Delta Q(z)$"])


def cmd_combine(a):
    if len(a.coef) != len(a.cubes):
        die("--coef needs one number per cube (%d cubes, %d coefficients)" % (len(a.cubes), len(a.coef)))
    first = read_cube(a.cubes[0])
    d = first.data * a.coef[0]
    for fn, cf in zip(a.cubes[1:], a.coef[1:]):
        c = read_cube(fn)
        same_grid(first, c, a.cubes[0], fn)
        d += cf * c.data
    res = Cube(first.comments, first.origin, first.n, first.vox, first.atoms, d)
    print("integral of the combination: %.6f" % (d.sum() * first.dv))
    write_cube(a.out, res, "linear combination " + " ".join("%+g*%s" % (cf, os.path.basename(fn))
                                                           for fn, cf in zip(a.cubes, a.coef)))


def cmd_psi2(a):
    re_c = read_cube(a.real)
    d = re_c.data ** 2
    if a.imag:
        im_c = read_cube(a.imag)
        same_grid(re_c, im_c, a.real, a.imag)
        d = d + im_c.data ** 2
    res = Cube(re_c.comments, re_c.origin, re_c.n, re_c.vox, re_c.atoms, d)
    print("norm of |psi|^2 in the cube: %.6f" % (d.sum() * re_c.dv))
    write_cube(a.out, res, "|psi|^2 from %s%s" % (os.path.basename(a.real),
                                                 " + " + os.path.basename(a.imag) if a.imag else ""),
               re_c.comments[0].strip())


def cmd_bader(a):
    import numpy as np
    first = read_cube(a.cubes[0])
    d = first.data.copy()
    for fn in a.cubes[1:]:
        c = read_cube(fn)
        same_grid(first, c, a.cubes[0], fn)
        d += c.data
    nneg = int((d < 0).sum())
    if nneg:
        print("%d grid points with negative density (min %.3e)%s"
              % (nneg, d.min(), "; clipped to 0" if a.clip_negative else ""))
        if a.clip_negative:
            d = np.clip(d, 0.0, None)
    atoms = first.atoms if a.keep_empty else [t for t in first.atoms if t[0] != 0]
    if len(atoms) != len(first.atoms):
        print("removed %d empty atoms (Z = 0) from the atom list" % (len(first.atoms) - len(atoms)))
    # bader uses the second column only for output; write the atomic number there
    atoms = [(t[0], float(t[0]), t[2], t[3], t[4]) for t in atoms]
    res = Cube(first.comments, first.origin, first.n, first.vox, atoms, d)
    print("total electrons in the cube: %.4f" % (d.sum() * first.dv))
    write_cube(a.out, res, "density for bader from %s" % " + ".join(os.path.basename(x) for x in a.cubes),
               "electrons/Bohr^3")
    print("next: bader %s   (charges in ACF.dat are VALENCE electrons of the OpenMX VPS; "
          "net charge = valence electrons - CHARGE column)" % a.out)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        prog="omx_cube_tools.py",
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", metavar="COMMAND")
    sub.required = True

    p = sub.add_parser("info", help="print header and statistics")
    p.add_argument("cube")
    p.set_defaults(func=cmd_info)

    p = sub.add_parser("integrate", help="integrate the volumetric data")
    p.add_argument("cube")
    p.add_argument("--axis", default="c", help="axis for --range (a/b/c, default c)")
    p.add_argument("--range", type=float, nargs=2, metavar=("LO", "HI"),
                   help="fractional range [LO, HI) along --axis")
    p.add_argument("--abs", action="store_true", help="integrate |data| instead of data")
    p.set_defaults(func=cmd_integrate)

    p = sub.add_parser("avg", help="planar average along one axis")
    p.add_argument("cube")
    p.add_argument("--axis", default="c", help="a/b/c (or 1/2/3), default c")
    p.add_argument("--unit", choices=["none", "eV"], default="none",
                   help="eV multiplies Hartree potentials by 27.2114 (default none)")
    p.add_argument("--macro", type=float, default=0.0, help="macroscopic-average window in Ang (0 = off)")
    p.add_argument("-o", "--out", default="", help="output .dat (default <cube>_avg_<axis>.dat)")
    p.add_argument("--plot", default="", help="figure file")
    p.add_argument("--show", action="store_true")
    p.set_defaults(func=cmd_avg)

    p = sub.add_parser("workfunc", help="work function from v0.cube and mu")
    p.add_argument("cube", help="System.Name.v0.cube (or v1.cube)")
    p.add_argument("--axis", default="c", help="surface normal a/b/c (default c)")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--mu", type=float, default=None, help="chemical potential in Hartree")
    g.add_argument("--outfile", default="", help="System.Name.out to read the chemical potential from")
    p.add_argument("--left", type=float, nargs=2, metavar=("LO", "HI"),
                   help="left vacuum window in Ang from the cube origin (default 5-15%% of the cell)")
    p.add_argument("--right", type=float, nargs=2, metavar=("LO", "HI"),
                   help="right vacuum window in Ang (default 85-95%% of the cell)")
    p.add_argument("--dat", default="", help="write z, V, V-mu (eV) to this file")
    p.add_argument("--plot", default="", help="figure file")
    p.add_argument("--show", action="store_true")
    p.set_defaults(func=cmd_workfunc)

    p = sub.add_parser("diff", help="rho(AB) - rho(A) - rho(B) ...")
    p.add_argument("cubes", nargs="+", help="AB cube first, then the cubes to subtract")
    p.add_argument("-o", "--out", default="", help="write the difference cube")
    p.add_argument("--axis", default="", help="also write plane-integrated drho(z) and dQ(z) along a/b/c")
    p.add_argument("--dat", default="", help="file for the planar data (with --axis)")
    p.add_argument("--plot", default="", help="figure of drho(z) and dQ(z) (with --axis)")
    p.add_argument("--show", action="store_true")
    p.set_defaults(func=cmd_diff)

    p = sub.add_parser("combine", help="linear combination of cubes")
    p.add_argument("cubes", nargs="+")
    p.add_argument("--coef", type=float, nargs="+", required=True, help="one coefficient per cube")
    p.add_argument("-o", "--out", required=True)
    p.set_defaults(func=cmd_combine)

    p = sub.add_parser("psi2", help="|psi|^2 from real (+ imaginary) MO cubes")
    p.add_argument("real", help="..._r.cube (or a real MO cube)")
    p.add_argument("imag", nargs="?", default="", help="..._i.cube (optional)")
    p.add_argument("-o", "--out", required=True)
    p.set_defaults(func=cmd_psi2)

    p = sub.add_parser("bader-prep", help="prepare a density cube for the bader program")
    p.add_argument("cubes", nargs="+", help="density cube(s); several are summed")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--keep-empty", action="store_true", help="keep empty atoms (Z = 0) in the atom list")
    p.add_argument("--clip-negative", action="store_true", help="set negative densities to zero")
    p.set_defaults(func=cmd_bader)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
