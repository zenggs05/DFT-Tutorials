#!/usr/bin/env python3
"""Free-energy profile/surface from a CP2K metadynamics HILLS file (1 or 2 CVs).

CP2K writes <PROJECT>-HILLS.metadynLog with one line per Gaussian:
    time[fs]  s_1..s_n  width_1..width_n  height[hartree]
(CV values and widths in CP2K internal units: bohr for distances, radian for
angles/torsions, dimensionless for coordination numbers).

Free energy estimate:
    F(s) = - sum_j h_j * prod_k exp(-0.5 * ((s_k - s_jk)/w_jk)^2)  + const
For well-tempered runs CP2K stores the heights already multiplied by
(T + DELTA_T)/DELTA_T (see src/metadynamics.F), so the same formula applies to
standard and well-tempered metadynamics.  Walls are not included in HILLS.
A width of 0 means "infinitely wide" in that direction (CP2K convention).

Precision: CP2K writes the HILLS file with format f13.5, i.e. heights in
hartree keep only 5 decimals (4.8E-4 Ha is written as 0.00048, late WT hills of
~3E-5 Ha lose ~15 %).  For quantitative work pass the CP2K restart file
(<PROJECT>-1.restart) instead: its &SPAWNED_HILLS_POS/_SCALE/_HEIGHT sections
hold the same hills at full precision (the file type is detected automatically;
the time column is then replaced by the hill index).

1 CV : F(s) is written to --out (optionally after 1/B, 2/B, ... of the hills).
2 CVs: F(s1, s2) is written to --out2d (gnuplot "splot" format) and the
       projection on CV --cv,  F(s) = -kT ln sum_other exp(-F(s1,s2)/kT),
       is written to --out (needs the simulation temperature -T).
More CVs: use CP2K's graph tool (graph.psmp -cp2k -file <restart>) instead.

Examples:
    python3 fes_from_hills.py nacl_wtmetad_dft-HILLS.metadynLog --bohr2ang
    python3 fes_from_hills.py HILLS.metadynLog --bohr2ang --grid 2.0 6.0 401 --blocks 4 -o fes.dat
    python3 fes_from_hills.py nacl_wt2d_xtb-HILLS.metadynLog --bohr2ang 1 -T 300 -o fes_d.dat
    python3 fes_from_hills.py nacl_wtmetad_dft-1.restart --bohr2ang --blocks 4
Put the HILLS file name first (before --bohr2ang, which takes optional values).
Needs numpy (matplotlib only for --png).
"""
import argparse
import sys

import numpy as np

BOHR_ANG = 0.529177210903
HARTREE = {"kJ/mol": 2625.4996394799, "kcal/mol": 627.5094740631, "eV": 27.211386245988}
KB = {"kJ/mol": 0.0083144626, "kcal/mol": 0.0019872043, "eV": 8.617333262e-5}


def read_hills(path):
    rows = []
    with open(path, errors="replace") as fh:
        for line in fh:
            p = line.split()
            if not p or p[0].startswith("#"):
                continue
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                continue
    if not rows:
        sys.exit(f"no data in {path}")
    ncol = {len(r) for r in rows}
    if len(ncol) != 1:
        sys.exit(f"{path}: lines with different numbers of columns {sorted(ncol)}")
    a = np.array(rows)
    ncv = (a.shape[1] - 2) // 2
    if ncv < 1 or a.shape[1] != 2 * ncv + 2:
        sys.exit(f"unexpected number of columns ({a.shape[1]}) in {path}")
    return a, ncv


def _restart_block(lines, name):
    """Rows of numbers inside &<name> ... &END <name> of a CP2K restart file."""
    rows, inside, buf = [], False, ""
    for raw in lines:
        t = raw.strip()
        u = t.upper()
        if not inside:
            if u.split() and u.split()[0] == "&" + name:
                inside = True
            continue
        if u.startswith("&END"):
            break
        if not t or t.startswith(("!", "#")):
            continue
        if t.endswith("\\"):
            buf += t[:-1] + " "
            continue
        rows.append([float(x.replace("D", "E").replace("d", "e")) for x in (buf + t).split()])
        buf = ""
    return rows


def read_restart(path):
    """Hills from &MOTION/&FREE_ENERGY/&METADYN/&SPAWNED_HILLS_* of a restart file."""
    with open(path, errors="replace") as fh:
        lines = fh.readlines()
    pos = _restart_block(lines, "SPAWNED_HILLS_POS")
    scale = _restart_block(lines, "SPAWNED_HILLS_SCALE")
    height = _restart_block(lines, "SPAWNED_HILLS_HEIGHT")
    if not pos:
        sys.exit(f"{path}: no &SPAWNED_HILLS_POS section (no hills spawned yet?)")
    if len(scale) != len(pos) or len(height) != len(pos):
        sys.exit(f"{path}: SPAWNED_HILLS_POS/SCALE/HEIGHT have different lengths")
    ncv = len(pos[0])
    rows = [[float(i + 1)] + p + s + [h[0]] for i, (p, s, h) in enumerate(zip(pos, scale, height))]
    return np.array(rows), ncv


def gauss_matrix(grid, centers, widths, period=None):
    """(n_hills, n_grid) matrix of exp(-0.5*((grid - c)/w)^2); w == 0 -> 1."""
    d = grid[None, :] - centers[:, None]
    if period:
        d = d - period * np.round(d / period)
    w = widths[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.exp(-0.5 * (d / np.where(w == 0.0, 1.0, w)) ** 2)
    g[np.broadcast_to(w == 0.0, g.shape)] = 1.0
    return g


def extrema(x, y):
    out = []
    for i in range(1, len(y) - 1):
        if y[i] < y[i - 1] and y[i] <= y[i + 1]:
            out.append(("min", x[i], y[i]))
        elif y[i] > y[i - 1] and y[i] >= y[i + 1]:
            out.append(("max", x[i], y[i]))
    return out


def make_grid(spec, s, w, npts):
    if spec:
        return np.linspace(spec[0], spec[1], int(spec[2]))
    wmax = w.max() if w.max() > 0 else 1.0
    return np.linspace(s.min() - 3 * wmax, s.max() + 3 * wmax, npts)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("hills", help="<PROJECT>-HILLS.metadynLog or <PROJECT>-1.restart (full precision)")
    ap.add_argument("--cv", type=int, default=1, help="2 CVs: CV to project on (1 or 2), default 1")
    ap.add_argument("--bohr2ang", nargs="*", type=int, default=None, metavar="K",
                    help="convert these CVs (1-based) from bohr to Angstrom; no value = all CVs")
    ap.add_argument("--periodic", nargs="+", type=float, default=None, metavar="P",
                    help="period of each CV in its (converted) unit, 0 = not periodic, e.g. 6.2831853")
    ap.add_argument("--grid", nargs=3, type=float, metavar=("MIN", "MAX", "N"),
                    help="grid of CV1 (default: range of centers +/- 3 widths)")
    ap.add_argument("--grid2", nargs=3, type=float, metavar=("MIN", "MAX", "N"), help="grid of CV2 (2 CVs only)")
    ap.add_argument("-T", "--temperature", type=float, default=None, help="simulation temperature in K (needed for 2 CVs)")
    ap.add_argument("--unit", choices=list(HARTREE), default="kJ/mol")
    ap.add_argument("--blocks", type=int, default=1, help="also write F after 1/B, 2/B, ... of the hills (convergence check)")
    ap.add_argument("-o", "--out", default="fes.dat", help="1D profile (or projection) output")
    ap.add_argument("--out2d", default="fes2d.dat", help="2D surface output (2 CVs only)")
    ap.add_argument("--png", help="plot the 1D profile/projection to this file (needs matplotlib)")
    args = ap.parse_args()

    with open(args.hills, errors="replace") as fh:
        is_restart = "&SPAWNED_HILLS_POS" in fh.read().upper()
    a, ncv = read_restart(args.hills) if is_restart else read_hills(args.hills)
    if ncv > 2:
        sys.exit(f"{ncv} CVs in file: use CP2K's graph tool or plumed sum_hills for > 2 CVs")
    time = a[:, 0]
    s = a[:, 1:1 + ncv].copy()
    w = a[:, 1 + ncv:1 + 2 * ncv].copy()
    h = a[:, -1] * HARTREE[args.unit]

    conv = [] if args.bohr2ang is None else (list(range(1, ncv + 1)) if len(args.bohr2ang) == 0 else args.bohr2ang)
    for k in conv:
        if not 1 <= k <= ncv:
            sys.exit(f"--bohr2ang {k}: file has only {ncv} CV(s)")
        s[:, k - 1] *= BOHR_ANG
        w[:, k - 1] *= BOHR_ANG
    unit_of = ["A" if (k + 1) in conv else "internal" for k in range(ncv)]
    period = [None] * ncv
    if args.periodic:
        if len(args.periodic) not in (1, ncv):
            sys.exit("--periodic needs one value per CV")
        vals = args.periodic * ncv if len(args.periodic) == 1 else args.periodic
        period = [p if p > 0 else None for p in vals]

    nb = max(1, args.blocks)
    cuts = [int(round(len(h) * b / nb)) for b in range(1, nb + 1)]
    if is_restart:
        print(f"# {len(h)} hills, {ncv} CV(s), read from restart file (full precision)")
    else:
        print(f"# {len(h)} hills, {ncv} CV(s), t = {time[0]:.1f} ... {time[-1]:.1f} fs"
              " (HILLS heights have only 5 decimals in hartree)")
    print(f"# first/last Gaussian height: {h[0]:.4f} / {h[-1]:.4f} {args.unit} -- should decay in WT-MetaD")

    if ncv == 1:
        x = make_grid(args.grid, s[:, 0], w[:, 0], 400)
        g = gauss_matrix(x, s[:, 0], w[:, 0], period[0])
        cols, labels = [], []
        for n in cuts:
            f = -(h[:n, None] * g[:n]).sum(axis=0)
            cols.append(f - f.min())
            labels.append(f"F_{n}hills")
        xlabel = f"CV({unit_of[0]})"
    else:
        if args.temperature is None:
            sys.exit("2 CVs: give the simulation temperature with -T (needed for the projection)")
        if args.cv not in (1, 2):
            sys.exit("--cv must be 1 or 2")
        kT = KB[args.unit] * args.temperature
        x1 = make_grid(args.grid, s[:, 0], w[:, 0], 151)
        x2 = make_grid(args.grid2, s[:, 1], w[:, 1], 151)
        g1 = gauss_matrix(x1, s[:, 0], w[:, 0], period[0])
        g2 = gauss_matrix(x2, s[:, 1], w[:, 1], period[1])
        cols, labels = [], []
        f2 = None
        for n in cuts:
            f2 = -((h[:n, None] * g1[:n]).T @ g2[:n])          # shape (len(x1), len(x2))
            f2 = f2 - f2.min()
            axis = 1 if args.cv == 1 else 0
            m = (-f2 / kT).max(axis=axis, keepdims=True)
            proj = -kT * (np.log(np.exp(-f2 / kT - m).sum(axis=axis)) + m.squeeze(axis))
            cols.append(proj - proj.min())
            labels.append(f"F_{n}hills")
        x = x1 if args.cv == 1 else x2
        xlabel = f"CV{args.cv}({unit_of[args.cv - 1]})"
        with open(args.out2d, "w") as fh:
            fh.write(f"# CV1({unit_of[0]})  CV2({unit_of[1]})  F[{args.unit}]  (all {len(h)} hills, min = 0)\n")
            for i, xi in enumerate(x1):
                for j, yj in enumerate(x2):
                    fh.write(f"{xi:14.6f} {yj:14.6f} {f2[i, j]:14.6f}\n")
                fh.write("\n")
        i, j = np.unravel_index(np.argmin(f2), f2.shape)
        print(f"# 2D surface written to {args.out2d}; global minimum at CV1 = {x1[i]:.3f}, CV2 = {x2[j]:.3f}")

    data = np.column_stack([x] + cols)
    np.savetxt(args.out, data, fmt="%14.6f", header=f"{xlabel} " + " ".join(labels) + f"   [{args.unit}]")
    print(f"# 1D profile written to {args.out}")
    for kind, xi, yi in extrema(x, cols[-1]):
        print(f"  local {kind}: CV = {xi:8.3f}   F = {yi:8.3f} {args.unit}")

    if args.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.8, 3.4))
        for c, lab in zip(cols, labels):
            ax.plot(x, c, label=lab)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f"F ({args.unit})")
        if nb > 1:
            ax.legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(args.png, dpi=200)
        print(f"# plot written to {args.png}")


if __name__ == "__main__":
    main()
