#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
neb_plot.py -- energy profile (MEP) of a CP2K RUN_TYPE BAND (NEB / CI-NEB) job.
               收集 CP2K BAND 计算各副本的能量，以累计位移为反应坐标画 MEP，
               标注正/逆反应能垒和反应能，并导出能量最高副本（过渡态近似）的 xyz。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Accepted sources (detected automatically)
-----------------------------------------
* main output  : summary blocks of &BAND/&PROGRAM_RUN_INFO
                   DISTANCES REP = ...   (bohr, 3N-dimensional distances)
                   ENERGIES [au] = ...   (hartree)
                 or, with PLOT_REL_ENERGY T (CP2K >= 2026.2), the table
                   REPLICA  ENERGY [au]  TEMPERATURE [K] ... / DISTANCE = ...
* <PROJECT>-*.ener (&BAND/&ENERGY print key): one line per band step,
                   step, E(1..N) [hartree], d(1..N-1) [bohr]
* <PROJECT>-FINAL-1.xyz (CP2K >= 2026.2): extended XYZ with
                   Replica=, Energy= and Ener_loc_max= in the titles
* the replica trajectories <PROJECT>-pos-Replica_nr_*-1.xyz (give all of them);
  the last frame of each is used, energy from the title "E = ...".
With xyz input the path length is computed from the coordinates (Angstrom);
otherwise the CP2K distances (bohr) are converted to Angstrom.
The per-replica logs <PROJECT>-BAND*.out belong to replica *groups* (several
replicas can share one group) - use them for error checking only.

Structures for --ts-xyz: taken from the xyz source, or from --xyz FILES
(a FINAL-1.xyz file or the pos-Replica_nr_* trajectories) when the energies
come from the main output or the .ener file.

Usage examples
--------------
    python neb_plot.py h_diffusion_cineb.out -o mep.png
    python neb_plot.py h_cu111_neb-1.ener --step -1 --units kJ/mol --csv mep.csv
    python neb_plot.py h_cu111_neb-FINAL-1.xyz --ts-xyz ts_from_neb.xyz --spline
    python neb_plot.py h_cu111_neb-pos-Replica_nr_*-1.xyz --ts-xyz ts_guess.xyz
    python neb_plot.py h_diffusion_cineb.out --xyz "h_cu111_neb-pos-Replica_nr_*-1.xyz" --ts-xyz ts.xyz

Dependencies: python >= 3.9, numpy, matplotlib (only for figures).
"""
import argparse
import glob
import os
import re
import sys

HA2EV = 27.211386245988
BOHR = 0.529177210903
NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
UNITS = {"eV": HA2EV, "kJ/mol": 2625.4996394799, "kcal/mol": 627.5094740631}


def die(msg):
    sys.exit("neb_plot.py: error: " + msg)


def need_numpy():
    try:
        import numpy as np
    except ImportError:
        die("numpy is required (pip install numpy)")
    return np


def fl(x):
    return float(x.replace("D", "E").replace("d", "e"))


# --------------------------------------------------------------------------- readers
def read_main_output(path):
    """List of (energies[Ha], distances[bohr]) for every band step found."""
    blocks = []
    key, vals, dist = None, [], None
    table, tab_e, tab_d = False, [], []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            # compact format
            m = re.match(r"\s*(DISTANCES REP|ENERGIES \[au\])\s*=(.*)", line)
            if m:
                if key == "DISTANCES REP":
                    dist = vals
                cur = [fl(x) for x in re.findall(NUM, m.group(2))]
                key, vals = m.group(1), cur
                continue
            if key and re.match(r"^\s+" + NUM + r"(\s+" + NUM + r")*\s*$", line):
                vals += [fl(x) for x in re.findall(NUM, line)]
                continue
            if key == "DISTANCES REP":
                dist = vals
            elif key == "ENERGIES [au]":
                blocks.append((vals, dist))
                dist = None
            key, vals = None, []
            # PLOT_REL_ENERGY table (2026.2+)
            if re.search(r"REPLICA\s+ENERGY \[au\]", line):
                table, tab_e, tab_d = True, [], []
                continue
            if table:
                m = re.match(r"\s*(\d+)\s+(" + NUM + r")\s+\(", line)
                if m:
                    tab_e.append(fl(m.group(2)))
                    continue
                m = re.match(r"\s*DISTANCE\s*=\s*(" + NUM + ")", line)
                if m:
                    tab_d.append(fl(m.group(1)))
                    continue
                if "v Nr." in line or "NUMBER OF LOCAL MAXIMA" in line:
                    if tab_e:
                        blocks.append((tab_e, tab_d))
                    table = False
    if key == "ENERGIES [au]":
        blocks.append((vals, dist))
    good = [(e, d) for e, d in blocks if d is not None and len(d) == len(e) - 1 and len(e) >= 2]
    if not good:
        die(f"{path}: no BAND summary blocks ('ENERGIES [au]' / 'DISTANCES REP' or the PLOT_REL_ENERGY table). "
            "Switch on &BAND/&PROGRAM_RUN_INFO (default level LOW in 2026.2, MEDIUM before).")
    return good


def read_ener(path):
    rows = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            t = line.split()
            try:
                v = [fl(x) for x in t]
            except ValueError:
                continue
            if len(v) < 4 or len(v) % 2:
                die(f"{path}: unexpected column count {len(v)} - expected step + N energies + (N-1) distances")
            n = len(v) // 2
            rows.append((v[1:1 + n], v[1 + n:]))
    if not rows:
        die(f"{path}: no data lines")
    return rows


def read_xyz_frames(path):
    """All frames of an xyz file: list of (title, symbols, coords[A])."""
    frames = []
    with open(path, "r", errors="replace") as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if not line.strip():
                continue
            try:
                na = int(line.split()[0])
            except ValueError:
                die(f"{path}: expected an atom count, got: {line.strip()[:60]}")
            title = fh.readline().rstrip("\n")
            atoms = [fh.readline() for _ in range(na)]
            if len(atoms) < na or not atoms[-1]:
                print(f"warning: {path}: incomplete last frame ignored")
                break
            syms, xyz = [], []
            for a in atoms:
                t = a.split()
                if len(t) < 4:
                    die(f"{path}: malformed atom line: {a.strip()[:60]}")
                syms.append(t[0])
                xyz.append([fl(t[1]), fl(t[2]), fl(t[3])])
            frames.append((title, syms, xyz))
    if not frames:
        die(f"{path}: no frames")
    return frames


def images_from_xyz(paths):
    """Return list of images (energy[Ha] or None, symbols, coords, title) and the lattice (or None)."""
    images, lattice = [], None
    if len(paths) == 1:
        frames = read_xyz_frames(paths[0])
        if any("Replica=" in f[0] for f in frames):
            for title, syms, xyz in frames:
                m = re.search(r"Replica=(\d+)", title)
                e = re.search(r"Energy=(" + NUM + ")", title)
                lat = re.search(r'Lattice="([^"]+)"', title)
                if lat and lattice is None:
                    lattice = [fl(x) for x in lat.group(1).split()]
                images.append((int(m.group(1)) if m else len(images) + 1,
                               fl(e.group(1)) if e else None, syms, xyz, title))
            images.sort(key=lambda t: t[0])
            return [i[1:] for i in images], lattice
        print(f"note: {paths[0]} is a single trajectory; treating every frame as an image")
        for title, syms, xyz in frames:
            e = re.search(r"\bE\s*=\s*(" + NUM + ")", title)
            images.append((fl(e.group(1)) if e else None, syms, xyz, title))
        return images, None

    def rep_no(p):
        m = re.search(r"Replica_nr_0*(\d+)", os.path.basename(p))
        return int(m.group(1)) if m else 10 ** 9

    for p in sorted(paths, key=rep_no):
        title, syms, xyz = read_xyz_frames(p)[-1]
        e = re.search(r"\bE\s*=\s*(" + NUM + ")", title)
        images.append((fl(e.group(1)) if e else None, syms, xyz, title))
    if any(rep_no(p) == 10 ** 9 for p in paths):
        print("warning: some xyz names lack 'Replica_nr_<n>'; images ordered by file name")
    return images, None


def expand(patterns):
    out = []
    for p in patterns:
        hits = sorted(glob.glob(p))
        if not hits and os.path.isfile(p):
            hits = [p]
        if not hits:
            die(f"no file matches {p}")
        out += hits
    return list(dict.fromkeys(out))


def path_from_coords(np, images, lattice, mic):
    X = [np.array(im[2], dtype=float) for im in images]
    if any(x.shape != X[0].shape for x in X):
        die("images have different numbers of atoms")
    d = []
    h = hinv = None
    if mic:
        if lattice is None:
            die("--mic needs a Lattice= entry (FINAL-1.xyz) in the xyz titles")
        h = np.array(lattice).reshape(3, 3)
        hinv = np.linalg.inv(h)
    for a, b in zip(X[:-1], X[1:]):
        diff = b - a
        if mic:
            f = diff @ hinv
            f -= np.round(f)
            diff = f @ h
        d.append(float(np.linalg.norm(diff)))
    return d


def natural_spline(np, x, y, xs):
    n = len(x)
    if n < 3:
        return np.interp(xs, x, y)
    h = np.diff(x)
    if (h <= 0).any():
        return np.interp(xs, x, y)
    A = np.zeros((n, n))
    rhs = np.zeros(n)
    A[0, 0] = A[-1, -1] = 1.0
    for i in range(1, n - 1):
        A[i, i - 1] = h[i - 1]
        A[i, i] = 2 * (h[i - 1] + h[i])
        A[i, i + 1] = h[i]
        rhs[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
    M = np.linalg.solve(A, rhs)
    idx = np.clip(np.searchsorted(x, xs) - 1, 0, n - 2)
    x0, x1, hh = x[idx], x[idx + 1], h[idx]
    return (M[idx] * (x1 - xs) ** 3 / (6 * hh) + M[idx + 1] * (xs - x0) ** 3 / (6 * hh)
            + (y[idx] / hh - M[idx] * hh / 6) * (x1 - xs) + (y[idx + 1] / hh - M[idx + 1] * hh / 6) * (xs - x0))


def write_xyz(path, syms, xyz, comment):
    with open(path, "w", newline="\n") as fh:
        fh.write(f"{len(syms)}\n{comment}\n")
        for s, (x, y, z) in zip(syms, xyz):
            fh.write(f"{s:<3s} {x:16.10f} {y:16.10f} {z:16.10f}\n")


# --------------------------------------------------------------------------- main
def build_parser():
    ap = argparse.ArgumentParser(
        prog="neb_plot.py",
        description="Plot the MEP of a CP2K BAND (NEB/CI-NEB) run and report barriers. "
                    "CP2K NEB 能量剖面、能垒与过渡态结构导出。",
        epilog="Examples:\n"
               "  python neb_plot.py h_diffusion_cineb.out -o mep.png\n"
               "  python neb_plot.py h_cu111_neb-FINAL-1.xyz --ts-xyz ts_from_neb.xyz --spline\n"
               "  python neb_plot.py h_cu111_neb-1.ener --xyz 'h_cu111_neb-pos-Replica_nr_*-1.xyz' --ts-xyz ts.xyz\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+",
                    help="main output, .ener file, FINAL-1.xyz, or all pos-Replica_nr_*-1.xyz files (globs allowed)")
    ap.add_argument("--step", type=int, default=-1,
                    help="band step to use from .out/.ener (1-based, negative counts from the end; default -1 = last)")
    ap.add_argument("--xyz", nargs="+", metavar="FILE",
                    help="structures for --ts-xyz when energies come from .out/.ener (FINAL-1.xyz or Replica files)")
    ap.add_argument("--ts-xyz", metavar="FILE", help="write the highest-energy intermediate image to this xyz file")
    ap.add_argument("--x", choices=["dist", "index"], default="dist",
                    help="x axis: cumulative path length in A (default) or image number")
    ap.add_argument("--mic", action="store_true",
                    help="minimum-image distances between images (needs Lattice= in FINAL-1.xyz)")
    ap.add_argument("--ref", choices=["first", "min"], default="first",
                    help="energy zero: first image (default) or lowest image")
    ap.add_argument("--units", choices=list(UNITS), default="eV", help="energy unit (default eV)")
    ap.add_argument("--spline", action="store_true", help="draw a natural cubic spline through the images")
    ap.add_argument("--csv", metavar="FILE", help="write image, s, E, E_rel to CSV")
    ap.add_argument("-o", "--output", default="mep.png", help="figure file (default mep.png; 'none' = no figure)")
    ap.add_argument("--title", default=None, help="figure title")
    ap.add_argument("--dpi", type=int, default=200, help="figure resolution (default 200)")
    ap.add_argument("--show", action="store_true", help="open an interactive window")
    return ap


def main():
    args = build_parser().parse_args()
    np = need_numpy()
    files = expand(args.inputs)
    xyz_in = all(f.lower().endswith(".xyz") for f in files)
    structs = None
    lattice = None
    if xyz_in:
        images, lattice = images_from_xyz(files)
        if any(im[0] is None for im in images):
            die("some images carry no energy in their xyz title")
        energies = [im[0] for im in images]
        dist = path_from_coords(np, images, lattice, args.mic) if len(images) > 1 else []
        structs = images
        src = f"xyz ({len(files)} file(s)), path length from coordinates"
    else:
        if len(files) != 1:
            die("give one main output or one .ener file (or only xyz files)")
        f = files[0]
        steps = read_ener(f) if f.lower().endswith(".ener") else read_main_output(f)
        try:
            e, d = steps[args.step - 1 if args.step > 0 else args.step]
        except IndexError:
            die(f"--step {args.step}: only {len(steps)} band step(s) found")
        energies = list(e)
        dist = [x * BOHR for x in d]
        src = f"{os.path.basename(f)}, band step {args.step if args.step > 0 else len(steps) + args.step + 1} of {len(steps)}"
        if args.xyz:
            structs, lattice = images_from_xyz(expand(args.xyz))
            if len(structs) != len(energies):
                die(f"--xyz gives {len(structs)} images but the energies list has {len(energies)}")
    n = len(energies)
    if n < 2:
        die("need at least two images")
    E = np.array(energies)
    s = np.concatenate([[0.0], np.cumsum(dist)]) if args.x == "dist" else np.arange(1, n + 1, dtype=float)
    e0 = E[0] if args.ref == "first" else E.min()
    fac = UNITS[args.units]
    rel = (E - e0) * fac
    imax = int(np.argmax(E[1:-1])) + 1 if n > 2 else int(np.argmax(E))
    maxima = [i for i in range(1, n - 1) if E[i] > E[i - 1] and E[i] > E[i + 1]]
    minima = [i for i in range(1, n - 1) if E[i] < E[i - 1] and E[i] < E[i + 1]]

    u = args.units
    print(f"# source: {src}")
    print(f"{'image':>5s} {'s (A)' if args.x == 'dist' else 'index':>10s} {'E (Ha)':>18s} {'E_rel (' + u + ')':>14s}")
    for i in range(n):
        tag = "  <- highest intermediate image" if i == imax else ""
        print(f"{i + 1:5d} {s[i]:10.4f} {E[i]:18.8f} {rel[i]:14.4f}{tag}")
    fw = (E[imax] - E[0]) * fac
    rv = (E[imax] - E[-1]) * fac
    de = (E[-1] - E[0]) * fac
    print(f"forward barrier  E(max) - E(1) = {fw:.4f} {u}")
    print(f"reverse barrier  E(max) - E(N) = {rv:.4f} {u}")
    print(f"reaction energy  E(N)   - E(1) = {de:.4f} {u}")
    if args.x == "dist":
        dd = np.diff(s)
        print(f"image spacing: min {dd.min():.3f}, max {dd.max():.3f}, mean {dd.mean():.3f} A (3N-dimensional)")
    if len(maxima) > 1:
        print(f"warning: {len(maxima)} local maxima (images {[i + 1 for i in maxima]}); "
              f"minima at {[i + 1 for i in minima]} - a stable intermediate? split the path into two NEBs")
    if E[imax] <= max(E[0], E[-1]):
        print("note: no intermediate image lies above both end points (barrierless on this path?)")
    print("note: the highest CI-NEB image approximates the TS; confirm with Dimer + frequency analysis")

    if args.ts_xyz:
        if structs is None:
            die("--ts-xyz needs structures: use an xyz source or --xyz FILES")
        _, syms, xyz, _ = structs[imax]
        write_xyz(args.ts_xyz, syms, xyz, f"image {imax + 1} of {n} (highest energy) E = {E[imax]:.10f} Ha; "
                                           "from neb_plot.py")
        print(f"wrote {args.ts_xyz} (image {imax + 1})")

    if args.csv:
        with open(args.csv, "w", newline="\n") as fh:
            fh.write(f"image,s_A,E_Ha,E_rel_{u.replace('/', '_per_')}\n")
            for i in range(n):
                fh.write(f"{i + 1},{s[i]:.6f},{E[i]:.10f},{rel[i]:.6f}\n")
        print(f"wrote {args.csv}")

    if args.output.lower() == "none" and not args.show:
        return
    try:
        import matplotlib
        if not args.show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        die("matplotlib is required for plotting (pip install matplotlib); use -o none")
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    if args.spline and n > 2:
        xs = np.linspace(s[0], s[-1], 300)
        ax.plot(xs, natural_spline(np, s, rel, xs), color="#2a78d6", lw=1.5)
        ax.plot(s, rel, "o", color="#2a78d6", ms=6, mec="white", mew=1.5)
    else:
        ax.plot(s, rel, "o-", color="#2a78d6", lw=1.5, ms=6, mec="white", mew=1.5)
    ax.plot([s[imax]], [rel[imax]], "o", color="#e34948", ms=8, mec="white", mew=1.5)
    ax.annotate(f"E$_a$ = {fw:.2f} {u}", (s[imax], rel[imax]), textcoords="offset points", xytext=(0, 8),
                ha="center", fontsize=9, color="#2f2f2d")
    ax.annotate(f"$\\Delta E$ = {de:.2f} {u}", (s[-1], rel[-1]), textcoords="offset points", xytext=(-4, -14),
                ha="right", fontsize=9, color="#2f2f2d")
    ax.axhline(0.0, color="#8a8a86", lw=0.8, ls=":")
    ax.set_xlabel(r"reaction coordinate ($\mathrm{\AA}$)" if args.x == "dist" else "image")
    ax.set_ylabel(f"E - E$_\\mathrm{{{'1' if args.ref == 'first' else 'min'}}}$ ({u})")
    ax.grid(alpha=0.25, lw=0.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if args.title:
        ax.set_title(args.title)
    fig.tight_layout()
    if args.output.lower() != "none":
        fig.savefig(args.output, dpi=args.dpi)
        print(f"wrote {args.output}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
