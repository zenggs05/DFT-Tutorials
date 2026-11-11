#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chgdiff.py -- charge-density difference and Bader helper for VASP CHGCAR-type files.
              差分电荷密度 Δρ = ρ(AB) - ρ(A) - ρ(B) 与 Bader 辅助（AECCAR0 + AECCAR2 求和）。

Purpose
-------
* ``diff``  : check that all CHGCAR files have the same FFT grid and cell, warn if the
              fragment atoms moved with respect to the combined system, compute
              rho(AB) - rho(A) - rho(B) - ... and write a CHGCAR-format file that VESTA
              (and vaspkit / bader) can read.  The planar-averaged difference along an
              axis, Delta-lambda(z) in e/Angstrom, and the cumulative charge transfer
              Delta-Q(z) are plotted and exported.
* ``sum``   : add CHGCAR-type files, typically AECCAR0 + AECCAR2 -> CHGCAR_sum (the
              reference density for ``bader CHGCAR -ref CHGCAR_sum``).
* ``planar``: planar average (e/Angstrom) of any CHGCAR-type file, e.g. an existing
              CHGCAR_diff.
* ``check`` : only compare grids/cells/atoms of several files.

For spin-polarised CHGCARs the first block (total density rho_up + rho_down) is
used by default; ``--block mag`` selects the magnetisation block (rho_up - rho_down).
Augmentation occupancies are not written to the output files (VESTA and bader do
not need them).

Usage examples
--------------
  python chgdiff.py diff AB/CHGCAR slab/CHGCAR mol/CHGCAR -o CHGCAR_diff
  python chgdiff.py diff AB/CHGCAR slab/CHGCAR mol/CHGCAR --axis z --dat cdd_z.dat
  python chgdiff.py sum AECCAR0 AECCAR2 -o CHGCAR_sum      # then: bader CHGCAR -ref CHGCAR_sum
  python chgdiff.py planar CHGCAR_diff --axis z -o planar.png
  python chgdiff.py check AB/CHGCAR slab/CHGCAR mol/CHGCAR

Dependencies
------------
  Python >= 3.9, numpy, matplotlib (only for plots).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from typing import Dict, List, Optional

import numpy as np

SERIES = ["#2a78d6", "#eb6834"]
INK = "#1f1f1f"
MUTED = "#8a8986"
FIX_EXP = re.compile(r"(\d)([+-]\d{2,3})(?=\s|$)")
AXES = {"x": 0, "a": 0, "1": 0, "y": 1, "b": 1, "2": 1, "z": 2, "c": 2, "3": 2}


def die(msg: str) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def check_file(path: str, what: str) -> None:
    if not os.path.isfile(path):
        die(f"{what} not found: {path}")
    if os.path.getsize(path) == 0:
        die(f"{what} is empty: {path}")


# --------------------------------------------------------------------------- reading
def _parse_floats(text: str) -> np.ndarray:
    try:
        return np.array(text.split(), dtype=float)
    except ValueError:
        return np.array(FIX_EXP.sub(r"\1E\2", text).split(), dtype=float)


def _read_n_values(fh, n: int, path: str) -> np.ndarray:
    first = fh.readline()
    if not first:
        die(f"{path}: grid data missing")
    ncol = max(len(first.split()), 1)
    nlines = int(math.ceil(n / ncol)) - 1
    vals = _parse_floats(first + "".join(fh.readline() for _ in range(nlines)))
    while vals.size < n:
        line = fh.readline()
        if not line:
            die(f"{path}: file truncated ({vals.size} of {n} grid values)")
        vals = np.concatenate([vals, _parse_floats(line)])
    return vals[:n]


def read_chgcar(path: str, block: str = "total") -> Dict:
    """Read structure + one data block ('total' or 'mag') of a CHGCAR-like file.

    Values are kept exactly as stored (rho * V_cell for CHGCAR/AECCAR/PARCHG) in the
    Fortran order of the file (x fastest).
    """
    check_file(path, "charge-density file")
    header: List[str] = []
    with open(path, "r", errors="replace") as fh:
        try:
            for _ in range(5):
                header.append(fh.readline())
            scale = [float(x) for x in header[1].split()[:3]]
            lat = np.array([[float(x) for x in header[i].split()[:3]] for i in (2, 3, 4)])
            line = fh.readline()
            header.append(line)
            toks = line.split()
            if toks and not re.match(r"^\d+$", toks[0]):
                symbols = [re.match(r"[A-Z][a-z]?", re.match(r"[A-Za-z]+", t).group(0).capitalize()).group(0)
                           for t in toks]
                line = fh.readline()
                header.append(line)
                counts = [int(x) for x in line.split()[:len(symbols)]]
            else:
                counts = [int(x) for x in toks]
                symbols = [f"X{i + 1}" for i in range(len(counts))]
            line = fh.readline()
            header.append(line)
            if line.strip()[:1] in ("s", "S"):
                line = fh.readline()
                header.append(line)
            cart = line.strip()[:1] in ("c", "C", "k", "K")
            nat = sum(counts)
            pos = []
            for _ in range(nat):
                line = fh.readline()
                header.append(line)
                pos.append([float(x) for x in line.split()[:3]])
            pos = np.array(pos)
        except (ValueError, IndexError, AttributeError) as exc:
            die(f"{path}: cannot parse the structure header ({exc})")
        if len(scale) == 3:
            fac = np.array(scale)
        elif scale[0] < 0:
            fac = np.full(3, (abs(scale[0]) / abs(np.linalg.det(lat))) ** (1 / 3))
        else:
            fac = np.full(3, scale[0])
        lat = lat * fac[None, :]
        if cart:
            pos = (pos * fac[None, :]) @ np.linalg.inv(lat)
        line = fh.readline()
        while line and not line.strip():
            line = fh.readline()
        try:
            ngrid = tuple(int(x) for x in line.split()[:3])
            assert len(ngrid) == 3
        except (ValueError, AssertionError):
            die(f"{path}: grid dimensions not found after the atomic positions")
        n = int(np.prod(ngrid))
        data = _read_n_values(fh, n, path)
        nblocks_seen = 1
        if block == "mag":
            data = None
            for line in fh:
                t = line.split()
                if len(t) == 3 and all(x.isdigit() for x in t) and tuple(map(int, t)) == ngrid:
                    data = _read_n_values(fh, n, path)
                    nblocks_seen = 2
                    break
            if data is None:
                die(f"{path}: no magnetisation block found (ISPIN = 1 run?)")
    symlist: List[str] = []
    for s, c in zip(symbols, counts):
        symlist += [s] * c
    return {"path": path, "header": header, "lattice": lat, "frac": pos % 1.0, "symbols": symlist,
            "grid": ngrid, "data": data, "volume": abs(np.linalg.det(lat)), "blocks": nblocks_seen}


# --------------------------------------------------------------------------- checks
def compare(ref: Dict, other: Dict, tol_cell: float, tol_pos: float, fragment: bool) -> bool:
    ok = True
    if ref["grid"] != other["grid"]:
        print(f"  [FAIL] grid {other['grid']} of {other['path']} != {ref['grid']} of {ref['path']}")
        ok = False
    dl = float(np.abs(ref["lattice"] - other["lattice"]).max())
    if dl > tol_cell:
        print(f"  [FAIL] lattice of {other['path']} differs by up to {dl:.5f} A")
        ok = False
    if not fragment:
        return ok
    # each fragment atom should sit on an atom of the same element in the reference
    lat = ref["lattice"]
    worst = 0.0
    missing = 0
    for f, s in zip(other["frac"], other["symbols"]):
        cand = [i for i, t in enumerate(ref["symbols"]) if t == s] or list(range(len(ref["symbols"])))
        d = ref["frac"][cand] - f
        d -= np.round(d)
        dist = np.linalg.norm(d @ lat, axis=1).min() if len(cand) else np.inf
        worst = max(worst, dist)
        if dist > tol_pos:
            missing += 1
    if missing:
        print(f"  [WARN] {missing} atom(s) of {other['path']} are > {tol_pos} A away from the matching atom "
              f"in {ref['path']} (max {worst:.3f} A): fragments must use the frozen geometry of the combined system")
    else:
        print(f"  [ OK ] atoms of {other['path']} coincide with {ref['path']} (max deviation {worst:.4f} A)")
    return ok


# --------------------------------------------------------------------------- writing / analysis
def write_chgcar(path: str, ref: Dict, values: np.ndarray, comment: Optional[str]) -> None:
    hdr = list(ref["header"])
    if comment:
        hdr[0] = comment.rstrip("\n") + "\n"
    with open(path, "w", newline="\n") as fh:
        fh.writelines(hdr)
        fh.write("\n")
        fh.write(" ".join(f"{g:5d}" for g in ref["grid"]) + "\n")
        n5 = values.size // 5 * 5
        if n5:
            np.savetxt(fh, values[:n5].reshape(-1, 5), fmt=" %.11E", delimiter="")
        if values.size > n5:
            fh.write("".join(f" {v:.11E}" for v in values[n5:]) + "\n")
    print(f"Written: {path}  (grid {ref['grid']}, CHGCAR format: values = rho * V_cell)")


def planar(ref: Dict, values: np.ndarray, axis: int):
    """Planar-averaged line density (e/Angstrom) and coordinate along the plane normal."""
    grid = ref["grid"]
    arr = values.reshape(grid, order="F")
    other = tuple(i for i in range(3) if i != axis)
    i, j = other
    area = np.linalg.norm(np.cross(ref["lattice"][i], ref["lattice"][j]))
    height = ref["volume"] / area
    n = grid[axis]
    dz = height / n
    lam = arr.mean(axis=other) / height  # (rho*V)/V * A  = e/Angstrom
    z = np.arange(n) * dz
    q = np.cumsum(lam) * dz
    return z, lam, q, height


def setup_matplotlib(show: bool):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_planar(args, z, lam, q, ref, axis, plt, title: str) -> None:
    fig, axs = plt.subplots(2, 1, figsize=tuple(args.figsize), sharex=True)
    axs[0].axhline(0, color=MUTED, lw=0.8)
    axs[0].plot(z, lam, color=SERIES[0], lw=1.4)
    axs[0].fill_between(z, 0, lam, where=lam >= 0, interpolate=True, color=SERIES[0], alpha=0.18, lw=0)
    axs[0].fill_between(z, 0, lam, where=lam < 0, interpolate=True, color=SERIES[1], alpha=0.18, lw=0)
    axs[0].set_ylabel("$\\Delta\\lambda$ (e/Å)")
    axs[1].axhline(0, color=MUTED, lw=0.8)
    axs[1].plot(z, q, color=INK, lw=1.4)
    axs[1].set_ylabel("$\\Delta Q(z)$ (e)")
    axs[1].set_xlabel(f"{'xyz'[axis]} (Å)")
    zat = np.sort(ref["frac"][:, axis]) * (z[1] - z[0]) * len(z)
    for ax in axs:
        ax.tick_params(direction="in", length=3)
        ymin = ax.get_ylim()[0]
        ax.plot(zat, np.full_like(zat, ymin), "|", color=INK, ms=7, alpha=0.6)
        ax.set_xlim(z[0], z[-1])
    axs[0].set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=args.dpi)
    print(f"Figure written: {args.plot}")
    if args.show:
        plt.show()


def write_planar_dat(path: str, z, lam, q, axis: int, src: str) -> None:
    with open(path, "w") as fh:
        fh.write(f"# planar average of {src} along {'xyz'[axis]}: line density (e/A) and cumulative charge (e)\n")
        fh.write("z(A) lambda(e/A) Q_cum(e)\n")
        np.savetxt(fh, np.column_stack([z, lam, q]), fmt="%14.6e")
    print(f"Data written: {path}")


# --------------------------------------------------------------------------- commands
def cmd_check(args) -> int:
    files = [read_chgcar(p, args.block) for p in args.files]
    print(f"Reference: {files[0]['path']}  grid {files[0]['grid']}  atoms {len(files[0]['symbols'])}")
    ok = all(compare(files[0], f, args.tol_cell, args.tol_pos, True) for f in files[1:])
    print("All grids and cells consistent." if ok else "Inconsistent files - do not subtract them.")
    return 0 if ok else 1


def cmd_diff(args) -> None:
    ab = read_chgcar(args.files[0], args.block)
    print(f"rho(AB): {ab['path']}  grid {ab['grid']}  atoms {len(ab['symbols'])}  V = {ab['volume']:.3f} A^3")
    delta = ab["data"].copy()
    ok = True
    for p in args.files[1:]:
        f = read_chgcar(p, args.block)
        ok &= compare(ab, f, args.tol_cell, args.tol_pos, True)
        if ab["grid"] == f["grid"]:
            delta -= f["data"]
    if not ok:
        die("grids or cells differ: use identical cells, ENCUT and NGX/NGY/NGZ (and NGXF/NGYF/NGZF) "
            "for all three calculations")
    npts = delta.size
    print(f"Sum of Delta-rho over the cell: {delta.sum() / npts:+.5f} e "
          "(should be ~0 when the electron counts add up)")
    print(f"Delta-rho range: {delta.min() / ab['volume']:+.5f} ... {delta.max() / ab['volume']:+.5f} e/A^3")
    write_chgcar(args.output, ab, delta, f"Delta rho = {' - '.join(args.files)}")
    axis = AXES[args.axis]
    z, lam, q, h = planar(ab, delta, axis)
    print(f"Planar average along {'xyz'[axis]}: max accumulation {lam.max():+.4f} e/A at z = {z[np.argmax(lam)]:.2f} A, "
          f"max depletion {lam.min():+.4f} e/A at z = {z[np.argmin(lam)]:.2f} A")
    if args.dat:
        write_planar_dat(args.dat, z, lam, q, axis, args.output)
    if args.plot and not args.no_plot:
        plt = setup_matplotlib(args.show)
        plot_planar(args, z, lam, q, ab, axis, plt, "Charge density difference (planar average)")


def cmd_sum(args) -> None:
    ref = read_chgcar(args.files[0], args.block)
    total = ref["data"].copy()
    for p in args.files[1:]:
        f = read_chgcar(p, args.block)
        if not compare(ref, f, args.tol_cell, args.tol_pos, False):
            die("grids/cells differ; files cannot be summed")
        total += f["data"]
    print(f"Total charge in the summed density: {total.sum() / total.size:.4f} e")
    write_chgcar(args.output, ref, total, f"Sum of {' + '.join(args.files)}")
    if os.path.basename(args.files[0]).upper().startswith("AECCAR"):
        print("Next: bader CHGCAR -ref " + args.output)


def cmd_planar(args) -> None:
    ref = read_chgcar(args.files[0], args.block)
    axis = AXES[args.axis]
    z, lam, q, h = planar(ref, ref["data"], axis)
    print(f"{ref['path']}: integral over the cell = {q[-1]:.5f} e; planar average written along {'xyz'[axis]}")
    if args.dat:
        write_planar_dat(args.dat, z, lam, q, axis, ref["path"])
    if not args.no_plot:
        plt = setup_matplotlib(args.show)
        plot_planar(args, z, lam, q, ref, axis, plt, f"Planar average of {os.path.basename(ref['path'])}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Charge-density difference / AECCAR sum / planar average for VASP CHGCAR files.\n"
                    "差分电荷（先检查网格、晶胞与原子位置一致）、AECCAR0+AECCAR2 求和、平面平均。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples / 示例:\n"
               "  python chgdiff.py diff AB/CHGCAR A/CHGCAR B/CHGCAR -o CHGCAR_diff --dat cdd.dat\n"
               "  python chgdiff.py sum AECCAR0 AECCAR2 -o CHGCAR_sum\n"
               "  python chgdiff.py planar CHGCAR_diff --axis z\n"
               "  python chgdiff.py check AB/CHGCAR A/CHGCAR B/CHGCAR\n")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp, out_default=None, nmin=1):
        sp.add_argument("files", nargs="+", help="CHGCAR-type files")
        sp.add_argument("--block", choices=["total", "mag"], default="total",
                        help="data block of spin-polarised files: total (default) or mag(netisation)")
        sp.add_argument("--tol-cell", type=float, default=1e-4, help="lattice tolerance in Angstrom (default 1e-4)")
        sp.add_argument("--tol-pos", type=float, default=0.02,
                        help="atom-position tolerance in Angstrom for fragment checks (default 0.02)")
        if out_default:
            sp.add_argument("-o", "--output", default=out_default, help=f"output file (default {out_default})")

    def plotting(sp, default_png):
        sp.add_argument("--axis", default="z", choices=sorted(AXES.keys()),
                        help="planar-average direction (default z = c)")
        sp.add_argument("--dat", default=None, help="write z, lambda(e/A), cumulative Q(e) to this file")
        sp.add_argument("--plot", default=default_png, help=f"figure file (default {default_png})")
        sp.add_argument("--no-plot", action="store_true")
        sp.add_argument("--figsize", type=float, nargs=2, default=[6.0, 4.6], metavar=("W", "H"))
        sp.add_argument("--dpi", type=int, default=300)
        sp.add_argument("--show", action="store_true", help="open an interactive window")

    d = sub.add_parser("diff", help="rho(AB) - rho(A) - rho(B) ... -> CHGCAR_diff (+ planar average)")
    common(d, "CHGCAR_diff")
    plotting(d, "chgdiff_planar.png")
    s = sub.add_parser("sum", help="sum of files, e.g. AECCAR0 + AECCAR2 -> CHGCAR_sum")
    common(s, "CHGCAR_sum")
    pl = sub.add_parser("planar", help="planar average of one CHGCAR-type file")
    common(pl)
    plotting(pl, "planar.png")
    c = sub.add_parser("check", help="only check grids/cells/atom positions")
    common(c)
    return p


def main(argv=None) -> None:
    try:
        sys.stdout.reconfigure(errors="backslashreplace")  # type: ignore[attr-defined]
    except (AttributeError, ValueError):
        pass
    args = build_parser().parse_args(argv)
    need = {"diff": 2, "sum": 2, "planar": 1, "check": 2}[args.cmd]
    if len(args.files) < need:
        die(f"'{args.cmd}' needs at least {need} file(s)")
    if args.cmd == "planar" and len(args.files) > 1:
        warn("planar uses only the first file")
    if args.cmd == "diff":
        cmd_diff(args)
    elif args.cmd == "sum":
        cmd_sum(args)
    elif args.cmd == "planar":
        cmd_planar(args)
    else:
        sys.exit(cmd_check(args))


if __name__ == "__main__":
    main()
