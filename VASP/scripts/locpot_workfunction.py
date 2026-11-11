#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
locpot_workfunction.py -- planar/macroscopic average of a VASP LOCPOT and work function.
                          LOCPOT 平面平均、宏观平均，自动识别真空平台并计算功函数。

Purpose
-------
* Read LOCPOT (or any CHGCAR-like volumetric file) and average the potential over
  the lattice planes perpendicular to a chosen axis (default: c / z).
* Optional macroscopic average (moving average over one or two periods, e.g. the
  interlayer spacing), useful for band alignment and bulk-like reference levels.
* Locate the vacuum region from the atomic positions and detect flat plateaus of
  the planar-averaged potential inside it.  With dipole correction
  (LDIPOL = .TRUE., IDIPOL = 3) two plateaus with different values are expected;
  both vacuum levels and both work functions Phi = V_vac - E_F are reported.
* E-fermi is read from OUTCAR (the last "E-fermi" line) unless --efermi is given;
  LDIPOL/IDIPOL and, for VASP >= 6.4.3 with LVACPOTAV = .TRUE., the line
  "vacuum level on the upper side and lower side of the slab" are echoed for
  comparison.
* Warns when the vacuum potential is not flat (missing dipole correction for an
  asymmetric slab, vacuum too thin, or LOCPOT not written with LVHAR = .TRUE.).

Usage examples
--------------
  python locpot_workfunction.py                          # ./LOCPOT, ./OUTCAR, axis z
  python locpot_workfunction.py -f LOCPOT --outcar OUTCAR -o wf.png --dat planar.dat
  python locpot_workfunction.py --efermi -1.234 --axis z --macro 2.27
  python locpot_workfunction.py --macro 3.2 2.9          # double macroscopic average (interface)

Dependencies
------------
  Python >= 3.9, numpy, matplotlib.

Notes
-----
  Write the electrostatic (Hartree + ionic) potential with LVHAR = .TRUE.; only the
  first data block of LOCPOT is used (this potential is spin independent).

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK = "#1f1f1f"
MUTED = "#8a8986"
FIX_EXP = re.compile(r"(\d)([+-]\d{2,3})(?=\s|$)")


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


# --------------------------------------------------------------------------- volumetric reader
def _read_n_values(fh, n: int, path: str) -> np.ndarray:
    first = fh.readline()
    if not first:
        die(f"{path}: no grid data after the grid dimensions")
    ncol = max(len(first.split()), 1)
    nlines = int(math.ceil(n / ncol)) - 1
    chunks = [first] + [fh.readline() for _ in range(nlines)]
    text = "".join(chunks)
    try:
        vals = np.array(text.split(), dtype=float)
    except ValueError:
        vals = np.array(FIX_EXP.sub(r"\1E\2", text).split(), dtype=float)
    while vals.size < n:  # irregular line lengths: keep reading
        line = fh.readline()
        if not line:
            die(f"{path}: file truncated ({vals.size} of {n} grid values)")
        vals = np.concatenate([vals, np.array(FIX_EXP.sub(r"\1E\2", line).split(), dtype=float)])
    return vals[:n]


def read_volumetric(path: str) -> Dict:
    """Read the structure and the first data block of LOCPOT/CHGCAR/ELFCAR files."""
    check_file(path, "volumetric file")
    with open(path, "r", errors="replace") as fh:
        try:
            comment = fh.readline()
            scale = [float(x) for x in fh.readline().split()[:3]]
            lat = np.array([[float(x) for x in fh.readline().split()[:3]] for _ in range(3)])
            line = fh.readline()
            toks = line.split()
            if toks and not re.match(r"^\d+$", toks[0]):
                symbols = [re.match(r"[A-Z][a-z]?", re.match(r"[A-Za-z]+", t).group(0).capitalize()).group(0)
                           for t in toks]
                counts = [int(x) for x in fh.readline().split()[:len(symbols)]]
            else:
                counts = [int(x) for x in toks]
                symbols = [f"X{i + 1}" for i in range(len(counts))]
            line = fh.readline()
            if line.strip()[:1] in ("s", "S"):
                line = fh.readline()
            cart = line.strip()[:1] in ("c", "C", "k", "K")
            nat = sum(counts)
            pos = np.array([[float(x) for x in fh.readline().split()[:3]] for _ in range(nat)])
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
        n = ngrid[0] * ngrid[1] * ngrid[2]
        vals = _read_n_values(fh, n, path)
    data = vals.reshape(ngrid, order="F")  # data[ix, iy, iz]
    symlist: List[str] = []
    for s, c in zip(symbols, counts):
        symlist += [s] * c
    return {"lattice": lat, "frac": pos % 1.0, "symbols": symlist, "grid": ngrid, "data": data,
            "comment": comment.strip()}


# --------------------------------------------------------------------------- analysis
def planar_average(data: np.ndarray, axis: int) -> np.ndarray:
    other = tuple(i for i in range(3) if i != axis)
    return data.mean(axis=other)


def plane_spacing(lat: np.ndarray, axis: int) -> float:
    """Length of the cell along the normal of the plane spanned by the other two vectors."""
    i, j = [k for k in range(3) if k != axis]
    area = np.linalg.norm(np.cross(lat[i], lat[j]))
    return abs(np.linalg.det(lat)) / area


def macro_average(v: np.ndarray, dz: float, periods: List[float]) -> np.ndarray:
    out = v.copy()
    n = len(v)
    for L in periods:
        w = max(1, int(round(L / dz)))
        kernel = np.zeros(n)
        idx = (np.arange(w) - w // 2) % n
        kernel[idx] = 1.0 / w
        out = np.real(np.fft.ifft(np.fft.fft(out) * np.fft.fft(kernel)))
    return out


def vacuum_region(frac_axis: np.ndarray) -> Tuple[float, float, float]:
    """Largest periodic gap between atomic planes: returns (top, bottom+1, width) in fractions."""
    z = np.sort(frac_axis % 1.0)
    if z.size == 0:
        return 0.0, 1.0, 1.0
    gaps = np.diff(np.concatenate([z, [z[0] + 1.0]]))
    k = int(np.argmax(gaps))
    top = z[k]
    bottom = z[(k + 1) % len(z)] + (1.0 if k + 1 >= len(z) else 0.0)
    return top, bottom, gaps[k]


def find_plateaus(z: np.ndarray, v: np.ndarray, start: int, stop: int, length: float,
                  grad_tol: float, min_len: float) -> List[Dict]:
    """Flat segments of v between grid indices start..stop (stop may exceed n: periodic)."""
    n = len(v)
    dz = length / n
    grad = (np.roll(v, -1) - np.roll(v, 1)) / (2 * dz)
    idx = np.arange(start, stop + 1)
    flat = np.abs(grad[idx % n]) < grad_tol
    segs: List[Dict] = []
    i = 0
    while i < len(idx):
        if flat[i]:
            j = i
            while j + 1 < len(idx) and flat[j + 1]:
                j += 1
            seg = idx[i:j + 1]
            if (len(seg) * dz) >= min_len:
                vals = v[seg % n]
                segs.append({"i0": int(seg[0]), "i1": int(seg[-1]), "z0": seg[0] * dz, "z1": seg[-1] * dz,
                             "len": len(seg) * dz, "mean": float(vals.mean()), "std": float(vals.std()),
                             "slope": float(np.polyfit(seg * dz, vals, 1)[0]) if len(seg) > 2 else 0.0})
            i = j + 1
        else:
            i += 1
    return segs


def read_outcar_info(path: Optional[str]) -> Dict:
    info: Dict = {}
    if not path or not os.path.isfile(path):
        return info
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if "E-fermi" in line and ":" in line:
                try:
                    info["efermi"] = float(line.split(":")[1].split()[0])
                except (ValueError, IndexError):
                    pass
            elif "Fermi energy:" in line and "efermi" not in info:
                try:
                    info["efermi_alt"] = float(line.split(":")[1].split()[0])
                except (ValueError, IndexError):
                    pass
            elif re.match(r"^\s*LDIPOL\s*=", line):
                info["ldipol"] = line.split("=")[1].split()[0]
            elif re.match(r"^\s*IDIPOL\s*=", line):
                info["idipol"] = line.split("=")[1].split()[0]
            elif "vacuum level on the upper side and lower side" in line:
                nums = re.findall(r"[-+]?\d+\.\d+", line)
                if len(nums) >= 2:
                    info["vac_outcar"] = (float(nums[0]), float(nums[1]))
            elif line.strip().startswith("dipolmoment"):
                nums = re.findall(r"[-+]?\d+\.\d+", line)
                if len(nums) >= 3:
                    info["dipole"] = tuple(float(x) for x in nums[:3])
    if "efermi" not in info and "efermi_alt" in info:
        info["efermi"] = info["efermi_alt"]
    return info


# --------------------------------------------------------------------------- plotting
def setup_matplotlib(show: bool):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def make_plot(args, z, v, vm, plateaus, efermi, atom_z, axis_name, plt) -> None:
    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    ax.plot(z, v, color=MUTED, lw=1.0, label="planar average", zorder=2)
    if vm is not None:
        ax.plot(z, vm, color=SERIES[0], lw=1.6, label="macroscopic average", zorder=3)
    L = z[-1] + (z[1] - z[0])
    for k, p in enumerate(plateaus):
        c = SERIES[1] if k == 0 else SERIES[2]
        z0, z1 = p["z0"] % L, p["z1"] % L
        spans = [(z0, z1)] if z1 >= z0 else [(z0, L), (0.0, z1)]
        for a, b in spans:
            ax.hlines(p["mean"], a, b, color=c, lw=3.0, zorder=4)
        ax.annotate(f"$V_\\mathrm{{vac}}$ = {p['mean']:.3f} eV", xy=(spans[0][0], p["mean"]),
                    xytext=(0, 6), textcoords="offset points", fontsize=8, color=INK)
    if efermi is not None:
        ax.axhline(efermi, color=INK, lw=0.9, ls="--", zorder=1, label=f"$E_\\mathrm{{F}}$ = {efermi:.3f} eV")
    ymin, ymax = float(np.min(v)), float(np.max(v))
    ax.plot(atom_z, np.full_like(atom_z, ymin), "|", color=INK, ms=8, alpha=0.6, zorder=5)
    ax.set_xlim(0, L)
    ax.set_ylim(ymin - 0.05 * (ymax - ymin), ymax + 0.15 * (ymax - ymin))
    ax.set_xlabel(f"{axis_name} (Å)")
    ax.set_ylabel("Potential energy (eV)")
    ax.tick_params(direction="in", length=3)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    if args.title:
        ax.set_title(args.title)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)
    print(f"Figure written: {args.output}")
    if args.show:
        plt.show()


# --------------------------------------------------------------------------- main
AXES = {"x": 0, "a": 0, "1": 0, "y": 1, "b": 1, "2": 1, "z": 2, "c": 2, "3": 2}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Planar/macroscopic average of a VASP LOCPOT, automatic vacuum-plateau detection and "
                    "work function (Phi = V_vac - E_F), incl. a check of dipole-corrected slabs.\n"
                    "LOCPOT 平面/宏观平均、真空平台识别、功函数计算以及偶极校正两侧真空能级检查。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples / 示例:\n"
               "  python locpot_workfunction.py -f LOCPOT --outcar OUTCAR\n"
               "  python locpot_workfunction.py --efermi -0.85 --macro 2.27 --dat planar.dat\n"
               "INCAR for the LOCPOT: LVHAR = .TRUE.; asymmetric slabs: LDIPOL = .TRUE., IDIPOL = 3")
    p.add_argument("-f", "--file", default="LOCPOT", help="LOCPOT file (default ./LOCPOT)")
    p.add_argument("--outcar", default=None, help="OUTCAR for E-fermi (default: OUTCAR next to LOCPOT)")
    p.add_argument("--efermi", type=float, default=None, help="Fermi energy in eV (overrides OUTCAR) / 费米能")
    p.add_argument("--axis", default="z", choices=sorted(AXES.keys()),
                   help="averaging direction (lattice vector index; default z = c) / 平均方向")
    p.add_argument("--macro", type=float, nargs="+", default=None, metavar="L",
                   help="macroscopic-average period(s) in Angstrom (one or two values) / 宏观平均周期")
    p.add_argument("--grad-tol", type=float, default=0.01,
                   help="max |dV/dz| (eV/Angstrom) for a point to count as flat (default 0.01)")
    p.add_argument("--min-len", type=float, default=1.5,
                   help="minimum plateau length in Angstrom (default 1.5)")
    p.add_argument("--figsize", type=float, nargs=2, default=[6.0, 3.8], metavar=("W", "H"))
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--title", default=None)
    p.add_argument("-o", "--output", default="workfunction.png", help="figure file (default workfunction.png)")
    p.add_argument("--dat", default=None, help="export z, planar average (and macro average) to a text file")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--show", action="store_true", help="open an interactive window")
    return p


def main(argv=None) -> None:
    try:
        sys.stdout.reconfigure(errors="backslashreplace")  # type: ignore[attr-defined]
    except (AttributeError, ValueError):
        pass
    args = build_parser().parse_args(argv)
    axis = AXES[args.axis]
    axis_name = "xyz"[axis] if args.axis in "xyz123" else "abc"[axis]
    vol = read_volumetric(args.file)
    lat, data = vol["lattice"], vol["data"]
    n = data.shape[axis]
    length = plane_spacing(lat, axis)
    dz = length / n
    z = np.arange(n) * dz
    v = planar_average(data, axis)
    vm = macro_average(v, dz, args.macro) if args.macro else None
    print(f"LOCPOT: {args.file}  grid {vol['grid']}  atoms {len(vol['symbols'])}  "
          f"averaging normal to lattice vector {'abc'[axis]}: length {length:.4f} A, {n} points ({dz:.4f} A)")

    folder = os.path.dirname(os.path.abspath(args.file))
    outcar = args.outcar or os.path.join(folder, "OUTCAR")
    info = read_outcar_info(outcar)
    if args.outcar and not os.path.isfile(args.outcar):
        die(f"OUTCAR not found: {args.outcar}")
    efermi = args.efermi if args.efermi is not None else info.get("efermi")
    if efermi is None:
        warn("E-fermi unknown (no OUTCAR and no --efermi): only vacuum levels are reported")
    else:
        print(f"E-fermi = {efermi:.4f} eV ({'command line' if args.efermi is not None else outcar})")
    if "ldipol" in info:
        print(f"OUTCAR: LDIPOL = {info['ldipol']}, IDIPOL = {info.get('idipol', '?')}"
              + (f", dipolmoment = {info['dipole']}" if "dipole" in info else ""))

    # vacuum region from the atomic positions
    fa = vol["frac"][:, axis]
    top, bottom, gap = vacuum_region(fa)
    vac_len = gap * length
    print(f"Slab region: {top * length:.3f} A (top) ... {bottom * length:.3f} A (next bottom, periodic); "
          f"vacuum thickness between atomic planes = {vac_len:.2f} A")
    if vac_len < 10.0:
        warn(f"vacuum thickness {vac_len:.1f} A is small (>= 15 A recommended for work functions); "
             "is this really a slab along this axis?")
    i_start = int(math.ceil(top * n))
    i_stop = int(math.floor(bottom * n))
    vv = v  # plateaus are detected on the planar average (the vacuum needs no macroscopic averaging)
    plateaus = find_plateaus(z, vv, i_start, i_stop, length, args.grad_tol, args.min_len)
    plateaus.sort(key=lambda p: -p["len"])
    plateaus = sorted(plateaus[:2], key=lambda p: p["i0"])  # upper side first (closest to the top surface)

    print("-" * 72)
    results = []
    if not plateaus:
        idx = np.arange(i_start, i_stop + 1) % n
        mid = idx[len(idx) // 4: 3 * len(idx) // 4] if len(idx) > 8 else idx
        slope = float(np.polyfit(np.arange(len(mid)) * dz, vv[mid], 1)[0]) if len(mid) > 2 else float("nan")
        vmax = float(vv[idx].max())
        print(f"No flat vacuum plateau found (|dV/dz| < {args.grad_tol} eV/A over >= {args.min_len} A).")
        print(f"Maximum potential in the vacuum region: {vmax:.4f} eV; slope in mid-vacuum: {slope:.4f} eV/A")
        if efermi is not None:
            print(f"Rough work function (max - E_F): {vmax - efermi:.4f} eV  (NOT reliable)")
        warn("vacuum potential is not flat: use LDIPOL = .TRUE. + IDIPOL = 3 for asymmetric slabs, a thicker "
             "vacuum, or relax --grad-tol / --min-len")
    else:
        names = ["upper side (above top surface)", "lower side (below bottom surface)"]
        if len(plateaus) == 1:
            names = ["vacuum"]
        for name, p in zip(names, plateaus):
            line = (f"{name:34s}: V_vac = {p['mean']:9.4f} eV  (z = {p['z0'] % length:.2f}-{p['z1'] % length:.2f} A, "
                    f"std {p['std'] * 1000:.1f} meV, slope {p['slope'] * 1000:.2f} meV/A)")
            print(line)
            if efermi is not None:
                print(f"{'':34s}  work function Phi = {p['mean'] - efermi:.4f} eV")
            results.append(p["mean"])
            if abs(p["slope"]) > 0.005:
                warn(f"plateau '{name}' still has a slope of {p['slope'] * 1000:.1f} meV/A")
        if len(plateaus) == 2:
            dv = plateaus[0]["mean"] - plateaus[1]["mean"]
            print(f"Difference of the two vacuum levels (upper - lower) = {dv:.4f} eV")
            if abs(dv) < 0.02:
                print("  -> both sides are equivalent within 20 meV (symmetric slab or no dipole).")
            else:
                print("  -> expected for an asymmetric slab with dipole correction (LDIPOL = .TRUE.); "
                      "for a symmetric slab this indicates an unconverged / asymmetric structure.")
                if info.get("ldipol", "").upper().startswith("F"):
                    warn("two different vacuum levels but LDIPOL = F in OUTCAR: the step comes from the "
                         "periodic boundary - rerun with LDIPOL = .TRUE., IDIPOL = 3")
        elif len(plateaus) == 1 and info.get("ldipol", "").upper().startswith("T"):
            print("Only one plateau found although LDIPOL = T: the dipole step may lie inside it "
                  "(try a smaller --min-len) or the slab is symmetric.")
    if "vac_outcar" in info:
        a, b = info["vac_outcar"]
        print(f"OUTCAR (LVACPOTAV, VASP >= 6.4.3): vacuum level upper = {a:.4f} eV, lower = {b:.4f} eV"
              + (f" -> Phi = {a - efermi:.4f} / {b - efermi:.4f} eV" if efermi is not None else ""))
    if vm is not None:
        zc = ((top + bottom - 1.0) / 2.0) % 1.0  # middle of the slab
        ic = int(round(zc * n)) % n
        print(f"Macroscopic average at the slab centre (z = {zc * length:.2f} A): {vm[ic]:.4f} eV")
    print("-" * 72)

    if args.dat:
        cols = [z, v] + ([vm] if vm is not None else [])
        with open(args.dat, "w") as fh:
            fh.write(f"# planar average of {args.file} along {axis_name}; E-fermi = "
                     f"{'unknown' if efermi is None else f'{efermi:.6f}'} eV\n")
            fh.write("z(A) V_planar(eV)" + (" V_macro(eV)" if vm is not None else "") + "\n")
            np.savetxt(fh, np.column_stack(cols), fmt="%14.6f")
        print(f"Data written: {args.dat}")
    if not args.no_plot:
        plt = setup_matplotlib(args.show)
        make_plot(args, z, v, vm, plateaus, efermi, np.sort(fa) * length, axis_name, plt)


if __name__ == "__main__":
    main()
