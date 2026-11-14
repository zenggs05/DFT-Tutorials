#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
optics_post.py -- post-processing of VASP dielectric / optical calculations.
                  VASP 介电与光学后处理：从 vasprun.xml 读取 eps(omega)（区分 density-density 与
                  current-current），计算吸收系数、折射率、消光系数、反射率、能量损失谱并画 Tauc 图；
                  从 OUTCAR / vasprun.xml 读取 LEPSILON / LCALCEPS 的静态介电张量与 Born 有效电荷。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Sub-commands
------------
optics   frequency-dependent dielectric function of a LOPTICS run (also GW/BSE blocks).
         Source: the <dielectricfunction> blocks of vasprun.xml (VASP 6 labels them
         comment="density-density" / "current-current"; unlabelled VASP 5 blocks are taken in
         the order density, current) or, as a fallback for a truncated vasprun.xml, the
         "frequency dependent IMAGINARY/REAL DIELECTRIC FUNCTION" blocks of OUTCAR.
         From eps = eps1 + i eps2 (isotropic average of xx, yy, zz by default):
             n     = sqrt((|eps| + eps1) / 2)            refractive index
             k     = sqrt((|eps| - eps1) / 2)            extinction coefficient
             alpha = 2 E k / (hbar c)                    absorption coefficient (cm^-1)
             R     = ((n-1)^2 + k^2) / ((n+1)^2 + k^2)   normal-incidence reflectivity
             L     = eps2 / (eps1^2 + eps2^2)            energy-loss function -Im(1/eps)
             A     = E L_z eps2 / (hbar c)               absorbance of a thin layer (--slab-length)
         Tauc plot (alpha E)^(1/r) with r = 1/2 (direct allowed) or 2 (indirect allowed);
         --fit EMIN EMAX extrapolates the linear part to the energy axis.
static   static dielectric tensors (every "... DIELECTRIC TENSOR ..." block of OUTCAR, e.g.
         eps_inf with/without local fields, the IONIC CONTRIBUTION of IBRION = 5-8) and the
         Born effective charges of LEPSILON / LCALCEPS runs; prints eps0 = eps_inf + eps_ion,
         averages and principal values, Z* per atom and the charge-neutrality sum rule.
         vasprun.xml is accepted too (varrays named epsilon*, array born_charges).
list     list the dielectric-function blocks and static tensors found in a file.

Usage examples
--------------
    python optics_post.py optics                              # ./vasprun.xml, density-density
    python optics_post.py optics run/vasprun.xml --kind current
    python optics_post.py optics --component zz --emax 8
    python optics_post.py optics --tauc direct --fit 3.0 3.6
    python optics_post.py optics --component inplane --slab-length 20   # 2D: absorbance
    python optics_post.py optics OUTCAR                       # fallback for a truncated vasprun.xml
    python optics_post.py static                              # ./OUTCAR of LEPSILON / LCALCEPS
    python optics_post.py static ionic/OUTCAR --json eps.json
    python optics_post.py list vasprun.xml

Outputs (optics): <prefix>.dat (E, wavelength, eps1, eps2, n, k, alpha, R, L[, A]),
<prefix>_tensor.dat (all 6 components), <prefix>.png and, with --tauc, <prefix>_tauc.dat/png.

Notes: alpha, n, k and R are bulk quantities; for slabs / 2D materials they depend on the
vacuum thickness (use --slab-length for the absorbance).  The independent-particle spectra
contain neither excitons nor phonon-assisted (indirect) transitions.

Dependencies: python >= 3.9, numpy, matplotlib (only for plots).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

__version__ = "1.0.0"

HBAR_C_EV_CM = 1.973269804e-5          # hbar*c in eV*cm
HC_EV_NM = 1239.84198                  # h*c in eV*nm
COMPONENTS = ["xx", "yy", "zz", "xy", "yz", "zx"]
FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?")
# light-mode categorical slots 1-3 and ink colours (validated reference palette)
SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"


def die(msg: str, code: int = 2):
    print(f"ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def floats(text: str) -> list[float]:
    return [float(x.replace("D", "E").replace("d", "e")) for x in FLOAT_RE.findall(text)]


def is_xml(path: Path) -> bool:
    if path.suffix.lower() == ".xml":
        return True
    try:
        with open(path, "rb") as fh:
            return fh.read(64).lstrip().startswith(b"<?xml")
    except OSError:
        return False


# ---------------------------------------------------------------------------
# readers
# ---------------------------------------------------------------------------
def read_vasprun(path: Path) -> dict:
    """Dielectric functions, static tensors, Born charges and atom names of vasprun.xml."""
    out = {"diel": [], "static": {}, "born": None, "atoms": None, "truncated": False}
    unlabelled = 0
    big = {"projected", "projected_kpoints_opt", "eigenvalues", "eigenvalues_kpoints_opt", "dos"}
    try:
        for _event, elem in ET.iterparse(str(path), events=("end",)):
            tag = elem.tag
            if tag == "dielectricfunction":
                comment = elem.attrib.get("comment", "").strip()
                if not comment:
                    label = ["density", "current"][min(unlabelled, 1)]
                    unlabelled += 1
                    comment = f"(unlabelled, assumed {label})"
                elif "density" in comment:
                    label = "density"
                elif "current" in comment:
                    label = "current"
                else:
                    label = comment
                tabs = {}
                for part in ("imag", "real"):
                    p = elem.find(part)
                    if p is None:
                        continue
                    rows = [floats(r.text or "") for r in p.iter("r")]
                    rows = [r for r in rows if len(r) >= 7]
                    if rows:
                        tabs[part] = np.array([r[:7] for r in rows])
                if "imag" in tabs and "real" in tabs and len(tabs["imag"]) == len(tabs["real"]):
                    out["diel"].append({"label": label, "comment": comment, "energy": tabs["imag"][:, 0],
                                        "re": tabs["real"][:, 1:7], "im": tabs["imag"][:, 1:7]})
                elem.clear()
            elif tag == "varray" and elem.attrib.get("name", "").startswith("epsilon"):
                rows = [floats(v.text or "") for v in elem.findall("v")]
                if len(rows) == 3 and all(len(r) >= 3 for r in rows):
                    out["static"][f"vasprun.xml varray '{elem.attrib['name']}'"] = np.array([r[:3] for r in rows])
            elif tag == "array" and elem.attrib.get("name") == "born_charges":
                born = []
                for s in elem.findall("set"):
                    rows = [floats(v.text or "") for v in s.findall("v")]
                    if len(rows) == 3 and all(len(r) >= 3 for r in rows):
                        born.append([r[:3] for r in rows])
                if born:
                    out["born"] = np.array(born)
            elif tag == "array" and elem.attrib.get("name") == "atoms":
                names = []
                for rc in elem.iter("rc"):
                    cs = rc.findall("c")
                    if cs:
                        names.append((cs[0].text or "X").strip())
                if names:
                    out["atoms"] = names
            elif tag in big:
                elem.clear()
    except ET.ParseError as exc:
        out["truncated"] = True
        warn(f"{path} is incomplete or malformed ({exc}); using the data read before that point")
    return out


def read_outcar_optics(path: Path) -> list[dict]:
    """'frequency dependent IMAGINARY/REAL DIELECTRIC FUNCTION' blocks of OUTCAR."""
    blocks = []          # (part, label, rows)
    count = {"IMAGINARY": 0, "REAL": 0}
    with open(path, errors="replace") as fh:
        lines = fh.readlines()
    i, nl = 0, len(lines)
    while i < nl:
        line = lines[i]
        i += 1
        if "frequency dependent" in line and "DIELECTRIC FUNCTION" in line:
            part = "IMAGINARY" if "IMAGINARY" in line else ("REAL" if "REAL" in line else None)
            if part is None:
                continue
            if "density-density" in line:
                label = "density"
            elif "current-current" in line:
                label = "current"
            else:
                label = ["density", "current"][min(count[part], 1)]
            count[part] += 1
            rows = []
            while i < nl:
                nxt = lines[i]
                vals = floats(nxt) if not re.search(r"[A-DF-Za-df-z]", nxt) else []
                if len(vals) == 7:
                    rows.append(vals)
                elif rows or "DIELECTRIC" in nxt:
                    break          # end of the data block (the line is examined again above)
                i += 1
            if rows:
                blocks.append((part, label, " ".join(line.split()), np.array(rows)))
    out = []
    for lab in ("density", "current"):
        im = [b for b in blocks if b[0] == "IMAGINARY" and b[1] == lab]
        re_ = [b for b in blocks if b[0] == "REAL" and b[1] == lab]
        if im and re_ and len(im[-1][3]) == len(re_[-1][3]):
            out.append({"label": lab, "comment": im[-1][2], "energy": im[-1][3][:, 0],
                        "re": re_[-1][3][:, 1:7], "im": im[-1][3][:, 1:7]})
    return out


def read_outcar_static(path: Path) -> dict:
    """Static dielectric tensors, Born charges and atom names of an OUTCAR."""
    out = {"static": {}, "born": None, "atoms": None, "nions": None}
    titles, per_type = [], None
    with open(path, errors="replace") as fh:
        lines = fh.readlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if out["nions"] is None and "NIONS" in line:
            m = re.search(r"NIONS\s*=\s*(\d+)", line)
            if m:
                out["nions"] = int(m.group(1))
        if "TITEL" in line:
            t = line.split("=", 1)[-1].split()
            if len(t) >= 2:
                titles.append(t[1].split("_")[0])
        if per_type is None and "ions per type" in line:
            per_type = [int(x) for x in line.split("=", 1)[1].split()]
        if "DIELECTRIC TENSOR" in line and "frequency dependent" not in line:
            head = " ".join(line.split())
            rows, j = [], i + 1
            while j < n and len(rows) < 3 and j < i + 8:
                s = lines[j].strip()
                if s and not s.startswith("-"):
                    v = floats(s) if not re.search(r"[A-DF-Za-df-z]", s) else []
                    if len(v) >= 3:
                        rows.append(v[:3])
                    else:
                        break
                j += 1
            if len(rows) == 3:
                out["static"][head] = np.array(rows)     # keep the last occurrence
            i = j
            continue
        if "BORN EFFECTIVE CHARGES" in line:
            born, cur, j = [], None, i + 1
            while j < n:
                s = lines[j].strip()
                if s.startswith("---"):
                    if born:
                        break
                elif re.match(r"^ion\s+\d+", s):
                    cur = []
                    born.append(cur)
                elif cur is not None and re.match(r"^[123]\s", s):
                    v = floats(s)
                    if len(v) >= 4:
                        cur.append(v[1:4])
                elif s and born:
                    break
                j += 1
            if born and all(len(b) == 3 for b in born):
                out["born"] = np.array(born)
            i = j
            continue
        i += 1
    if titles and per_type and len(titles) >= len(per_type):
        names = []
        for el, cnt in zip(titles[:len(per_type)], per_type):
            names += [el] * cnt
        out["atoms"] = names
    return out


# ---------------------------------------------------------------------------
# optics
# ---------------------------------------------------------------------------
def optical_constants(energy, e1, e2):
    mod = np.sqrt(e1 ** 2 + e2 ** 2)
    n = np.sqrt(np.clip((mod + e1) / 2.0, 0.0, None))
    k = np.sqrt(np.clip((mod - e1) / 2.0, 0.0, None))
    alpha = 2.0 * energy * k / HBAR_C_EV_CM
    refl = ((n - 1.0) ** 2 + k ** 2) / ((n + 1.0) ** 2 + k ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        loss = np.where(mod > 0, e2 / mod ** 2, 0.0)
    return n, k, alpha, refl, loss


def select_component(re6, im6, comp):
    if comp == "avg":
        return re6[:, :3].mean(axis=1), im6[:, :3].mean(axis=1)
    if comp == "inplane":
        return re6[:, :2].mean(axis=1), im6[:, :2].mean(axis=1)
    i = COMPONENTS.index(comp)
    return re6[:, i], im6[:, i]


def load_dielectric(path: Path) -> list[dict]:
    if is_xml(path):
        blocks = read_vasprun(path)["diel"]
        if not blocks:
            die(f"no complete <dielectricfunction> block in {path}: was LOPTICS = .TRUE. set and did VASP "
                "finish? (for a truncated vasprun.xml try the OUTCAR)")
    else:
        blocks = read_outcar_optics(path)
        if not blocks:
            die(f"no 'frequency dependent IMAGINARY/REAL DIELECTRIC FUNCTION' blocks in {path}")
    return blocks


def pick_block(blocks, kind: str, block: str | None):
    if block:
        hit = [b for b in blocks if block.lower() in b["comment"].lower()]
        if not hit:
            die(f"no block whose comment contains '{block}' (use the 'list' command)")
        return hit[-1]
    hit = [b for b in blocks if b["label"] == kind]
    if not hit:
        avail = ", ".join(f"'{b['comment']}'" for b in blocks)
        die(f"dielectric function '{kind}' not found; available blocks: {avail} (use --block)")
    return hit[-1]


def style(ax, ylabel):
    ax.set_ylabel(ylabel, color=INK)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2)


def cmd_optics(args) -> int:
    path = Path(args.file)
    if not path.is_file():
        die(f"file '{path}' not found")
    blk = pick_block(load_dielectric(path), args.kind, args.block)
    energy, re6, im6 = blk["energy"], blk["re"], blk["im"]
    if len(energy) < 3:
        die("fewer than 3 frequency points in the dielectric function")
    e1, e2 = select_component(re6, im6, args.component)
    n, k, alpha, refl, loss = optical_constants(energy, e1, e2)
    wavelength = np.full_like(energy, np.inf)
    wavelength[energy > 0] = HC_EV_NM / energy[energy > 0]

    print(f"file          : {path}")
    print(f"block         : {blk['comment']}  ({len(energy)} points, E = {energy[0]:.3f} ... {energy[-1]:.3f} eV)")
    print(f"component     : {args.component}")
    i0 = int(np.argmin(np.abs(energy)))
    diag = "  ".join(f"{c}={re6[i0, j]:.3f}" for j, c in enumerate(COMPONENTS[:3]))
    print(f"eps1(E->0)    : {diag}   (compare with the LEPSILON value without local fields)")
    print(f"n(E->0)       : {n[i0]:.3f}")
    offd = np.max(np.abs(np.concatenate([re6[:, 3:], im6[:, 3:]])))
    if offd > 1e-3 * max(np.max(np.abs(re6[:, :3])), 1e-12):
        print(f"off-diagonal  : max |eps_ij| = {offd:.3g} (non-cubic / low symmetry: inspect {args.prefix}_tensor.dat)")
    if np.min(im6[:, :3]) < -1e-3:
        warn(f"eps2 becomes negative (min {np.min(im6[:, :3]):.3g}): numerical artefact, check CSHIFT/NEDOS")
    onset = energy[np.argmax(alpha >= args.onset)] if np.any(alpha >= args.onset) else None
    if onset is not None:
        print(f"onset         : alpha >= {args.onset:.0e} cm^-1 from E = {onset:.3f} eV (crude absorption-edge marker)")
    if blk["label"] == "density" and np.any(e1 < 0) and energy[np.argmax(e1 < 0)] < 0.5:
        print("note          : eps1 < 0 at low energy - metal? Interband (density-density) spectra lack the Drude term.")

    cols = [energy, wavelength, e1, e2, n, k, alpha, refl, loss]
    head = "E(eV) lambda(nm) eps1 eps2 n k alpha(cm^-1) R loss(-Im(1/eps))"
    if args.slab_length:
        absorb = energy * (args.slab_length * 1e-8) * e2 / HBAR_C_EV_CM
        cols.append(absorb)
        head += f" absorbance(Lz={args.slab_length:g}A)"
    head += f"   [{blk['comment']}, component={args.component}]"
    np.savetxt(f"{args.prefix}.dat", np.column_stack(cols), fmt="%14.6e", header=head)
    np.savetxt(f"{args.prefix}_tensor.dat", np.column_stack([energy, re6, im6]), fmt="%14.6e",
               header="E(eV) " + " ".join(f"re_{c}" for c in COMPONENTS) + " " + " ".join(f"im_{c}" for c in COMPONENTS))
    print(f"wrote {args.prefix}.dat, {args.prefix}_tensor.dat")

    tauc = None
    if args.tauc:
        power = 2.0 if args.tauc == "direct" else 0.5
        y = (alpha * energy) ** power
        tauc = {"power": power, "y": y, "eg": None}
        if args.fit:
            lo, hi = sorted(args.fit)
            w = (energy >= lo) & (energy <= hi)
            if w.sum() < 3:
                die("the --fit window contains fewer than 3 points")
            slope, icpt = np.polyfit(energy[w], y[w], 1)
            if slope <= 0:
                warn("the fitted slope is not positive - choose the linear rising part of the Tauc curve")
            else:
                tauc.update(eg=-icpt / slope, slope=slope, icpt=icpt, window=(lo, hi))
                print(f"Tauc ({args.tauc}): extrapolated gap = {tauc['eg']:.3f} eV (fit {lo:g}-{hi:g} eV)")
        if args.tauc == "indirect":
            warn("independent-particle spectra contain only direct (vertical) transitions; an 'indirect' "
                 "Tauc plot of them has no physical meaning")
        np.savetxt(f"{args.prefix}_tauc.dat", np.column_stack([energy, y]), fmt="%14.6e",
                   header=f"E(eV) (alpha*E)^{power:g}  [alpha in cm^-1, E in eV]")
        print(f"wrote {args.prefix}_tauc.dat")

    if not args.no_plot:
        plot_optics(args, blk, energy, re6, im6, e1, e2, n, k, alpha, refl, loss, tauc)
    return 0


def plot_optics(args, blk, energy, re6, im6, e1, e2, n, k, alpha, refl, loss, tauc):
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m = (energy >= args.emin) & (energy <= args.emax)
    if m.sum() < 2:
        die(f"no data between --emin {args.emin} and --emax {args.emax} eV")
    x = energy[m]
    fig, ax = plt.subplots(2, 2, figsize=(9.6, 7.0), sharex=True)
    lw = 2.0
    # (a) dielectric function: diagonal components or the selected one
    a = ax[0, 0]
    if args.component == "avg":
        for j, c in enumerate(COMPONENTS[:3]):
            a.plot(x, im6[m, j], color=SERIES[j], lw=lw, label=rf"$\varepsilon_2^{{{c}}}$")
        for j, c in enumerate(COMPONENTS[:3]):
            a.plot(x, re6[m, j], color=SERIES[j], lw=lw, ls="--", label=rf"$\varepsilon_1^{{{c}}}$")
    else:
        a.plot(x, e2[m], color=SERIES[0], lw=lw, label=r"$\varepsilon_2$")
        a.plot(x, e1[m], color=SERIES[1], lw=lw, label=r"$\varepsilon_1$")
    a.axhline(0, color=INK2, lw=0.8)
    style(a, "dielectric function")
    a.legend(frameon=False, fontsize=8, ncol=2)
    a.set_title(f"{blk['label']} ({args.component})", color=INK, fontsize=10, loc="left")
    # (b) n and k
    a = ax[0, 1]
    a.plot(x, n[m], color=SERIES[0], lw=lw, label="n")
    a.plot(x, k[m], color=SERIES[1], lw=lw, label="k")
    style(a, "n, k")
    a.legend(frameon=False, fontsize=8)
    # (c) absorption coefficient
    a = ax[1, 0]
    a.semilogy(x, np.clip(alpha[m], 1e2, None), color=SERIES[0], lw=lw)
    style(a, r"$\alpha$ (cm$^{-1}$)")
    a.set_ylim(bottom=1e2)
    # (d) reflectivity and loss function on one scale (L normalised to its maximum, no second axis)
    a = ax[1, 1]
    a.plot(x, refl[m], color=SERIES[0], lw=lw, label="R")
    lmax = max(float(np.max(loss[m])), 1e-12)
    a.plot(x, loss[m] / lmax, color=SERIES[1], lw=lw, label=f"L / {lmax:.3g}")
    style(a, "R, L (normalised)")
    a.legend(frameon=False, fontsize=8)
    for a in ax[1]:
        a.set_xlabel("photon energy (eV)", color=INK)
        a.set_xlim(args.emin, args.emax)
    fig.tight_layout()
    fig.savefig(f"{args.prefix}.png", dpi=args.dpi)
    print(f"wrote {args.prefix}.png")

    if tauc is not None:
        fig2, a = plt.subplots(figsize=(5.2, 4.0))
        a.plot(x, tauc["y"][m], color=SERIES[0], lw=lw, label="Tauc curve")
        if tauc.get("eg") is not None:
            xs = np.linspace(tauc["eg"], tauc["window"][1], 50)
            a.plot(xs, tauc["slope"] * xs + tauc["icpt"], color=SERIES[1], lw=lw, ls="--",
                   label=f"fit: Eg = {tauc['eg']:.2f} eV")
        style(a, rf"$(\alpha E)^{{{tauc['power']:g}}}$")
        a.set_xlabel("photon energy (eV)", color=INK)
        a.set_xlim(args.emin, args.emax)
        a.set_ylim(bottom=0)
        a.legend(frameon=False, fontsize=8)
        fig2.tight_layout()
        fig2.savefig(f"{args.prefix}_tauc.png", dpi=args.dpi)
        print(f"wrote {args.prefix}_tauc.png")
    if args.show:
        plt.show()
    plt.close("all")


# ---------------------------------------------------------------------------
# static tensors and Born charges
# ---------------------------------------------------------------------------
def fmt_tensor(t, indent="    "):
    return "\n".join(indent + " ".join(f"{v:12.6f}" for v in row) for row in t)


def tensor_summary(t):
    sym = 0.5 * (t + t.T)
    ev = np.linalg.eigvalsh(sym)
    return f"average (trace/3) = {np.trace(t) / 3:.4f}; principal values = " + ", ".join(f"{v:.4f}" for v in ev)


def choose_eps(static: dict):
    """(eps_inf key, eps_ion key) among the OUTCAR headings."""
    keys = list(static)
    ion = [k for k in keys if "IONIC" in k.upper()]
    elec = [k for k in keys if "MACROSCOPIC" in k.upper() and "IONIC" not in k.upper()]
    pref = [k for k in elec if "including local field effects" in k and "RPA" not in k]
    if pref:
        e = pref[-1]
    elif elec:
        e = elec[-1]
    else:
        e = next((k for k in keys if k.endswith("'epsilon'")), None)      # vasprun.xml
    if not ion:
        ion = [k for k in keys if k.endswith("'epsilon_ion'")]
    return e, (ion[-1] if ion else None)


def cmd_static(args) -> int:
    path = Path(args.file)
    if not path.is_file():
        die(f"file '{path}' not found")
    data = read_vasprun(path) if is_xml(path) else read_outcar_static(path)
    static, born = data["static"], data["born"]
    if not static and born is None:
        die(f"no static dielectric tensor and no Born charges in {path} (LEPSILON / LCALCEPS run finished?)")
    report = {"file": str(path), "tensors": {}, "born": None}
    print(f"file: {path}\n")
    for key, t in static.items():
        print(key)
        print(fmt_tensor(t))
        print("    " + tensor_summary(t) + "\n")
        report["tensors"][key] = t.tolist()
    ekey, ikey = choose_eps(static)
    if ekey:
        print(f"eps_inf  <- {ekey}")
        report["eps_inf"] = static[ekey].tolist()
        if ikey:
            eps0 = static[ekey] + static[ikey]
            print(f"eps_ion  <- {ikey}")
            print("eps_0 = eps_inf + eps_ion:")
            print(fmt_tensor(eps0))
            print("    " + tensor_summary(eps0))
            report["eps_0"] = eps0.tolist()
        print()
    if born is not None:
        names = data.get("atoms")
        if not names or len(names) != len(born):
            names = [f"atom{i + 1}" for i in range(len(born))]
        print(f"Born effective charges Z* (e), {len(born)} atoms; rows as printed by VASP")
        for i, (nm, z) in enumerate(zip(names, born), 1):
            print(f"  {i:4d} {nm:<3s}  Z*_iso = {np.trace(z) / 3:8.4f}")
            print(fmt_tensor(z, indent="         "))
        s = born.sum(axis=0)
        dev = float(np.max(np.abs(s)))
        print("\n  sum over atoms (charge neutrality, should be ~0):")
        print(fmt_tensor(s, indent="         "))
        print(f"  max |sum| = {dev:.4f}" + ("   -> large: check the k-point convergence" if dev > args.asr_tol else ""))
        report["born"] = {"atoms": names, "Z": born.tolist(), "sum": s.tolist()}
        print("\n  For phonopy NAC write BORN with: phonopy-vasp-born [--pa ...] > BORN")
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=1)
        print(f"\nwrote {args.json}")
    return 0


def cmd_list(args) -> int:
    path = Path(args.file)
    if not path.is_file():
        die(f"file '{path}' not found")
    if is_xml(path):
        data = read_vasprun(path)
        blocks = data["diel"]
    else:
        blocks = read_outcar_optics(path)
        data = read_outcar_static(path)
    print(f"file: {path}")
    print(f"dielectric-function blocks: {len(blocks)}")
    for b in blocks:
        print(f"  kind={b['label']:<8s} points={len(b['energy']):<6d} comment: {b['comment']}")
    print(f"static tensors: {len(data['static'])}")
    for key in data["static"]:
        print(f"  {key}")
    print(f"Born charges: {'yes (' + str(len(data['born'])) + ' atoms)' if data['born'] is not None else 'no'}")
    return 0


# ---------------------------------------------------------------------------
# command line
# ---------------------------------------------------------------------------
def build_parser():
    ap = argparse.ArgumentParser(
        prog="optics_post.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=("VASP dielectric/optical post-processing: eps(omega) -> n, k, alpha, R, loss, Tauc plot; "
                     "static dielectric tensors and Born effective charges.\n"
                     "VASP 介电与光学后处理：光学常数、Tauc 图、静态介电张量与 Born 有效电荷。"),
        epilog=("Examples / 示例:\n"
                "  python optics_post.py optics\n"
                "  python optics_post.py optics vasprun.xml --kind current --component zz\n"
                "  python optics_post.py optics --tauc direct --fit 3.0 3.6\n"
                "  python optics_post.py static OUTCAR --json eps.json\n"
                "  python optics_post.py list vasprun.xml\n"
                "Use 'optics_post.py COMMAND -h' for the options of each command.\n"
                "Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License"))
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="command", metavar="COMMAND")

    po = sub.add_parser("optics", help="eps(omega), optical constants and Tauc plot of a LOPTICS run",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    po.add_argument("file", nargs="?", default="vasprun.xml", help="vasprun.xml (or OUTCAR as fallback)")
    po.add_argument("--kind", choices=("density", "current"), default="density",
                    help="density-density or current-current response")
    po.add_argument("--block", help="select another <dielectricfunction> block by a substring of its "
                                    "comment (see 'list'); overrides --kind")
    po.add_argument("--component", default="avg", choices=["avg", "inplane"] + COMPONENTS,
                    help="avg = (xx+yy+zz)/3, inplane = (xx+yy)/2, or one Cartesian component")
    po.add_argument("--emin", type=float, default=0.0, help="lower energy of the plots (eV)")
    po.add_argument("--emax", type=float, default=10.0, help="upper energy of the plots (eV)")
    po.add_argument("--tauc", choices=("direct", "indirect"),
                    help="Tauc plot (alpha E)^2 (direct) or (alpha E)^(1/2) (indirect)")
    po.add_argument("--fit", nargs=2, type=float, metavar=("EMIN", "EMAX"),
                    help="linear fit window of the Tauc curve (eV); prints the extrapolated gap")
    po.add_argument("--slab-length", type=float, metavar="LZ",
                    help="cell length normal to a slab/2D layer in Angstrom: adds the thin-layer absorbance "
                         "A = E*Lz*eps2/(hbar c) (weak-absorption, normal-incidence limit)")
    po.add_argument("--onset", type=float, default=1e4, help="alpha threshold (cm^-1) for the reported onset")
    po.add_argument("-o", "--prefix", default="optics", help="prefix of the output files")
    po.add_argument("--no-plot", action="store_true", help="data files only")
    po.add_argument("--show", action="store_true", help="show the figures interactively")
    po.add_argument("--dpi", type=int, default=200, help="figure resolution")

    ps = sub.add_parser("static", help="static dielectric tensors and Born charges (LEPSILON / LCALCEPS)",
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ps.add_argument("file", nargs="?", default="OUTCAR", help="OUTCAR (or vasprun.xml)")
    ps.add_argument("--asr-tol", type=float, default=0.1,
                    help="flag max |sum of Z*| above this value (charge-neutrality sum rule)")
    ps.add_argument("--json", help="write tensors and Born charges to this JSON file")

    pl = sub.add_parser("list", help="list the dielectric blocks / tensors found in a file")
    pl.add_argument("file", nargs="?", default="vasprun.xml", help="vasprun.xml or OUTCAR")
    return ap


def main(argv=None) -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except (AttributeError, ValueError):
        pass
    ap = build_parser()
    args = ap.parse_args(argv)
    if args.command is None:
        ap.print_help()
        return 1
    try:
        if args.command == "optics":
            return cmd_optics(args)
        if args.command == "static":
            return cmd_static(args)
        return cmd_list(args)
    except OSError as exc:
        die(str(exc))
    return 1


if __name__ == "__main__":
    sys.exit(main())
