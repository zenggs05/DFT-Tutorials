#!/usr/bin/env python3
"""plot_band_dos_matplotlib.py - publication-style band structure + (P)DOS with matplotlib.

Reads VASP vasprun.xml files with pymatgen and draws a two-panel figure
(bands on the left, DOS on the right, shared energy axis). Handles
spin-polarised runs (spin up solid, spin down dashed / mirrored DOS),
optional element-projected "fat bands" and element/orbital PDOS.

Typical use (directory layout of this example):
    python plot_band_dos_matplotlib.py --band band --dos dos
    python plot_band_dos_matplotlib.py --band band --dos dos --project Ga,As \
        --pdos-orbitals --ylim -6 4 --formats pdf png --out gaas_band_dos
    python plot_band_dos_matplotlib.py --band band --zero vbm      # bands only

Energy reference (--zero):
    efermi  E-fermi of the DOS run (or --efermi value); default
    vbm     valence-band maximum of the band run (semiconductors)
The E-fermi printed by a non-SCF line-mode run is NOT reliable; this script
never uses it unless you pass it explicitly with --efermi.

Requirements: pymatgen (2024+), numpy, matplotlib.
The band run must be a line-mode run (KPOINTS with labels); the DOS run
should use LORBIT = 11 for projections. For HSE zero-weight runs use
sumo-bandplot or VASPKIT 252 instead.
"""
import argparse
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

try:
    from pymatgen.io.vasp import Vasprun
    from pymatgen.electronic_structure.core import Spin
except ImportError:  # pragma: no cover
    sys.exit("pymatgen is required: pip install pymatgen")

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_STYLE = os.path.join(HERE, "dft_paper.mplstyle")
GAMMA_CHAR = "\u0393"  # Greek capital gamma (escape keeps the file ASCII-only)
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--band", help="directory of the line-mode band run (vasprun.xml + KPOINTS)")
    p.add_argument("--dos", help="directory of the DOS run (vasprun.xml)")
    p.add_argument("--efermi", type=float, default=None,
                   help="reference Fermi energy in eV (default: E-fermi of the DOS run)")
    p.add_argument("--zero", choices=["efermi", "vbm"], default="efermi",
                   help="energy zero: Fermi level or valence-band maximum")
    p.add_argument("--ylim", type=float, nargs=2, default=(-6.0, 4.0), metavar=("EMIN", "EMAX"))
    p.add_argument("--project", default=None,
                   help="comma-separated elements for fat bands, e.g. 'Ga,As' (needs LORBIT=11)")
    p.add_argument("--fat-scale", type=float, default=12.0, help="marker size scale for fat bands")
    p.add_argument("--pdos-orbitals", action="store_true",
                   help="split element PDOS into s/p/d/f instead of one curve per element")
    p.add_argument("--smear", type=float, default=0.0,
                   help="extra Gaussian smearing of the DOS in eV (0 = none)")
    p.add_argument("--style", default=DEFAULT_STYLE, help="matplotlib style file")
    p.add_argument("--width", type=float, default=3.35, help="figure width in inches")
    p.add_argument("--height", type=float, default=2.6, help="figure height in inches")
    p.add_argument("--out", default="band_dos", help="output file name without extension")
    p.add_argument("--formats", nargs="+", default=["pdf", "png"],
                   help="output formats, e.g. pdf svg png")
    a = p.parse_args()
    if not a.band and not a.dos:
        p.error("give at least --band or --dos")
    return a


def pretty_label(lab):
    """Convert KPOINTS labels such as GAMMA, \\Gamma, G to a Greek letter."""
    if lab is None:
        return ""
    s = lab.strip().replace("$", "")
    if s.upper() in ("GAMMA", "\\GAMMA", "G", "GM", GAMMA_CHAR):
        return GAMMA_CHAR
    if s.startswith("\\"):
        s = s[1:]
    s = s.replace("_", "")
    return s


def read_bands(band_dir, efermi, need_projections):
    vr = Vasprun(os.path.join(band_dir, "vasprun.xml"), parse_dos=False,
                 parse_projected_eigen=need_projections, parse_potcar_file=False)
    kfile = os.path.join(band_dir, "KPOINTS")
    bs = vr.get_band_structure(kpoints_filename=kfile, line_mode=True,
                               efermi=efermi if efermi is not None else "smart")
    return vr, bs


def band_ticks(bs):
    """Collect (distance, label) pairs, merging discontinuities such as U|K."""
    ticks = []
    for i, kp in enumerate(bs.kpoints):
        if kp.label is None:
            continue
        x, lab = bs.distance[i], pretty_label(kp.label)
        if ticks and abs(ticks[-1][0] - x) < 1e-6:
            if ticks[-1][1] != lab:
                ticks[-1] = (x, ticks[-1][1] + "|" + lab)
            continue
        ticks.append((x, lab))
    return ticks


def plot_bands(ax, bs, zero, project, fat_scale):
    dist = np.asarray(bs.distance)
    spins = [Spin.up] + ([Spin.down] if bs.is_spin_polarized else [])
    for spin in spins:
        energies = np.asarray(bs.bands[spin]) - zero          # (nbands, nkpts)
        ls = "-" if spin == Spin.up else "--"
        for br in bs.branches:
            s, e = br["start_index"], br["end_index"] + 1
            for ib in range(energies.shape[0]):
                ax.plot(dist[s:e], energies[ib, s:e], color="0.15", lw=0.8, ls=ls, zorder=2)
        if project:
            # projections[spin] has shape (nbands, nkpts, norbitals, nions)
            proj = np.asarray(bs.projections[spin])
            symbols = [site.specie.symbol for site in bs.structure]
            for ic, el in enumerate(project):
                idx = [i for i, sym in enumerate(symbols) if sym == el]
                if not idx:
                    print(f"WARNING: element {el} not in structure {set(symbols)}")
                    continue
                w = proj[:, :, :, idx].sum(axis=(2, 3))         # (nbands, nkpts)
                xx = np.broadcast_to(dist, energies.shape)
                ax.scatter(xx.ravel(), energies.ravel(), s=fat_scale * w.ravel(),
                           color=OKABE_ITO[ic % len(OKABE_ITO)], alpha=0.6, lw=0,
                           label=el if spin == Spin.up else None, zorder=3)
    ticks = band_ticks(bs)
    ax.set_xticks([t[0] for t in ticks])
    ax.set_xticklabels([t[1] for t in ticks])
    for x, _ in ticks:
        ax.axvline(x, color="0.6", lw=0.5, zorder=1)
    ax.axhline(0.0, color="0.4", lw=0.5, ls=":", zorder=1)
    ax.set_xlim(dist[0], dist[-1])
    ax.tick_params(axis="x", which="both", top=False, direction="out")
    if project:
        ax.legend(loc="upper right", markerscale=1.5, handletextpad=0.2)


def smear(dos_obj, spin, sigma):
    if sigma > 0:
        return np.asarray(dos_obj.get_smeared_densities(sigma)[spin])
    return np.asarray(dos_obj.densities[spin])


def plot_dos(ax, vr, zero, orbital_split, sigma, horizontal):
    cdos = vr.complete_dos
    energy = np.asarray(cdos.energies) - zero
    spin_pol = Spin.down in cdos.densities
    curves = [("Total", cdos, "0.2")]
    ic = 0
    for el, el_dos in cdos.get_element_dos().items():
        if orbital_split:
            for orb, orb_dos in cdos.get_element_spd_dos(el).items():
                curves.append((f"{el.symbol}-{orb.name}", orb_dos, OKABE_ITO[ic % 7]))
                ic += 1
        else:
            curves.append((el.symbol, el_dos, OKABE_ITO[ic % 7]))
            ic += 1
    for label, d, color in curves:
        lw = 0.9 if label == "Total" else 0.8
        up = smear(d, Spin.up, sigma)
        kw = dict(color=color, lw=lw)
        if horizontal:
            ax.plot(up, energy, label=label, **kw)
        else:
            ax.plot(energy, up, label=label, **kw)
        if spin_pol:
            dn = -smear(d, Spin.down, sigma)
            if horizontal:
                ax.plot(dn, energy, **kw)
            else:
                ax.plot(energy, dn, **kw)
    if horizontal:
        ax.axhline(0.0, color="0.4", lw=0.5, ls=":")
        if spin_pol:
            ax.axvline(0.0, color="0.6", lw=0.4)
        ax.set_xlabel("DOS (states/eV)")
        ax.set_xticks([])
    else:
        ax.axvline(0.0, color="0.4", lw=0.5, ls=":")
        ax.set_ylabel("DOS (states/eV)")
    ax.legend(loc="best", handlelength=1.0)


def main():
    a = parse_args()
    if a.style and os.path.isfile(a.style):
        plt.style.use(a.style)
    project = [s.strip() for s in a.project.split(",")] if a.project else None

    vr_dos = None
    efermi = a.efermi
    if a.dos:
        vr_dos = Vasprun(os.path.join(a.dos, "vasprun.xml"), parse_potcar_file=False)
        if efermi is None:
            efermi = vr_dos.efermi
    if efermi is None and a.zero == "efermi":
        sys.exit("No reference E-fermi: give --dos <dir> or --efermi <value> (from the SCF OUTCAR)")

    bs = None
    if a.band:
        _, bs = read_bands(a.band, efermi, project is not None)

    if a.zero == "vbm" and bs is None:
        sys.exit("--zero vbm needs --band")
    zero = efermi if efermi is not None else bs.efermi
    if a.zero == "vbm":
        if bs.is_metal():
            print("WARNING: the band structure looks metallic; using E-fermi as zero")
        else:
            zero = bs.get_vbm()["energy"]
    if bs is not None and not bs.is_metal():
        gap = bs.get_band_gap()
        print(f"band gap = {gap['energy']:.3f} eV  direct = {gap['direct']}  transition = {gap['transition']}")
    print(f"energy zero = {zero:.4f} eV ({a.zero})")

    ylabel = r"$E - E_\mathrm{F}$ (eV)" if a.zero == "efermi" else r"$E - E_\mathrm{VBM}$ (eV)"
    if bs is not None and vr_dos is not None:
        fig, (axb, axd) = plt.subplots(1, 2, sharey=True, figsize=(a.width, a.height),
                                       gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05})
        plot_bands(axb, bs, zero, project, a.fat_scale)
        plot_dos(axd, vr_dos, zero, a.pdos_orbitals, a.smear, horizontal=True)
        axb.set_ylabel(ylabel)
        axb.set_ylim(*a.ylim)
        axd.tick_params(axis="y", which="both", left=True, right=True, labelleft=False)
    elif bs is not None:
        fig, axb = plt.subplots(figsize=(a.width * 0.8, a.height))
        plot_bands(axb, bs, zero, project, a.fat_scale)
        axb.set_ylabel(ylabel)
        axb.set_ylim(*a.ylim)
    else:
        fig, axd = plt.subplots(figsize=(a.width, a.height * 0.8))
        plot_dos(axd, vr_dos, zero, a.pdos_orbitals, a.smear, horizontal=False)
        axd.set_xlabel(ylabel)
        axd.set_xlim(*a.ylim)

    for fmt in a.formats:
        fn = f"{a.out}.{fmt}"
        fig.savefig(fn)
        print("written", fn)


if __name__ == "__main__":
    main()
