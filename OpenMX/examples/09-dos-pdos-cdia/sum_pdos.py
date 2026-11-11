#!/usr/bin/env python3
"""Sum OpenMX DosMain PDOS files by element (or species) and angular momentum.

DosMain writes, for every selected atom N:
    PREFIX.PDOS.<Method>.atomN          (sum over all orbitals of atom N)
    PREFIX.PDOS.<Method>.atomN.s1       (s, all radial functions summed)
    PREFIX.PDOS.<Method>.atomN.p1..p3   (px, py, pz)
    PREFIX.PDOS.<Method>.atomN.d1..d5   (d3z^2-r^2, dx^2-y^2, dxy, dxz, dyz)
    PREFIX.PDOS.<Method>.atomN.f1..f7
with <Method> = Tetrahedron or Gaussian. Columns:
    non-spin / non-collinear : E(eV)  DOS  integrated
    collinear spin-polarized : E(eV)  DOS_up  -DOS_down  int_up  int_down
(In the OpenMX 4.0 DosMain source the down-spin DOS column is written with a
minus sign; this script keeps that convention in its output.)

Usage:
    python3 sum_pdos.py cdia_dos --dat Cdia_dos.dat
    python3 sum_pdos.py nio_afm --dat NiO_AFM_pdos.dat --by atom --resolve-m --plot
Requires numpy; matplotlib only for --plot.
"""
import argparse
import glob
import os
import re
from collections import OrderedDict

import numpy as np

TRAPZ = getattr(np, "trapezoid", None) or getattr(np, "trapz")

M_NAMES = {
    "s": ["s"],
    "p": ["px", "py", "pz"],
    "d": ["d3z2-r2", "dx2-y2", "dxy", "dxz", "dyz"],
    "f": ["f5z3-3zr2", "f5xz2-xr2", "f5yz2-yr2", "fzx2-zy2", "fxyz", "fx3-3xy2", "f3yx2-y3"],
}


def strip(line):
    return line.split("#", 1)[0].strip()


def block(lines, name):
    out, inside = [], False
    for raw in lines:
        s = strip(raw)
        if s.lower().startswith("<" + name.lower()):
            inside = True
            continue
        if inside and s.lower().startswith(name.lower() + ">"):
            break
        if inside and s:
            out.append(s.split())
    return out


def atom_groups(datfile, by):
    """Map atom index -> element or species name using the OpenMX input file."""
    lines = open(datfile).read().splitlines()
    pao = {p[0]: p[1] for p in block(lines, "Definition.of.Atomic.Species")}
    groups = {}
    for p in block(lines, "Atoms.SpeciesAndCoordinates"):
        idx, spe = int(p[0]), p[1]
        if by == "species":
            groups[idx] = spe
        elif by == "atom":
            m = re.match(r"([A-Z][a-z]?)", pao.get(spe, spe))
            groups[idx] = "%s%d" % (m.group(1) if m else spe, idx)
        else:
            m = re.match(r"([A-Z][a-z]?)", pao.get(spe, spe))
            groups[idx] = m.group(1) if m else spe
    return groups


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="System.Name used for the DOS run")
    ap.add_argument("--dat", default="", help="OpenMX input file to map atoms to elements/species")
    ap.add_argument("--by", choices=["element", "species", "atom"], default="element",
                    help="group atoms by element (default), species name, or keep atoms separate")
    ap.add_argument("--method", default="", help="Tetrahedron or Gaussian (auto-detected)")
    ap.add_argument("--resolve-m", action="store_true", help="also write m-resolved sums (e.g. d orbitals)")
    ap.add_argument("--plot", action="store_true", help="save PREFIX_pdos.png")
    ap.add_argument("--emin", type=float, default=-15.0)
    ap.add_argument("--emax", type=float, default=10.0)
    args = ap.parse_args()

    methods = [args.method] if args.method else ["Tetrahedron", "Gaussian"]
    pat = re.compile(r"\.atom(\d+)\.([spdfg])(\d+)$")
    files, method = [], None
    for m in methods:
        files = [f for f in glob.glob("%s.PDOS.%s.atom*.*" % (args.prefix, m)) if pat.search(f)]
        if files:
            method = m
            break
    if not files:
        raise SystemExit("no %s.PDOS.*.atomN.<l><m> files found; run DosMain (PDOS) first" % args.prefix)

    groups = atom_groups(args.dat, args.by) if args.dat else {}
    if args.by == "atom" and not args.dat:
        groups = {}
    sums, msums, energy, spin = OrderedDict(), OrderedDict(), None, False
    for f in sorted(files, key=lambda x: (int(pat.search(x).group(1)), x)):
        atom, l, m = pat.search(f).groups()
        atom, m = int(atom), int(m)
        d = np.loadtxt(f, ndmin=2)
        if energy is None:
            energy, spin = d[:, 0], d.shape[1] >= 5
        dos = d[:, 1:3] if spin else d[:, 1:2]
        g = groups.get(atom, "atom%d" % atom)
        sums[(g, l)] = sums.get((g, l), 0.0) + dos
        if args.resolve_m:
            name = M_NAMES.get(l, [l + str(i) for i in range(1, 10)])[m - 1]
            msums[(g, name)] = msums.get((g, name), 0.0) + dos

    ef_idx = energy <= 0.0
    print("method: %s, %d files, spin-polarized: %s" % (method, len(files), spin))
    header = "E(eV) DOS_up -DOS_down" if spin else "E(eV) DOS"
    for table, tag in ((sums, "l"), (msums, "m")):
        for (g, l), dos in table.items():
            out = "%s.pdos_%s_%s.dat" % (args.prefix, g, l)
            np.savetxt(out, np.column_stack([energy, dos]), fmt="%12.6f", header=header)
            # electrons below E_F: integrate from the lower end of Dos.Erange to 0 eV
            n_occ = [TRAPZ(np.abs(dos[ef_idx, i]), energy[ef_idx]) for i in range(dos.shape[1])]
            occ = " + ".join("%.3f" % x for x in n_occ)
            print("%-28s -> %-34s N(E<E_F) = %s" % ("%s %s" % (g, l), out, occ))

    tot_file = "%s.DOS.%s" % (args.prefix, method)
    if os.path.exists(tot_file):
        t = np.loadtxt(tot_file, ndmin=2)
        tot = t[:, 1:3] if spin else t[:, 1:2]
        if len(tot) == len(energy):
            psum = sum(v for (g, l), v in sums.items())
            err = np.max(np.abs(psum - tot))
            print("max |sum(PDOS) - DOS| = %.4f states/eV (small only if PDOS was done for ALL atoms)" % err)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4))
        for (g, l), dos in sums.items():
            line, = ax.plot(energy, dos[:, 0], lw=1.0, label="%s-%s" % (g, l))
            if spin:
                ax.plot(energy, dos[:, 1], lw=1.0, color=line.get_color())
        ax.axvline(0.0, color="gray", ls="--", lw=0.8)
        if spin:
            ax.axhline(0.0, color="k", lw=0.5)
        ax.set_xlim(args.emin, args.emax)
        ax.set_xlabel("E - E$_F$ (eV)")
        ax.set_ylabel("PDOS (states/eV)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(args.prefix + "_pdos.png", dpi=200)
        print("figure written to %s_pdos.png" % args.prefix)


if __name__ == "__main__":
    main()
