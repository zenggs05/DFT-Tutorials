#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
omx_dos_auto.py -- drive OpenMX DosMain non-interactively and sum the PDOS.
                   非交互地运行 DosMain（四面体/高斯展宽、DOS/PDOS、原子列表），再把
                   *.PDOS.* 按元素（或物种/原子）与 s/p/d/f 求和，处理自旋列，输出数据与图。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

What it does
------------
1. Reads the header of System.Name.Dos.val (mode, Kgrid, Nspin, NonCol,
   atomnum, <WhatSpecies) and feeds DosMain exactly the answers it asks for
   (questions checked against OpenMX 4.0 DosMain.c, function input_main):
       mode 1, Kgrid != 1 1 1 : "Tetrahedron(1), Gaussian Broadening(2)" [+ width]
       mode 1, Kgrid  = 1 1 1 : Gaussian chosen automatically -> width only
       mode 5 (DC/Krylov)     : Gaussian chosen automatically -> width only
       mode 6 (NEGF), 7 (DosGauss.fileout): no method question
       then "Dos(1) or PDos(2)" and, if atomnum > 1, the atom list.
   DosMain is run in the directory of the .Dos.val file; it names its output
   after the part of the file name BEFORE THE FIRST DOT, e.g.
   si.Dos.val -> si.DOS.Tetrahedron, si.PDOS.Tetrahedron.atom1.p2 ...
2. Sums PDOS files  <base>.PDOS.<Method>.atomN.<l><m>  (l = s,p,d,f,g;
   m = 1..2l+1, radial multiplicity already summed by DosMain) by element,
   species or atom and by l (optionally m).  Column layout written by DosMain:
       3 columns : E(eV)  DOS  integrated            (non-spin)
       5 columns : E(eV)  DOS_up  -DOS_down  int_up  int_down
                   (collinear spin, and non-collinear runs where DosMain projects
                    on up/down spinors; the down DOS is written NEGATIVE)
   Energies are relative to the chemical potential (E_F = 0).
3. Writes <base>_pdos_<group>_<l>.dat (or one CSV with --csv) and a figure.

Usage examples
--------------
    python omx_dos_auto.py si                       # tetrahedron, total DOS + PDOS of all atoms
    python omx_dos_auto.py nio_afm --method gauss --width 0.1 --dat NiO_AFM_pdos.dat
    python omx_dos_auto.py mos2 --atoms "1-3" --by atom --resolve-m --emin -8 --emax 5
    python omx_dos_auto.py si --no-run --dat Si.dat          # only sum existing PDOS files
    python omx_dos_auto.py si --dry-run                      # print the DosMain answers only
    DOSMAIN=/opt/openmx4.0/bin/DosMain python omx_dos_auto.py si

Element names: taken from Definition.of.Atomic.Species (VPS name, e.g. Ni_PBE19
-> Ni; the empty-atom VPS "E" -> E) and Atoms.SpeciesAndCoordinates of --dat
(or <base>.dat / <base>.out if they exist).  Without them the species index
from <WhatSpecies of the .Dos.val file is used.

Dependencies: Python >= 3.9, numpy, matplotlib (plot); DosMain from OpenMX
(only when the PDOS/DOS files still have to be generated).
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from collections import OrderedDict

L_ORDER = "spdfg"
M_NAMES = {
    "s": ["s"],
    "p": ["px", "py", "pz"],
    "d": ["d3z2-r2", "dx2-y2", "dxy", "dxz", "dyz"],
    "f": ["f1", "f2", "f3", "f4", "f5", "f6", "f7"],
}


def die(msg):
    sys.exit("omx_dos_auto.py: error: " + msg)


# ----------------------------------------------------------------------------
# .Dos.val header
# ----------------------------------------------------------------------------
def read_dosval_header(fname):
    hdr = {"mode": 1, "noncol": 0, "nspin": 0, "kgrid": None, "atomnum": None, "whatspecies": []}
    try:
        f = open(fname)
    except OSError as exc:
        die("cannot read %s (%s)" % (fname, exc))
    with f:
        lines = iter(f)
        for line in lines:
            p = line.split()
            if not p:
                continue
            key = p[0].lower()
            try:
                if key == "mode":
                    hdr["mode"] = int(p[1])
                elif key == "noncol":
                    hdr["noncol"] = int(p[1])
                elif key == "nspin":
                    hdr["nspin"] = int(p[1])
                elif key == "kgrid":
                    hdr["kgrid"] = [int(v) for v in p[1:4]]
                elif key == "atomnum":
                    hdr["atomnum"] = int(p[1])
                elif key == "<whatspecies":
                    vals = []
                    for l2 in lines:
                        if l2.strip().lower().startswith("whatspecies>"):
                            break
                        vals += [int(v) for v in l2.split()]
                    hdr["whatspecies"] = vals
                elif key.startswith("<eigenvalues") or key.startswith("<an"):
                    break
            except (IndexError, ValueError):
                die("malformed header line in %s: %s" % (fname, line.strip()))
    if hdr["atomnum"] is None:
        die("%s has no 'atomnum' line - is it a .Dos.val file?" % fname)
    return hdr


def parse_atoms(text, natom):
    """'1 2 5-8' or '1,2,5-8' -> sorted unique list within 1..natom."""
    out = []
    for tok in re.split(r"[,\s]+", text.strip()):
        if not tok:
            continue
        m = re.fullmatch(r"(\d+)-(\d+)", tok)
        try:
            vals = range(int(m.group(1)), int(m.group(2)) + 1) if m else [int(tok)]
        except ValueError:
            die("bad atom specification: %s" % tok)
        out.extend(vals)
    out = sorted(set(out))
    bad = [a for a in out if a < 1 or a > natom]
    if bad:
        die("atom index out of range 1..%d: %s" % (natom, bad))
    return out


def dosmain_answers(hdr, todo, method, width, atoms):
    """Lines DosMain reads from stdin (see input_main in DosMain.c)."""
    ans = []
    mode, kg = hdr["mode"], hdr["kgrid"] or [0, 0, 0]
    chosen = method
    if mode == 5:
        ans.append("%g" % width)
        chosen = "gauss"
    elif mode in (6, 7):
        chosen = "negf" if mode == 6 else "gauss-file"
    elif kg == [1, 1, 1]:
        ans.append("%g" % width)
        chosen = "gauss"
    else:
        if method == "gauss":
            ans += ["2", "%g" % width]
        else:
            ans.append("1")
    ans.append(str(todo))
    if todo == 2 and hdr["atomnum"] > 1:
        line = " ".join(str(a) for a in atoms)
        if len(line) > 990:
            die("the atom list is longer than DosMain's 1000-character input buffer; "
                "split it into several runs with --atoms")
        ans.append(line)
    return ans, chosen


def run_dosmain(exe, valfile, vecfile, answers, logname):
    wd = os.path.dirname(os.path.abspath(valfile))
    cmd = [exe, os.path.basename(valfile), os.path.basename(vecfile)]
    print("running: (cd %s && %s) < answers [%s]" % (wd, " ".join(cmd), " | ".join(answers)))
    try:
        proc = subprocess.run(cmd, input="\n".join(answers) + "\n", cwd=wd,
                              capture_output=True, text=True)
    except OSError as exc:
        die("cannot start DosMain '%s' (%s); set --dosmain or $DOSMAIN" % (exe, exc))
    with open(os.path.join(wd, logname), "w") as f:
        f.write(proc.stdout)
        f.write(proc.stderr)
    if proc.returncode != 0:
        tail = "\n".join((proc.stdout + proc.stderr).strip().splitlines()[-8:])
        die("DosMain exited with code %d; last lines:\n%s\n(full log: %s)"
            % (proc.returncode, tail, os.path.join(wd, logname)))
    return proc.stdout


# ----------------------------------------------------------------------------
# atom -> group mapping
# ----------------------------------------------------------------------------
def _block(lines, name):
    out, inside = [], False
    for raw in lines:
        s = raw.split("#", 1)[0].strip()
        low = s.lower()
        if low.startswith("<" + name.lower()):
            if low == "<" + name.lower():
                inside = True
                out = []
            continue
        if inside and low.startswith(name.lower() + ">"):
            return out
        if inside and s:
            out.append(s.split())
    return out


def element_of(spe, pao, vps):
    """Element symbol: VPS name first (Ni_PBE19 -> Ni), empty atoms (VPS 'E') -> 'E'."""
    if vps == "E":
        return "E"
    for cand in (vps, pao, spe):
        m = re.match(r"([A-Z][a-z]?)", cand or "")
        if m:
            return m.group(1)
    return spe


def atom_groups(source, by):
    try:
        lines = open(source, errors="replace").read().splitlines()
    except OSError as exc:
        die("cannot read %s (%s)" % (source, exc))
    spec = {}
    for p in _block(lines, "Definition.of.Atomic.Species"):
        spec[p[0]] = (p[1] if len(p) > 1 else "", p[2] if len(p) > 2 else "")
    coords = _block(lines, "Atoms.SpeciesAndCoordinates")
    groups = {}
    for p in coords:
        try:
            idx = int(p[0])
        except (ValueError, IndexError):
            continue
        spe = p[1]
        el = element_of(spe, *spec.get(spe, ("", "")))
        groups[idx] = spe if by == "species" else ("%s%d" % (el, idx) if by == "atom" else el)
    return groups


def find_group_source(args, base, wd):
    if args.dat:
        return args.dat
    for cand in (os.path.join(wd, base + ".dat"), os.path.join(wd, base + ".out")):
        if os.path.exists(cand):
            return cand
    return ""


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        prog="omx_dos_auto.py",
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="System.Name (path allowed); expects <prefix>.Dos.val/.Dos.vec")
    ap.add_argument("--dosmain", default=os.environ.get("DOSMAIN", "DosMain"),
                    help="DosMain executable (default $DOSMAIN or 'DosMain' in PATH)")
    ap.add_argument("--method", choices=["tetra", "gauss"], default="tetra",
                    help="broadening when DosMain asks (mode 1, Kgrid != 1 1 1); default tetra")
    ap.add_argument("--width", type=float, default=0.1,
                    help="Gaussian width a in eV for exp(-(E/a)^2) (default 0.1)")
    ap.add_argument("--what", choices=["both", "dos", "pdos"], default="both",
                    help="run DosMain for total DOS, PDOS or both (default both)")
    ap.add_argument("--atoms", default="", help='atoms for PDOS, e.g. "1 2 5-8" (default: all)')
    ap.add_argument("--no-run", action="store_true", help="do not call DosMain; use existing files")
    ap.add_argument("--dry-run", action="store_true", help="print the DosMain answers and exit")
    ap.add_argument("--dat", default="", help="OpenMX input file used to map atoms to elements")
    ap.add_argument("--by", choices=["element", "species", "atom"], default="element",
                    help="group PDOS by element (default), species name or atom")
    ap.add_argument("--resolve-m", action="store_true", help="also sum m-resolved components (px, dxy ...)")
    ap.add_argument("--emin", type=float, default=-10.0, help="plot range, eV (default -10)")
    ap.add_argument("--emax", type=float, default=5.0, help="plot range, eV (default 5)")
    ap.add_argument("--csv", default="", help="write all summed curves to one CSV instead of .dat files")
    ap.add_argument("-o", "--out", default="", help="figure name (default <base>_pdos.png)")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--show", action="store_true", help="interactive window instead of Agg")
    args = ap.parse_args()

    import numpy as np

    trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    valfile, vecfile = args.prefix + ".Dos.val", args.prefix + ".Dos.vec"
    wd = os.path.dirname(os.path.abspath(valfile))
    base = os.path.basename(args.prefix).split(".", 1)[0]
    if "." in os.path.basename(args.prefix):
        print("note: System.Name contains '.', DosMain will name its output '%s.*'" % base)

    hdr = None
    if os.path.exists(valfile):
        hdr = read_dosval_header(valfile)
        print("%s: mode=%d NonCol=%d Nspin=%d Kgrid=%s atomnum=%d"
              % (valfile, hdr["mode"], hdr["noncol"], hdr["nspin"], hdr["kgrid"], hdr["atomnum"]))
    elif not args.no_run:
        die("%s not found (run OpenMX with Dos.fileout on, or use --no-run)" % valfile)

    if not args.no_run:
        atoms = parse_atoms(args.atoms, hdr["atomnum"]) if args.atoms else list(range(1, hdr["atomnum"] + 1))
        todo = {"dos": [1], "pdos": [2], "both": [1, 2]}[args.what]
        jobs = []
        for t in todo:
            ans, chosen = dosmain_answers(hdr, t, args.method, args.width, atoms)
            jobs.append((t, ans, chosen))
            print("DosMain answers for %s (%s): %s" % ("DOS" if t == 1 else "PDOS", chosen, " | ".join(ans)))
        if args.dry_run:
            return
        if not os.path.exists(vecfile):
            die("%s not found" % vecfile)
        if shutil.which(args.dosmain) is None and not os.path.exists(args.dosmain):
            die("DosMain '%s' not found; compile it (cd openmx4.0/source && make DosMain) and "
                "put it in PATH, or pass --dosmain / set $DOSMAIN" % args.dosmain)
        for t, ans, _ in jobs:
            out = run_dosmain(args.dosmain, valfile, vecfile, ans,
                              "%s.dosmain_%s.log" % (base, "dos" if t == 1 else "pdos"))
            made = re.findall(r"make\s+(\S+)", out)
            if made:
                print("  DosMain wrote: %s%s" % (", ".join(made[:4]), " ..." if len(made) > 4 else ""))
    elif args.dry_run:
        die("--dry-run needs the .Dos.val file and no --no-run")

    # ---------------- collect PDOS files ----------------
    pat = re.compile(r"\.PDOS\.(\w+)\.atom(\d+)\.([spdfg])(\d+)$")
    files = [f for f in glob.glob(os.path.join(wd, base + ".PDOS.*.atom*.*")) if pat.search(f)]
    methods = sorted({pat.search(f).group(1) for f in files})
    if len(methods) > 1:
        newest = max(methods, key=lambda m: max(os.path.getmtime(f) for f in files if ".PDOS.%s." % m in f))
        print("note: PDOS files of several methods found %s; using the newest: %s" % (methods, newest))
        methods = [newest]
    files = [f for f in files if methods and ".PDOS.%s." % methods[0] in f]
    total_files = sorted(glob.glob(os.path.join(wd, base + ".DOS.*")))
    if not files and not total_files:
        die("no %s.PDOS.* or %s.DOS.* files in %s" % (base, base, wd))

    src = find_group_source(args, base, wd)
    groups = atom_groups(src, args.by) if src else {}
    if src and not groups:
        print("warning: no Atoms.SpeciesAndCoordinates found in %s" % src)
    if not groups and hdr and hdr["whatspecies"] and args.by != "atom":
        groups = {i + 1: "spe%d" % s for i, s in enumerate(hdr["whatspecies"])}
        print("note: atoms grouped by species index from <WhatSpecies (give --dat for element names)")

    sums, msums, energy, ncol = OrderedDict(), OrderedDict(), None, None
    for f in sorted(files, key=lambda x: (int(pat.search(x).group(2)),
                                          L_ORDER.index(pat.search(x).group(3)),
                                          int(pat.search(x).group(4)))):
        _, atom, l, m = pat.search(f).groups()
        atom, m = int(atom), int(m)
        try:
            d = np.loadtxt(f, ndmin=2)
        except ValueError as exc:
            die("cannot parse %s (%s)" % (f, exc))
        if d.shape[1] not in (3, 5):
            die("%s has %d columns; expected 3 (non-spin) or 5 (spin)" % (f, d.shape[1]))
        if energy is None:
            energy, ncol = d[:, 0], d.shape[1]
        elif d.shape[1] != ncol or len(d) != len(energy):
            die("%s has a different shape than the other PDOS files" % f)
        dos = d[:, 1:3] if ncol == 5 else d[:, 1:2]
        g = groups.get(atom, "atom%d" % atom)
        sums[(g, l)] = sums.get((g, l), 0.0) + dos
        if args.resolve_m:
            names = M_NAMES.get(l, [])
            name = names[m - 1] if m - 1 < len(names) else "%s%d" % (l, m)
            msums[(g, name)] = msums.get((g, name), 0.0) + dos

    spin = ncol == 5
    tot = None
    if total_files:
        tf = [f for f in total_files if not methods or f.endswith("." + methods[0])] or total_files
        t = np.loadtxt(tf[0], ndmin=2)
        tot = (t[:, 0], t[:, 1:3] if t.shape[1] >= 5 else t[:, 1:2])
        print("total DOS: %s" % tf[0])
        if energy is None:
            energy, spin = tot[0], t.shape[1] >= 5

    if files:
        print("PDOS: %d files, method %s, %s" % (len(files), methods[0],
                                                 "spin-resolved (down column negative)" if spin else "non-spin"))
        below = energy <= 0.0
        rows = []
        for (g, l), dos in list(sums.items()) + list(msums.items()):
            n_occ = [trapz(np.abs(dos[below, i]), energy[below]) for i in range(dos.shape[1])]
            rows.append((g, l, n_occ))
            print("  %-10s %-8s N(E<E_F) = %s" % (g, l, " + ".join("%.3f" % v for v in n_occ)))
        print("  (integrals start at the lower end of Dos.Erange; compare with valence electrons)")
        if tot is not None and len(tot[0]) == len(energy):
            psum = sum(sums.values())
            print("  max |sum(PDOS) - DOS| = %.4f states/eV (small only if all atoms were projected)"
                  % float(np.max(np.abs(psum - tot[1]))))

        header = "E(eV) DOS_up -DOS_down" if spin else "E(eV) DOS"
        if args.csv:
            cols, names = [energy], ["E_eV"]
            for (g, l), dos in list(sums.items()) + list(msums.items()):
                if spin:
                    cols += [dos[:, 0], dos[:, 1]]
                    names += ["%s_%s_up" % (g, l), "%s_%s_dn" % (g, l)]
                else:
                    cols.append(dos[:, 0])
                    names.append("%s_%s" % (g, l))
            np.savetxt(args.csv, np.column_stack(cols), delimiter=",", fmt="%.6f",
                       header=",".join(names), comments="")
            print("summed PDOS written to %s" % args.csv)
        else:
            for (g, l), dos in list(sums.items()) + list(msums.items()):
                out = os.path.join(wd, "%s_pdos_%s_%s.dat" % (base, g, l))
                np.savetxt(out, np.column_stack([energy, dos]), fmt="%12.6f", header=header)
            print("summed PDOS written to %s_pdos_<group>_<l>.dat" % os.path.join(wd, base))

    if args.no_plot:
        return
    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    if tot is not None:
        ax.fill_between(tot[0], tot[1][:, 0], color="0.85", lw=0, label="total")
        if tot[1].shape[1] == 2:
            ax.fill_between(tot[0], tot[1][:, 1], color="0.85", lw=0)
    for (g, l), dos in sums.items():
        line, = ax.plot(energy, dos[:, 0], lw=1.1, label="%s-%s" % (g, l))
        if spin:
            ax.plot(energy, dos[:, 1], lw=1.1, color=line.get_color())
    ax.axvline(0.0, color="0.4", ls="--", lw=0.8)
    if spin:
        ax.axhline(0.0, color="k", lw=0.5)
    ax.set_xlim(args.emin, args.emax)
    sel = (energy >= args.emin) & (energy <= args.emax)
    curves = [d[sel] for d in sums.values()]
    if tot is not None:
        tsel = (tot[0] >= args.emin) & (tot[0] <= args.emax)
        curves.append(tot[1][tsel])
    if curves:
        top = max(float(np.max(c)) for c in curves)
        bot = min(float(np.min(c)) for c in curves)
        pad = 0.05 * max(top - bot, 1e-6)
        ax.set_ylim(bot - pad if spin else 0.0, top + pad)
    ax.set_xlabel(r"$E - E_\mathrm{F}$ (eV)")
    ax.set_ylabel("DOS (states/eV)")
    ax.legend(fontsize=8, frameon=False, ncol=2)
    fig.tight_layout()
    if args.show:
        plt.show()
    else:
        out = args.out or os.path.join(wd, base + "_pdos.png")
        fig.savefig(out, dpi=300)
        print("figure written to %s" % out)


if __name__ == "__main__":
    main()
