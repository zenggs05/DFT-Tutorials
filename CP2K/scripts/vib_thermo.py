#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
vib_thermo.py -- frequencies, IR intensities and thermochemistry from a CP2K
                 RUN_TYPE VIBRATIONAL_ANALYSIS job.
                 解析 CP2K 振动分析输出与 -VIBRATIONS-1.mol（molden），列出频率和 IR 强度，
                 检查虚频，计算 ZPE 以及谐振子近似下的 U、H、S、G（吸附物模式或理想气体分子模式）。

Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

Inputs (one or both, detected by content)
-----------------------------------------
* main output: the LAST "NORMAL MODES - CARTESIAN DISPLACEMENTS" table
    VIB|Frequency (cm^-1)   ...   (negative value = imaginary mode)
    VIB|IR int (KM/Mole)    ...   (only with INTENSITIES T)
    VIB|Red.Masses (a.u.)   ...
  plus the electronic energy "VIB| Minimum Structure - Energy and Forces" /
  "Total energy:" (hartree) and CP2K's own THERMOCHEMICAL DATA block (if
  printed) for comparison.
* <PROJECT>-VIBRATIONS-1.mol (&PRINT/&MOLDEN_VIB): [FREQ], [FR-COORD]
  (bohr) and [INT].  The [INT] numbers are not in km/mol (relative only);
  the geometry is needed for --mode gas.
* --xyz FILE: geometry in Angstrom for --mode gas when no molden file is given.

Models
------
  --mode adsorbate (default): every mode is a vibration (harmonic limit),
      no translation/rotation - the usual treatment of adsorbates and of
      partial-Hessian runs (&MODE_SELECTIVE / fixed atoms):
        ZPE = sum h nu / 2,  U_vib = ZPE + sum h nu / (e^x - 1),
        S_vib = k sum [x/(e^x - 1) - ln(1 - e^-x)],  F_vib = U_vib - T S_vib
  --mode gas: ideal-gas RRHO molecule: + translation (Sackur-Tetrode at
      pressure p), rigid rotor with symmetry number --sigma (linear molecules
      detected from the moments of inertia), electronic S = R ln(2S+1)
      (--mult); H = U + RT, G = H - TS.
  Imaginary modes are listed and always excluded from the sums (CP2K's
  built-in thermochemistry includes them and uses symmetry number 1).
  --floor F raises real modes below F cm^-1 to F; --qrrho uses Grimme's
  quasi-RRHO entropy (Chem. Eur. J. 18, 9955 (2012); w0 = --qrrho-w0,
  B_av = 1e-44 kg m^2) for the vibrational entropy.

Usage examples
--------------
    python vib_thermo.py h2o_vib.out h2o_vib-VIBRATIONS-1.mol --mode gas --sigma 2 -T 298.15
    python vib_thermo.py co_on_pt_partial_hessian.out -T 300 500 --floor 50
    python vib_thermo.py ts_freq_check.out --ts
    python vib_thermo.py co_pt111_vib-VIBRATIONS-1.mol --qrrho -T 298.15 --csv modes.csv

Dependencies: python >= 3.9, numpy.
"""
import argparse
import math
import os
import re
import sys

# CODATA 2018
H = 6.62607015e-34
KB = 1.380649e-23
C_CM = 2.99792458e10
NA = 6.02214076e23
R = KB * NA
AMU = 1.66053906660e-27
BOHR_M = 0.529177210903e-10
HA2KJ = 2625.4996394799
KJ2EV = 1.0 / 96.48533212
NUM = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"

MASSES = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867,
    "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95, "Tc": 98.0,
    "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33, "La": 138.91,
    "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 145.0, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25,
    "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08,
    "Au": 196.97, "Hg": 200.59, "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209.0, "At": 210.0,
    "Rn": 222.0, "D": 2.014,
}


def die(msg):
    sys.exit("vib_thermo.py: error: " + msg)


def need_numpy():
    try:
        import numpy as np
    except ImportError:
        die("numpy is required (pip install numpy)")
    return np


def fl(x):
    return float(x.replace("D", "E").replace("d", "e"))


# --------------------------------------------------------------------------- readers
def is_molden(path):
    with open(path, "r", errors="replace") as fh:
        head = fh.read(4096)
    return "[Molden Format]" in head or "[FREQ]" in head or "[FR-COORD]" in head


def read_output(path):
    """Frequencies, IR (km/mol), reduced masses of the last table; E0; CP2K thermo block."""
    tables, cur = [], None
    e0, pending_min = None, False
    thermo = {}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if "NORMAL MODES - CARTESIAN DISPLACEMENTS" in s:
                cur = {"freq": [], "ir": [], "rmass": []}
                tables.append(cur)
                continue
            if cur is not None:
                for key, tag in (("freq", "VIB|Frequency (cm^-1)"), ("ir", "VIB|IR int (KM/Mole)"),
                                 ("rmass", "VIB|Red.Masses (a.u.)")):
                    if s.startswith(tag):
                        cur[key] += [fl(x) for x in re.findall(NUM, s[len(tag):])]
            if "Minimum Structure - Energy and Forces" in s:
                pending_min = True
                continue
            if pending_min and "Total energy:" in s:
                m = re.search(r"Total energy:\s*(" + NUM + ")", s)
                if m:
                    e0 = fl(m.group(1))
                pending_min = False
            for key, pat in (("T", r"Temperature \[K\]:\s*(" + NUM + ")"),
                             ("p", r"Pressure \[Pa\]:\s*(" + NUM + ")"),
                             ("zpe", r"Zero-point correction \[kJ/mol\]:\s*(" + NUM + ")"),
                             ("S", r"Entropy \[kJ/\(mol K\)\]:\s*(" + NUM + ")"),
                             ("H-U", r"Enthalpy correction \(H-U\) \[kJ/mol\]:\s*(" + NUM + ")"),
                             ("G", r"Gibbs energy correction \[kJ/mol\]:\s*(" + NUM + ")")):
                if s.startswith("VIB|"):
                    m = re.search(pat, s)
                    if m:
                        thermo[key] = fl(m.group(1))
    tables = [t for t in tables if t["freq"]]
    if not tables:
        die(f"{path}: no 'VIB|Frequency (cm^-1)' lines - not a finished VIBRATIONAL_ANALYSIS output?")
    t = tables[-1]
    n = len(t["freq"])
    if t["ir"] and len(t["ir"]) != n:
        print("warning: IR intensity count differs from the frequency count; IR ignored")
        t["ir"] = []
    if t["rmass"] and len(t["rmass"]) != n:
        t["rmass"] = []
    return t, e0, thermo


def read_molden(path):
    sec = None
    freq, ints, syms, coords = [], [], [], []
    atoms_au = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if s.startswith("["):
                sec = s.split("]")[0].upper() + "]"
                unit_au = "AU" in s.upper()
                continue
            t = s.split()
            if sec == "[FREQ]":
                freq.append(fl(t[0]))
            elif sec == "[INT]":
                ints.append(fl(t[0]))
            elif sec == "[FR-COORD]" and len(t) >= 4:
                syms.append(t[0])
                coords.append([fl(t[1]), fl(t[2]), fl(t[3])])
            elif sec == "[ATOMS]" and len(t) >= 6:
                f = 1.0 if unit_au else 1.0 / 0.529177210903
                atoms_au.append((t[0], [fl(t[3]) * f, fl(t[4]) * f, fl(t[5]) * f]))
    if not freq:
        die(f"{path}: no [FREQ] block")
    if not coords and atoms_au:
        syms = [a[0] for a in atoms_au]
        coords = [a[1] for a in atoms_au]
    if ints and len(ints) != len(freq):
        ints = []
    return freq, ints, syms, coords          # coords in bohr


def read_xyz(path):
    with open(path, "r", errors="replace") as fh:
        lines = fh.readlines()
    try:
        na = int(lines[0].split()[0])
    except (IndexError, ValueError):
        die(f"{path}: not an xyz file")
    syms, xyz = [], []
    for ln in lines[2:2 + na]:
        t = ln.split()
        if len(t) < 4:
            die(f"{path}: malformed atom line: {ln.strip()[:60]}")
        syms.append(t[0])
        xyz.append([fl(t[1]) / 0.529177210903, fl(t[2]) / 0.529177210903, fl(t[3]) / 0.529177210903])
    if len(syms) != na:
        die(f"{path}: expected {na} atoms, found {len(syms)}")
    return syms, xyz                         # bohr


def element(sym):
    m = re.match(r"([A-Z][a-z]?)", sym.strip().capitalize() if len(sym) <= 2 else sym)
    return m.group(1) if m else sym


# --------------------------------------------------------------------------- thermo
def vib_terms(nu, T, qrrho=False, w0=100.0):
    """Per-mode ZPE, thermal U (incl. ZPE), S for real frequencies nu (cm^-1); J/mol, J/mol/K."""
    zpe = u = s = 0.0
    for v in nu:
        e = H * C_CM * v                       # J
        x = e / (KB * T)
        zpe += 0.5 * e * NA
        u += (0.5 * e + e / math.expm1(x)) * NA
        s_ho = R * (x / math.expm1(x) - math.log(-math.expm1(-x)))
        if qrrho:
            mu = H / (8 * math.pi ** 2 * C_CM * v)
            bav = 1.0e-44
            mup = mu * bav / (mu + bav)
            s_fr = R * (0.5 + math.log(math.sqrt(8 * math.pi ** 3 * mup * KB * T / H ** 2)))
            w = 1.0 / (1.0 + (w0 / v) ** 4)
            s += w * s_ho + (1 - w) * s_fr
        else:
            s += s_ho
    return zpe, u, s


def inertia(np, syms, coords_bohr, masses):
    m = np.array(masses) * AMU
    r = np.array(coords_bohr) * BOHR_M
    com = (m[:, None] * r).sum(axis=0) / m.sum()
    r = r - com
    Imat = np.zeros((3, 3))
    for mi, ri in zip(m, r):
        Imat += mi * (np.dot(ri, ri) * np.eye(3) - np.outer(ri, ri))
    return np.sort(np.linalg.eigvalsh(Imat)), m.sum()


# --------------------------------------------------------------------------- main
def build_parser():
    ap = argparse.ArgumentParser(
        prog="vib_thermo.py",
        description="Frequencies, IR intensities, imaginary-mode check and harmonic thermochemistry "
                    "(ZPE, U, H, S, G) from a CP2K VIBRATIONAL_ANALYSIS run. CP2K 振动分析与热化学。",
        epilog="Examples:\n"
               "  python vib_thermo.py h2o_vib.out h2o_vib-VIBRATIONS-1.mol --mode gas --sigma 2\n"
               "  python vib_thermo.py co_on_pt_partial_hessian.out -T 300 500 --floor 50\n"
               "  python vib_thermo.py ts_freq_check.out --ts\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="CP2K main output and/or <PROJECT>-VIBRATIONS-1.mol")
    ap.add_argument("-T", "--temperature", nargs="+", type=float, default=[298.15], metavar="K",
                    help="temperature(s) in K (default 298.15)")
    ap.add_argument("-p", "--pressure", type=float, default=101325.0, help="pressure in Pa for --mode gas (default 101325)")
    ap.add_argument("--mode", choices=["adsorbate", "gas"], default="adsorbate",
                    help="adsorbate: vibrations only (default); gas: ideal-gas molecule (trans + rot + vib)")
    ap.add_argument("--sigma", type=int, default=1, help="rotational symmetry number for --mode gas (default 1)")
    ap.add_argument("--mult", type=int, default=1, help="spin multiplicity 2S+1 for the electronic entropy (default 1)")
    ap.add_argument("--linear", choices=["auto", "yes", "no"], default="auto",
                    help="linear molecule? (default auto, from the moments of inertia)")
    ap.add_argument("--xyz", metavar="FILE", help="geometry (Angstrom) for --mode gas if no molden file is given")
    ap.add_argument("--mass", nargs="+", metavar="EL=M", default=[], help="override atomic masses, e.g. H=2.014")
    ap.add_argument("--floor", type=float, default=0.0, metavar="CM",
                    help="raise real frequencies below CM cm^-1 to CM (default 0 = off)")
    ap.add_argument("--qrrho", action="store_true", help="Grimme quasi-RRHO vibrational entropy")
    ap.add_argument("--qrrho-w0", type=float, default=100.0, help="quasi-RRHO cut-off w0 in cm^-1 (default 100)")
    ap.add_argument("--soft", type=float, default=50.0,
                    help="flag real modes below this value (cm^-1) as soft (default 50)")
    ap.add_argument("--ts", action="store_true", help="transition state: expect exactly one imaginary mode")
    ap.add_argument("--e0", type=float, metavar="HA",
                    help="electronic energy in hartree (default: 'Minimum Structure' total energy of the output)")
    ap.add_argument("--csv", metavar="FILE", help="write the mode list to CSV")
    return ap


def main():
    args = build_parser().parse_args()
    np = need_numpy()
    out = mol = None
    for f in args.files:
        if not os.path.isfile(f):
            die(f"file not found: {f}")
        if is_molden(f):
            mol = f
        else:
            out = f
    if any(t <= 0 for t in args.temperature):
        die("temperatures must be positive")
    if args.sigma < 1 or args.mult < 1:
        die("--sigma and --mult must be >= 1")

    masses = dict(MASSES)
    for item in args.mass:
        try:
            el, val = item.split("=")
            masses[el] = float(val)
        except ValueError:
            die(f"--mass: '{item}' is not of the form EL=MASS")

    freq, ir, rmass, ir_unit, e0, thermo = None, [], [], "", None, {}
    syms, coords = [], []
    if out:
        tab, e0, thermo = read_output(out)
        freq, ir, rmass = tab["freq"], tab["ir"], tab["rmass"]
        ir_unit = "km/mol" if ir else ""
    if mol:
        mfreq, mint, syms, coords = read_molden(mol)
        if freq is None:
            freq = mfreq
            ir = mint
            ir_unit = "rel. (molden)" if mint else ""
        elif len(mfreq) != len(freq) or max(abs(a - b) for a, b in zip(mfreq, freq)) > 0.5:
            print("warning: frequencies in the molden file differ from the output table; using the output")
    if args.xyz:
        syms, coords = read_xyz(args.xyz)
    if args.e0 is not None:
        e0 = args.e0

    # mode table -------------------------------------------------------------
    imag = [(i, v) for i, v in enumerate(freq) if v < 0]
    soft = [(i, v) for i, v in enumerate(freq) if 0 <= v < args.soft]
    print(f"source: {', '.join(x for x in (out, mol, args.xyz) if x)}")
    print(f"{len(freq)} modes, {len(imag)} imaginary, {len(soft)} soft (< {args.soft:g} cm^-1)")
    head = f"{'mode':>5s} {'freq (cm^-1)':>14s}"
    if ir:
        head += f" {'IR (' + ir_unit + ')':>16s}"
    if rmass:
        head += f" {'red. mass':>10s}"
    print(head)
    for i, v in enumerate(freq):
        row = f"{i + 1:5d} {abs(v):12.2f}{'i' if v < 0 else ' '} "
        if ir:
            row += f" {ir[i]:16.4f}"
        if rmass:
            row += f" {rmass[i]:10.4f}"
        if v < 0:
            row += "   <- imaginary"
        elif v < args.soft:
            row += "   <- soft"
        print(row)
    if args.ts:
        big = [v for _, v in imag if v < -args.soft]
        if len(imag) == 1:
            print(f"TS check: OK - exactly one imaginary mode ({abs(imag[0][1]):.1f}i cm^-1)")
        else:
            print(f"TS check: FAILED - {len(imag)} imaginary modes ({len(big)} larger than {args.soft:g}i)")
    elif imag:
        print("warning: imaginary mode(s) at a supposed minimum - tighten the optimisation/SCF or displace "
              "along the mode and re-optimise; they are excluded below")

    if args.csv:
        with open(args.csv, "w", newline="\n") as fh:
            fh.write("mode,freq_cm-1,imaginary" + (",ir" if ir else "") + (",red_mass" if rmass else "") + "\n")
            for i, v in enumerate(freq):
                fh.write(f"{i + 1},{abs(v):.4f},{int(v < 0)}" + (f",{ir[i]:.6f}" if ir else "")
                         + (f",{rmass[i]:.6f}" if rmass else "") + "\n")
        print(f"wrote {args.csv}")

    real = [v for v in freq if v > 0]
    if args.floor > 0:
        nlow = sum(1 for v in real if v < args.floor)
        real = [max(v, args.floor) for v in real]
        if nlow:
            print(f"--floor: {nlow} mode(s) raised to {args.floor:g} cm^-1")
    if not real:
        die("no real frequencies left")

    # gas-phase pieces ---------------------------------------------------------
    gas = args.mode == "gas"
    if gas:
        if not syms:
            die("--mode gas needs the geometry: give the molden file or --xyz FILE")
        try:
            mlist = [masses[element(s)] if element(s) in masses else masses[s] for s in syms]
        except KeyError as exc:
            die(f"no mass for element {exc}; use --mass EL=M")
        I, mtot = inertia(np, syms, coords, mlist)
        linear = I[0] < 1e-3 * I[2] if args.linear == "auto" else args.linear == "yes"
        nat = len(syms)
        nexp = 3 * nat - (5 if linear else 6)
        if len(freq) != nexp:
            print(f"warning: {len(freq)} modes but a {'linear' if linear else 'non-linear'} molecule with "
                  f"{nat} atoms has {nexp}; partial Hessian or fixed atoms? --mode gas may not be appropriate")
        print(f"gas phase: M = {mtot / AMU:.4f} amu, linear = {linear}, sigma = {args.sigma}, "
              f"mult = {args.mult}, p = {args.pressure:g} Pa; I = " + " ".join(f"{x / AMU / 1e-20:.4f}" for x in I)
              + " amu A^2")

    print()
    print("All energies in kJ/mol (eV in brackets); imaginary modes excluded"
          + ("; quasi-RRHO entropy" if args.qrrho else ""))
    for T in args.temperature:
        zpe, uvib, svib = vib_terms(real, T, args.qrrho, args.qrrho_w0)
        u, s = uvib, svib
        lines = [("ZPE", zpe), ("U_vib (incl. ZPE)", uvib), ("T*S_vib", T * svib)]
        if gas:
            kt = KB * T
            qt = (2 * math.pi * mtot * kt / H ** 2) ** 1.5 * kt / args.pressure
            s_tr = R * (math.log(qt) + 2.5)
            u_tr = 1.5 * R * T
            if linear:
                qr = 8 * math.pi ** 2 * I[2] * kt / (args.sigma * H ** 2)
                s_rot, u_rot = R * (math.log(qr) + 1.0), R * T
            else:
                qr = (math.sqrt(math.pi * I[0] * I[1] * I[2]) / args.sigma
                      * (8 * math.pi ** 2 * kt / H ** 2) ** 1.5)
                s_rot, u_rot = R * (math.log(qr) + 1.5), 1.5 * R * T
            s_el = R * math.log(args.mult)
            u = uvib + u_tr + u_rot
            s = svib + s_tr + s_rot + s_el
            hcorr = u + R * T
            lines += [("U_trans + U_rot", u_tr + u_rot), ("T*S_trans", T * s_tr), ("T*S_rot", T * s_rot),
                      ("T*S_elec", T * s_el), ("U correction", u), ("H correction (U + RT)", hcorr),
                      ("T*S total", T * s), ("G correction (H - TS)", hcorr - T * s)]
            gcorr = hcorr - T * s
        else:
            gcorr = uvib - T * svib
            lines += [("F_vib = U_vib - T*S_vib", gcorr)]
        print(f"--- T = {T:.2f} K")
        for name, val in lines:
            print(f"  {name:26s} {val / 1000:14.4f}   [{val / 1000 * KJ2EV:10.5f} eV]")
        print(f"  {'S total (J/mol/K)':26s} {s:14.4f}")
        if e0 is not None:
            print(f"  E0 = {e0:.8f} Ha;  E0 + ZPE = {e0 + zpe / 1000 / HA2KJ:.8f} Ha;  "
                  f"E0 + {'G' if gas else 'F'}corr = {e0 + gcorr / 1000 / HA2KJ:.8f} Ha")
    if thermo:
        print()
        print("CP2K's own THERMOCHEMICAL DATA block (symmetry number 1, imaginary modes included):")
        for k, lab in (("T", "Temperature [K]"), ("p", "Pressure [Pa]"), ("zpe", "ZPE [kJ/mol]"),
                       ("S", "Entropy [kJ/(mol K)]"), ("H-U", "Enthalpy correction (H-U) [kJ/mol]"),
                       ("G", "Gibbs energy correction [kJ/mol]")):
            if k in thermo:
                print(f"  {lab:36s} {thermo[k]:16.6f}")
    if not gas and syms and len(freq) >= 3 * len(syms) - 6 and len(syms) > 1 and not thermo:
        print("hint: for an isolated molecule use --mode gas --sigma N to include translation and rotation")


if __name__ == "__main__":
    main()
