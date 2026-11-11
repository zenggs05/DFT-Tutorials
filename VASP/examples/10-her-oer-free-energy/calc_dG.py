#!/usr/bin/env python3
"""Free-energy helper for the computational hydrogen electrode (CHE) model.

Sub-commands
  thermo   ZPE, int Cp dT (Cv for adsorbates), -TS and G-correction from the
           frequencies of a VASP IBRION=5/6 OUTCAR.
           * adsorbate (default): harmonic approximation, all 3N_free modes,
             low frequencies raised to --cutoff (default 50 cm^-1).
           * --gas: ideal-gas model through ase.thermochemistry.IdealGasThermo
             (needs ASE and the structure file of the molecule).
  her      Delta G_H* = Delta E_H* + correction for the hydrogen evolution reaction.
  oer      4-electron OER/ORR free-energy steps, limiting potentials,
           overpotentials and (optionally) a free-energy diagram.

Energies can be given either as numbers (eV) or as a VASP run directory /
OUTCAR file; in the latter case the last "energy(sigma->0)" is used.

Examples
  python calc_dG.py thermo freq_OH/OUTCAR -T 298.15
  python calc_dG.py thermo H2_freq/OUTCAR --gas --geometry linear \
         --symmetry 2 --spin 0 --structure H2_freq/POSCAR
  python calc_dG.py her --slab slab --slabH H_fcc --h2 H2_gas
  python calc_dG.py her --slab -95.1 --slabH -98.7 --h2 -6.77 --corr 0.24
  python calc_dG.py oer oer_input_demo.json --U 1.23 --plot oer.png

Requires numpy; matplotlib (plots) and ase (--gas) are optional.
"""
import argparse
import json
import os
import re
import sys

import numpy as np

KB = 8.617333262e-5          # Boltzmann constant, eV/K
CM2EV = 1.239841984e-4       # 1 cm^-1 in eV
E_H2O_FORMATION = 4.92       # 2 H2O -> O2 + 2 H2 at 298 K (2 x 2.46 eV), CHE convention
EQ_POT = 1.23                # equilibrium potential of O2/H2O, V vs RHE

ENERGY_PAT = re.compile(r"energy\(sigma->0\)\s*=\s*([-+0-9.Ee]+)")
FREQ_PAT = re.compile(r"^\s*\d+\s+f(/i)?\s*=\s*[-0-9.Ee+]+\s+THz\s+[-0-9.Ee+]+\s+2PiTHz"
                      r"\s+([-0-9.Ee+]+)\s+cm-1")


# --------------------------------------------------------------------------
# input helpers
# --------------------------------------------------------------------------
def outcar_path(p):
    return os.path.join(p, "OUTCAR") if os.path.isdir(p) else p


def energy_of(value):
    """A float, or a run directory / OUTCAR -> last energy(sigma->0)."""
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(value)
    except ValueError:
        pass
    f = outcar_path(value)
    if not os.path.isfile(f):
        sys.exit("ERROR: '%s' is neither a number nor an existing run/OUTCAR" % value)
    e = None
    with open(f, errors="replace") as fh:
        for line in fh:
            m = ENERGY_PAT.search(line)
            if m:
                e = float(m.group(1))
    if e is None:
        sys.exit("ERROR: no energy(sigma->0) in %s" % f)
    return e


def read_frequencies(path):
    """Return (real_cm, imag_cm) from the FIRST frequency block of an OUTCAR.

    VASP prints the modes twice when selective dynamics / IBRION=5 is used
    ("Eigenvectors and eigenvalues of the dynamical matrix" and
    "Eigenvectors after division by SQRT(mass)"); only the first block is read.
    """
    f = outcar_path(path)
    real, imag = [], []
    inside = False
    with open(f, errors="replace") as fh:
        for line in fh:
            if "Eigenvectors and eigenvalues of the dynamical matrix" in line:
                if inside or real or imag:
                    break
                inside = True
                continue
            if inside and "Eigenvectors after division by SQRT(mass)" in line:
                break
            if inside:
                m = FREQ_PAT.match(line)
                if m:
                    (imag if m.group(1) else real).append(float(m.group(2)))
    if not real and not imag:
        sys.exit("ERROR: no vibrational frequencies found in %s "
                 "(IBRION=5/6 run finished?)" % f)
    return np.array(real), np.array(imag)


# --------------------------------------------------------------------------
# harmonic thermochemistry (adsorbates)
# --------------------------------------------------------------------------
def harmonic(freqs_cm, T):
    """ZPE, U_vib(T)-ZPE, S_vib (eV/K) for a list of real frequencies in cm^-1."""
    eps = np.asarray(freqs_cm) * CM2EV
    zpe = 0.5 * eps.sum()
    if T <= 0:
        return zpe, 0.0, 0.0
    x = eps / (KB * T)
    du = np.sum(eps / np.expm1(x))
    s = KB * np.sum(x / np.expm1(x) - np.log(-np.expm1(-x)))
    return zpe, du, s


def cmd_thermo(a):
    real, imag = read_frequencies(a.outcar)
    print("OUTCAR: %s" % outcar_path(a.outcar))
    print("  real modes: %d   imaginary modes: %d" % (len(real), len(imag)))
    if len(imag):
        print("  imaginary (cm^-1): " + " ".join("%.1f" % v for v in imag))
    T = a.temperature
    if not a.gas:
        used = np.maximum(real, a.cutoff) if a.cutoff > 0 else real
        nraise = int(np.sum(real < a.cutoff)) if a.cutoff > 0 else 0
        if len(imag):
            print("  WARNING: imaginary modes are ignored. A minimum should have none;"
                  " small ones (< ~50 cm^-1) usually mean loose convergence.")
        zpe, du, s = harmonic(used, T)
        g = zpe + du - T * s
        print("  model: harmonic adsorbate, %d modes, %d raised to %.0f cm^-1"
              % (len(used), nraise, a.cutoff))
        print("  T = %.2f K" % T)
        print("  ZPE                 = %8.4f eV" % zpe)
        print("  int Cv dT (0->T)    = %8.4f eV" % du)
        print("  T*S                 = %8.4f eV" % (T * s))
        print("  G_corr = ZPE + int Cv dT - TS = %8.4f eV" % g)
        print("  (add G_corr to the DFT energy of the adsorbate structure)")
        return
    # ideal gas via ASE
    try:
        from ase.io import read
        from ase.thermochemistry import IdealGasThermo
    except ImportError:
        sys.exit("ERROR: --gas needs ASE (pip install ase)")
    if not a.structure:
        sys.exit("ERROR: --gas needs --structure (POSCAR/CONTCAR of the molecule)")
    atoms = read(a.structure)
    atoms.pbc = False          # IdealGasThermo expects an isolated molecule
    vib = real * CM2EV
    thermo = IdealGasThermo(vib_energies=vib, geometry=a.geometry, potentialenergy=0.0,
                            atoms=atoms, symmetrynumber=a.symmetry, spin=a.spin)
    g = thermo.get_gibbs_energy(temperature=T, pressure=a.pressure, verbose=False)
    h = thermo.get_enthalpy(temperature=T, verbose=False)
    s = thermo.get_entropy(temperature=T, pressure=a.pressure, verbose=False)
    zpe = thermo.get_ZPE_correction()
    print("  model: ideal gas (%s, sigma=%d, S=%.1f), p = %.1f Pa, T = %.2f K"
          % (a.geometry, a.symmetry, a.spin, a.pressure, T))
    print("  ZPE                 = %8.4f eV" % zpe)
    print("  int Cp dT (0->T)    = %8.4f eV" % (h - zpe))
    print("  T*S                 = %8.4f eV" % (T * s))
    print("  G_corr = ZPE + int Cp dT - TS = %8.4f eV" % g)


# --------------------------------------------------------------------------
# HER
# --------------------------------------------------------------------------
def cmd_her(a):
    es, esh, eh2 = energy_of(a.slab), energy_of(a.slabH), energy_of(a.h2)
    de = (esh - es - 0.5 * a.nH * eh2) / a.nH
    if a.gcorr_h is not None and a.gcorr_h2 is not None:
        corr = a.gcorr_h - 0.5 * a.gcorr_h2
        src = "computed: G_corr(H*) - 1/2 G_corr(H2)"
    else:
        corr = a.corr
        src = "constant (Norskov et al. 2005 use 0.24 eV)"
    dg = de + corr
    print("E(slab)      = %14.6f eV" % es)
    print("E(slab+nH)   = %14.6f eV   (n = %g)" % (esh, a.nH))
    print("E(H2)        = %14.6f eV" % eh2)
    print("Delta E_H*   = %8.3f eV per H" % de)
    print("correction   = %8.3f eV   [%s]" % (corr, src))
    print("Delta G_H*   = %8.3f eV per H   (ideal HER catalyst: close to 0)" % dg)


# --------------------------------------------------------------------------
# OER / ORR
# --------------------------------------------------------------------------
def cmd_oer(a):
    with open(a.json) as fh:
        d = json.load(fh)
    E = {k: energy_of(v) for k, v in d["E"].items()}
    C = d.get("G_corr", {})
    need = ["slab", "OH", "O", "OOH", "H2O", "H2"]
    miss = [k for k in need if k not in E]
    if miss:
        sys.exit("ERROR: missing energies in json: %s" % miss)
    G = {k: E[k] + float(C.get(k, 0.0)) for k in need}
    h2 = 0.5 * G["H2"]
    dG1 = G["OH"] - G["slab"] - G["H2O"] + h2
    dG2 = G["O"] - G["OH"] + h2
    dG3 = G["OOH"] - G["O"] - G["H2O"] + h2
    dG4 = E_H2O_FORMATION - (dG1 + dG2 + dG3)
    steps = np.array([dG1, dG2, dG3, dG4])
    # adsorption free energies relative to H2O/H2 (common descriptors)
    g_oh = G["OH"] - G["slab"] - (G["H2O"] - h2)
    g_o = G["O"] - G["slab"] - (G["H2O"] - G["H2"])
    g_ooh = G["OOH"] - G["slab"] - (2 * G["H2O"] - 3 * h2)

    names = ["* + H2O -> OH* + (H+ + e-)", "OH* -> O* + (H+ + e-)",
             "O* + H2O -> OOH* + (H+ + e-)", "OOH* -> * + O2 + (H+ + e-)"]
    print("Descriptors (vs H2O/H2, U = 0 V_RHE):")
    print("  dG_OH* = %6.3f   dG_O* = %6.3f   dG_OOH* = %6.3f eV" % (g_oh, g_o, g_ooh))
    print("  dG_OOH* - dG_OH* = %6.3f eV   (scaling relation: ~3.2 eV on many surfaces)"
          % (g_ooh - g_oh))
    print("OER steps at U = %.2f V_RHE:" % a.U)
    for i, (n, v) in enumerate(zip(names, steps), 1):
        print("  dG%d = %6.3f eV   %s" % (i, v - a.U, n))
    ul_oer = steps.max()
    ul_orr = steps.min()
    print("OER: limiting potential U_L = %.3f V, overpotential eta = %.3f V (step %d)"
          % (ul_oer, ul_oer - EQ_POT, int(steps.argmax()) + 1))
    print("ORR: limiting potential U_L = %.3f V, overpotential eta = %.3f V (step %d, reverse)"
          % (ul_orr, EQ_POT - ul_orr, int(steps.argmin()) + 1))

    if a.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not available - no plot")
            return
        labels = ["*", "OH*", "O*", "OOH*", "O$_2$"]
        fig, ax = plt.subplots(figsize=(5.5, 3.8))
        for U, ls, col in ((0.0, "-", "black"), (EQ_POT, "--", "tab:blue"),
                           (ul_oer, ":", "tab:red")):
            levels = np.concatenate([[0.0], np.cumsum(steps - U)])
            for i, g in enumerate(levels):
                ax.plot([i - 0.3, i + 0.3], [g, g], ls, color=col, lw=2,
                        label=("U = %.2f V" % U) if i == 0 else None)
                if i:
                    ax.plot([i - 0.7, i - 0.3], [levels[i - 1], g], ls, lw=0.8,
                            color=col, alpha=0.6)
        ax.set_xticks(range(5))
        ax.set_xticklabels(labels)
        ax.set_ylabel("Free energy (eV)")
        ax.set_title("OER free-energy diagram (CHE)")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=200)
        print("diagram written to", a.plot)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("thermo", help="ZPE / TS / G correction from an IBRION=5 OUTCAR")
    t.add_argument("outcar", help="OUTCAR or run directory")
    t.add_argument("-T", "--temperature", type=float, default=298.15)
    t.add_argument("--cutoff", type=float, default=50.0,
                   help="raise real frequencies below this value (cm^-1); 0 = off")
    t.add_argument("--gas", action="store_true", help="ideal-gas model (needs ASE)")
    t.add_argument("--geometry", default="nonlinear",
                   choices=["monatomic", "linear", "nonlinear"])
    t.add_argument("--symmetry", type=int, default=1, help="rotational symmetry number")
    t.add_argument("--spin", type=float, default=0.0, help="total spin S (O2: 1)")
    t.add_argument("--pressure", type=float, default=101325.0, help="Pa")
    t.add_argument("--structure", default=None, help="POSCAR/CONTCAR of the molecule")
    t.set_defaults(func=cmd_thermo)

    h = sub.add_parser("her", help="Delta G_H* for HER")
    h.add_argument("--slab", required=True)
    h.add_argument("--slabH", required=True)
    h.add_argument("--h2", required=True)
    h.add_argument("--nH", type=float, default=1.0, help="number of H atoms added")
    h.add_argument("--corr", type=float, default=0.24,
                   help="Delta(ZPE) - T Delta S correction if not computed (eV)")
    h.add_argument("--gcorr-h", type=float, default=None, dest="gcorr_h",
                   help="G_corr of adsorbed H from 'thermo' (eV)")
    h.add_argument("--gcorr-h2", type=float, default=None, dest="gcorr_h2",
                   help="G_corr of gas H2 from 'thermo --gas' (eV)")
    h.set_defaults(func=cmd_her)

    o = sub.add_parser("oer", help="4e OER/ORR steps from a json file")
    o.add_argument("json")
    o.add_argument("--U", type=float, default=0.0, help="electrode potential, V vs RHE")
    o.add_argument("--plot", default=None, help="png file for the free-energy diagram")
    o.set_defaults(func=cmd_oer)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
