#!/usr/bin/env python3
"""Harmonic-limit vibrational thermochemistry from a CP2K frequency run.

Reads the frequencies of a RUN_TYPE VIBRATIONAL_ANALYSIS job either from the
main output ("VIB|Frequency (cm^-1)" lines) or from the Molden file
(<PROJECT>-VIBRATIONS-1.mol, [FREQ] block) and prints, for temperature T:

    ZPE   = sum_i h*nu_i / 2
    U_vib = sum_i [ h*nu_i/2 + h*nu_i / (exp(x_i) - 1) ]          x_i = h*nu_i / kT
    S_vib = k * sum_i [ x_i / (exp(x_i) - 1) - ln(1 - exp(-x_i)) ]
    F_vib = U_vib - T*S_vib

This is the usual treatment of adsorbates (all modes are vibrations, no
translational/rotational partition functions).  Imaginary modes (printed by
CP2K as negative numbers) are listed and always excluded from the sums; for
a transition state exactly one of them is expected.

Examples:
    python3 harmonic_thermo.py co_on_pt_partial_hessian.out
    python3 harmonic_thermo.py co_pt111_vib-VIBRATIONS-1.mol -T 500 --floor 50
Only the Python standard library is needed.
"""
import argparse
import math
import re
import sys

HC_EV_CM = 1.239841984e-4      # h*c in eV*cm  (E[eV] = HC_EV_CM * nu[cm^-1])
KB_EV = 8.617333262e-5         # Boltzmann constant in eV/K
EV_KJMOL = 96.48533212         # 1 eV in kJ/mol
NUM = r"[-+]?\d+\.\d*(?:[EeDd][-+]?\d+)?"


def freqs_from_output(path):
    """Frequencies of the LAST normal-mode table in a CP2K output."""
    blocks, cur = [], []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if "NORMAL MODES - CARTESIAN DISPLACEMENTS" in line:
                if cur:
                    blocks.append(cur)
                cur = []
            if line.strip().startswith("VIB|Frequency (cm^-1)"):
                cur += [float(x.replace("D", "E")) for x in re.findall(NUM, line.split(")", 1)[1])]
    if cur:
        blocks.append(cur)
    return blocks[-1] if blocks else []


def freqs_from_molden(path):
    out, inside = [], False
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("["):
                inside = s.upper().startswith("[FREQ]")
                continue
            if inside and s:
                out.append(float(s.split()[0]))
    return out


def cp2k_thermo_lines(path):
    """Echo CP2K's own THERMOCHEMISTRY summary, if present (gas-phase RRHO)."""
    keys = ("Zero-point correction", "Entropy [kJ/(mol K)]", "Enthalpy correction",
            "Gibbs energy correction", "Temperature [K]", "Symmetry number")
    lines = []
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                if line.strip().startswith("VIB|") and any(k in line for k in keys):
                    lines.append(line.rstrip())
    except OSError:
        pass
    return lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", help="CP2K output (*.out) or Molden vibration file (*.mol)")
    ap.add_argument("-T", "--temperature", type=float, default=298.15, help="temperature in K (default 298.15)")
    ap.add_argument("--floor", type=float, default=0.0,
                    help="raise real frequencies below this value (cm^-1) to it, e.g. 50; default: off")
    args = ap.parse_args()

    if args.file.lower().endswith((".mol", ".molden")):
        nu = freqs_from_molden(args.file)
    else:
        nu = freqs_from_output(args.file)
    if not nu:
        sys.exit("no frequencies found (wrong file, or the run did not finish)")

    imag = [f for f in nu if f < 0.0]
    real = [f for f in nu if f >= 0.0]
    print(f"# {len(nu)} modes read from {args.file}")
    print("# frequencies (cm^-1): " + " ".join(f"{f:.1f}" for f in nu))
    if imag:
        print(f"# {len(imag)} imaginary mode(s) (excluded): " + " ".join(f"{f:.1f}i" for f in (-x for x in imag)))
    else:
        print("# no imaginary modes")
    if args.floor > 0.0:
        n_low = sum(1 for f in real if f < args.floor)
        real = [max(f, args.floor) for f in real]
        print(f"# {n_low} real mode(s) below {args.floor:g} cm^-1 raised to {args.floor:g} cm^-1")

    T = args.temperature
    kT = KB_EV * T
    zpe = u = s = 0.0
    for f in real:
        e = HC_EV_CM * f
        zpe += 0.5 * e
        u += 0.5 * e
        if f < 1e-6:
            continue
        x = e / kT
        u += e / math.expm1(x)
        s += KB_EV * (x / math.expm1(x) - math.log(-math.expm1(-x)))
    if not real:
        sys.exit("no real modes left")
    ts = T * s
    print(f"# T = {T:.2f} K, harmonic limit, {len(real)} real modes")
    print(f"ZPE        = {zpe:10.4f} eV  = {zpe * EV_KJMOL:9.3f} kJ/mol")
    print(f"U_vib(T)   = {u:10.4f} eV  = {u * EV_KJMOL:9.3f} kJ/mol   (includes ZPE)")
    print(f"S_vib(T)   = {s * 1000:10.4f} meV/K = {s * EV_KJMOL * 1000:8.3f} J/(mol K)")
    print(f"T*S_vib    = {ts:10.4f} eV")
    print(f"F_vib(T)   = {u - ts:10.4f} eV  = {(u - ts) * EV_KJMOL:9.3f} kJ/mol   (U_vib - T*S_vib)")

    cp2k = cp2k_thermo_lines(args.file)
    if cp2k:
        print("# CP2K THERMOCHEMISTRY block (gas-phase RRHO, symmetry number 1):")
        for line in cp2k:
            print("#   " + line.strip())


if __name__ == "__main__":
    main()
