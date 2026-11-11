#!/usr/bin/env python3
"""Extract nearest (J1) and next-nearest (J2) Ni-Ni exchange constants of NiO
from three collinear total energies (FM, AFM-I, AFM-II), and give a rough
mean-field estimate of the Neel temperature.

Heisenberg model used here (classical unit vectors, every pair counted ONCE):

    E = E0 - sum_{<ij>} J_ij e_i . e_j          (J > 0 : ferromagnetic)

With 12 nearest and 6 next-nearest Ni neighbours, the energy per Ni is
    FM     : E0 - 6 J1 - 3 J2
    AFM-I  : E0 + 2 J1 - 3 J2
    AFM-II : E0        + 3 J2
so that
    J1 = (E_AFM1 - E_FM) / 8
    J2 = (E_AFM2 - E_FM) / 6 - J1
J_ij absorbs the spin length (S^2); compare with literature only after
converting conventions (factor 2 for double counting, S(S+1) vs S^2, sign).

Usage:
    python exchange_J.py --fm FM_R/OUTCAR --afm2 AFM2/OUTCAR --afm1 AFM1/OUTCAR \
                         [--fm-check FM_T/OUTCAR] [--nni 2]
"""
import argparse
import re
import sys

KB_MEV = 0.08617333262  # Boltzmann constant in meV/K


def last_energy(outcar):
    """Return the last 'energy(sigma->0)' value (eV) found in an OUTCAR."""
    pattern = re.compile(r"energy\(sigma->0\)\s*=\s*(-?\d+\.\d+)")
    energy = None
    try:
        with open(outcar) as fh:
            for line in fh:
                m = pattern.search(line)
                if m:
                    energy = float(m.group(1))
    except OSError as err:
        sys.exit(f"cannot read {outcar}: {err}")
    if energy is None:
        sys.exit(f"no 'energy(sigma->0)' line in {outcar} - did the run finish?")
    return energy


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fm", required=True, help="OUTCAR of the FM run")
    p.add_argument("--afm1", required=True, help="OUTCAR of the AFM-I run")
    p.add_argument("--afm2", required=True, help="OUTCAR of the AFM-II run")
    p.add_argument("--fm-check", help="optional OUTCAR of FM in the other cell")
    p.add_argument("--nni", type=int, default=2,
                   help="number of Ni atoms per cell (default 2, same for all runs)")
    a = p.parse_args()

    n = a.nni
    e_fm = last_energy(a.fm) / n * 1000.0    # meV per Ni
    e_a1 = last_energy(a.afm1) / n * 1000.0
    e_a2 = last_energy(a.afm2) / n * 1000.0

    print("Energies relative to FM (meV per Ni):")
    print(f"  FM     : {0.0:10.2f}")
    print(f"  AFM-I  : {e_a1 - e_fm:10.2f}")
    print(f"  AFM-II : {e_a2 - e_fm:10.2f}")
    if a.fm_check:
        e_chk = last_energy(a.fm_check) / n * 1000.0
        print(f"  FM (other cell) - FM : {e_chk - e_fm:8.2f}  "
              "(should be ~0; otherwise k-points/ENCUT are not converged)")

    j1 = (e_a1 - e_fm) / 8.0
    j2 = (e_a2 - e_fm) / 6.0 - j1
    print("\nExchange constants (meV, pairs counted once, J>0 = FM):")
    print(f"  J1 = {j1:8.3f}")
    print(f"  J2 = {j2:8.3f}")

    # Mean-field estimate for the AFM-II order: k_B T_N = (1/3) |sum_j J_0j s_j|
    # with s_j = +-1 the relative orientation in AFM-II  ->  k_B T_N = -2 J2
    tn = -2.0 * j2 / KB_MEV
    print(f"\nMean-field T_N (AFM-II) ~ {tn:7.0f} K "
          "(mean-field usually OVERestimates; use Monte Carlo for numbers)")


if __name__ == "__main__":
    main()
