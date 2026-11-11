#!/usr/bin/env python3
"""Read the elastic tensor written by VASP (IBRION=6, ISIF>=3) and derive
polycrystalline moduli.

What it does
  1. parses the LAST block "TOTAL ELASTIC MODULI (kBar)" in OUTCAR
     (or "SYMMETRIZED ELASTIC MODULI (kBar)" = clamped-ion with --clamped);
  2. converts kBar -> GPa (divide by 10);
  3. reorders VASP's XX YY ZZ XY YZ ZX into Voigt order xx yy zz yz zx xy;
  4. prints C_ij, its eigenvalues (all > 0 <=> Born stable), and the
     Voigt / Reuss / Hill bulk and shear moduli, Young's modulus,
     Poisson ratio, Pugh ratio B/G, universal anisotropy A^U;
  5. optionally (--poscar CONTCAR) the density, sound velocities and the
     Debye temperature (Anderson 1963).

Usage
  python voigt_vrh.py                      # reads ./OUTCAR
  python voigt_vrh.py -f run/OUTCAR --poscar run/CONTCAR
  python voigt_vrh.py --clamped            # clamped-ion (rigid-ion) tensor

Only numpy is required.  Tested with VASP 6.x OUTCAR layout.
"""
import argparse
import re
import sys

import numpy as np

# VASP row/column order: XX YY ZZ XY YZ ZX  ->  Voigt: xx yy zz yz zx xy
VASP_TO_VOIGT = [0, 1, 2, 4, 5, 3]

# atomic masses (g/mol) for the density / Debye temperature estimate
MASS = {
    "H": 1.008, "He": 4.003, "Li": 6.94, "Be": 9.012, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180, "Na": 22.990,
    "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Sc": 44.956,
    "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845,
    "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723,
    "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906,
    "Mo": 95.95, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87,
    "Cd": 112.41, "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Te": 127.60,
    "I": 126.90, "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Os": 190.23,
    "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Tl": 204.38,
    "Pb": 207.2, "Bi": 208.98,
}


def read_moduli(outcar, header):
    """Return the last 6x6 block (kBar, VASP order) following `header`."""
    with open(outcar, errors="replace") as fh:
        lines = fh.readlines()
    starts = [i for i, l in enumerate(lines) if header in l]
    if not starts:
        sys.exit(f"'{header}' not found in {outcar}. "
                 "Was the run IBRION=6 with ISIF>=3, and did it finish?")
    rows = []
    for line in lines[starts[-1] + 1:]:
        m = re.match(r"\s*(XX|YY|ZZ|XY|YZ|ZX)\s+(.*)", line)
        if m:
            rows.append([float(x) for x in m.group(2).split()[:6]])
            if len(rows) == 6:
                break
    if len(rows) != 6:
        sys.exit("could not read 6 rows of the elastic tensor")
    return np.array(rows)


def to_voigt_gpa(c_vasp_kbar):
    c = c_vasp_kbar[np.ix_(VASP_TO_VOIGT, VASP_TO_VOIGT)] / 10.0
    return 0.5 * (c + c.T)  # remove tiny numerical asymmetry


def vrh(c):
    s = np.linalg.inv(c)
    bv = (c[0, 0] + c[1, 1] + c[2, 2] + 2 * (c[0, 1] + c[1, 2] + c[0, 2])) / 9
    gv = (c[0, 0] + c[1, 1] + c[2, 2] - (c[0, 1] + c[1, 2] + c[0, 2])
          + 3 * (c[3, 3] + c[4, 4] + c[5, 5])) / 15
    br = 1 / (s[0, 0] + s[1, 1] + s[2, 2] + 2 * (s[0, 1] + s[1, 2] + s[0, 2]))
    gr = 15 / (4 * (s[0, 0] + s[1, 1] + s[2, 2])
               - 4 * (s[0, 1] + s[1, 2] + s[0, 2])
               + 3 * (s[3, 3] + s[4, 4] + s[5, 5]))
    return bv, br, gv, gr


def read_poscar(path):
    """Minimal VASP5 POSCAR/CONTCAR reader: returns volume (A^3) and symbols."""
    with open(path) as fh:
        lines = [l.split() for l in fh.readlines()]
    scale = float(lines[1][0])
    lat = np.array([[float(x) for x in lines[i][:3]] for i in (2, 3, 4)])
    vol = abs(np.linalg.det(lat))
    vol = vol * scale ** 3 if scale > 0 else abs(scale)
    species = lines[5]
    counts = [int(x) for x in lines[6]]
    symbols = [s for s, n in zip(species, counts) for _ in range(n)]
    return vol, symbols


def debye(b, g, vol, symbols):
    """Anderson's estimate. b, g in GPa, vol in A^3."""
    mass_g = sum(MASS[s.split("_")[0]] for s in symbols) / 6.02214076e23
    rho = mass_g / (vol * 1e-24) * 1e3                 # kg/m^3
    vl = np.sqrt((b + 4 * g / 3) * 1e9 / rho)          # m/s
    vt = np.sqrt(g * 1e9 / rho)
    vm = ((2 / vt ** 3 + 1 / vl ** 3) / 3) ** (-1 / 3)
    n = len(symbols) / (vol * 1e-30)                   # atoms per m^3
    h, kb = 6.62607015e-34, 1.380649e-23
    theta = h / kb * (3 * n / (4 * np.pi)) ** (1 / 3) * vm
    return rho, vl, vt, vm, theta


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--outcar", default="OUTCAR")
    ap.add_argument("--clamped", action="store_true",
                    help="use SYMMETRIZED ELASTIC MODULI (clamped/rigid ions)")
    ap.add_argument("--poscar", help="CONTCAR/POSCAR for density and Debye T")
    args = ap.parse_args()

    header = ("SYMMETRIZED ELASTIC MODULI (kBar)" if args.clamped
              else "TOTAL ELASTIC MODULI (kBar)")
    c = to_voigt_gpa(read_moduli(args.outcar, header))
    np.set_printoptions(precision=2, suppress=True, linewidth=120)
    print(f"# source: {header}  (converted to GPa, Voigt order 1..6 = xx yy zz yz zx xy)")
    print(c)

    eig = np.linalg.eigvalsh(c)
    print("\neigenvalues of C (GPa):", eig)
    print("Born stability (all eigenvalues > 0):", "YES" if np.all(eig > 0) else "NO")
    if not np.all(eig > 0):
        print("C is not positive definite -> skip averages")
        return

    bv, br, gv, gr = vrh(c)
    b, g = (bv + br) / 2, (gv + gr) / 2
    e = 9 * b * g / (3 * b + g)
    nu = (3 * b - 2 * g) / (2 * (3 * b + g))
    au = 5 * gv / gr + bv / br - 6
    print(f"\n{'':>18}{'Voigt':>10}{'Reuss':>10}{'Hill':>10}")
    print(f"{'B (GPa)':>18}{bv:10.2f}{br:10.2f}{b:10.2f}")
    print(f"{'G (GPa)':>18}{gv:10.2f}{gr:10.2f}{g:10.2f}")
    print(f"Young's modulus E (Hill) = {e:.2f} GPa")
    print(f"Poisson ratio nu (Hill)  = {nu:.3f}")
    print(f"Pugh ratio B/G           = {b / g:.2f}  (> 1.75: ductile tendency)")
    print(f"Universal anisotropy A^U = {au:.3f}  (0 for isotropic)")
    if abs(c[0, 0] - c[1, 1]) < 1e-3 * c[0, 0] and abs(c[0, 0] - c[2, 2]) < 1e-3 * c[0, 0]:
        print(f"Cubic-like: Cauchy pressure C12-C44 = {c[0, 1] - c[3, 3]:.2f} GPa, "
              f"Zener A = 2C44/(C11-C12) = {2 * c[3, 3] / (c[0, 0] - c[0, 1]):.3f}")

    if args.poscar:
        vol, symbols = read_poscar(args.poscar)
        rho, vl, vt, vm, theta = debye(b, g, vol, symbols)
        print(f"\ndensity = {rho / 1000:.4f} g/cm^3")
        print(f"v_l = {vl:.0f} m/s, v_t = {vt:.0f} m/s, v_m = {vm:.0f} m/s")
        print(f"Debye temperature (Anderson) = {theta:.1f} K")


if __name__ == "__main__":
    main()
