#!/usr/bin/env python3
"""Planar-averaged electrostatic potential and work function from a VASP LOCPOT.

* reads LOCPOT (VASP >= 5 format, first data block only) and E-fermi from OUTCAR
* averages the potential over the plane spanned by the two other lattice vectors
* locates the vacuum gap (largest atom-free region along the chosen axis),
  detects the step produced by the dipole correction (if any) and averages the
  potential on both sides of it -> upper / lower vacuum levels
* prints  Phi = V_vac - E_F  for both sides, the field left in the vacuum
  (a slope means: missing / unconverged dipole correction or too little vacuum)
* optional --vbm/--cbm: ionisation potential / electron affinity (semiconductors)
* writes planar_avg.dat (z, V, V - E_F) and planar_avg.png

Use LVHAR = .TRUE. (ionic + Hartree potential) for work functions.

Usage:
    python plot_locpot_planar.py LOCPOT --outcar OUTCAR [--axis 3] [--margin 3.0]
"""
import argparse
import re
import sys

import numpy as np


def read_locpot(path):
    with open(path) as fh:
        title = fh.readline()
        scale = float(fh.readline().split()[0])
        cell = np.array([[float(x) for x in fh.readline().split()[:3]] for _ in range(3)])
        if scale < 0:
            scale = (-scale / abs(np.linalg.det(cell))) ** (1.0 / 3.0)
        cell *= scale
        line = fh.readline().split()
        if line[0].isdigit():               # VASP 4 style: no element line
            counts = [int(x) for x in line]
        else:
            counts = [int(x) for x in fh.readline().split()]
        nat = sum(counts)
        mode = fh.readline().strip()
        if mode[:1] in "sS":                # selective dynamics (rare in LOCPOT)
            mode = fh.readline().strip()
        frac = np.array([[float(x) for x in fh.readline().split()[:3]] for _ in range(nat)])
        if mode[:1] in "cCkK":
            frac = np.linalg.solve(cell.T, (frac * scale).T).T
        line = fh.readline()
        while not line.strip():
            line = fh.readline()
        ng = [int(x) for x in line.split()]
        ntot = ng[0] * ng[1] * ng[2]
        data = []
        while len(data) < ntot:
            data.extend(float(x) for x in fh.readline().split())
    pot = np.array(data[:ntot]).reshape(ng, order="F")     # x runs fastest
    return title.strip(), cell, frac % 1.0, pot


def read_efermi(outcar):
    ef = None
    with open(outcar) as fh:
        for line in fh:
            m = re.search(r"E-fermi\s*:\s*(-?\d+\.\d+)", line)
            if m:
                ef = float(m.group(1))
    if ef is None:
        sys.exit(f"no E-fermi in {outcar}")
    return ef


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("locpot", nargs="?", default="LOCPOT")
    p.add_argument("--outcar", default="OUTCAR", help="OUTCAR for E-fermi")
    p.add_argument("--efermi", type=float, help="override the Fermi energy (eV)")
    p.add_argument("--axis", type=int, default=3, choices=(1, 2, 3),
                   help="lattice vector normal to the surface (IDIPOL), default 3")
    p.add_argument("--margin", type=float, default=3.0,
                   help="distance (A) from the outermost atoms excluded from the average")
    p.add_argument("--vbm", type=float, help="valence-band maximum (eV) -> IP")
    p.add_argument("--cbm", type=float, help="conduction-band minimum (eV) -> EA")
    p.add_argument("--no-plot", action="store_true")
    a = p.parse_args()

    title, cell, frac, pot = read_locpot(a.locpot)
    ax = a.axis - 1
    others = tuple(i for i in range(3) if i != ax)
    vz = pot.mean(axis=others)
    n = vz.size
    # length of the axis projected on the surface normal
    normal = np.cross(cell[others[0]], cell[others[1]])
    normal /= np.linalg.norm(normal)
    L = abs(np.dot(cell[ax], normal))
    z = np.arange(n) * L / n
    ef = a.efermi if a.efermi is not None else read_efermi(a.outcar)

    # ---- largest atom-free gap along the axis (periodic) ----
    s = np.sort(frac[:, ax])
    gaps = np.diff(np.append(s, s[0] + 1.0))
    k = int(np.argmax(gaps))
    g0, g1 = s[k], s[k] + gaps[k]                  # fractional, g1 may exceed 1
    m = a.margin / L
    if g1 - g0 <= 2 * m:
        sys.exit(f"vacuum gap {gaps[k] * L:.2f} A is smaller than 2 x margin")
    fz = np.arange(n) / n
    idx = [i for i in range(n) for f in (fz[i], fz[i] + 1.0) if g0 + m <= f <= g1 - m]
    idx = sorted(set(idx), key=lambda i: (fz[i] - g0) % 1.0)   # order from upper side
    seg = vz[idx]

    # ---- dipole step: largest jump between consecutive grid points ----
    jumps = np.abs(np.diff(seg))
    j = int(np.argmax(jumps)) if jumps.size else 0
    dz = L / n
    has_step = jumps.size > 0 and jumps[j] > 10 * np.median(jumps) and jumps[j] > 0.05
    pad = max(1, int(round(1.0 / dz)))             # stay 1 A away from the step
    if has_step:
        upper, lower = seg[: max(j - pad + 1, 1)], seg[j + 1 + pad:]
        if lower.size == 0:
            lower = seg[j + 1:]
    else:
        upper = lower = seg

    def stats(v):
        slope = np.polyfit(np.arange(v.size) * dz, v, 1)[0] if v.size > 2 else 0.0
        return v.mean(), v.std(), slope

    vu, su, fu = stats(upper)
    vl, sl, fl = stats(lower)
    print(f"file            : {a.locpot}  ({title})")
    print(f"axis            : a{a.axis}, length along normal {L:.3f} A, {n} grid points")
    print(f"vacuum gap      : {gaps[k] * L:.2f} A (atoms excluded within {a.margin} A)")
    print(f"dipole step     : {'yes, %.3f eV' % (seg[j + 1] - seg[j]) if has_step else 'not found'}")
    print(f"E_F             : {ef:.4f} eV")
    print(f"V_vac upper side: {vu:.4f} eV (std {su:.4f}, slope {fu:+.4f} eV/A)")
    print(f"V_vac lower side: {vl:.4f} eV (std {sl:.4f}, slope {fl:+.4f} eV/A)")
    print(f"Phi upper       : {vu - ef:.4f} eV")
    print(f"Phi lower       : {vl - ef:.4f} eV")
    print(f"Delta Phi       : {(vu - vl):+.4f} eV (upper - lower)")
    if max(abs(fu), abs(fl)) > 0.01:
        print("WARNING: potential is not flat in the vacuum (|slope| > 0.01 eV/A): "
              "check LDIPOL/IDIPOL/DIPOL, SCF convergence and vacuum thickness")
    for label, e in (("IP (V_vac - VBM)", a.vbm), ("EA (V_vac - CBM)", a.cbm)):
        if e is not None:
            print(f"{label}: upper {vu - e:.4f} eV, lower {vl - e:.4f} eV")

    np.savetxt("planar_avg.dat", np.column_stack([z, vz, vz - ef]),
               header="z(A)  V(eV)  V-E_F(eV)", fmt="%12.5f")
    print("wrote planar_avg.dat")
    if not a.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axp = plt.subplots(figsize=(6, 3.5))
        axp.plot(z, vz - ef, lw=1.2)
        axp.axhline(vu - ef, ls="--", lw=0.8, color="C1", label=f"upper: {vu - ef:.2f} eV")
        axp.axhline(vl - ef, ls=":", lw=0.8, color="C2", label=f"lower: {vl - ef:.2f} eV")
        axp.axhline(0.0, color="grey", lw=0.6)
        axp.set_xlabel(f"position along a{a.axis} (A)")
        axp.set_ylabel("V - E$_F$ (eV)")
        axp.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig("planar_avg.png", dpi=200)
        print("wrote planar_avg.png")


if __name__ == "__main__":
    main()
