#!/usr/bin/env python3
"""Optical constants from a VASP LOPTICS run (vasprun.xml).

    python plot_optics.py                          # ./vasprun.xml, density-density
    python plot_optics.py -f run/vasprun.xml --kind current
    python plot_optics.py --tauc direct --fit 3.0 3.6
    python plot_optics.py --component zz           # one Cartesian component

From eps = eps1 + i eps2 (isotropic average of xx, yy, zz by default):
    n     = sqrt((|eps| + eps1)/2)          refractive index
    k     = sqrt((|eps| - eps1)/2)          extinction coefficient
    alpha = 2 E k / (hbar c)                absorption coefficient (cm^-1)
    R     = ((n-1)^2 + k^2)/((n+1)^2 + k^2) normal-incidence reflectivity
    L     = eps2/(eps1^2 + eps2^2)          energy-loss function -Im(1/eps)
Tauc: (alpha E)^(1/r) versus E with r = 1/2 (direct allowed) or 2 (indirect
allowed); --fit Emin Emax extrapolates the linear part to the energy axis.

Outputs optics.dat (columns documented in its header) and optics.png
(and tauc.png).  Requires numpy and matplotlib.

Note: alpha, n, k, R are only meaningful for bulk (3D) cells.  For slabs /
2D materials the values depend on the vacuum thickness.
"""
import argparse
import sys
import xml.etree.ElementTree as ET

import numpy as np

HBAR_C_EV_CM = 1.973269804e-5  # hbar*c in eV*cm
COMPONENTS = ["xx", "yy", "zz", "xy", "yz", "zx"]


def read_dielectric(path, kind):
    """Return (energy, eps_real[N,6], eps_imag[N,6], label)."""
    blocks = []
    for _, elem in ET.iterparse(path, events=("end",)):
        if elem.tag == "dielectricfunction":
            blocks.append(elem)
    if not blocks:
        sys.exit("no <dielectricfunction> found: was LOPTICS = .TRUE. set and did VASP finish?")
    labelled = []
    unlabelled = 0
    for b in blocks:
        c = b.attrib.get("comment", "")
        if "density" in c:
            lab = "density"
        elif "current" in c:
            lab = "current"
        elif c == "":
            # VASP 5.x: first unlabelled block = density, second = current
            lab = ["density", "current"][min(unlabelled, 1)]
            unlabelled += 1
        else:
            lab = c
        labelled.append((lab, b))
    chosen = [b for lab, b in labelled if lab == kind]
    if not chosen:
        found = ", ".join(lab for lab, _ in labelled)
        sys.exit(f"dielectric function of kind '{kind}' not found (available: {found})")
    b = chosen[-1]

    def table(tag):
        rows = [[float(x) for x in r.text.split()] for r in b.find(tag).iter("r")]
        return np.array(rows)

    im, re = table("imag"), table("real")
    return im[:, 0], re[:, 1:7], im[:, 1:7], kind


def optical_constants(energy, e1, e2):
    mod = np.sqrt(e1 ** 2 + e2 ** 2)
    n = np.sqrt(np.clip((mod + e1) / 2, 0, None))
    k = np.sqrt(np.clip((mod - e1) / 2, 0, None))
    alpha = 2 * energy * k / HBAR_C_EV_CM
    refl = ((n - 1) ** 2 + k ** 2) / ((n + 1) ** 2 + k ** 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        loss = np.where(mod > 0, e2 / mod ** 2, 0.0)
    return n, k, alpha, refl, loss


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-f", "--file", default="vasprun.xml")
    ap.add_argument("--kind", default="density", choices=("density", "current"),
                    help="density-density (default) or current-current response")
    ap.add_argument("--component", default="avg", choices=["avg"] + COMPONENTS[:3],
                    help="isotropic average or a single diagonal component")
    ap.add_argument("--emax", type=float, default=10.0, help="upper energy for plots (eV)")
    ap.add_argument("--tauc", choices=("direct", "indirect"), help="also write tauc.png")
    ap.add_argument("--fit", nargs=2, type=float, metavar=("EMIN", "EMAX"),
                    help="linear Tauc fit window (eV)")
    args = ap.parse_args()

    energy, re6, im6, label = read_dielectric(args.file, args.kind)
    if args.component == "avg":
        e1, e2 = re6[:, :3].mean(axis=1), im6[:, :3].mean(axis=1)
    else:
        i = COMPONENTS.index(args.component)
        e1, e2 = re6[:, i], im6[:, i]
    n, k, alpha, refl, loss = optical_constants(energy, e1, e2)

    print(f"read {len(energy)} frequency points ({label}-{label}), "
          f"E = {energy[0]:.3f} ... {energy[-1]:.3f} eV")
    print(f"eps1(omega->0) = {e1[0]:.3f}   (compare with LEPSILON, no local fields)")
    print(f"n(omega->0)    = {n[0]:.3f}")

    header = ("E(eV) eps1 eps2 n k alpha(cm^-1) R loss   "
              f"[{label}-{label}, component={args.component}]")
    np.savetxt("optics.dat", np.column_stack([energy, e1, e2, n, k, alpha, refl, loss]),
               fmt="%12.6e", header=header)
    print("wrote optics.dat")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m = energy <= args.emax
    fig, ax = plt.subplots(2, 2, figsize=(9, 6.5))
    ax[0, 0].plot(energy[m], e1[m], label=r"$\varepsilon_1$")
    ax[0, 0].plot(energy[m], e2[m], label=r"$\varepsilon_2$")
    ax[0, 0].axhline(0, color="0.6", lw=0.5)
    ax[0, 0].legend()
    ax[0, 0].set_ylabel("dielectric function")
    ax[0, 1].plot(energy[m], n[m], label="n")
    ax[0, 1].plot(energy[m], k[m], label="k")
    ax[0, 1].legend()
    ax[0, 1].set_ylabel("n, k")
    ax[1, 0].semilogy(energy[m], np.clip(alpha[m], 1e2, None))
    ax[1, 0].set_ylabel(r"$\alpha$ (cm$^{-1}$)")
    ax[1, 1].plot(energy[m], refl[m], label="R")
    ax[1, 1].plot(energy[m], loss[m] / max(loss[m].max(), 1e-12), label="loss (scaled)")
    ax[1, 1].legend()
    for a in ax.flat:
        a.set_xlabel("photon energy (eV)")
        a.set_xlim(0, args.emax)
    fig.tight_layout()
    fig.savefig("optics.png", dpi=200)
    print("wrote optics.png")

    if args.tauc:
        power = 2.0 if args.tauc == "direct" else 0.5
        y = (alpha * energy) ** power
        fig2, ax2 = plt.subplots(figsize=(5, 4))
        ax2.plot(energy[m], y[m])
        ax2.set_xlabel("photon energy (eV)")
        ax2.set_ylabel(r"$(\alpha E)^{%g}$" % power)
        if args.fit:
            w = (energy >= args.fit[0]) & (energy <= args.fit[1])
            if w.sum() < 3:
                sys.exit("fit window contains fewer than 3 points")
            slope, icpt = np.polyfit(energy[w], y[w], 1)
            eg = -icpt / slope
            x = np.linspace(eg, args.fit[1], 50)
            ax2.plot(x, slope * x + icpt, "--", label=f"fit: Eg = {eg:.2f} eV")
            ax2.legend()
            print(f"Tauc ({args.tauc}) extrapolated gap = {eg:.3f} eV")
        ax2.set_xlim(0, args.emax)
        ax2.set_ylim(bottom=0)
        fig2.tight_layout()
        fig2.savefig("tauc.png", dpi=200)
        print("wrote tauc.png")


if __name__ == "__main__":
    main()
