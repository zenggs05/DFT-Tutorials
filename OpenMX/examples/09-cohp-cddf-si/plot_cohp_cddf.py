#!/usr/bin/env python3
"""Plot OpenMX COHP/COOP files and CDDF optical spectra.

COHP (post-processing code "cohp", OpenMX >= 4.0, manual Sec. 24), with
PREFIX = System.Name (the binary input of cohp is PREFIX.cohp):
    PREFIX_i<pair>_s<spin>_tetra.cohp           (summed over all orbitals of both atoms)
    PREFIX_i<pair>_o<a>_s<spin>_tetra.cohp      (orbital a of atom 1, summed over atom 2)
    PREFIX_i<pair>_o<a>_o<b>_s<spin>_tetra.cohp (fully orbital resolved, a/b start at 0)
    (.coop files have the same format; "_gauss" replaces "_tetra" for Gaussian broadening)
    columns: E (eV, E_F = 0), COHP (dimensionless), integrated COHP up to E
    Sign convention of the manual: negative COHP = bonding, positive = antibonding;
    positive COOP = bonding. The plot shows -COHP so that bonding points to the right.

CDDF (manual Sec. 62): PREFIX.df_re, PREFIX.df_im, PREFIX.cd_re, ... with
    columns: energy(eV), xx, xy, xz, yx, yy, yz, zx, zy, zz, (xx+yy+zz)/3

Usage:
    python3 plot_cohp_cddf.py cohp si_cohp_i1_s0_tetra.cohp si_cohp_i3_s0_tetra.cohp
    python3 plot_cohp_cddf.py optics si_cddf
Requires numpy and matplotlib.
"""
import argparse
import os

import numpy as np


def plot_cohp(files, emin, emax, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(8, 4.5), sharey=True)
    for f in files:
        d = np.loadtxt(f, comments="#", ndmin=2)
        e, c, ic = d[:, 0], d[:, 1], d[:, 2]
        icohp_ef = np.interp(0.0, e, ic)
        if f.endswith(".coop"):
            print("%-45s ICOOP(E_F) = %8.4f   (bonding if positive)" % (os.path.basename(f), icohp_ef))
        else:
            print("%-45s ICOHP(E_F) = %8.4f   (bonding if negative)" % (os.path.basename(f), icohp_ef))
        label = os.path.basename(f).replace(".cohp", "").replace(".coop", "")
        sign = -1.0 if f.endswith(".cohp") else 1.0
        ax[0].plot(sign * c, e, lw=1.0, label=label)
        ax[1].plot(sign * ic, e, lw=1.0)
    for a in ax:
        a.axhline(0.0, color="gray", ls="--", lw=0.8)
        a.axvline(0.0, color="k", lw=0.5)
    ax[0].set_ylim(emin, emax)
    ax[0].set_ylabel("E - E$_F$ (eV)")
    ax[0].set_xlabel("-COHP (COOP)   bonding ->")
    ax[1].set_xlabel("-ICOHP (ICOOP)")
    ax[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print("figure written to", out)


def plot_optics(prefix, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    re_ = np.loadtxt(prefix + ".df_re", comments="#", ndmin=2)
    im_ = np.loadtxt(prefix + ".df_im", comments="#", ndmin=2)
    e, eps1, eps2 = re_[:, 0], re_[:, -1], im_[:, -1]            # last column: (xx+yy+zz)/3
    print("energy grid: %.4f - %.4f eV, %d points" % (e[0], e[-1], len(e)))
    print("eps1 at the lowest energy (%.3f eV) = %.3f  (electronic, no local-field effects)" % (e[0], eps1[0]))
    # local maxima of eps2 above 1 eV
    pk = [i for i in range(1, len(e) - 1)
          if eps2[i] > eps2[i - 1] and eps2[i] >= eps2[i + 1] and e[i] > 1.0 and eps2[i] > 0.2 * eps2.max()]
    print("eps2 peaks near: " + ", ".join("%.2f eV" % e[i] for i in pk[:5]))
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.8))
    ax[0].plot(e, eps1, label=r"$\varepsilon_1$")
    ax[0].plot(e, eps2, label=r"$\varepsilon_2$")
    ax[0].axhline(0.0, color="k", lw=0.5)
    ax[0].set_xlabel("photon energy (eV)")
    ax[0].set_ylabel("dielectric function (trace/3)")
    ax[0].legend()
    if os.path.exists(prefix + ".absorption"):
        ab = np.loadtxt(prefix + ".absorption", comments="#", ndmin=2)
        ax[1].plot(ab[:, 0], ab[:, -1])
        ax[1].set_ylabel("absorption (trace/3, units: see file header)")
    else:
        cd = np.loadtxt(prefix + ".cd_re", comments="#", ndmin=2)
        ax[1].plot(cd[:, 0], cd[:, -1])
        ax[1].set_ylabel(r"Re $\sigma$ (S/m)")
    ax[1].set_xlabel("photon energy (eV)")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print("figure written to", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("cohp", help="plot one or more .cohp/.coop files")
    c.add_argument("files", nargs="+")
    c.add_argument("--emin", type=float, default=-14.0)
    c.add_argument("--emax", type=float, default=6.0)
    c.add_argument("--out", default="cohp.png")
    o = sub.add_parser("optics", help="plot PREFIX.df_re / df_im (and absorption)")
    o.add_argument("prefix")
    o.add_argument("--out", default="optics.png")
    args = ap.parse_args()
    if args.cmd == "cohp":
        plot_cohp(args.files, args.emin, args.emax, args.out)
    else:
        plot_optics(args.prefix, args.out)


if __name__ == "__main__":
    main()
