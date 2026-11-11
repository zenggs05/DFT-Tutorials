#!/usr/bin/env python3
"""Read an OpenMX System.Name.Band file, find VBM/CBM/gap and plot the bands.

Format of System.Name.Band (written by Band_DFT_kpath.c, read by bandgnu13.c):
    line 1 : n_eigen  spin_flag  chemical_potential(Hartree)
             spin_flag = 0 (non-spin or non-collinear), 1 (collinear spin-polarized)
    line 2 : reciprocal lattice vectors rtv (9 numbers, Bohr^-1, rows b1 b2 b3)
    line 3 : number of k-path segments
    next   : one line per segment: n_k  k1 k2 k3  k1' k2' k3'  label1 label2
    then   : for each segment, each k point, each spin:
             "n_eigen k1 k2 k3" (fractional, in units of rtv) followed by n_eigen
             eigenvalues in Hartree.
Energies are converted to eV and shifted so that the chemical potential is 0,
exactly as bandgnu13 does; the x axis is the path length in Bohr^-1.

Usage:
    python3 band_parse.py cdia.Band
    python3 band_parse.py si.Band --emin -13 --emax 6 --plot si_band.png
Requires numpy; matplotlib only for --plot.
"""
import argparse

import numpy as np

HA2EV = 27.2113845          # same factor as bandgnu13.c


def read_band(fname):
    with open(fname) as f:
        tokens = f.read().split()
    pos = 0

    def nxt():
        nonlocal pos
        pos += 1
        return tokens[pos - 1]

    neig, spin_flag, chemp = int(nxt()), int(nxt()), float(nxt())
    rtv = np.array([float(nxt()) for _ in range(9)]).reshape(3, 3)
    nseg = int(nxt())
    segs = []
    for _ in range(nseg):
        nk = int(nxt())
        k0 = [float(nxt()) for _ in range(3)]
        k1 = [float(nxt()) for _ in range(3)]
        segs.append((nk, k0, k1, nxt(), nxt()))
    nspin = spin_flag + 1 if spin_flag in (0, 1) else 1
    kfrac, eig = [], [[] for _ in range(nspin)]
    for nk, _, _, _, _ in segs:
        for _ in range(nk):
            for s in range(nspin):
                n = int(nxt())
                k = [float(nxt()) for _ in range(3)]
                e = [float(nxt()) for _ in range(n)]
                if s == 0:
                    kfrac.append(k)
                eig[s].append(e[:neig])
    nmin = min(len(e) for es in eig for e in es)
    eig = np.array([[e[:nmin] for e in es] for es in eig])          # [spin, k, band]
    eig = (eig - chemp) * HA2EV
    kfrac = np.array(kfrac)
    # path length: zero-length jump at segment boundaries (as in bandgnu13)
    kcart = kfrac @ rtv
    x = np.zeros(len(kfrac))
    i = 0
    ticks, point_label = [], {}
    for nk, _, _, l0, l1 in segs:
        for j in range(nk):
            if j == 0:
                x[i] = x[i - 1] if i > 0 else 0.0
            else:
                x[i] = x[i - 1] + np.linalg.norm(kcart[i] - kcart[i - 1])
            i += 1
        ticks.append((x[i - nk], l0))
        ticks.append((x[i - 1], l1))
        point_label[i - nk] = pretty(l0)
        point_label[i - 1] = pretty(l1)
    return dict(chemp=chemp, rtv=rtv, segs=segs, x=x, kfrac=kfrac, eig=eig,
                ticks=merge_ticks(ticks), point_label=point_label)


def pretty(label):
    return "G" if label.lower() in ("g", "gamma", "gm") else label


def tex(label):
    """Matplotlib label: G -> Gamma (mathtext keeps this script ASCII-only)."""
    return "|".join(r"$\Gamma$" if p == "G" else p for p in label.split("|"))


def merge_ticks(ticks):
    """Merge labels at the same x; 'K' + 'U' at a jump becomes 'K|U'."""
    out = []
    for xv, lab in ticks:
        lab = pretty(lab)
        if out and abs(out[-1][0] - xv) < 1e-8:
            if out[-1][1] != lab:
                out[-1] = (xv, out[-1][1] + "|" + lab)
        else:
            out.append((xv, lab))
    return out


def nearest_label(data, idx):
    """Label of a k point: high-symmetry name if it is a segment end, else its coordinates."""
    kf = "k=(%.3f %.3f %.3f)" % tuple(data["kfrac"][idx])
    if idx in data["point_label"]:
        return "%s %s" % (data["point_label"][idx], kf)
    return kf


def gap_report(data):
    eig = data["eig"]
    nspin = eig.shape[0]
    below = eig <= 0.0
    # a band is metallic if it has states on both sides of E_F along the path
    crossing = [(s, b) for s in range(nspin) for b in range(eig.shape[2])
                if below[s, :, b].any() and (~below[s, :, b]).any()]
    if crossing:
        print("bands crossing E_F (metal or semimetal along this path):",
              ", ".join("spin%d/band%d" % (s, b + 1) for s, b in crossing[:6]))
    occ = np.where(below, eig, -np.inf)
    emp = np.where(~below, eig, np.inf)
    s_v, k_v, b_v = np.unravel_index(np.argmax(occ), occ.shape)
    s_c, k_c, b_c = np.unravel_index(np.argmin(emp), emp.shape)
    vbm, cbm = occ[s_v, k_v, b_v], emp[s_c, k_c, b_c]
    print("E_F (chemical potential) = %.4f eV (set to 0)" % (data["chemp"] * HA2EV))
    print("VBM = %8.4f eV at %s (spin %d, band %d)" % (vbm, nearest_label(data, k_v), s_v, b_v + 1))
    print("CBM = %8.4f eV at %s (spin %d, band %d)" % (cbm, nearest_label(data, k_c), s_c, b_c + 1))
    if not crossing:
        kind = "direct" if k_v == k_c else "indirect"
        print("band gap along the path = %.4f eV (%s)" % (cbm - vbm, kind))
        direct = np.min(np.min(emp, axis=2) - np.max(occ, axis=2))
        print("smallest direct gap along the path = %.4f eV" % direct)
    print("note: VBM/CBM are searched only on the k points of the path")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bandfile", help="System.Name.Band")
    ap.add_argument("--emin", type=float, default=-15.0)
    ap.add_argument("--emax", type=float, default=10.0)
    ap.add_argument("--plot", default="", help="output image, e.g. band.png")
    args = ap.parse_args()

    data = read_band(args.bandfile)
    print("%d segments, %d k points, %d spin channel(s), %d bands"
          % (len(data["segs"]), len(data["x"]), data["eig"].shape[0], data["eig"].shape[2]))
    print("ticks:", "  ".join("%s=%.4f" % (lab, xv) for xv, lab in data["ticks"]))
    gap_report(data)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4.5))
        colors = ["tab:blue", "tab:red"]
        for s in range(data["eig"].shape[0]):
            for b in range(data["eig"].shape[2]):
                ax.plot(data["x"], data["eig"][s, :, b], color=colors[s], lw=1.0,
                        label=("spin up" if s == 0 else "spin down") if (b == 0 and data["eig"].shape[0] == 2) else None)
        for xv, _ in data["ticks"]:
            ax.axvline(xv, color="k", lw=0.5)
        ax.axhline(0.0, color="gray", ls="--", lw=0.8)
        ax.set_xticks([xv for xv, _ in data["ticks"]])
        ax.set_xticklabels([tex(lab) for _, lab in data["ticks"]])
        ax.set_xlim(data["x"][0], data["x"][-1])
        ax.set_ylim(args.emin, args.emax)
        ax.set_ylabel("E - E$_F$ (eV)")
        if data["eig"].shape[0] == 2:
            ax.legend(loc="best")
        fig.tight_layout()
        fig.savefig(args.plot, dpi=200)
        print("figure written to", args.plot)


if __name__ == "__main__":
    main()
