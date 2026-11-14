#!/usr/bin/env python3
"""Compare Wannier-interpolated bands with the original OpenMX bands.

Both files must use the OpenMX .Band layout (System.Name.Band,
System.Name.Wannier_Band from MLWF, System.Name.CWF.Band from CWF) and the
same Band.kpath:
    line 1 : n_eigen  spin_flag  chemical_potential(Hartree)
    line 2 : reciprocal lattice vectors (9 numbers)
    line 3 : number of k-path segments, then one line per segment
    then   : per k point and spin: "n_eigen k1 k2 k3" + n_eigen eigenvalues (Hartree)

For every interpolated eigenvalue inside [--emin, --emax] (eV, relative to the
chemical potential of each file) the closest original eigenvalue at the same k
point is searched and the deviation is reported (max / mean / rms).

Usage:
    python3 band_compare.py si_mlwf_sp3.Band si_mlwf_sp3.Wannier_Band --emin -14 --emax 1
    python3 band_compare.py gaas_cwf_ao.Band gaas_cwf_ao.CWF.Band --emax 1 --plot cmp.png
Requires numpy; matplotlib only for --plot.
"""
import argparse
import sys

import numpy as np

HA2EV = 27.2113845      # same factor as bandgnu13.c


def read_band(fname):
    with open(fname) as f:
        tok = f.read().split()
    pos = [0]

    def nxt():
        pos[0] += 1
        return tok[pos[0] - 1]

    neig, spin_flag, chemp = int(nxt()), int(nxt()), float(nxt())
    rtv = np.array([float(nxt()) for _ in range(9)]).reshape(3, 3)
    nseg = int(nxt())
    segs = []
    for _ in range(nseg):
        nk = int(nxt())
        k0 = [float(nxt()) for _ in range(3)]
        k1 = [float(nxt()) for _ in range(3)]
        segs.append((nk, k0, k1, nxt(), nxt()))
    nspin = 2 if spin_flag == 1 else 1
    kfrac, eig = [], [[] for _ in range(nspin)]
    for nk, _, _, _, _ in segs:
        for _ in range(nk):
            for s in range(nspin):
                n = int(nxt())
                k = [float(nxt()) for _ in range(3)]
                e = [float(nxt()) for _ in range(n)]
                if s == 0:
                    kfrac.append(k)
                eig[s].append((np.array(e[:neig]) - chemp) * HA2EV)
    kfrac = np.array(kfrac)
    kcart = kfrac @ rtv
    x = np.zeros(len(kfrac))
    i = 0
    ticks = []
    for nk, _, _, l0, l1 in segs:
        for j in range(nk):
            if i > 0:
                x[i] = x[i - 1] + (0.0 if j == 0 else np.linalg.norm(kcart[i] - kcart[i - 1]))
            i += 1
        ticks.append((x[i - nk], l0))
        ticks.append((x[i - 1], l1))
    return dict(chemp=chemp, kfrac=kfrac, eig=eig, x=x, ticks=ticks, nspin=nspin)


def compare(ref, wan, emin, emax):
    if len(ref["kfrac"]) != len(wan["kfrac"]):
        sys.exit("different number of k points: %d vs %d (use the same Band.kpath)"
                 % (len(ref["kfrac"]), len(wan["kfrac"])))
    if np.abs(ref["kfrac"] - wan["kfrac"]).max() > 1e-6:
        print("warning: k points differ between the two files")
    dev, worst = [], (0.0, None)
    for s in range(min(ref["nspin"], wan["nspin"])):
        for ik in range(len(ref["kfrac"])):
            e_ref = ref["eig"][s][ik]
            for ew in wan["eig"][s][ik]:
                if emin <= ew <= emax:
                    d = np.abs(e_ref - ew).min()
                    dev.append(d)
                    if d > worst[0]:
                        worst = (d, (s, ik, ew))
    return np.array(dev), worst


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("ref", help="original band file (System.Name.Band)")
    ap.add_argument("wan", help="interpolated band file (.Wannier_Band or .CWF.Band)")
    ap.add_argument("--emin", type=float, default=-100.0, help="lower energy bound (eV)")
    ap.add_argument("--emax", type=float, default=1.0, help="upper energy bound (eV)")
    ap.add_argument("--plot", help="write a comparison figure (png/pdf)")
    ap.add_argument("--ylim", type=float, nargs=2, default=None, help="y range of the plot")
    a = ap.parse_args()

    ref, wan = read_band(a.ref), read_band(a.wan)
    print("chemical potential (eV): ref %.4f  wannier %.4f"
          % (ref["chemp"] * HA2EV, wan["chemp"] * HA2EV))
    dev, worst = compare(ref, wan, a.emin, a.emax)
    if dev.size == 0:
        sys.exit("no interpolated eigenvalue inside [%.2f, %.2f] eV" % (a.emin, a.emax))
    print("window [%.2f, %.2f] eV: %d eigenvalues compared" % (a.emin, a.emax, dev.size))
    print("  max |dE| = %.4f eV   mean = %.4f eV   rms = %.4f eV"
          % (dev.max(), dev.mean(), np.sqrt((dev ** 2).mean())))
    s, ik, ew = worst[1]
    print("  worst point: spin %d, k #%d (%.4f %.4f %.4f), E_wannier = %.3f eV"
          % (s, ik + 1, *ref["kfrac"][ik], ew))

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 4))
        for s in range(ref["nspin"]):
            e = np.array([v for v in ref["eig"][s]])
            ax.plot(ref["x"], e, color="0.3", lw=1.0)
            ew_all = wan["eig"][s]
            for ik, row in enumerate(ew_all):
                ax.plot([wan["x"][ik]] * len(row), row, "o", ms=2.0, color="tab:red")
        for xv, lab in ref["ticks"]:
            ax.axvline(xv, color="0.8", lw=0.5)
        ax.set_xticks([t[0] for t in ref["ticks"]])
        ax.set_xticklabels([r"$\Gamma$" if t[1].upper() in ("G", "GAMMA") else t[1]
                            for t in ref["ticks"]])
        ax.axhline(0.0, ls="--", color="0.5", lw=0.6)
        ax.set_xlim(ref["x"][0], ref["x"][-1])
        if a.ylim:
            ax.set_ylim(*a.ylim)
        ax.set_ylabel("E - E_F (eV)")
        ax.set_title("lines: original   dots: Wannier", fontsize=9)
        fig.tight_layout()
        fig.savefig(a.plot, dpi=200)
        print("figure written to", a.plot)


if __name__ == "__main__":
    main()
