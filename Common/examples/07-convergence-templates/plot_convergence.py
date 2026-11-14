#!/usr/bin/env python3
"""Plot and judge a convergence scan from any of the four codes.

Input: one or more whitespace tables, '#' lines are ignored.
  column 1 : scanned parameter (ENCUT, N1*N2*N3, CUTOFF, scf.energycutoff, ...)
  column 2 : total energy of the cell / molecule
  column 3 : number of atoms (optional; otherwise use --natoms)
Extra columns are ignored, so encut.dat, kpoints.dat, cutoff.dat ... can be used directly.

The reference is the LAST (most accurate) point unless --ref is given.
A parameter is "converged" when it and every larger value stay within
--thresh (meV/atom) of the reference.

Examples
  python plot_convergence.py encut_scan/encut.dat --xlabel "ENCUT (eV)"
  python plot_convergence.py kpoints.dat --xlabel "N1*N2*N3" --logx
  python plot_convergence.py cp2k_cutoff_scan/cutoff.dat --unit hartree --xlabel "CUTOFF (Ry)"
  python plot_convergence.py a.dat b.dat --labels Si Al --out conv.png
"""
import argparse
import sys

HARTREE_EV = 27.211386245988      # CODATA 2018
UNITS = {"ev": 1.0, "hartree": HARTREE_EV, "ry": HARTREE_EV / 2.0,
         "kcal/mol": 0.0433641043, "kj/mol": 0.0103642688}


def read_table(path, natoms_default):
    xs, es, ns = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            try:
                x, e = float(parts[0]), float(parts[1])
            except (ValueError, IndexError):
                continue
            n = natoms_default
            if natoms_default is None:
                try:
                    n = int(float(parts[2]))
                except (ValueError, IndexError):
                    n = 1
            xs.append(x)
            es.append(e)
            ns.append(n)
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    return [xs[i] for i in order], [es[i] for i in order], [ns[i] for i in order]


def analyse(xs, es, ns, unit, ref, thresh):
    conv = UNITS[unit]
    if ref is None:
        iref = len(xs) - 1
    else:
        iref = min(range(len(xs)), key=lambda i: abs(xs[i] - ref))
    eref = es[iref] / ns[iref]
    de = [1000.0 * conv * (e / n - eref) for e, n in zip(es, ns)]   # meV/atom
    first = None
    for i in range(len(xs)):
        if all(abs(d) < thresh for d in de[i:]):
            first = xs[i]
            break
    return de, xs[iref], first


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+")
    p.add_argument("--unit", default="ev", choices=sorted(UNITS), help="unit of column 2")
    p.add_argument("--natoms", type=int, default=None, help="atoms per cell if column 3 is absent")
    p.add_argument("--ref", type=float, default=None, help="x value used as reference (default: largest x)")
    p.add_argument("--thresh", type=float, default=1.0, help="tolerance in meV/atom (default 1.0)")
    p.add_argument("--xlabel", default="parameter")
    p.add_argument("--labels", nargs="*", default=None)
    p.add_argument("--logx", action="store_true")
    p.add_argument("--out", default="convergence.png")
    p.add_argument("--no-plot", action="store_true", help="only print the table")
    a = p.parse_args()

    results = []
    for k, f in enumerate(a.files):
        xs, es, ns = read_table(f, a.natoms)
        if len(xs) < 2:
            sys.exit(f"ERROR: {f} needs at least two numeric rows")
        de, xref, first = analyse(xs, es, ns, a.unit, a.ref, a.thresh)
        label = a.labels[k] if a.labels and k < len(a.labels) else f
        results.append((label, xs, de))
        print(f"== {label}  (reference x = {xref:g})")
        print(f"{'x':>12s} {'dE (meV/atom)':>16s}")
        for x, d in zip(xs, de):
            flag = "" if abs(d) < a.thresh else "  *"
            print(f"{x:12g} {d:16.3f}{flag}")
        if first is None:
            print(f"not converged within {a.thresh} meV/atom: extend the scan")
        else:
            print(f"converged (|dE| < {a.thresh} meV/atom from here on): x = {first:g}")
        print("note: the reference is itself a finite value; check that the last two points agree.")

    if a.no_plot:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.8, 3.4), constrained_layout=True)
    for label, xs, de in results:
        ax.plot(xs, [abs(d) if d != 0 else float("nan") for d in de], "o-", label=label)
    ax.axhline(a.thresh, ls="--", lw=1, color="0.4")
    ax.set_yscale("log")
    if a.logx:
        ax.set_xscale("log")
    ax.set_xlabel(a.xlabel)
    ax.set_ylabel("|E - E_ref| (meV/atom)")
    if len(results) > 1:
        ax.legend(frameon=False)
    fig.savefig(a.out, dpi=200)
    print(f"figure written: {a.out}")


if __name__ == "__main__":
    main()
