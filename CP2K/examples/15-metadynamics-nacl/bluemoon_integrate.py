#!/usr/bin/env python3
"""Blue-moon (constrained MD) free-energy profile for ONE distance constraint.

Each window directory is named window_<target in Angstrom> (e.g. window_3.20)
and contains bluemoon.LagrangeMultLog written by constrained_md_distance.inp:

    Shake  Lagrangian Multipliers:        <lambda>
    Rattle Lagrangian Multipliers:        <lambda>

For a single interatomic distance r the mass-metric factor is constant and its
curvature term vanishes, so the mean force is simply

    dA/dr = - < lambda_SHAKE >          (hartree/bohr, CP2K sign convention)

and A(r) follows by trapezoidal integration over the windows.  A(r) here is
-kT ln P(r) and therefore contains the entropic r^2 (Jacobian) term; the
"PMF" w(r) = -kT ln g(r) is obtained as w(r) = A(r) + 2 kT ln r (option --pmf).
Constraints on angles, combined CVs, several constraints, fixed atoms etc. need
the full mass-metric correction of Sprik & Ciccotti, J. Chem. Phys. 109, 7737
(1998), https://doi.org/10.1063/1.477419 -- not implemented here.
Sign convention as in the official CP2K exercise "Free Energy Profile of NaCl
Dissociation" (A = - integral of <lambda>).

Example:
    python3 bluemoon_integrate.py window_*/ --skip 2000 --pmf -T 300 -o bluemoon_fes.dat
Needs numpy (matplotlib only for --png).
"""
import argparse
import math
import os
import re
import sys

import numpy as np

BOHR_ANG = 0.529177210903
HARTREE_KJMOL = 2625.4996394799
KB_KJMOL = 0.0083144626


def read_lambdas(path):
    vals = []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("Shake"):
                nums = re.findall(r"[-+]?\d+\.\d+(?:[EeDd][-+]?\d+)?", line.split(":", 1)[1])
                if len(nums) != 1:
                    sys.exit(f"{path}: expected exactly one multiplier per line (one constraint)")
                vals.append(float(nums[0].replace("D", "E")))
    return np.array(vals)


def block_error(x, nblocks=5):
    n = len(x) // nblocks
    if n < 2:
        return float("nan")
    means = [x[i * n:(i + 1) * n].mean() for i in range(nblocks)]
    return float(np.std(means, ddof=1) / math.sqrt(nblocks))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("windows", nargs="+", help="window_<r> directories")
    ap.add_argument("--file", default="bluemoon.LagrangeMultLog", help="multiplier file name in each window")
    ap.add_argument("--skip", type=int, default=2000, help="discard the first N SHAKE records (equilibration), default 2000")
    ap.add_argument("-T", "--temperature", type=float, default=300.0)
    ap.add_argument("--pmf", action="store_true", help="also print w(r) = A(r) + 2kT ln r")
    ap.add_argument("-o", "--out", default="bluemoon_fes.dat")
    ap.add_argument("--png", help="plot to this file (needs matplotlib)")
    args = ap.parse_args()

    rows = []
    for d in args.windows:
        d = d.rstrip("/\\")
        m = re.search(r"window_([-+]?\d+(?:\.\d+)?)$", os.path.basename(d))
        if not m:
            print(f"# skip {d}: name is not window_<r>", file=sys.stderr)
            continue
        f = os.path.join(d, args.file)
        if not os.path.isfile(f):
            print(f"# skip {d}: {args.file} not found", file=sys.stderr)
            continue
        lam = read_lambdas(f)
        if len(lam) <= args.skip + 10:
            print(f"# skip {d}: only {len(lam)} records", file=sys.stderr)
            continue
        lam = lam[args.skip:]
        # mean force in kJ/mol/A
        mf = -lam.mean() * HARTREE_KJMOL / BOHR_ANG
        err = block_error(-lam) * HARTREE_KJMOL / BOHR_ANG
        rows.append((float(m.group(1)), mf, err, len(lam)))
    if len(rows) < 2:
        sys.exit("need at least two usable windows")
    rows.sort()
    r = np.array([x[0] for x in rows])
    mf = np.array([x[1] for x in rows])
    err = np.array([x[2] for x in rows])
    a = np.concatenate([[0.0], np.cumsum(0.5 * (mf[1:] + mf[:-1]) * np.diff(r))])
    kT = KB_KJMOL * args.temperature
    w = a + 2.0 * kT * np.log(r)
    a -= a.min()
    w -= w.min()

    hdr = "r(A)  dA/dr(kJ/mol/A)  err  A(r)(kJ/mol)" + ("  w(r)=A+2kTln(r)" if args.pmf else "") + "  nsamples"
    cols = [r, mf, err, a] + ([w] if args.pmf else []) + [np.array([x[3] for x in rows], float)]
    np.savetxt(args.out, np.column_stack(cols), fmt="%12.5f", header=hdr)
    print("#" + hdr)
    for row in np.column_stack(cols):
        print(" ".join(f"{v:12.4f}" for v in row))
    print(f"# written {args.out}; A(r) set to zero at its minimum")

    if args.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.8, 3.4))
        ax.plot(r, a, "o-", label="A(r)")
        if args.pmf:
            ax.plot(r, w, "s--", label="w(r) = A + 2kT ln r")
        ax.set_xlabel("d(Na-Cl) (A)")
        ax.set_ylabel("free energy (kJ/mol)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(args.png, dpi=200)
        print(f"# plot written to {args.png}")


if __name__ == "__main__":
    main()
