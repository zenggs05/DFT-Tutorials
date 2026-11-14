#!/usr/bin/env python3
"""Arithmetic on VASP volumetric files (CHGCAR / AECCAR / LOCPOT) with pymatgen.

Two common jobs:

  sum  : CHGCAR_sum = AECCAR0 + AECCAR2  (reference density for Bader).
         Replacement for chgsum.pl when that script prints "Atoms in file1: 0".
  diff : charge density difference  d_rho = rho(AB) - rho(A) - rho(B)
         written as CHGCAR_diff (open in VESTA) plus the plane-integrated
         profile d_lambda(z) in e/A, saved to diff_z.dat (and diff_z.png).

All files must come from runs with the SAME cell, ENCUT and FFT grid
(set NGXF/NGYF/NGZF explicitly); sub-systems must keep the atomic positions
of the combined system and must not be relaxed.

Usage:
  python pmg_chgcar_ops.py sum  AECCAR0 AECCAR2           [-o CHGCAR_sum]
  python pmg_chgcar_ops.py diff AB/CHGCAR A/CHGCAR B/CHGCAR [-o CHGCAR_diff] [--axis 2]
"""
import argparse

import numpy as np
from pymatgen.io.vasp.outputs import Chgcar


def check_grids(objs, names):
    dims = {tuple(o.dim) for o in objs}
    if len(dims) != 1:
        info = ", ".join(f"{n}: {tuple(o.dim)}" for o, n in zip(objs, names))
        raise SystemExit("FFT grids differ (" + info + "). Re-run with identical "
                         "NGXF/NGYF/NGZF, ENCUT and cell.")


def plane_profile(chg, axis):
    """Plane-integrated density along `axis` in e/A.

    CHGCAR stores rho*V_cell, so rho = data/V.  Averaging over the plane and
    multiplying by the plane area A = V/h (h = cell height along axis) gives
    the line density d_lambda(z) = <data>_plane / h.
    """
    lat = chg.structure.lattice
    vol = lat.volume
    other = [i for i in range(3) if i != axis]
    area = np.linalg.norm(np.cross(lat.matrix[other[0]], lat.matrix[other[1]]))
    height = vol / area
    z = np.array(chg.get_axis_grid(axis))
    prof = chg.get_average_along_axis(axis) / height
    return z, prof


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["sum", "diff"])
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--out")
    ap.add_argument("--axis", type=int, default=2, help="0/1/2 = a/b/c for the profile")
    a = ap.parse_args()

    objs = [Chgcar.from_file(f) for f in a.files]
    check_grids(objs, a.files)

    if a.mode == "sum":
        if len(objs) != 2:
            raise SystemExit("sum needs exactly two files: AECCAR0 AECCAR2")
        out = a.out or "CHGCAR_sum"
        total = objs[0] + objs[1]
        total.write_file(out)
        nel = total.data["total"].sum() / total.ngridpts
        print(f"written {out}; integrated electrons = {nel:.3f} (core + valence)")
        print(f"next: bader CHGCAR -ref {out}")
        return

    if len(objs) < 3:
        raise SystemExit("diff needs AB A B [C ...]")
    out = a.out or "CHGCAR_diff"
    diff = objs[0]
    for o in objs[1:]:
        diff = diff - o
    diff.write_file(out)
    # Sum of d_rho over the cell should be ~0 (charge conservation check)
    net = diff.data["total"].sum() / diff.ngridpts
    print(f"written {out}; net charge of d_rho = {net:.4f} e (should be close to 0)")

    z, prof = plane_profile(diff, a.axis)
    np.savetxt("diff_z.dat", np.column_stack([z, prof]),
               header="coordinate(A)  plane-integrated d_rho (e/A)")
    print("written diff_z.dat")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.plot(z, prof)
        ax.axhline(0, c="gray", lw=0.5)
        ax.set_xlabel("position along axis %d (A)" % a.axis)
        ax.set_ylabel(r"$\Delta\lambda$ (e/A)")
        fig.savefig("diff_z.png", dpi=300, bbox_inches="tight")
        print("written diff_z.png")
    except ImportError:
        pass


if __name__ == "__main__":
    main()
