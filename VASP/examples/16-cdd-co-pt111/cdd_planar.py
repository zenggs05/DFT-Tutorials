#!/usr/bin/env python3
"""Charge density difference  d_rho = rho(AB) - rho(A) - rho(B) - ...  from CHGCAR files.

* checks that lattice and FFT grid of all files are identical
* writes the difference in CHGCAR format (readable by VESTA)
* plane-averages along one lattice vector and integrates the charge displacement
  dQ(z) = sum over planes up to z (in electrons)

Only the first data block (total density) of each file is used.
Requires numpy; matplotlib is optional (for cdd_planar.png).

Example
-------
    python cdd_planar.py AB/CHGCAR slab/CHGCAR co/CHGCAR --axis z --out CHGCAR_cdd
Output: CHGCAR_cdd, cdd_planar.dat (z, <d_rho>, d_lambda, dQ), cdd_planar.png
"""
import argparse
import sys

import numpy as np


def read_chgcar(path):
    """Return (header_lines, lattice(3x3, A), grid(3,), data[NGX,NGY,NGZ]) of the first block."""
    with open(path) as fh:
        lines = fh.readlines()
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
    if scale < 0:                                   # negative scale = target volume
        lat *= (abs(scale) / abs(np.linalg.det(lat))) ** (1.0 / 3.0)
    else:
        lat *= scale
    counts = [int(x) for x in lines[6].split()]
    i = 7
    if lines[i].strip()[:1] in "sS":
        i += 1
    i += 1 + sum(counts)                            # coordinate mode line + atoms
    header = lines[:i]
    while not lines[i].split():                     # blank line(s)
        i += 1
    grid = [int(x) for x in lines[i].split()[:3]]
    ntot = grid[0] * grid[1] * grid[2]
    vals = []
    i += 1
    while len(vals) < ntot:
        vals.extend(lines[i].split())
        i += 1
    data = np.array(vals[:ntot], dtype=float).reshape(grid[2], grid[1], grid[0]).transpose(2, 1, 0)
    return header, lat, np.array(grid), data


def write_chgcar(path, header, grid, data):
    flat = data.ravel(order="F")                    # x fastest, z slowest
    with open(path, "w", newline="\n") as fh:
        fh.writelines(header)
        fh.write("\n %5d %5d %5d\n" % tuple(grid))
        n5 = len(flat) // 5 * 5
        np.savetxt(fh, flat[:n5].reshape(-1, 5), fmt=" %17.10E", delimiter="")
        if len(flat) > n5:
            fh.write("".join(" %17.10E" % v for v in flat[n5:]) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="CHGCAR of AB first, then the fragments to subtract")
    p.add_argument("--axis", default="z", choices=["x", "y", "z", "a", "b", "c"],
                   help="lattice vector along which to plane-average (x=a, y=b, z=c)")
    p.add_argument("--out", default="CHGCAR_cdd")
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args()
    if len(args.files) < 2:
        sys.exit("need at least two files")

    header, lat, grid, total = read_chgcar(args.files[0])
    nel = [total.sum() / total.size]
    for f in args.files[1:]:
        _, lat2, grid2, d = read_chgcar(f)
        if not np.array_equal(grid, grid2):
            sys.exit("grid mismatch: %s %s vs %s %s -> use identical ENCUT/PREC or set NGXF/NGYF/NGZF"
                     % (args.files[0], grid, f, grid2))
        if np.abs(lat - lat2).max() > 1e-5:
            sys.exit("lattice mismatch between %s and %s -> fragments must use the same cell" % (args.files[0], f))
        nel.append(d.sum() / d.size)
        total = total - d
    for f, n in zip(args.files, nel):
        print("%-30s electrons = %10.4f" % (f, n))
    print("integral of d_rho (should be ~0): %.5f e" % (total.sum() / total.size))

    write_chgcar(args.out, header, grid, total)
    print("written", args.out, "(open in VESTA, isosurface e.g. +/-0.002 e/A^3)")

    ax = {"x": 0, "a": 0, "y": 1, "b": 1, "z": 2, "c": 2}[args.axis]
    vol = abs(np.linalg.det(lat))
    others = [k for k in range(3) if k != ax]
    area = np.linalg.norm(np.cross(lat[others[0]], lat[others[1]]))
    height = vol / area                              # thickness of the cell along the plane normal
    n = grid[ax]
    mean = total.mean(axis=tuple(others))            # plane average of (rho*V) per layer
    zz = np.arange(n) * height / n
    drho = mean / vol                                # e/A^3
    dlam = mean / height                             # e/A  (charge per unit length)
    dq = np.cumsum(mean) / n                         # electrons
    np.savetxt("cdd_planar.dat", np.column_stack([zz, drho, dlam, dq]),
               fmt="%12.5f %14.6e %14.6e %12.6f",
               header="z(A)  <d_rho>(e/A^3)  d_lambda(e/A)  dQ(z)(e)")
    print("written cdd_planar.dat; max |dQ| = %.4f e" % np.abs(dq).max())

    if args.no_plot:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available - skip plot")
        return
    fig, ax1 = plt.subplots(figsize=(6, 3.5))
    ax1.plot(zz, dlam, color="tab:blue", label=r"$\Delta\lambda(z)$")
    ax1.axhline(0, color="k", lw=0.5)
    ax1.set_xlabel("position along %s (A)" % args.axis)
    ax1.set_ylabel(r"$\Delta\lambda$ (e/A)", color="tab:blue")
    ax2 = ax1.twinx()
    ax2.plot(zz, dq, color="tab:red", label=r"$\Delta Q(z)$")
    ax2.set_ylabel(r"$\Delta Q$ (e)", color="tab:red")
    fig.tight_layout()
    fig.savefig("cdd_planar.png", dpi=200)
    print("written cdd_planar.png")


if __name__ == "__main__":
    main()
