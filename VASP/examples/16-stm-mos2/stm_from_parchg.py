#!/usr/bin/env python3
"""Tersoff-Hamann STM images from a VASP PARCHG (works with any VASP version).

constant height : rho(x, y, z_top + height)            (--mode height --height 3.0)
constant current: z(x, y) where rho = iso, scanning down from the vacuum (--mode current --iso 1e-4)

Assumptions: slab with the surface normal along z and the third lattice vector along z;
the image is taken above the topmost atom. Only the first data block is used
(for ISPIN = 2 this is the total partial density). Requires numpy + matplotlib.

Examples
--------
    python stm_from_parchg.py PARCHG --mode height  --height 3.0 --repeat 6 -o stm_h.png
    python stm_from_parchg.py PARCHG --mode current --iso 1e-4   --repeat 6 -o stm_c.png
"""
import argparse
import sys

import numpy as np


def read_parchg(path):
    """Return lattice (3x3, A), cartesian positions (N,3), data[NGX,NGY,NGZ] (= rho*V)."""
    with open(path) as fh:
        lines = fh.readlines()
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
    lat = lat * ((abs(scale) / abs(np.linalg.det(lat))) ** (1 / 3) if scale < 0 else scale)
    counts = [int(x) for x in lines[6].split()]
    i = 7
    if lines[i].strip()[:1] in "sS":
        i += 1
    cart = lines[i].strip()[:1] in "cCkK"
    i += 1
    nat = sum(counts)
    xyz = np.array([[float(x) for x in lines[i + k].split()[:3]] for k in range(nat)])
    pos = xyz * (scale if scale > 0 else 1.0) if cart else xyz @ lat
    i += nat
    while not lines[i].split():
        i += 1
    grid = [int(x) for x in lines[i].split()[:3]]
    ntot = grid[0] * grid[1] * grid[2]
    vals, i = [], i + 1
    while len(vals) < ntot:
        vals.extend(lines[i].split())
        i += 1
    data = np.array(vals[:ntot], dtype=float).reshape(grid[2], grid[1], grid[0]).transpose(2, 1, 0)
    return lat, pos, data


def layer(rho, zfrac):
    """Linear interpolation of rho between two z layers (periodic)."""
    nz = rho.shape[2]
    x = (zfrac % 1.0) * nz
    k0 = int(np.floor(x)) % nz
    w = x - np.floor(x)
    return (1 - w) * rho[:, :, k0] + w * rho[:, :, (k0 + 1) % nz]


def constant_current(rho, z_top, c_len, iso, zmax_above=8.0):
    """Height map (A above the top atom) of the iso-surface rho = iso, scanning downwards."""
    nz = rho.shape[2]
    dz = c_len / nz
    k_start = int(np.ceil((z_top + zmax_above) / dz))
    k_end = int(np.floor((z_top - 1.0) / dz))
    height = np.full(rho.shape[:2], np.nan)
    prev = rho[:, :, k_start % nz]
    for k in range(k_start - 1, k_end - 1, -1):
        cur = rho[:, :, k % nz]
        hit = np.isnan(height) & (cur >= iso) & (prev < iso)
        if hit.any():
            frac = (iso - prev[hit]) / (cur[hit] - prev[hit])      # between k+1 (frac=0) and k (frac=1)
            height[hit] = (k + 1 - frac) * dz - z_top
        prev = cur
    if np.isnan(height).all():
        sys.exit("iso value never reached - lower --iso or check the PARCHG")
    return np.where(np.isnan(height), np.nanmin(height), height)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("parchg", nargs="?", default="PARCHG")
    p.add_argument("--mode", choices=["height", "current"], default="height")
    p.add_argument("--height", type=float, default=3.0, help="tip height above the topmost atom (A)")
    p.add_argument("--iso", type=float, default=1e-4, help="iso value (e/A^3) for constant current")
    p.add_argument("--repeat", type=int, default=5, help="lateral repetition of the cell")
    p.add_argument("--cmap", default="afmhot")
    p.add_argument("-o", "--output", default="stm.png")
    args = p.parse_args()

    lat, pos, data = read_parchg(args.parchg)
    if abs(lat[2, 0]) > 1e-6 or abs(lat[2, 1]) > 1e-6 or abs(lat[0, 2]) > 1e-6 or abs(lat[1, 2]) > 1e-6:
        sys.exit("this script expects a and b in the xy plane and c along z")
    vol = abs(np.linalg.det(lat))
    rho = data / vol                                   # e/A^3
    c_len = lat[2, 2]
    z_top = (pos[:, 2] % c_len).max()
    print("grid %s, topmost atom at z = %.3f A, cell height %.3f A" % (data.shape, z_top, c_len))

    if args.mode == "height":
        img = layer(rho, (z_top + args.height) / c_len)
        label = "rho at %.1f A above top atom (e/A^3)" % args.height
        print("rho range on the plane: %.3e ... %.3e e/A^3" % (img.min(), img.max()))
    else:
        img = constant_current(rho, z_top, c_len, args.iso)
        label = "tip height for rho = %.1e e/A^3 (A)" % args.iso
        print("corrugation: %.3f A" % (img.max() - img.min()))

    r = args.repeat
    tile = np.tile(img, (r, r))
    nx, ny = img.shape
    i = np.arange(r * nx + 1) / nx
    j = np.arange(r * ny + 1) / ny
    ii, jj = np.meshgrid(i, j, indexing="ij")
    X = ii * lat[0, 0] + jj * lat[1, 0]
    Y = ii * lat[0, 1] + jj * lat[1, 1]
    np.savetxt(args.output.rsplit(".", 1)[0] + ".dat", img, fmt="%.6e",
               header="%s; rows = a index, columns = b index" % label)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))
    m = ax.pcolormesh(X, Y, tile, cmap=args.cmap, shading="flat")
    ax.set_aspect("equal")
    ax.set_xlabel("x (A)")
    ax.set_ylabel("y (A)")
    fig.colorbar(m, ax=ax, label=label)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print("written", args.output)


if __name__ == "__main__":
    main()
