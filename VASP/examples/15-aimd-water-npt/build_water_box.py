#!/usr/bin/env python3
"""Build a cubic box of water molecules for VASP (POSCAR, O first then H).

    python build_water_box.py                 # 32 H2O, rho = 0.997 g/cm^3
    python build_water_box.py --cells 3 --density 1.0 --seed 7 -o POSCAR_108

Molecules sit on an fcc lattice (4 per cubic sub-cell, --cells^3 sub-cells,
so 32 molecules for --cells 2) with random orientations; geometry
O-H = 0.97 A, H-O-H = 104.5 deg.  Orientations are re-drawn until every
intermolecular H...H / O...H distance exceeds a safety threshold.
Only numpy is required.  The structure is a starting guess: equilibrate it
with MD before collecting statistics.
"""
import argparse

import numpy as np

NA = 6.02214076e23
M_H2O = 18.015  # g/mol
R_OH, ANGLE = 0.97, np.deg2rad(104.5)


def random_rotation(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([
        [a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
        [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
        [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d],
    ])


def min_image_dist(x, y, box):
    d = x[:, None, :] - y[None, :, :]
    d -= box * np.round(d / box)
    return np.linalg.norm(d, axis=-1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells", type=int, default=2, help="fcc sub-cells per edge (n_mol = 4 n^3)")
    ap.add_argument("--density", type=float, default=0.997, help="g/cm^3")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--dmin", type=float, default=1.6, help="min intermolecular distance (A)")
    ap.add_argument("-o", "--output", default="POSCAR")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    nmol = 4 * args.cells ** 3
    vol = nmol * M_H2O / NA / args.density * 1e24       # A^3
    box = vol ** (1 / 3)
    sub = box / args.cells
    fcc = np.array([[0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    centres = np.array([(f + [i, j, k]) * sub for i in range(args.cells)
                        for j in range(args.cells) for k in range(args.cells) for f in fcc])
    centres += 0.25 * sub                                # move away from the cell origin
    local_h = R_OH * np.array([[np.sin(ANGLE / 2), 0, np.cos(ANGLE / 2)],
                               [-np.sin(ANGLE / 2), 0, np.cos(ANGLE / 2)]])

    for attempt in range(1000):
        hyd = np.vstack([c + local_h @ random_rotation(rng).T for c in centres])
        owner = np.repeat(np.arange(nmol), 2)
        d_hh = min_image_dist(hyd, hyd, box)
        d_oh = min_image_dist(centres, hyd, box)
        same_hh = owner[:, None] == owner[None, :]
        same_oh = np.arange(nmol)[:, None] == owner[None, :]
        if d_hh[~same_hh].min() > args.dmin and d_oh[~same_oh].min() > args.dmin:
            break
    else:
        raise SystemExit("could not place molecules; lower --dmin or --density")

    frac = np.vstack([centres, hyd]) / box % 1.0
    lines = [f"(H2O){nmol} rho={args.density} g/cm3 box={box:.4f} A seed={args.seed}",
             "   1.0",
             f"   {box:14.8f}  {0:14.8f}  {0:14.8f}",
             f"   {0:14.8f}  {box:14.8f}  {0:14.8f}",
             f"   {0:14.8f}  {0:14.8f}  {box:14.8f}",
             "   O    H",
             f"   {nmol}   {2 * nmol}",
             "Direct"]
    lines += [f"  {p[0]:.8f}  {p[1]:.8f}  {p[2]:.8f}" for p in frac]
    with open(args.output, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {args.output}: {nmol} molecules, box {box:.3f} A, "
          f"min O-O {min_image_dist(centres, centres, box)[~np.eye(nmol, dtype=bool)].min():.2f} A "
          f"(attempt {attempt + 1})")


if __name__ == "__main__":
    main()
