#!/usr/bin/env python3
"""Search commensurate supercells of two hexagonal 2D layers and build a
van der Waals heterostructure (bottom layer fixed, top layer rotated and
strained to fit).

  python build_heterostructure.py --list                 # candidates only
  python build_heterostructure.py --pair graphene-MoS2 --pick 0
  python build_heterostructure.py --pair MoS2-WS2         # 1x1, tiny strain
  python build_heterostructure.py --pair graphene-hBN --tol 0.02 --pick 1

Ideas (a mini version of the Zur-McGill / ZSL lattice matching):
  * a supercell vector of a hexagonal lattice is v = i*a1 + j*a2; with the
    120-degree convention |v| = a*sqrt(i^2 - i*j + j^2) and the supercell
    matrix [[i, j], [-j, i-j]] contains n = i^2 - i*j + j^2 cells;
  * two layers match when |v1| ~ |v2|; strain = |v1|/|v2| - 1 is put on
    the TOP layer; the twist angle is the angle between v1 and v2;
  * candidates are sorted by the total number of atoms.

Default lattice constants (approximate experimental values; replace them
by values optimized with YOUR functional + vdW correction):
  graphene 2.46 A, h-BN 2.50 A, MoS2 3.16 A, WS2 3.15 A.
The interlayer distance (default 3.4 A) is a starting guess only - relax
the structure with a vdW-corrected functional.
"""
import argparse
import math

import numpy as np
from ase.build import graphene, make_supercell, mx2

LAYERS = {
    "graphene": lambda: graphene("C2", a=2.46, vacuum=5.0),
    "hBN": lambda: graphene("BN", a=2.50, vacuum=5.0),
    "MoS2": lambda: mx2("MoS2", "2H", a=3.16, thickness=3.17, vacuum=5.0),
    "WS2": lambda: mx2("WS2", "2H", a=3.15, thickness=3.14, vacuum=5.0),
}


def hex_vectors(a, nmax):
    """All supercell vectors (i, j) of a 120-degree hexagonal lattice."""
    a1 = np.array([a, 0.0])
    a2 = np.array([-a / 2.0, a * math.sqrt(3.0) / 2.0])
    out = []
    for i in range(0, nmax + 1):
        for j in range(-nmax, nmax + 1):
            if i == 0 and j <= 0:
                continue
            v = i * a1 + j * a2
            n = i * i - i * j + j * j
            ang = math.degrees(math.atan2(v[1], v[0]))
            out.append((float(np.linalg.norm(v)), ang, i, j, n))
    return out


def candidates(lay1, lay2, nmax, tol, natoms_max):
    a1 = lay1.cell.lengths()[0]
    a2 = lay2.cell.lengths()[0]
    v1s, v2s = hex_vectors(a1, nmax), hex_vectors(a2, nmax)
    seen, res = set(), []
    for L1, ang1, i1, j1, n1 in v1s:
        for L2, ang2, i2, j2, n2 in v2s:
            strain = L1 / L2 - 1.0
            if abs(strain) > tol:
                continue
            twist = (ang1 - ang2) % 60.0
            twist = min(twist, 60.0 - twist)
            nat = n1 * len(lay1) + n2 * len(lay2)
            if nat > natoms_max:
                continue
            key = (n1, n2, round(twist, 2))
            if key in seen:
                continue
            seen.add(key)
            res.append(dict(n1=n1, n2=n2, ij1=(i1, j1), ij2=(i2, j2),
                            L=L1, strain=strain, twist=twist,
                            rot=ang1 - ang2, natoms=nat))
    res.sort(key=lambda c: (c["natoms"], abs(c["strain"])))
    return res


def supercell(atoms, ij):
    i, j = ij
    P = [[i, j, 0], [-j, i - j, 0], [0, 0, 1]]
    return make_supercell(atoms, P)


def build(lay1, lay2, cand, dist, vacuum):
    bot = supercell(lay1, cand["ij1"])
    top = supercell(lay2, cand["ij2"])
    # rotate the top layer so that its first supercell vector is parallel
    # to the one of the bottom layer
    top.rotate(cand["rot"], "z", rotate_cell=True)
    # strain the top layer onto the bottom in-plane cell
    cell = top.cell[:].copy()
    cell[0], cell[1] = bot.cell[0], bot.cell[1]
    top.set_cell(cell, scale_atoms=True)
    # stack: top layer 'dist' above the highest atom of the bottom layer
    zb = bot.positions[:, 2]
    zt = top.positions[:, 2]
    top.positions[:, 2] += zb.max() + dist - zt.min()
    het = bot.copy()
    het.extend(top)
    het.center(vacuum=vacuum / 2.0, axis=2)
    het.pbc = True
    return het, len(bot)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pair", default="graphene-MoS2",
                   help="bottom-top, from: %s" % ", ".join(LAYERS))
    p.add_argument("--nmax", type=int, default=7)
    p.add_argument("--tol", type=float, default=0.03,
                   help="max |strain| on the top layer (0.03 = 3%%)")
    p.add_argument("--natoms-max", type=int, default=150)
    p.add_argument("--list", action="store_true")
    p.add_argument("--pick", type=int, default=0,
                   help="index in the candidate table")
    p.add_argument("--dist", type=float, default=3.4,
                   help="initial interlayer distance in A")
    p.add_argument("--vacuum", type=float, default=20.0)
    args = p.parse_args()

    name1, name2 = args.pair.split("-")
    lay1, lay2 = LAYERS[name1](), LAYERS[name2]()
    cands = candidates(lay1, lay2, args.nmax, args.tol, args.natoms_max)
    print("%s (bottom, a=%.3f) / %s (top, a=%.3f): %d candidates"
          % (name1, lay1.cell.lengths()[0], name2, lay2.cell.lengths()[0],
             len(cands)))
    print("%4s %10s %10s %8s %9s %8s %7s"
          % ("idx", "bottom", "top", "L(A)", "strain%", "twist", "atoms"))
    for k, c in enumerate(cands[:15]):
        print("%4d %10s %10s %8.3f %9.2f %8.2f %7d"
              % (k, "%dx(%d,%d)" % (c["n1"], *c["ij1"]),
                 "%dx(%d,%d)" % (c["n2"], *c["ij2"]), c["L"],
                 100 * c["strain"], c["twist"], c["natoms"]))
    if args.list or not cands:
        return
    cand = cands[args.pick]
    het, nbot = build(lay1, lay2, cand, args.dist, args.vacuum)
    fn = "POSCAR_%s_%s" % (name1, name2)
    from ase.io import write
    with open(fn, "w", newline="\n") as f:
        write(f, het, format="vasp", direct=True, sort=True)
    z = het.positions[:, 2]
    print("\nBuilt candidate %d: %s, %d atoms (%d bottom + %d top)"
          % (args.pick, het.get_chemical_formula(), len(het), nbot,
             len(het) - nbot))
    print("in-plane a = b = %.4f A, gamma = %.1f deg, c = %.2f A, "
          "vacuum gap = %.2f A" % (het.cell.lengths()[0],
                                   het.cell.angles()[2], het.cell[2, 2],
                                   het.cell[2, 2] - (z.max() - z.min())))
    print("top-layer strain %.2f %%, twist %.2f deg -> %s"
          % (100 * cand["strain"], cand["twist"], fn))


if __name__ == "__main__":
    main()
