#!/usr/bin/env python3
"""Build a slab with ASE and print OpenMX blocks in the ESM convention.

OpenMX ESM requires the FIRST unit vector (the a-axis) to be the surface
normal, parallel to a Cartesian axis (ESM.direction x|y|z, default x), with
b and c orthogonal to it.  ASE builds slabs with the normal along z, so this
script rotates/permutes the coordinates accordingly, centres the slab, and can
add layers of empty atoms (species E, E.vps) on both sides.

Examples
    # 6-layer Al(111), 50 Ang cell along the normal, 2 empty layers per side
    python make_slab_omx.py --surface fcc111 --element Al --a 4.05 \
        --layers 6 --length 50 --empty 2 > al111_blocks.txt

    # 4-layer Pt(111), normal along z (ESM.direction z)
    python make_slab_omx.py --surface fcc111 --element Pt --a 3.92 \
        --layers 4 --length 40 --direction z

The output contains Atoms.Number, the coordinate block (Ang) and the
Atoms.UnitVectors block; paste them into your .dat file.  Initial up/down
electrons are zval/2 each; give --zval for elements not in the small table.
Requires: numpy, ase.
"""
import argparse
import sys

import numpy as np

# valence electrons of the DFT_DATA19 *_PBE19 pseudopotentials used most often for slabs
ZVAL = {"Al": 3.0, "Pt": 16.0, "Au": 17.0, "Ag": 17.0, "Pd": 16.0, "Ir": 15.0,
        "Mg": 8.0, "Si": 4.0, "C": 4.0, "Ti": 12.0, "W": 12.0, "Mo": 14.0}


def build(args):
    from ase import build as ab
    func = getattr(ab, args.surface, None)
    if func is None:
        sys.exit(f"error: ase.build has no surface function '{args.surface}'")
    nlay = args.layers + 2 * args.empty
    kw = dict(size=(args.nx, args.ny, nlay), a=args.a, vacuum=None, orthogonal=False)
    if args.c is not None:
        kw["c"] = args.c
    try:
        slab = func(args.element, **kw)
    except TypeError:
        kw.pop("orthogonal")
        slab = func(args.element, **kw)
    return slab


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--surface", default="fcc111", help="ase.build function: fcc111, fcc100, bcc110, hcp0001 ...")
    ap.add_argument("--element", required=True)
    ap.add_argument("--a", type=float, required=True, help="bulk lattice constant (Ang)")
    ap.add_argument("--c", type=float, default=None, help="c lattice constant for hcp")
    ap.add_argument("--layers", type=int, default=6, help="number of real atomic layers")
    ap.add_argument("--empty", type=int, default=0, help="empty-atom layers added on EACH side")
    ap.add_argument("--nx", type=int, default=1)
    ap.add_argument("--ny", type=int, default=1)
    ap.add_argument("--length", type=float, default=40.0, help="cell length along the normal (Ang)")
    ap.add_argument("--direction", choices=["x", "y", "z"], default="x", help="value of ESM.direction")
    ap.add_argument("--zval", type=float, default=None, help="valence electrons of the VPS")
    ap.add_argument("--empty-name", default="E", help="species name of empty atoms")
    args = ap.parse_args()

    zval = args.zval if args.zval is not None else ZVAL.get(args.element)
    if zval is None:
        sys.exit(f"error: give --zval for {args.element}")

    slab = build(args)
    cell = slab.get_cell().array
    pos = slab.get_positions()
    tags = slab.get_tags()            # ASE: layer 1 = top ... layer N = bottom
    nlay = args.layers + 2 * args.empty
    is_empty = (tags <= args.empty) | (tags > nlay - args.empty)

    # in-plane vectors of the ASE cell lie in the xy plane, normal along z
    v1, v2 = cell[0].copy(), cell[1].copy()
    v1[2] = v2[2] = 0.0
    z = pos[:, 2]
    zmid = 0.5 * (z.max() + z.min())
    pos = pos.copy()
    pos[:, 2] = z - zmid + 0.5 * args.length   # centre the slab inside [0, length]
    if pos[:, 2].min() < 0.0 or pos[:, 2].max() > args.length:
        sys.exit("error: slab does not fit into --length")

    # cyclic permutation keeps the cell right-handed: (x, y, z) -> new frame
    perm = {"x": (2, 0, 1), "y": (1, 2, 0), "z": (0, 1, 2)}[args.direction]

    def rot(v):
        out = np.zeros(3)
        for new_axis, old_axis in enumerate(perm):
            out[new_axis] = v[old_axis]
        return out

    normal = rot(np.array([0.0, 0.0, args.length]))
    b, c = rot(v1), rot(v2)
    newpos = np.array([rot(p) for p in pos])
    order = np.argsort(newpos[:, "xyz".index(args.direction)])

    gap_lo = pos[:, 2].min()
    gap_hi = args.length - pos[:, 2].max()
    print(f"# {args.surface} {args.element} a={args.a} layers={args.layers} empty/side={args.empty} "
          f"ESM.direction={args.direction}")
    print(f"# distance of outermost atoms to the ESM boundaries: {gap_lo:.2f} / {gap_hi:.2f} Ang")
    print(f"Atoms.Number                  {len(slab)}")
    print("Atoms.SpeciesAndCoordinates.Unit   Ang")
    print("<Atoms.SpeciesAndCoordinates")
    for i, j in enumerate(order, 1):
        p = newpos[j]
        if is_empty[j]:
            name, up = args.empty_name, 0.0
        else:
            name, up = args.element, 0.5 * zval
        print(f"  {i:4d}  {name:<3s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}   {up:5.2f} {up:5.2f}")
    print("Atoms.SpeciesAndCoordinates>")
    print("Atoms.UnitVectors.Unit        Ang")
    print("<Atoms.UnitVectors")
    for v in (normal, b, c):
        print(f"  {v[0]:12.6f} {v[1]:12.6f} {v[2]:12.6f}")
    print("Atoms.UnitVectors>")


if __name__ == "__main__":
    main()
