#!/usr/bin/env python3
"""Change the vacuum (length of a3 along z) of a slab / 2D POSCAR while keeping
the Cartesian geometry of the atoms, and re-centre the layer at c/2.

Assumes a1, a2 lie in the xy plane and a3 is parallel to z (true for the
POSCARs in this repository). Selective-dynamics flags are kept.

Usage: python set_vacuum.py POSCAR --c 25.0 > POSCAR_c25
"""
import argparse

import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("poscar")
    p.add_argument("--c", type=float, required=True, help="new length of a3 (A)")
    a = p.parse_args()

    lines = open(a.poscar).read().splitlines()
    scale = float(lines[1].split()[0])
    if scale <= 0:
        raise SystemExit("negative/zero scaling factor not supported here")
    cell = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * scale
    if abs(cell[2, 0]) > 1e-6 or abs(cell[2, 1]) > 1e-6 or abs(cell[0, 2]) > 1e-6 or abs(cell[1, 2]) > 1e-6:
        raise SystemExit("a3 must be parallel to z and a1, a2 in-plane")
    species, counts = lines[5].split(), [int(x) for x in lines[6].split()]
    nat = sum(counts)
    i = 7
    selective = lines[i].strip().lower().startswith("s")
    if selective:
        i += 1
    direct = lines[i].strip().lower()[0] in "d"
    i += 1
    pos, flags = [], []
    for line in lines[i:i + nat]:
        parts = line.split()
        pos.append([float(x) for x in parts[:3]])
        flags.append(" ".join(parts[3:6]) if selective else "")
    pos = np.array(pos)
    cart = pos @ cell if direct else pos * scale
    z = cart[:, 2]
    thickness = z.max() - z.min()
    if thickness >= a.c:
        raise SystemExit(f"layer thickness {thickness:.2f} A >= requested c")
    cart[:, 2] = z - (z.max() + z.min()) / 2.0 + a.c / 2.0
    new = cell.copy()
    new[2] = [0.0, 0.0, a.c]
    frac = cart @ np.linalg.inv(new)

    out = [lines[0] + f" | c = {a.c:.3f} A, vacuum = {a.c - thickness:.3f} A", "1.0"]
    out += ["  %14.10f %14.10f %14.10f" % tuple(v) for v in new]
    out += [lines[5], lines[6]]
    if selective:
        out.append("Selective dynamics")
    out.append("Direct")
    for f, fl in zip(frac, flags):
        out.append(("  %14.10f %14.10f %14.10f" % tuple(f)) + ("  " + fl if fl else ""))
    print("\n".join(out))


if __name__ == "__main__":
    main()
