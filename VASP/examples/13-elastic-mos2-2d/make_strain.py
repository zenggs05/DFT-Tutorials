#!/usr/bin/env python3
"""Apply a homogeneous in-plane strain to a 2D POSCAR (vacuum along z).

    python make_strain.py POSCAR_in MODE DELTA POSCAR_out

MODE  : xx | yy | biaxial | shear
DELTA : strain magnitude, e.g. 0.010 for 1 %

Strain tensors (symmetric, the z row/column is never touched):
    xx      : [[d, 0], [0, 0]]
    yy      : [[0, 0], [0, d]]
    biaxial : [[d, 0], [0, d]]
    shear   : [[0, d/2], [d/2, 0]]   (engineering shear strain gamma = d)

New lattice (row vectors): A' = A (I + eps).  Fractional coordinates are kept,
so atoms move affinely with the cell; VASP then relaxes them (ISIF=2).
"""
import sys

import numpy as np


def read_poscar(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)])
    if scale < 0:  # negative scale = target volume
        lat *= (abs(scale) / abs(np.linalg.det(lat))) ** (1 / 3)
    else:
        lat *= scale
    return lines, lat


def strain_tensor(mode, d):
    e = np.zeros((3, 3))
    if mode == "xx":
        e[0, 0] = d
    elif mode == "yy":
        e[1, 1] = d
    elif mode == "biaxial":
        e[0, 0] = e[1, 1] = d
    elif mode == "shear":
        e[0, 1] = e[1, 0] = d / 2
    else:
        sys.exit(f"unknown mode '{mode}' (use xx, yy, biaxial, shear)")
    return e


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    src, mode, delta, dst = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
    lines, lat = read_poscar(src)
    new = lat @ (np.eye(3) + strain_tensor(mode, delta))
    out = [f"{lines[0].strip()} | strain {mode} {delta:+.4f}", "   1.0"]
    out += ["  " + "  ".join(f"{x:16.10f}" for x in row) for row in new]
    out += lines[5:]  # species, counts, (Selective), Direct/Cartesian, coords
    if any(l.strip().lower().startswith(("c", "k")) for l in lines[7:9]):
        sys.exit("please provide fractional (Direct) coordinates")
    with open(dst, "w", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
