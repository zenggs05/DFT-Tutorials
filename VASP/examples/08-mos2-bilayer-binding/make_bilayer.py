#!/usr/bin/env python3
"""Write POSCARs for a rigid interlayer-distance scan of 2H-stacked bilayer MoS2
(and the monolayer reference in the same cell).

Geometry: 1H layers with lattice constant a and S-S height t. In 2H (AA')
stacking the upper layer is rotated by 180 deg: its Mo sits above the S of the
lower layer and vice versa. d is the vertical Mo-Mo distance between the two
layers (bulk 2H-MoS2: d = c_bulk/2 ~ 6.15 A).

Usage:
    python make_bilayer.py --d 6.2 > POSCAR             # bilayer
    python make_bilayer.py --mono > POSCAR              # monolayer reference
Options: --a 3.19 --t 3.13 --c 30.0
"""
import argparse
import math


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--a", type=float, default=3.19, help="in-plane lattice constant (A)")
    p.add_argument("--t", type=float, default=3.13, help="S-S height of one layer (A)")
    p.add_argument("--c", type=float, default=30.0, help="cell length along z (A)")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--d", type=float, help="Mo-Mo interlayer distance (A)")
    g.add_argument("--mono", action="store_true", help="monolayer reference")
    a = p.parse_args()

    A1 = (a.a, 0.0, 0.0)
    A2 = (-a.a / 2.0, a.a * math.sqrt(3.0) / 2.0, 0.0)
    A3 = (0.0, 0.0, a.c)
    P1, P2 = (1.0 / 3.0, 2.0 / 3.0), (2.0 / 3.0, 1.0 / 3.0)

    mo, s = [], []
    if a.mono:
        z0 = 0.5 * a.c
        mo.append((P1, z0))
        s += [(P2, z0 - a.t / 2), (P2, z0 + a.t / 2)]
        title = f"MoS2 monolayer reference a={a.a} t={a.t} c={a.c}"
    else:
        if a.d <= a.t:
            raise SystemExit("d must be larger than the layer height t")
        z1 = 0.5 * a.c - a.d / 2.0          # lower layer centre
        z2 = 0.5 * a.c + a.d / 2.0          # upper layer centre
        mo += [(P1, z1), (P2, z2)]          # 2H: upper Mo above lower S
        s += [(P2, z1 - a.t / 2), (P2, z1 + a.t / 2),
              (P1, z2 - a.t / 2), (P1, z2 + a.t / 2)]
        title = f"MoS2 bilayer 2H (AA') d(Mo-Mo)={a.d:.3f} a={a.a} t={a.t} c={a.c}"

    print(title)
    print("1.0")
    for v in (A1, A2, A3):
        print("  %14.10f %14.10f %14.10f" % v)
    print("Mo S")
    print(f"{len(mo)} {len(s)}")
    print("Direct")
    for (x, y), z in mo + s:
        print("  %14.10f %14.10f %14.10f" % (x, y, z / a.c))


if __name__ == "__main__":
    main()
