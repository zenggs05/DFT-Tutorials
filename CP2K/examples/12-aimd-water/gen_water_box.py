#!/usr/bin/env python3
"""Generate a 32-molecule water box (fcc arrangement of O, random orientations).

The result is only a reasonable *starting* geometry for AIMD: no overlaps,
correct density (0.997 g/cm^3), but no hydrogen-bond network yet.
Equilibrate (NVT) before collecting any statistics.

Usage: python gen_water_box.py [-o water32.xyz] [--seed 2026]
Only the Python standard library is needed.
"""
import argparse
import math
import random

ROH = 0.9572           # Angstrom, experimental gas-phase geometry
THETA = math.radians(104.52)
NA = 6.02214076e23


def box_length(nmol, density=0.997):
    mass = nmol * 18.01528 / NA          # g
    return (mass / density * 1e24) ** (1.0 / 3.0)   # Angstrom


def random_rotation(rng):
    # uniform random rotation from a random unit quaternion
    u1, u2, u3 = rng.random(), rng.random(), rng.random()
    q = (math.sqrt(1 - u1) * math.sin(2 * math.pi * u2),
         math.sqrt(1 - u1) * math.cos(2 * math.pi * u2),
         math.sqrt(u1) * math.sin(2 * math.pi * u3),
         math.sqrt(u1) * math.cos(2 * math.pi * u3))
    x, y, z, w = q
    return [[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]]


def mic(d, L):
    return d - L * round(d / L)


def dist(a, b, L):
    return math.sqrt(sum(mic(a[i] - b[i], L) ** 2 for i in range(3)))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("-o", "--output", default="water32.xyz")
    p.add_argument("--seed", type=int, default=2026)
    args = p.parse_args()
    rng = random.Random(args.seed)

    ncell = 2                      # 2x2x2 fcc cells x 4 sites = 32 molecules
    L = box_length(4 * ncell ** 3)
    a = L / ncell
    basis = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    oxygens = [((i + bx + 0.25) * a, (j + by + 0.25) * a, (k + bz + 0.25) * a)
               for i in range(ncell) for j in range(ncell) for k in range(ncell)
               for (bx, by, bz) in basis]
    # H positions in the molecular frame (O at origin, molecule in xz plane)
    h_local = [(ROH * math.sin(THETA / 2), 0.0, ROH * math.cos(THETA / 2)),
               (-ROH * math.sin(THETA / 2), 0.0, ROH * math.cos(THETA / 2))]
    atoms = []
    hydrogens = []
    for o in oxygens:
        for _ in range(1000):
            R = random_rotation(rng)
            hs = [tuple(o[r] + sum(R[r][c] * h[c] for c in range(3)) for r in range(3))
                  for h in h_local]
            ok = all(dist(h, h2, L) > 1.8 for h in hs for h2 in hydrogens)
            ok = ok and all(dist(h, o2, L) > 2.0 for h in hs for o2 in oxygens if o2 != o)
            if ok:
                break
        else:
            raise SystemExit("could not place molecule without clashes")
        hydrogens.extend(hs)
        atoms.append(("O", o))
        atoms.extend(("H", h) for h in hs)
    with open(args.output, "w", newline="\n") as f:
        f.write("%d\n" % len(atoms))
        f.write("32 H2O, cubic box L = %.4f Angstrom (density 0.997 g/cm3), seed %d\n"
                % (L, args.seed))
        for s, (x, y, z) in atoms:
            f.write("%-2s %12.6f %12.6f %12.6f\n" % (s, x % L, y % L, z % L))
    print("wrote %s  (%d atoms, L = %.4f A)" % (args.output, len(atoms), L))


if __name__ == "__main__":
    main()
