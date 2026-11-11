#!/usr/bin/env python3
"""Build an explicit KPOINTS file for hybrid-functional (HSE/PBE0) or meta-GGA band
structures with the zero-weight method:

    part 1: the weighted irreducible k-points of a regular mesh, copied from IBZKPT
    part 2: k-points along a band path with weight 0

The band path is read from a VASP line-mode file (same format as a line-mode KPOINTS):
    comment
    N                     <- points per segment (end points included)
    Line-mode
    Reciprocal
    kx ky kz LABEL        <- start of segment 1
    kx ky kz LABEL        <- end of segment 1
    (blank line optional)
    ...

Usage:
    python make_zero_weight_kpoints.py --ibzkpt ../0-pbe-scf/IBZKPT --kpath KPATH --out KPOINTS
    python make_zero_weight_kpoints.py --ibzkpt IBZKPT --kpath KPATH --npts 8   # override N

Only the standard library is needed.
"""
import argparse
import sys


def read_ibzkpt(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    nk = int(lines[1].split()[0])
    mode = lines[2].strip()[:1].lower()
    if mode in ("c", "k"):
        sys.exit("IBZKPT uses Cartesian coordinates; this script expects reciprocal (fractional) ones")
    kpts = []
    for line in lines[3:3 + nk]:
        tok = line.split()
        kpts.append((float(tok[0]), float(tok[1]), float(tok[2]), float(tok[3])))
    if len(kpts) != nk:
        sys.exit("IBZKPT: expected %d k-points, found %d" % (nk, len(kpts)))
    return kpts


def read_kpath(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    npts = int(lines[1].split()[0])
    if not lines[2].strip()[:1].lower() == "l":
        sys.exit("%s: third line must start with 'L' (line mode)" % path)
    if lines[3].strip()[:1].lower() in ("c", "k"):
        sys.exit("%s: Cartesian path coordinates are not supported, use Reciprocal" % path)
    pts = []
    for line in lines[4:]:
        tok = line.split()
        if len(tok) < 3:
            continue
        label = tok[3] if len(tok) > 3 else ""
        pts.append(((float(tok[0]), float(tok[1]), float(tok[2])), label))
    if len(pts) % 2:
        sys.exit("%s: the path must contain pairs of points (start/end of each segment)" % path)
    return npts, [(pts[i], pts[i + 1]) for i in range(0, len(pts), 2)]


def main():
    ap = argparse.ArgumentParser(description="IBZKPT + zero-weight band path -> explicit KPOINTS")
    ap.add_argument("--ibzkpt", default="IBZKPT", help="IBZKPT of a calculation with the SCF mesh")
    ap.add_argument("--kpath", default="KPATH", help="line-mode file with the band path")
    ap.add_argument("--out", default="KPOINTS", help="output file (default KPOINTS)")
    ap.add_argument("--npts", type=int, default=None, help="override points per segment")
    args = ap.parse_args()

    ibz = read_ibzkpt(args.ibzkpt)
    npts, segments = read_kpath(args.kpath)
    if args.npts:
        npts = args.npts
    if npts < 2:
        sys.exit("need at least 2 points per segment")

    path = []
    for (k0, l0), (k1, l1) in segments:
        for i in range(npts):
            t = i / (npts - 1.0)
            k = tuple(a + t * (b - a) for a, b in zip(k0, k1))
            label = l0 if i == 0 else (l1 if i == npts - 1 else "")
            path.append((k, label))

    with open(args.out, "w", newline="\n") as fh:
        fh.write("Explicit k-points: %d weighted (IBZKPT) + %d zero-weight (band path)\n" % (len(ibz), len(path)))
        fh.write("%d\n" % (len(ibz) + len(path)))
        fh.write("Reciprocal lattice\n")
        for kx, ky, kz, w in ibz:
            fh.write("%16.10f %16.10f %16.10f %10.4f\n" % (kx, ky, kz, w))
        for (kx, ky, kz), label in path:
            line = "%16.10f %16.10f %16.10f %10.4f %s" % (kx, ky, kz, 0.0, label)
            fh.write(line.rstrip() + "\n")
    print("wrote %s: %d weighted + %d zero-weight k-points (%d segments x %d points)"
          % (args.out, len(ibz), len(path), len(segments), npts))
    print("cost note: every zero-weight point is diagonalized in every SCF step of the hybrid run")


if __name__ == "__main__":
    main()
