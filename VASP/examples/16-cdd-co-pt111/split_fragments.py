#!/usr/bin/env python3
"""Cut a fragment out of a POSCAR/CONTCAR without touching cell or coordinates.

Used for charge-density-difference (CDD) calculations: the fragments A and B
must have exactly the same lattice and atomic positions as the combined system AB.
Coordinates are copied as text, so no rounding happens. Selective-dynamics flags
are dropped (the fragments are single points).

Examples
--------
    python split_fragments.py POSCAR --keep Pt   --out slab/POSCAR
    python split_fragments.py POSCAR --keep C O  --out co/POSCAR
    python split_fragments.py POSCAR --indices 17-18 --out co/POSCAR   (1-based)
"""
import argparse
import os
import sys


def parse_indices(text):
    out = set()
    for part in text.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-", 1)
            out.update(range(int(a), int(b) + 1))
        elif part:
            out.add(int(part))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("poscar")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--keep", nargs="+", help="element symbols to keep")
    g.add_argument("--indices", help="1-based atom indices to keep, e.g. 17-18 or 1,2,5-8")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    with open(args.poscar) as fh:
        lines = fh.read().splitlines()
    names = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    if len(names) != len(counts):
        sys.exit("line 6 (elements) and line 7 (counts) do not match - VASP5 format required")
    i = 7
    if lines[i].strip()[:1] in "sS":          # Selective dynamics
        i += 1
    mode = lines[i].strip()
    i += 1
    natoms = sum(counts)
    coords = lines[i:i + natoms]
    species = [n for n, c in zip(names, counts) for _ in range(c)]

    if args.keep:
        keep = [k for k in range(natoms) if species[k] in args.keep]
    else:
        sel = parse_indices(args.indices)
        keep = [k for k in range(natoms) if (k + 1) in sel]
    if not keep:
        sys.exit("nothing selected")

    new_names, new_counts, new_coords = [], [], []
    for k in keep:
        if not new_names or new_names[-1] != species[k]:
            new_names.append(species[k])
            new_counts.append(0)
        new_counts[-1] += 1
        new_coords.append("  " + "  ".join(coords[k].split()[:3]))
    if len(set(new_names)) != len(new_names):
        sys.exit("selected atoms of the same element are not contiguous; reorder the POSCAR first")

    out = [lines[0].strip() + " | fragment: " + " ".join(new_names), lines[1]]
    out += lines[2:5]
    out.append("  " + "  ".join(new_names))
    out.append("  " + "  ".join(str(c) for c in new_counts))
    out.append(mode)
    out += new_coords
    d = os.path.dirname(args.out)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(args.out, "w", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")
    print("wrote %s: %s" % (args.out, " ".join("%s%d" % nc for nc in zip(new_names, new_counts))))


if __name__ == "__main__":
    main()
