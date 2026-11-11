#!/usr/bin/env python3
"""Prepare a POSCAR for a VASP frequency run (IBRION = 5) with selective dynamics.

Only the selected atoms are allowed to move (T T T); every other atom is
frozen (F F F).  Typical use: take the CONTCAR of a relaxed adsorption
structure and free only the adsorbate atoms, so that VASP builds the Hessian
block of the adsorbate only (3 x N_free modes).

Atoms can be selected by element symbol, by 1-based index (ranges allowed),
or by a minimum height (Cartesian z in Angstrom).  The coordinates are copied
verbatim from the input file; only the selective-dynamics flags change.

Examples
    python set_freq_selective.py ../3_co_top/CONTCAR -e C O -o POSCAR
    python set_freq_selective.py CONTCAR -i 17-18 -o POSCAR
    python set_freq_selective.py CONTCAR --zmin 15.0 -o POSCAR
    python set_freq_selective.py CONTCAR -i 37 --zmin 16.0

Criteria are combined with OR: an atom is freed if it matches any of them.

Only the Python standard library is needed (VASP 5/6 POSCAR format).
"""
import argparse
import sys


def parse_indices(tokens):
    """Turn ['1-3', '7'] into {1, 2, 3, 7}."""
    out = set()
    for tok in tokens:
        for part in tok.split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                a, b = part.split("-", 1)
                out.update(range(int(a), int(b) + 1))
            else:
                out.add(int(part))
    return out


def read_poscar(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh]
    comment = lines[0]
    scale = float(lines[1].split()[0])
    cell = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    tok5 = lines[5].split()
    if all(t.lstrip("-").isdigit() for t in tok5):
        sys.exit("ERROR: POSCAR without element line (VASP4 format) is not supported.")
    symbols_line = tok5
    counts = [int(x) for x in lines[6].split()]
    i = 7
    selective = lines[i].strip()[:1].lower() == "s"
    if selective:
        i += 1
    coord_type = lines[i].strip()
    i += 1
    natoms = sum(counts)
    coords = []
    for k in range(natoms):
        parts = lines[i + k].split()
        coords.append(parts[:3])
    symbols = []
    for s, n in zip(symbols_line, counts):
        symbols += [s] * n
    return dict(comment=comment, scale=scale, cell=cell, symbols_line=symbols_line,
                counts=counts, coord_type=coord_type, coords=coords, symbols=symbols,
                header=lines[:7])


def cart_z(p, k):
    """Cartesian z (Angstrom) of atom k, needed only for --zmin."""
    x = [float(v) for v in p["coords"][k]]
    if p["coord_type"][:1].lower() in ("c", "k"):
        return x[2] * p["scale"]
    cell = p["cell"]
    return sum(x[j] * cell[j][2] for j in range(3)) * p["scale"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("poscar", help="input POSCAR/CONTCAR")
    ap.add_argument("-e", "--elements", nargs="+", default=[],
                    help="free all atoms of these elements, e.g. -e C O")
    ap.add_argument("-i", "--indices", nargs="+", default=[],
                    help="free these 1-based atom indices, e.g. -i 17-18 25")
    ap.add_argument("--zmin", type=float, default=None,
                    help="free atoms with Cartesian z >= ZMIN (Angstrom)")
    ap.add_argument("-o", "--output", default="POSCAR_freq",
                    help="output file (default POSCAR_freq)")
    args = ap.parse_args()

    p = read_poscar(args.poscar)
    nat = len(p["symbols"])
    idx = parse_indices(args.indices)
    if any(k < 1 or k > nat for k in idx):
        sys.exit("ERROR: index out of range 1..%d" % nat)

    free = []
    for k in range(nat):
        f = (p["symbols"][k] in args.elements) or ((k + 1) in idx)
        if args.zmin is not None and cart_z(p, k) >= args.zmin:
            f = True
        free.append(f)
    nfree = sum(free)
    if nfree == 0:
        sys.exit("ERROR: no atom selected - nothing would move.")

    out = list(p["header"])
    out.append("Selective dynamics")
    out.append(p["coord_type"])
    for k in range(nat):
        flag = "T T T" if free[k] else "F F F"
        out.append("  %s %s %s   %s" % (p["coords"][k][0], p["coords"][k][1],
                                       p["coords"][k][2], flag))
    with open(args.output, "w", newline="\n") as fh:
        fh.write("\n".join(out) + "\n")

    print("wrote %s: %d atoms, %d free atoms = %d degrees of freedom"
          % (args.output, nat, nfree, 3 * nfree))
    print("IBRION=5, NFREE=2 -> %d displaced structures, %d vibrational modes in OUTCAR"
          % (2 * 3 * nfree, 3 * nfree))
    print("free atoms (1-based):", " ".join("%d%s" % (k + 1, p["symbols"][k])
                                             for k in range(nat) if free[k]))


if __name__ == "__main__":
    main()
