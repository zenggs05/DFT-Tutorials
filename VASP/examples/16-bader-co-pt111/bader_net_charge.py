#!/usr/bin/env python3
"""Convert Bader populations (ACF.dat) into net atomic charges.

net charge = ZVAL(element, read from POTCAR) - Bader population (column CHARGE of ACF.dat)

Only the Python standard library is used.

Examples
--------
    python bader_net_charge.py
    python bader_net_charge.py --acf ACF.dat --poscar POSCAR --potcar POTCAR --group 17-18
    python bader_net_charge.py --group 17-18 --group 13,14,15,16 --csv net_charges.csv
"""
import argparse
import re
import sys


def read_poscar_species(path):
    """Return the list of element symbols, one per atom (VASP5 POSCAR/CONTCAR)."""
    with open(path) as fh:
        lines = [fh.readline() for _ in range(7)]
    names = lines[5].split()
    counts = lines[6].split()
    if not names or not all(c.isdigit() for c in counts):
        sys.exit("POSCAR line 6 must contain element symbols (VASP5 format) and line 7 the counts")
    species = []
    for name, cnt in zip(names, counts):
        # labels such as 'Fe_pv' or 'Fe/abc' -> 'Fe'
        species.extend([re.split(r"[_/]", name)[0]] * int(cnt))
    return names, [int(c) for c in counts], species


def read_potcar_zval(path):
    """Return (titles, zvals) in the order of the POTCAR blocks."""
    titles, zvals = [], []
    with open(path, errors="replace") as fh:
        for line in fh:
            if "TITEL" in line:
                titles.append(line.split("=", 1)[1].strip())
            m = re.search(r"ZVAL\s*=\s*([-+0-9.Ee]+)", line)
            if m:
                zvals.append(float(m.group(1)))
    if not zvals:
        sys.exit("no ZVAL found in %s" % path)
    return titles, zvals


def read_acf(path):
    """Return list of (index, population) and the trailer dictionary of ACF.dat."""
    rows, trailer = [], {}
    with open(path) as fh:
        for line in fh:
            tok = line.split()
            if len(tok) >= 5 and tok[0].isdigit():
                rows.append((int(tok[0]), float(tok[4])))
            elif ":" in line:
                key, val = line.split(":", 1)
                try:
                    trailer[key.strip()] = float(val.split()[0])
                except (ValueError, IndexError):
                    pass
    return rows, trailer


def parse_group(text, natoms):
    """'17-18' or '1,3,5-7' (1-based) -> sorted list of indices."""
    idx = set()
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            idx.update(range(int(a), int(b) + 1))
        else:
            idx.add(int(part))
    bad = [i for i in idx if i < 1 or i > natoms]
    if bad:
        sys.exit("group %s contains invalid atom indices %s" % (text, bad))
    return sorted(idx)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--acf", default="ACF.dat")
    p.add_argument("--poscar", default="POSCAR", help="POSCAR/CONTCAR with the same atom order as the run")
    p.add_argument("--potcar", default="POTCAR")
    p.add_argument("--group", action="append", default=[],
                   help="1-based atom indices to sum, e.g. 17-18 or 1,2,5-8 (can be repeated)")
    p.add_argument("--csv", help="optional CSV output file")
    args = p.parse_args()

    names, counts, species = read_poscar_species(args.poscar)
    titles, zvals = read_potcar_zval(args.potcar)
    if len(zvals) != len(names):
        sys.exit("POTCAR has %d blocks but POSCAR lists %d species - check the element order!"
                 % (len(zvals), len(names)))
    for name, title in zip(names, titles):
        tt = title.split()
        if len(tt) > 1 and re.split(r"[_/]", name)[0] not in tt[1]:
            print("WARNING: POSCAR species %s vs POTCAR '%s' - check the order" % (name, title))

    zval_atom = []
    for z, cnt in zip(zvals, counts):
        zval_atom.extend([z] * cnt)

    rows, trailer = read_acf(args.acf)
    if len(rows) != len(species):
        sys.exit("ACF.dat has %d atoms, POSCAR has %d" % (len(rows), len(species)))

    print("%5s %4s %12s %8s %10s" % ("atom", "el", "Bader_pop", "ZVAL", "net_q(e)"))
    out = []
    for (i, pop), el, z in zip(rows, species, zval_atom):
        q = z - pop
        out.append((i, el, pop, z, q))
        print("%5d %4s %12.4f %8.3f %+10.4f" % (i, el, pop, z, q))

    total_pop = sum(r[2] for r in out)
    print("-" * 44)
    print("sum of Bader populations : %10.4f" % total_pop)
    print("sum of ZVAL (= NELECT)   : %10.4f" % sum(zval_atom))
    if "NUMBER OF ELECTRONS" in trailer:
        print("ACF.dat NUMBER OF ELECTRONS: %10.4f" % trailer["NUMBER OF ELECTRONS"])
    if abs(total_pop - sum(zval_atom)) > 0.05:
        print("WARNING: populations do not add up to NELECT (charged cell, grid too coarse or vacuum option?)")

    for g in args.group:
        idx = parse_group(g, len(out))
        qsum = sum(out[i - 1][4] for i in idx)
        label = "".join(out[i - 1][1] for i in idx) if len(idx) <= 6 else "%d atoms" % len(idx)
        print("group %-12s (%s): net charge = %+8.4f e" % (g, label, qsum))

    if args.csv:
        with open(args.csv, "w") as fh:
            fh.write("atom,element,bader_population,zval,net_charge\n")
            for r in out:
                fh.write("%d,%s,%.5f,%.3f,%.5f\n" % r)
        print("written", args.csv)


if __name__ == "__main__":
    main()
