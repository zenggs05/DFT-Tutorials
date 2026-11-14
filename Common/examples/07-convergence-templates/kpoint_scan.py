#!/usr/bin/env python3
"""k-point convergence scan for VASP (Monkhorst-Pack / Gamma meshes or KSPACING).

The script only uses the Python standard library. It never starts VASP itself
unless you pass --submit, in which case the given command (e.g. "sbatch job.sh")
is executed once inside every new directory.

Examples
--------
  # Gamma-centred meshes scaled to the cell (n subdivisions along the LONGEST b_i,
  # the other directions scaled down to keep a similar k-point density)
  python kpoint_scan.py prepare --template encut_scan --encut 520 --mesh 4 6 8 10 12 14

  # KSPACING values (1/Angstrom, VASP convention including 2*pi); duplicates are skipped
  python kpoint_scan.py prepare --template encut_scan --encut 520 --kspacing 0.5 0.4 0.3 0.25 0.2 0.15

  # collect E0 and write kpoints.dat (columns: label  E0  NIONS  NKPTS  N1xN2xN3)
  python kpoint_scan.py collect --outdir kscan

The template directory must contain INCAR.template (with @ENCUT@), POSCAR and POTCAR.
"""
import argparse
import math
import pathlib
import re
import shutil
import subprocess
import sys


def read_poscar_lattice(path):
    """Return the three lattice vectors (Angstrom) of a VASP POSCAR."""
    lines = pathlib.Path(path).read_text().splitlines()
    scale = [float(x) for x in lines[1].split()]
    vecs = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    if len(scale) == 1 and scale[0] < 0:          # negative value = target volume
        vol = abs(det3(vecs))
        s = (abs(scale[0]) / vol) ** (1.0 / 3.0)
        return [[s * x for x in v] for v in vecs]
    if len(scale) == 3:
        return [[scale[j] * v[j] for j in range(3)] for v in vecs]
    return [[scale[0] * x for x in v] for v in vecs]


def det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))


def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def recip_lengths(lat):
    """|b_i| WITHOUT the 2*pi factor (1/Angstrom)."""
    vol = det3(lat)
    b = [cross(lat[1], lat[2]), cross(lat[2], lat[0]), cross(lat[0], lat[1])]
    return [math.sqrt(sum(x * x for x in v)) / abs(vol) for v in b]


def mesh_from_kspacing(blen, ksp):
    # VASP wiki: N_i = max(1, ceiling(|b_i| * 2*pi / KSPACING))
    return [max(1, math.ceil(b * 2.0 * math.pi / ksp - 1e-8)) for b in blen]


def mesh_from_n(blen, n):
    # n subdivisions along the longest reciprocal vector, others scaled to keep
    # roughly the same k-point density
    bmax = max(blen)
    return [max(1, int(round(n * b / bmax))) for b in blen]


def prepare(args):
    tdir = pathlib.Path(args.template)
    for f in ("INCAR.template", "POSCAR", "POTCAR"):
        if not (tdir / f).is_file():
            sys.exit(f"ERROR: {tdir / f} not found (POTCAR must be created by you)")
    incar_t = (tdir / "INCAR.template").read_text()
    incar_t = incar_t.replace("@ENCUT@", str(args.encut))
    blen = recip_lengths(read_poscar_lattice(tdir / "POSCAR"))
    out = pathlib.Path(args.outdir)
    out.mkdir(exist_ok=True)
    seen = set()
    made = []
    if args.kspacing:
        items = [("ksp_%.3f" % k, mesh_from_kspacing(blen, k), k) for k in args.kspacing]
    else:
        items = [("k_%02d" % n, mesh_from_n(blen, n), None) for n in args.mesh]
    for label, mesh, ksp in items:
        key = tuple(mesh)
        if key in seen:
            print(f"skip {label}: mesh {mesh} already in the scan")
            continue
        seen.add(key)
        d = out / label
        d.mkdir(exist_ok=True)
        shutil.copy(tdir / "POSCAR", d / "POSCAR")
        shutil.copy(tdir / "POTCAR", d / "POTCAR")
        incar = incar_t
        if ksp is not None:
            incar = re.sub(r"(?m)^\s*#?\s*KSPACING.*$", "", incar)
            incar += f"\nKSPACING = {ksp}\nKGAMMA   = .TRUE.\n"
            (d / "KPOINTS").unlink(missing_ok=True)
        else:
            style = "Gamma" if args.gamma else "Monkhorst-Pack"
            (d / "KPOINTS").write_text(
                f"{label}\n0\n{style}\n {mesh[0]} {mesh[1]} {mesh[2]}\n 0 0 0\n")
        (d / "INCAR").write_text(incar)
        (d / "mesh.txt").write_text("%d %d %d\n" % tuple(mesh))
        made.append(d)
        print(f"{label:10s} mesh {mesh[0]:3d} x {mesh[1]:3d} x {mesh[2]:3d}")
    if args.submit:
        for d in made:
            print(f"[{d}] {args.submit}")
            subprocess.run(args.submit, shell=True, cwd=d, check=False)


def last_float(pattern, text, field=-1, first=False):
    hits = [l for l in text.splitlines() if pattern in l]
    if not hits:
        return None
    try:
        return float(hits[0 if first else -1].split()[field])
    except (ValueError, IndexError):
        return None


def collect(args):
    out = pathlib.Path(args.outdir)
    rows = []
    for d in sorted(p for p in out.iterdir() if p.is_dir()):
        oc = d / "OUTCAR"
        if not oc.is_file():
            print(f"WARNING: {oc} missing", file=sys.stderr)
            continue
        txt = oc.read_text(errors="replace")
        if "General timing" not in txt:
            print(f"WARNING: {oc} did not finish normally", file=sys.stderr)
        e0 = last_float("energy  without entropy", txt)
        nions = last_float("NIONS", txt, first=True)
        nk = last_float("NKPTS", txt, 3, first=True)  # "k-points NKPTS = 29 ..."
        mesh = (d / "mesh.txt").read_text().split() if (d / "mesh.txt").is_file() else ["?"] * 3
        if e0 is None:
            continue
        ntot = 1
        for m in mesh:
            ntot *= int(m) if m.isdigit() else 1
        rows.append((ntot, d.name, e0, int(nions or 1), int(nk or 0), "x".join(mesh)))
    rows.sort()
    with open(args.datfile, "w") as fh:
        fh.write("# x=N1*N2*N3  E0(eV)  NIONS  NKPTS(irr)  mesh  label\n")
        for ntot, label, e0, n, nk, mesh in rows:
            fh.write(f"{ntot} {e0:.8f} {n} {nk} {mesh} {label}\n")
    print(f"{len(rows)} points written to {args.datfile}")
    print("plot: python plot_convergence.py %s --xlabel 'N1*N2*N3'" % args.datfile)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("prepare", help="create one directory per mesh")
    a.add_argument("--template", default="encut_scan", help="dir with INCAR.template, POSCAR, POTCAR")
    a.add_argument("--encut", type=float, required=True, help="converged ENCUT (eV) from the ENCUT scan")
    g = a.add_mutually_exclusive_group(required=True)
    g.add_argument("--mesh", type=int, nargs="+", help="subdivisions along the longest b_i")
    g.add_argument("--kspacing", type=float, nargs="+", help="KSPACING values (1/Angstrom)")
    a.add_argument("--mp", dest="gamma", action="store_false", help="Monkhorst-Pack instead of Gamma-centred")
    a.add_argument("--outdir", default="kscan")
    a.add_argument("--submit", default=None, help='command run in every new dir, e.g. "sbatch ../../job_vasp.sh"')
    c = sub.add_parser("collect", help="parse OUTCARs and write a table")
    c.add_argument("--outdir", default="kscan")
    c.add_argument("--datfile", default="kpoints.dat")
    args = p.parse_args()
    prepare(args) if args.cmd == "prepare" else collect(args)


if __name__ == "__main__":
    main()
