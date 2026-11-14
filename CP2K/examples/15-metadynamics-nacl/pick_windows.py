#!/usr/bin/env python3
"""Pick starting structures for blue-moon windows from a (cheap) trajectory.

Typical use (chapter 15, "xTB pre-sampling -> DFT refinement"): a GFN1-xTB
well-tempered metadynamics run has visited the whole range of d(Na-Cl).  For
every target distance this script takes the frame whose distance is closest to
the target and writes it to window_<target>/start.xyz.  The DFT windows can then
be run independently (in parallel) with

    CHAIN=0 bash run_bluemoon.sh            # or one job per window

CONSTRAINT_INIT T in constrained_md_distance.inp moves the two atoms exactly
onto the target before the first MD step, so a small mismatch is harmless.

The distance uses the minimum-image convention for an orthorhombic box.
The window directory names use the target strings exactly as given, so pass the
same strings as TARGETS in run_bluemoon.sh (e.g. 2.60, not 2.6).

Example:
    python3 pick_windows.py nacl_wtmetad_xtb-pos-1.xyz --box 12.42 \\
        --targets 2.60 2.80 3.00 3.20 3.40 3.60 3.80 4.00 4.20 4.40 4.60 4.80 5.00 5.20
Only the Python standard library is needed.
"""
import argparse
import math
import os
import sys


def frames(path):
    """Yield (comment, lines) for every frame of a multi-frame xyz file."""
    with open(path, errors="replace") as fh:
        while True:
            head = fh.readline()
            if not head:
                return
            if not head.strip():
                continue
            try:
                nat = int(head.split()[0])
            except ValueError:
                sys.exit(f"{path}: cannot read the atom count from line '{head.strip()}'")
            comment = fh.readline()
            atoms = [fh.readline() for _ in range(nat)]
            if len(atoms) != nat or not atoms[-1]:
                print(f"# warning: incomplete last frame in {path} ignored", file=sys.stderr)
                return
            yield comment, atoms


def mic_distance(a, b, box):
    d2 = 0.0
    for k in range(3):
        d = b[k] - a[k]
        if box[k] > 0:
            d -= box[k] * round(d / box[k])
        d2 += d * d
    return math.sqrt(d2)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("traj", help="multi-frame xyz trajectory (e.g. <PROJECT>-pos-1.xyz)")
    ap.add_argument("--targets", nargs="+", required=True, help="target distances in Angstrom")
    ap.add_argument("--atoms", nargs=2, type=int, default=[1, 2], help="1-based atom indices, default 1 2")
    ap.add_argument("--box", nargs="+", type=float, required=True,
                    help="orthorhombic box a [b c] in Angstrom (0 = not periodic)")
    ap.add_argument("--skip", type=int, default=0, help="ignore the first N frames")
    ap.add_argument("--tol", type=float, default=0.10, help="warn if the best frame is further than this (A)")
    ap.add_argument("--outdir", default=".", help="where to create the window_<target> directories")
    ap.add_argument("--dry-run", action="store_true", help="only print the selection")
    args = ap.parse_args()

    box = args.box * 3 if len(args.box) == 1 else args.box
    if len(box) != 3:
        sys.exit("--box needs 1 or 3 values")
    i1, i2 = args.atoms[0] - 1, args.atoms[1] - 1
    targets = [(t, float(t)) for t in args.targets]
    best = {t: (float("inf"), None, None, None) for t, _ in targets}

    nframes = 0
    for n, (comment, atoms) in enumerate(frames(args.traj)):
        nframes += 1
        if n < args.skip:
            continue
        try:
            p1 = [float(x) for x in atoms[i1].split()[1:4]]
            p2 = [float(x) for x in atoms[i2].split()[1:4]]
        except (IndexError, ValueError):
            sys.exit(f"frame {n + 1}: cannot read atoms {args.atoms}")
        d = mic_distance(p1, p2, box)
        for t, tv in targets:
            if abs(d - tv) < best[t][0]:
                best[t] = (abs(d - tv), n + 1, d, (comment, atoms))
    if nframes == 0:
        sys.exit(f"no frames in {args.traj}")
    print(f"# {nframes} frames read from {args.traj}")
    print("# target(A)  frame   d(A)    |d - target|")
    for t, tv in targets:
        dev, frame, d, data = best[t]
        if frame is None:
            print(f"  {t:>8}   no frame available (check --skip)")
            continue
        flag = "   <-- far from target, sample more or check the CV range" if dev > args.tol else ""
        print(f"  {t:>8}  {frame:6d}  {d:7.3f}  {dev:7.3f}{flag}")
        if args.dry_run:
            continue
        wdir = os.path.join(args.outdir, f"window_{t}")
        os.makedirs(wdir, exist_ok=True)
        comment, atoms = data
        with open(os.path.join(wdir, "start.xyz"), "w") as fh:
            fh.write(f"{len(atoms)}\n")
            fh.write(f"picked from {os.path.basename(args.traj)} frame {frame}, d = {d:.4f} A, target {t} A\n")
            fh.writelines(atoms)
    if not args.dry_run:
        print("# wrote window_<target>/start.xyz; run them with CHAIN=0 bash run_bluemoon.sh")


if __name__ == "__main__":
    main()
