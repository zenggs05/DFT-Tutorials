#!/usr/bin/env python3
"""Displace a structure along one normal mode taken from a CP2K Molden file.

Typical uses
  * remove a spurious imaginary mode of a "minimum": push along the mode and
    re-optimize (try both signs, keep the lower-energy result);
  * "poor man's IRC" for a transition state: push +/- along the imaginary mode
    and optimize both structures; they should end in the initial and final state.

The atomic coordinates are read from an xyz file in Angstrom (same atom order
as the frequency run), the mode from the [FR-NORM-COORD] block of
<PROJECT>-VIBRATIONS-1.mol.  The mode is scaled so that the atom that moves
most is displaced by --amp Angstrom.

Example:
    python3 displace_mode.py ts_dimer.xyz ts_freq-VIBRATIONS-1.mol --mode 1 --amp 0.10
    -> writes ts_dimer_mode1_plus.xyz and ts_dimer_mode1_minus.xyz
Only the Python standard library is needed.
"""
import argparse
import math
import os


def read_xyz(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    sym, xyz = [], []
    for line in lines[2:2 + n]:
        p = line.split()
        sym.append(p[0])
        xyz.append([float(v) for v in p[1:4]])
    return lines[1], sym, xyz


def read_mode(path, k):
    """Return (frequency, list of 3-vectors) of vibration k (1-based)."""
    freqs, modes, cur, section = [], {}, None, None
    with open(path, errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("["):
                section = s.upper()
                cur = None
                continue
            if not s:
                continue
            if section == "[FREQ]":
                freqs.append(float(s.split()[0]))
            elif section == "[FR-NORM-COORD]":
                if s.lower().startswith("vibration"):
                    cur = int(s.split()[1])
                    modes[cur] = []
                elif cur is not None:
                    modes[cur].append([float(v) for v in s.split()[:3]])
    if k not in modes:
        raise SystemExit(f"vibration {k} not found in {path} (found {len(modes)} modes)")
    return freqs[k - 1] if k - 1 < len(freqs) else float("nan"), modes[k]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xyz", help="structure used for the frequency run (Angstrom)")
    ap.add_argument("molden", help="<PROJECT>-VIBRATIONS-1.mol")
    ap.add_argument("--mode", type=int, default=1, help="mode index, 1 = lowest (default 1)")
    ap.add_argument("--amp", type=float, default=0.10, help="largest atomic displacement in Angstrom (default 0.10)")
    args = ap.parse_args()

    title, sym, xyz = read_xyz(args.xyz)
    freq, vec = read_mode(args.molden, args.mode)
    if len(vec) != len(xyz):
        raise SystemExit(f"atom count mismatch: xyz has {len(xyz)}, mode has {len(vec)}")
    vmax = max(math.sqrt(sum(c * c for c in v)) for v in vec)
    if vmax == 0.0:
        raise SystemExit("mode vector is zero")
    label = f"{freq:.1f}i" if freq < 0 else f"{freq:.1f}"
    stem = os.path.splitext(os.path.basename(args.xyz))[0]
    for sign, tag in ((+1, "plus"), (-1, "minus")):
        out = f"{stem}_mode{args.mode}_{tag}.xyz"
        with open(out, "w", newline="\n") as fh:
            fh.write(f"{len(sym)}\n{title.strip()} | displaced {tag} along mode {args.mode} ({label} cm^-1), amp {args.amp} A\n")
            for s, r, v in zip(sym, xyz, vec):
                d = [sign * args.amp * c / vmax for c in v]
                fh.write(f"{s:<2s} {r[0] + d[0]:14.6f} {r[1] + d[1]:14.6f} {r[2] + d[2]:14.6f}\n")
        print(f"wrote {out}")
    print(f"mode {args.mode}: {label} cm^-1, max displacement {args.amp} A")


if __name__ == "__main__":
    main()
