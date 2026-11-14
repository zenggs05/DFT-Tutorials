#!/usr/bin/env python3
"""Convert DFT outputs into extended-XYZ training data for NEP or MACE.

Reads every frame (or the last one) of files that ASE understands and that
carry energies and forces: vasprun.xml, OUTCAR, Gaussian .log, extxyz, ...
(for CP2K/ABACUS trajectories use dpdata or the scripts in Scripts/nep/).

--format nep   energy=  virial="9 values, eV"  Properties=...:forces:R:3     (GPUMD nep)
--format mace  REF_energy=  REF_stress="9 values, eV/A^3"  ...:REF_forces:R:3  (mace_run_train)

Sign conventions: ASE stress is positive under tension; virial = -stress * volume,
positive under compression, as expected by NEP (https://gpumd.org/nep/input_files/train_test_xyz.html).

Examples
  python dft_to_extxyz.py run*/vasprun.xml --format nep -o train.xyz
  python dft_to_extxyz.py aimd/vasprun.xml --every 10 --format mace -o mace_train.xyz
  python dft_to_extxyz.py atoms/Mo/OUTCAR atoms/S/OUTCAR --last --format mace --config-type IsolatedAtom -o e0.xyz
  python dft_to_extxyz.py h2o_opt.log --format nep --box 15 -o mol.xyz     # molecule in a 15 A box
"""
import argparse
import sys

import numpy as np
from ase.io import read


def fmt_vec(v):
    return " ".join(f"{x:.10f}" for x in np.asarray(v).reshape(-1))


def frame_lines(atoms, fmt, energy_kind, config_type, box):
    if box and not atoms.pbc.any():
        atoms = atoms.copy()
        atoms.set_cell([box, box, box])
        atoms.center()
        atoms.pbc = True
    try:
        e = atoms.get_potential_energy(force_consistent=(energy_kind == "free"))
    except Exception:
        e = atoms.get_potential_energy()
    f = atoms.get_forces()
    stress = None
    if atoms.pbc.all():
        try:
            stress = atoms.get_stress(voigt=False)
        except Exception:
            stress = None
    periodic = atoms.pbc.any()
    if fmt == "nep" and not periodic:
        raise ValueError("NEP needs a cell: use --box L for molecules")
    ekey, fkey = ("energy", "forces") if fmt == "nep" else ("REF_energy", "REF_forces")
    head = []
    if periodic:
        head.append(f'Lattice="{fmt_vec(atoms.cell.array)}"')
    head.append(f"Properties=species:S:1:pos:R:3:{fkey}:R:3")
    head.append(f"{ekey}={e:.10f}")
    if stress is not None:
        if fmt == "nep":
            head.append(f'virial="{fmt_vec(-stress * atoms.get_volume())}"')
        else:
            head.append(f'REF_stress="{fmt_vec(stress)}"')
    if config_type:
        head.append(f"config_type={config_type}")
    pbc = " ".join("T" if x else "F" for x in atoms.pbc)
    head.append(f'pbc="{pbc}"')
    lines = [str(len(atoms)), " ".join(head)]
    for s, r, fr in zip(atoms.get_chemical_symbols(), atoms.positions, f):
        lines.append(f"{s:2s} {r[0]:16.10f} {r[1]:16.10f} {r[2]:16.10f} {fr[0]:16.10f} {fr[1]:16.10f} {fr[2]:16.10f}")
    return lines, e / len(atoms)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+")
    p.add_argument("--format", required=True, choices=["nep", "mace"])
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--last", action="store_true", help="only the last frame of each file")
    p.add_argument("--every", type=int, default=1, help="take every N-th frame (decorrelate AIMD)")
    p.add_argument("--energy", default="free", choices=["free", "sigma0"],
                   help="free = force-consistent free energy (default, matches the forces); sigma0 = E(sigma->0)")
    p.add_argument("--config-type", default=None, help="e.g. IsolatedAtom for MACE E0 frames")
    p.add_argument("--box", type=float, default=None, help="cubic box (A) for non-periodic molecules")
    p.add_argument("--max-eatom", type=float, default=None, help="drop frames above this energy (eV/atom)")
    a = p.parse_args()

    n_out = n_skip = 0
    with open(a.output, "w") as fh:
        for path in a.files:
            index = "-1" if a.last else f"::{a.every}"
            try:
                frames = read(path, index=index)
            except Exception as err:
                print(f"WARNING: cannot read {path}: {err}", file=sys.stderr)
                continue
            frames = frames if isinstance(frames, list) else [frames]
            for at in frames:
                try:
                    lines, eat = frame_lines(at, a.format, a.energy, a.config_type, a.box)
                except Exception as err:
                    print(f"WARNING: {path}: frame skipped ({err})", file=sys.stderr)
                    n_skip += 1
                    continue
                if a.max_eatom is not None and eat > a.max_eatom:
                    n_skip += 1
                    continue
                fh.write("\n".join(lines) + "\n")
                n_out += 1
    print(f"{n_out} frames written to {a.output} ({n_skip} skipped)")
    print("check: python -c \"from ase.io import read; print(len(read('%s', ':')))\"" % a.output)


if __name__ == "__main__":
    main()
