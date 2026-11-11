#!/usr/bin/env python3
"""Collect labelled frames (energy, forces, virial) from VASP / CP2K / Gaussian /
OpenMX outputs with dpdata and export them for machine-learning potentials.

Supported input formats (dpdata >= 0.2.x, tested with 1.1):
  vasp/outcar   OUTCAR of relaxations / AIMD
  vasp/xml      vasprun.xml
  cp2k/output   single-point or geo-opt main output (needs the force print)
  cp2k/aimd_output   directory of a CP2K MD run (pos/frc/ener files + output)
  gaussian/log  Gaussian log with Force keyword
  openmx/md     OpenMX MD (System.Name.dat + System.Name.md)
Output formats used here:
  deepmd/npy    DeePMD-kit training directory
  extxyz        extended xyz (MACE, NequIP, GPUMD ... can read it)

Usage:
  python dpdata_convert.py run1/OUTCAR run2/OUTCAR --fmt vasp/outcar --out dataset
  python dpdata_convert.py md/OUTCAR --fmt vasp/outcar --stride 10 --test-frac 0.1
"""
import argparse
import os

import dpdata
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--fmt", default="vasp/outcar")
    ap.add_argument("--out", default="dataset")
    ap.add_argument("--stride", type=int, default=1, help="keep every n-th frame")
    ap.add_argument("--test-frac", type=float, default=0.0, help="fraction kept for testing")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    ms = dpdata.MultiSystems()
    for p in a.inputs:
        ls = dpdata.LabeledSystem(p, fmt=a.fmt)
        if a.stride > 1:
            ls = ls.sub_system(list(range(0, len(ls), a.stride)))
        e = ls["energies"]
        print(f"{p}: {len(ls)} frames, {ls.get_natoms()} atoms, formula {ls.formula}, "
              f"E range {e.min():.3f} .. {e.max():.3f} eV")
        ms.append(ls)

    rng = np.random.default_rng(a.seed)
    train, test = dpdata.MultiSystems(), dpdata.MultiSystems()
    for sysname, s in ms.systems.items():
        idx = rng.permutation(len(s))
        ntest = int(round(a.test_frac * len(s)))
        if ntest:
            test.append(s.sub_system(sorted(idx[:ntest])))
        train.append(s.sub_system(sorted(idx[ntest:])))

    os.makedirs(a.out, exist_ok=True)
    train.to("deepmd/npy", os.path.join(a.out, "train_npy"))
    train.to("extxyz", os.path.join(a.out, "train.xyz"))
    if a.test_frac > 0:
        test.to("deepmd/npy", os.path.join(a.out, "test_npy"))
        test.to("extxyz", os.path.join(a.out, "test.xyz"))
    print(f"train frames: {train.get_nframes()}  test frames: {test.get_nframes()}  -> {a.out}/")


if __name__ == "__main__":
    main()
