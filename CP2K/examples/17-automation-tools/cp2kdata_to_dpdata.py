#!/usr/bin/env python3
"""Parse CP2K outputs with cp2kdata and (optionally) convert them to a DeePMD data set.

Chapter 17 of the CP2K tutorial.

Single-point (ENERGY_FORCE) outputs:
    python cp2kdata_to_dpdata.py ef  run1/run.out run2/run.out ... [--deepmd outdir]
MD run directory (needs <project>-pos-1.xyz, -frc-1.xyz, -1.ener and the output file):
    python cp2kdata_to_dpdata.py md  md_dir --output md.out [--deepmd outdir]

Requirements: python >= 3.9, numpy, cp2kdata (pip install cp2kdata); dpdata for --deepmd.

Things to know (cp2kdata 0.7.x):
  * the CP2K output must be written with PRINT_LEVEL MEDIUM (LOW is rejected),
  * forces must be printed (&FORCE_EVAL / &PRINT / &FORCES ON),
  * energies/forces from Cp2kOutput are in CP2K units (hartree, hartree/bohr);
    the dpdata plugin converts them to eV and eV/A,
  * single points whose SCF did not converge are silently returned as empty systems
    by the dpdata plugin -- check the "frames" count printed below,
  * the energy parser looks for "ENERGY| Total FORCE_EVAL ( QS )", i.e. Quickstep
    outputs; QM/MM or FIST outputs are not covered.
"""
import argparse
import sys

import numpy as np

HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANG = 0.529177210903


def summarise_ef(files):
    from cp2kdata import Cp2kOutput

    for f in files:
        out = Cp2kOutput(f)
        e = out.get_energies_list()
        frc = out.get_atomic_forces_list()
        print(f"== {f}")
        print(f"   CP2K version : {out.get_version_string()}")
        print(f"   run type     : {out.get_run_type()}")
        print(f"   atoms        : {out.get_num_atoms()}  symbols: {' '.join(out.get_chemical_symbols())}")
        if e is not None:
            print(f"   energy       : {e[-1]:.10f} Ha = {e[-1] * HARTREE_TO_EV:.6f} eV")
        if frc is not None:
            fmax = np.abs(frc[-1]).max() * HARTREE_TO_EV / BOHR_TO_ANG
            print(f"   max |F_i|    : {fmax:.4f} eV/A")
        else:
            print("   forces       : not found (add &FORCE_EVAL/&PRINT/&FORCES ON)")


def to_deepmd(inputs, mode, outdir, output_name=None):
    import dpdata

    system = None
    if mode == "ef":
        for f in inputs:
            s = dpdata.LabeledSystem(f, fmt="cp2kdata/e_f")
            print(f"{f}: {s.get_nframes()} frame(s)")
            if s.get_nframes() == 0:
                print("   -> skipped (SCF not converged or nothing parsed)", file=sys.stderr)
                continue
            system = s if system is None else system + s
    else:
        system = dpdata.LabeledSystem(inputs[0], cp2k_output_name=output_name, fmt="cp2kdata/md")
        print(f"{inputs[0]}: {system.get_nframes()} frame(s)")
    if system is None or system.get_nframes() == 0:
        sys.exit("nothing to write")
    system.to_deepmd_npy(outdir)
    print(f"wrote {system.get_nframes()} frames to {outdir}/ (deepmd/npy)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=("ef", "md"), help="ef: ENERGY_FORCE outputs, md: MD directory")
    ap.add_argument("inputs", nargs="+", help="output files (ef) or one MD directory (md)")
    ap.add_argument("--output", help="name of the CP2K output file inside the MD directory")
    ap.add_argument("--deepmd", metavar="OUTDIR", help="also write a deepmd/npy data set")
    args = ap.parse_args()

    if args.mode == "ef":
        summarise_ef(args.inputs)
    elif args.output is None:
        sys.exit("md mode needs --output <name of the CP2K output file>")
    if args.deepmd:
        to_deepmd(args.inputs, args.mode, args.deepmd, args.output)


if __name__ == "__main__":
    main()
