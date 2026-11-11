#!/usr/bin/env python3
"""Exchange coupling constant J from a high-spin (HS) and a broken-symmetry (BS)
Gaussian 16 calculation.

Usage examples
--------------
  # both states in one Link1 log: HS = job step 1, BS = job step 3
  python calc_j_yamaguchi.py --hs cu2_dimer_bs_guess_fragment.log:1 \
                             --bs cu2_dimer_bs_guess_fragment.log:3 --smax 1

  # two separate logs (the last step with an SCF energy is used)
  python calc_j_yamaguchi.py --hs hs.log --bs bs.log --smax 1

Conventions
-----------
  H = -2J S_A.S_B (positive J = ferromagnetic, negative J = antiferromagnetic).
  Values for the H = -J S_A.S_B convention are 2x larger and are printed too.
  Yamaguchi:  J = (E_BS - E_HS) / (<S^2>_HS - <S^2>_BS)
  Noodleman:  J = (E_BS - E_HS) / Smax^2               (weak-overlap limit)
  Ruiz/Bencini: J = (E_BS - E_HS) / (Smax (Smax + 1))  (strong-overlap limit)
  Smax = S_A + S_B, e.g. Smax = 1 for two Cu(II) (S = 1/2) centres.
  The Noodleman/Ruiz forms and the projected gap assume S_A = S_B.

<S^2> is taken from "S**2 before annihilation" (the unprojected value).
Only the Python standard library is used (Python >= 3.8).
"""
import argparse
import re
import sys

HARTREE_TO_CM = 219474.6313705
HARTREE_TO_KCAL = 627.509474

RE_SCF = re.compile(r"SCF Done:\s+E\((\S+)\)\s+=\s+(-?\d+\.\d+)")
RE_S2 = re.compile(r"S\*\*2 before annihilation\s+(-?\d+\.\d+),\s+after\s+(-?\d+\.\d+)")
RE_END = re.compile(r"(Normal|Error) termination of Gaussian")


def read_steps(path):
    """Split a (possibly multi-step) log into job steps and return, for each
    step, the last SCF energy and the last <S**2> before/after annihilation."""
    steps, cur = [], {"E": None, "method": None, "S2": None, "S2A": None}
    with open(path, errors="replace") as fh:
        for line in fh:
            m = RE_SCF.search(line)
            if m:
                cur["method"], cur["E"] = m.group(1), float(m.group(2))
                cur["S2"] = cur["S2A"] = None   # S**2 lines follow SCF Done
                continue
            m = RE_S2.search(line)
            if m:
                cur["S2"], cur["S2A"] = float(m.group(1)), float(m.group(2))
                continue
            if RE_END.search(line):
                steps.append(cur)
                cur = {"E": None, "method": None, "S2": None, "S2A": None}
    if cur["E"] is not None:          # unfinished last step
        steps.append(cur)
    return steps


def pick(spec):
    """spec = 'file.log' or 'file.log:N' (N = 1-based job step)."""
    path, step = spec, None
    m = re.match(r"^(.*):(\d+)$", spec)
    if m:
        path, step = m.group(1), int(m.group(2))
    steps = read_steps(path)
    if not steps:
        sys.exit(f"no job steps found in {path}")
    if step is None:
        cand = [s for s in steps if s["E"] is not None]
        if not cand:
            sys.exit(f"no SCF energy found in {path}")
        data = cand[-1]
    else:
        if step < 1 or step > len(steps):
            sys.exit(f"{path}: step {step} requested, but only {len(steps)} step(s) found")
        data = steps[step - 1]
    if data["E"] is None or data["S2"] is None:
        sys.exit(f"{spec}: no SCF energy or no <S**2> in this step "
                 f"(was it an unrestricted calculation?)")
    return data


def main():
    ap = argparse.ArgumentParser(description="J from HS and BS Gaussian logs "
                                 "(Yamaguchi / Noodleman / Ruiz formulas).")
    ap.add_argument("--hs", required=True, help="high-spin log[:step]")
    ap.add_argument("--bs", required=True, help="broken-symmetry log[:step]")
    ap.add_argument("--smax", type=float, required=True,
                    help="Smax = S_A + S_B of the HS state (1 for two S=1/2 centres)")
    a = ap.parse_args()

    hs, bs = pick(a.hs), pick(a.bs)
    smax = a.smax
    dE = bs["E"] - hs["E"]                       # Hartree
    dS2 = hs["S2"] - bs["S2"]
    if abs(dS2) < 1e-6:
        sys.exit("<S^2>_HS equals <S^2>_BS: the BS calculation probably "
                 "collapsed to the HS (or closed-shell) solution")

    print(f"HS : E({hs['method']}) = {hs['E']:.8f} Eh   <S**2> = {hs['S2']:.4f} "
          f"(ideal {smax * (smax + 1):.4f})")
    print(f"BS : E({bs['method']}) = {bs['E']:.8f} Eh   <S**2> = {bs['S2']:.4f} "
          f"(weak-overlap ideal {smax:.4f})")
    print(f"E_BS - E_HS = {dE:.8f} Eh = {dE * HARTREE_TO_CM:.1f} cm-1 "
          f"= {dE * HARTREE_TO_KCAL:.3f} kcal/mol")
    print()
    rows = [("Yamaguchi", dE / dS2),
            ("Noodleman (weak overlap)", dE / smax ** 2),
            ("Ruiz/Bencini (strong overlap)", dE / (smax * (smax + 1)))]
    print(f"{'formula':32s}{'J (H=-2J S.S)':>16s}{'J (H=-J S.S)':>16s}   [cm-1]")
    for name, j in rows:
        print(f"{name:32s}{j * HARTREE_TO_CM:16.1f}{2 * j * HARTREE_TO_CM:16.1f}")
    jy = dE / dS2
    lowest_s = smax % 1.0       # 0 for integer Smax, 1/2 otherwise
    gap = -jy * (lowest_s * (lowest_s + 1) - smax * (smax + 1))
    print()
    print(f"Yamaguchi-projected E(S={lowest_s:g}) - E(S={smax:g}) = "
          f"{gap * HARTREE_TO_CM:.1f} cm-1 ({gap * HARTREE_TO_KCAL:.3f} kcal/mol)")
    print("J < 0: antiferromagnetic (low-spin ground state); J > 0: ferromagnetic.")


if __name__ == "__main__":
    main()
