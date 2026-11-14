"""Shared helper: build an ASE calculator for a machine-learned potential.

Supported kinds (each needs its own package):
  mace-mp   universal MACE-MP foundation model    pip install mace-torch
            (model "medium" = MACE-MP-0 medium; mace-torch >= 0.3.10 uses
             "medium-mpa-0" when no model is given - we pass "medium" explicitly)
  mace-off  MACE-OFF organic foundation model     pip install mace-torch
  mace      your own trained MACE .model file     pip install mace-torch
  deepmd    frozen DeePMD-kit model (.pb / .pth)  pip install deepmd-kit
  nep       NEP nep.txt through calorine (CPU)    pip install calorine
  emt       ASE EMT toy potential (Cu, Ag, Au, Ni, Pd, Pt, Al, H, C, N, O) for testing scripts
"""


def make_calculator(kind, model=None, device="cpu", dtype="float64", d3=False):
    kind = kind.lower()
    if kind == "mace-mp":
        from mace.calculators import mace_mp
        return mace_mp(model=model or "medium", device=device, default_dtype=dtype, dispersion=d3)
    if kind == "mace-off":
        from mace.calculators import mace_off
        return mace_off(model=model or "medium", device=device, default_dtype=dtype)
    if kind == "mace":
        if not model:
            raise SystemExit("--model path/to/your.model is required for kind 'mace'")
        from mace.calculators import MACECalculator
        return MACECalculator(model_paths=model, device=device, default_dtype=dtype)
    if kind == "deepmd":
        if not model:
            raise SystemExit("--model frozen_model.pb/.pth is required for kind 'deepmd'")
        from deepmd.calculator import DP
        return DP(model=model)
    if kind == "nep":
        from calorine.calculators import CPUNEP
        return CPUNEP(model or "nep.txt")
    if kind == "emt":
        from ase.calculators.emt import EMT
        return EMT()
    raise SystemExit(f"unknown calculator kind: {kind}")


KINDS = ["mace-mp", "mace-off", "mace", "deepmd", "nep", "emt"]
