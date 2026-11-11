#!/usr/bin/env python3
"""Li-Fe-P-O phase diagram, E_above_hull and formation energy of LiFePO4.

Data come from the (new) Materials Project API.  The legacy API was shut
down in 2025; old 16-character keys no longer work.  Get the 32-character
key from https://materialsproject.org/dashboard and export it:
    export MP_API_KEY=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

Usage:
  python pmg_phase_diagram_LiFePO4.py                  # query MP, save entries
  python pmg_phase_diagram_LiFePO4.py --offline        # reuse LiFePO_entries.json
  python pmg_phase_diagram_LiFePO4.py --my-energy -191.3 --my-formula LiFePO4
      # add your own calculation (total energy in eV for that formula unit
      # cell, computed with MP-compatible settings, e.g. MPRelaxSet + GGA+U)
"""
import argparse
import os

from monty.serialization import dumpfn, loadfn
from pymatgen.analysis.phase_diagram import PDEntry, PDPlotter, PhaseDiagram
from pymatgen.core import Composition

CACHE = "LiFePO_entries.json"


def get_entries(offline):
    if offline:
        return loadfn(CACHE)
    from pymatgen.ext.matproj import MPRester   # no mp-api dependency needed
    key = os.environ.get("MP_API_KEY") or os.environ.get("PMG_MAPI_KEY")
    if not key:
        raise SystemExit("set MP_API_KEY (new 32-character key from the MP dashboard)")
    with MPRester(key) as mpr:
        # GGA/GGA+U entries with MaterialsProject2020Compatibility corrections
        entries = mpr.get_entries_in_chemsys(["Li", "Fe", "P", "O"])
        docs = mpr.summary_search(formula="LiFePO4",
                                  _fields=["material_id", "energy_above_hull",
                                           "formation_energy_per_atom", "band_gap",
                                           "symmetry"])
    for d in docs:
        print("MP summary:", d.get("material_id"), "E_hull =", d.get("energy_above_hull"),
              "eV/atom  Ef =", d.get("formation_energy_per_atom"), "eV/atom  gap =",
              d.get("band_gap"), "eV")
    dumpfn(entries, CACHE)
    print(f"{len(entries)} entries cached in {CACHE}")
    return entries


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--my-energy", type=float, help="total energy (eV) of your own cell")
    ap.add_argument("--my-formula", default="LiFePO4", help="composition of that cell")
    a = ap.parse_args()

    entries = list(get_entries(a.offline))
    if a.my_energy is not None:
        # NOTE: your energy must be on the same footing (functional, U values,
        # POTCARs, corrections) as the MP entries, otherwise E_hull is meaningless.
        entries.append(PDEntry(Composition(a.my_formula), a.my_energy, name="my_calc"))

    pd = PhaseDiagram(entries)
    target = Composition("LiFePO4").reduced_composition
    cands = [e for e in entries if e.composition.reduced_composition == target]
    print(f"{'entry_id':>22} {'E_hull(eV/at)':>14} {'E_form(eV/at)':>14}")
    for e in sorted(cands, key=pd.get_e_above_hull):
        eid = getattr(e, "entry_id", None) or getattr(e, "name", "")
        print(f"{str(eid):>22} {pd.get_e_above_hull(e):14.4f} {pd.get_form_energy_per_atom(e):14.4f}")

    best = min(cands, key=pd.get_e_above_hull)
    decomp, ehull = pd.get_decomp_and_e_above_hull(best)
    print("lowest LiFePO4 decomposes to:",
          {k.composition.reduced_formula: round(v, 3) for k, v in decomp.items()},
          f"(E_hull = {ehull:.4f} eV/atom)")
    print("stable phases:", sorted({e.composition.reduced_formula for e in pd.stable_entries}))

    # A quaternary diagram is hard to read in 2D; as an example plot the
    # Li-Fe-O ternary sub-system (PDPlotter handles 2-4 elements).
    sub = PhaseDiagram([e for e in entries
                        if set(e.composition.chemical_system.split("-")) <= {"Li", "Fe", "O"}])
    fig = PDPlotter(sub, backend="matplotlib", show_unstable=0.05).get_plot()
    fig.figure.savefig("LiFeO_pd.png", dpi=300, bbox_inches="tight")
    print("written LiFeO_pd.png (Li-Fe-O ternary sub-diagram)")


if __name__ == "__main__":
    main()
