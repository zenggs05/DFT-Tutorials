#!/usr/bin/env python3
"""Parse a finished (or interrupted) VASP run with pymatgen.

Reads vasprun.xml (+ OUTCAR / OSZICAR if present) in a directory and prints
energy, convergence, band gap, VBM/CBM, magnetisation; optionally plots the
band structure / DOS and the planar-averaged LOCPOT.

Usage:
  python pmg_parse_vasprun.py RUN_DIR                  # summary only
  python pmg_parse_vasprun.py RUN_DIR --band           # + band.png (line-mode run)
  python pmg_parse_vasprun.py RUN_DIR --dos            # + dos.png (needs LORBIT=11)
  python pmg_parse_vasprun.py RUN_DIR --locpot         # + planar average of LOCPOT
"""
import argparse
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from pymatgen.io.vasp.outputs import Locpot, Oszicar, Outcar, Vasprun  # noqa: E402


def f4(x):
    """Format a float that may be None (e.g. efermi of a truncated vasprun.xml)."""
    return "n/a" if x is None else f"{x:.4f}"


def find(d, name):
    for cand in (name, name + ".gz"):
        p = os.path.join(d, cand)
        if os.path.exists(p):
            return p
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir", nargs="?", default=".")
    ap.add_argument("--band", action="store_true")
    ap.add_argument("--dos", action="store_true")
    ap.add_argument("--locpot", action="store_true")
    a = ap.parse_args()

    xml = find(a.rundir, "vasprun.xml")
    if xml is None:
        raise SystemExit("no vasprun.xml in " + a.rundir)

    # exception_on_bad_xml=False: an interrupted job still yields the
    # ionic steps that were written completely (a warning is issued).
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        vr = Vasprun(xml, parse_potcar_file=False, exception_on_bad_xml=False,
                     parse_projected_eigen=a.band)
    for w in wlist:
        if "XML" in str(w.message) or "incomplete" in str(w.message).lower():
            print("WARNING:", w.message)

    s = vr.final_structure
    print(f"formula            : {s.composition.reduced_formula}  ({len(s)} atoms)")
    print(f"run type           : {vr.run_type}   ISPIN={vr.parameters.get('ISPIN')}")
    print(f"ionic steps        : {len(vr.ionic_steps)}  NSW={vr.parameters.get('NSW')}")
    print(f"converged (e/ion)  : {vr.converged_electronic} / {vr.converged_ionic}")
    print(f"final energy       : {vr.final_energy:.6f} eV  ({vr.final_energy / len(s):.6f} eV/atom)")
    print(f"E-fermi            : {f4(vr.efermi)} eV")
    # Occupation-based estimate; for metals it can report a spurious small
    # "gap" -- use bs.is_metal() / bs.get_band_gap() (--band) to decide.
    try:
        gap, cbm, vbm, direct = vr.eigenvalue_band_properties
        print(f"gap (eigenvalues)  : {f4(gap)} eV  direct={direct}  VBM={f4(vbm)}  CBM={f4(cbm)}")
    except Exception as exc:  # truncated files may lack eigenvalues
        print("gap (eigenvalues)  : not available:", exc)
    print(f"lattice abc (A)    : {s.lattice.abc}")

    outcar = find(a.rundir, "OUTCAR")
    if outcar:
        oc = Outcar(outcar)
        if oc.total_mag is not None:
            print(f"total magnetisation: {oc.total_mag:.3f} muB")
        with open(outcar, errors="ignore") as f:
            print("OUTCAR 'reached required accuracy':", "reached required accuracy" in f.read())
    oszicar = find(a.rundir, "OSZICAR")
    if oszicar:
        osz = Oszicar(oszicar)
        # '******' overflow values are returned as nan with a warning (pymatgen >= 2026)
        print(f"OSZICAR final energy: {osz.final_energy:.6f} eV ({len(osz.ionic_steps)} steps)")

    if a.band:
        from pymatgen.electronic_structure.plotter import BSPlotter
        # line_mode needs the line-mode KPOINTS (with labels) of this run; with
        # KPOINTS_OPT runs (VASP >= 6.3) the KPOINTS_OPT file must be present.
        kp = find(a.rundir, "KPOINTS")
        try:
            bs = vr.get_band_structure(kpoints_filename=kp, line_mode=True, efermi="smart")
        except Exception as exc:
            raise SystemExit(f"cannot build band structure: {exc}")
        bg = bs.get_band_gap()
        print(f"band structure gap : {bg['energy']:.4f} eV  direct={bg['direct']}  {bg['transition']}")
        print(f"is_metal           : {bs.is_metal()}")
        ax = BSPlotter(bs).get_plot(ylim=(-6, 6), vbm_cbm_marker=True)
        ax.figure.savefig("band.png", dpi=300, bbox_inches="tight")
        print("written band.png")

    if a.dos:
        from pymatgen.electronic_structure.plotter import DosPlotter
        cdos = vr.complete_dos
        dp = DosPlotter(sigma=0.05)            # extra Gaussian broadening (eV)
        dp.add_dos("Total", cdos)
        dp.add_dos_dict(cdos.get_element_dos())  # element-projected (needs LORBIT)
        ax = dp.get_plot(xlim=(-8, 6))
        ax.figure.savefig("dos.png", dpi=300, bbox_inches="tight")
        print("written dos.png")

    if a.locpot:
        lp = find(a.rundir, "LOCPOT")
        if lp is None:
            print("no LOCPOT (set LVHAR = .TRUE. in INCAR)")
        else:
            loc = Locpot.from_file(lp)
            z = loc.get_axis_grid(2)
            v = loc.get_average_along_axis(2)
            print(f"vacuum level (max planar avg): {v.max():.3f} eV ; "
                  f"work function = {v.max() - vr.efermi:.3f} eV")
            fig, ax = plt.subplots(figsize=(4, 3))
            ax.plot(z, v)
            ax.axhline(vr.efermi, ls="--", c="gray")
            ax.set_xlabel("z (A)")
            ax.set_ylabel("planar average V (eV)")
            fig.savefig("locpot_z.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
