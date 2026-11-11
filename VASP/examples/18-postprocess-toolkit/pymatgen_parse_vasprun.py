#!/usr/bin/env python3
"""Quick look at a finished VASP run with pymatgen.

* convergence, final energy, Fermi energy, band gap (from vasprun.xml)
* magnetisation and timing (from OUTCAR, optional)
* band-structure plot for line-mode runs (--band) and element-projected DOS (--dos)
* optional JSON summary for scripts / databases

Requires pymatgen (>= 2024; tested against the API of pymatgen 2026.10, which installs
pymatgen-core) and matplotlib.

Examples
--------
    python pymatgen_parse_vasprun.py                    # current directory
    python pymatgen_parse_vasprun.py path/to/run --dos --json summary.json
    python pymatgen_parse_vasprun.py band_run --band --kpoints band_run/KPOINTS
"""
import argparse
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")

from pymatgen.io.vasp import Outcar, Vasprun  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("directory", nargs="?", default=".")
    p.add_argument("--band", action="store_true", help="plot the band structure (line-mode KPOINTS run)")
    p.add_argument("--kpoints", help="KPOINTS file with labels (default: <directory>/KPOINTS)")
    p.add_argument("--dos", action="store_true", help="plot total + element-projected DOS")
    p.add_argument("--json", help="write a JSON summary to this file")
    args = p.parse_args()

    xml = os.path.join(args.directory, "vasprun.xml")
    if not os.path.isfile(xml):
        sys.exit("no vasprun.xml in %s" % args.directory)

    # parse_potcar_file=False: do not look for the POTCAR (it is not in the repository)
    vr = Vasprun(xml, parse_potcar_file=False, parse_dos=True, parse_eigen=True)
    st = vr.final_structure
    summary = {
        "formula": st.composition.reduced_formula,
        "natoms": len(st),
        "lattice_abc": [round(x, 5) for x in st.lattice.abc],
        "volume": round(st.volume, 4),
        "converged_electronic": vr.converged_electronic,
        "converged_ionic": vr.converged_ionic,
        "final_energy_eV": vr.final_energy,
        "energy_per_atom_eV": vr.final_energy / len(st),
        "efermi_eV": vr.efermi,
        "incar": {k: vr.incar[k] for k in ("ENCUT", "ISMEAR", "SIGMA", "ISPIN", "IBRION", "ISIF") if k in vr.incar},
    }
    gap, cbm, vbm, direct = vr.eigenvalue_band_properties
    summary.update({"band_gap_eV": gap, "cbm_eV": cbm, "vbm_eV": vbm, "gap_is_direct": direct})

    outcar_path = os.path.join(args.directory, "OUTCAR")
    if os.path.isfile(outcar_path):
        oc = Outcar(outcar_path)
        summary["total_magnetization"] = oc.total_mag
        summary["elapsed_time_s"] = oc.run_stats.get("Elapsed time (sec)")
        if oc.magnetization:
            summary["site_magnetization"] = [round(m.get("tot", 0.0), 3) for m in oc.magnetization]

    for k, v in summary.items():
        print("%-22s %s" % (k, v))
    if not (vr.converged_electronic and vr.converged_ionic):
        print("WARNING: the run is not converged - do not use these numbers")

    if args.band:
        from pymatgen.electronic_structure.plotter import BSPlotter

        kp = args.kpoints or os.path.join(args.directory, "KPOINTS")
        bs = vr.get_band_structure(kpoints_filename=kp, line_mode=True)
        info = bs.get_band_gap()
        print("band gap from band structure: %.4f eV (direct: %s, transition: %s)"
              % (info["energy"], info["direct"], info["transition"]))
        ax = BSPlotter(bs).get_plot(zero_to_efermi=True, ylim=(-6, 6))
        ax.figure.savefig(os.path.join(args.directory, "band_pymatgen.png"), dpi=200, bbox_inches="tight")
        print("written band_pymatgen.png")

    if args.dos:
        from pymatgen.electronic_structure.plotter import DosPlotter

        cdos = vr.complete_dos
        plotter = DosPlotter(sigma=0.05)
        plotter.add_dos("Total", cdos)
        plotter.add_dos_dict(cdos.get_element_dos())
        ax = plotter.get_plot(xlim=(-8, 6))
        ax.figure.savefig(os.path.join(args.directory, "dos_pymatgen.png"), dpi=200, bbox_inches="tight")
        print("written dos_pymatgen.png")

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(summary, fh, indent=2, default=str)
        print("written", args.json)


if __name__ == "__main__":
    main()
