#!/usr/bin/env python3
"""pyprocar_fatband.py - fat bands, PDOS and a combined band+DOS figure with PyProcar 6.x.

PyProcar 6 changed the API: functions take a code name and a *directory*
(code="vasp", dirname="band") instead of a PROCAR path. Old calls such as
pyprocar.bandsplot("PROCAR", outcar="OUTCAR") no longer work.

Files needed in each directory (VASP): PROCAR (LORBIT = 11), OUTCAR, KPOINTS,
POSCAR, vasprun.xml. Always pass the Fermi energy of the SCF/DOS run with
--fermi: the value written by a non-SCF line-mode run can be off.

Examples (layout of this example directory):
    python pyprocar_fatband.py --fermi 3.05                      # all plots
    python pyprocar_fatband.py --fermi 3.05 --mode parametric --atoms 1 --orbitals 1 2 3
    python pyprocar_fatband.py --fermi 3.05 --spins 0            # spin-up only (ISPIN = 2)
    python pyprocar_fatband.py --fermi 3.05 --fermi3d fs_dir     # 3D Fermi surface (metals)

Indices are 0-based. VASP orbital order with LORBIT = 11:
    0 s | 1 py 2 pz 3 px | 4 dxy 5 dyz 6 dz2 7 dxz 8 dx2-y2 | 9-15 f
Requirements: pip install pyprocar  (6.x; the 3D Fermi surface needs pyvista).
"""
import argparse

import matplotlib

matplotlib.use("Agg")

try:
    import pyprocar
except ImportError as exc:  # pragma: no cover
    raise SystemExit("pyprocar is required: pip install pyprocar") from exc


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--band", default="band", help="directory of the line-mode band run")
    p.add_argument("--dos", default="dos", help="directory of the DOS run")
    p.add_argument("--fermi", type=float, required=True,
                   help="E-fermi of the SCF or DOS run (grep E-fermi OUTCAR)")
    p.add_argument("--mode", default="parametric",
                   choices=["plain", "parametric", "scatter", "overlay_species", "overlay_orbitals"],
                   help="band plot mode")
    p.add_argument("--atoms", type=int, nargs="+", default=None, help="0-based atom indices to project on")
    p.add_argument("--orbitals", type=int, nargs="+", default=None, help="0-based orbital indices")
    p.add_argument("--spins", type=int, nargs="+", default=None,
                   help="spin channels: 0 = up, 1 = down for ISPIN = 2")
    p.add_argument("--elimit", type=float, nargs=2, default=(-6.0, 4.0), metavar=("EMIN", "EMAX"))
    p.add_argument("--out", default="gaas", help="prefix of the output PNG files")
    p.add_argument("--fermi3d", default=None, metavar="DIR",
                   help="directory of a uniform-mesh run for a 3D Fermi surface (optional)")
    return p.parse_args()


def main():
    a = parse_args()
    elimit = list(a.elimit)

    # 1) plain bands
    pyprocar.bandsplot(code="vasp", dirname=a.band, mode="plain", fermi=a.fermi,
                       elimit=elimit, savefig=f"{a.out}_bands_plain.png", show=False)

    # 2) projected bands: colour (parametric) or marker size (scatter) = weight
    pyprocar.bandsplot(code="vasp", dirname=a.band, mode=a.mode, fermi=a.fermi,
                       atoms=a.atoms, orbitals=a.orbitals, spins=a.spins,
                       elimit=elimit, savefig=f"{a.out}_bands_{a.mode}.png", show=False)

    # 3) PDOS stacked by species
    pyprocar.dosplot(code="vasp", dirname=a.dos, mode="stack_species", fermi=a.fermi,
                     elimit=elimit, savefig=f"{a.out}_dos_species.png", show=False)

    # 4) bands + DOS in one figure (the two settings dicts are passed to
    #    bandsplot/dosplot; both need 'mode' and 'dirname')
    bands_settings = {"dirname": a.band, "mode": "plain", "fermi": a.fermi}
    dos_settings = {"dirname": a.dos, "mode": "plain", "fermi": a.fermi}
    pyprocar.bandsdosplot(code="vasp", bands_settings=bands_settings, dos_settings=dos_settings,
                          elimit=elimit, savefig=f"{a.out}_bands_dos.png", show=False)

    # 5) optional 3D Fermi surface (metals; needs a dense uniform mesh run)
    if a.fermi3d:
        fs = pyprocar.FermiHandler(code="vasp", dirname=a.fermi3d, fermi=a.fermi)
        fs.plot_fermi_surface(mode="plain", show=False, save_2d=f"{a.out}_fermi3d.png")

    print("done: see", a.out + "_*.png")


if __name__ == "__main__":
    main()
