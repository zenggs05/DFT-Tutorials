#!/bin/bash
# sumo_plot.sh - publication-ready band/DOS figures and effective masses with sumo 3.x
# Usage:  bash sumo_plot.sh            (run in Common/examples/05-band-dos-plotting
#                                       after job_band_dos.slurm has finished)
# Needs:  pip install sumo   (sumo >= 3.0, pulls in pymatgen)
# Output: band/*.pdf, dos/*.pdf, band/band.dat, band/sumo-bandstats.log ...
#
# Energy zero: sumo puts 0 eV at the VBM for semiconductors (E-fermi for metals);
# VASPKIT and py4vasp use E-fermi by default. Say which one you use in captions.
set -euo pipefail

for exe in sumo-bandplot sumo-dosplot sumo-bandstats; do
    command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: $exe not found (pip install sumo)"; exit 1; }
done
[ -s band/vasprun.xml ] || { echo "ERROR: band/vasprun.xml missing"; exit 1; }
[ -s dos/vasprun.xml ]  || { echo "ERROR: dos/vasprun.xml missing"; exit 1; }

# --- 0) (before the calculation) k-path generation, for reference only ---------
# sumo-kgen -p POSCAR --pymatgen          # Setyawan-Curtarolo path -> KPOINTS_band
# sumo-kgen -p POSCAR --seekpath          # SeeK-path (Hinuma) convention
# sumo-kgen --hybrid --split 30 --folders # HSE: needs IBZKPT in the directory
# If sumo writes POSCAR_prim, the band run MUST use that primitive cell.

cd band
# --- 1) plain band structure, band edges marked ------------------------------
sumo-bandplot --band-edges --ymin -6 --ymax 4 --width 3.35 --height 2.6 --format pdf

# --- 2) orbital-projected bands (RGB mode: at most three projections) ---------
sumo-bandplot --project Ga.s,Ga.p,As.p --mode rgb --ymin -6 --ymax 4 \
              --prefix proj --format pdf

# --- 3) bands + PDOS side by side (DOS from the dense-mesh run) --------------
sumo-bandplot --dos ../dos/vasprun.xml --elements Ga.s.p.d,As.s.p \
              --ymin -6 --ymax 4 --prefix banddos --format pdf

# --- 4) effective masses from the line-mode band structure ------------------
sumo-bandstats 2>&1 | tee bandstats.txt   # sumo logs to stderr; also writes sumo-bandstats.log
# (default: parabolic fit; add -n for the nonparabolic model, -s N for N sample points)
cd ..

# --- 5) DOS alone: element/orbital resolved, Gaussian broadening 0.05 eV -----
cd dos
sumo-dosplot --elements Ga.s.p.d,As.s.p --gaussian 0.05 --xmin -6 --xmax 4 \
             --prefix gaas --format pdf
cd ..

ls -1 band/*.pdf dos/*.pdf
echo "Raw data for Origin/Excel: band/band.dat, band/proj_band.dat, dos/gaas_total_dos.dat, dos/gaas_<El>_dos.dat"
