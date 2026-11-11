#!/bin/bash
# Build (and optionally submit) all strained MoS2 cells for the 2D
# energy-strain method.
#   bash run_strain.sh          # only create directories
#   bash run_strain.sh submit   # create and sbatch every directory
# Directory layout: strain_<mode>/<delta>/{INCAR,KPOINTS,POSCAR,POTCAR,job.slurm}
set -e

MODES="xx yy biaxial shear"
STRAINS="-0.015 -0.010 -0.005 0.000 0.005 0.010 0.015"

# POTCAR is not part of the repository: build it from your own PAW library.
if [ ! -s POTCAR ]; then
    PP=${VASP_PP_PATH:?please set VASP_PP_PATH (directory containing potpaw_PBE)}
    cat "$PP/potpaw_PBE/Mo_pv/POTCAR" "$PP/potpaw_PBE/S/POTCAR" > POTCAR
fi
grep -q "Mo" POTCAR || { echo "POTCAR does not look right"; exit 1; }

for m in $MODES; do
    for s in $STRAINS; do
        d="strain_${m}/${s}"
        mkdir -p "$d"
        cp INCAR KPOINTS POTCAR job.slurm "$d/"
        python3 make_strain.py POSCAR "$m" "$s" "$d/POSCAR"
        if [ "${1:-}" = "submit" ]; then
            (cd "$d" && sbatch job.slurm)
        fi
    done
done
echo "created $(ls -d strain_*/* | wc -l) strained cells"
