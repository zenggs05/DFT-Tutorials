#!/bin/bash
# Create three sub-directories (FM / AFM / NM) from the same POSCAR/KPOINTS,
# differing only in the spin settings, and submit them with job.slurm.
# Usage:  bash run_fm_afm_nm.sh            (prepare + submit)
#         bash run_fm_afm_nm.sh collect    (print energies and moments)
set -e

if [ "$1" == "collect" ]; then
    printf "%-5s %16s %14s\n" "conf" "E0(sigma->0)/eV" "mag(total)/muB"
    for c in FM AFM NM; do
        [ -f $c/OUTCAR ] || { echo "$c: no OUTCAR"; continue; }
        e=$(grep "energy(sigma->0)" $c/OUTCAR | tail -1 | awk '{print $NF}')
        m=$(grep "number of electron" $c/OUTCAR | tail -1 | awk '{if ($NF ~ /^-?[0-9.]+$/) print $NF; else print 0}')
        printf "%-5s %16s %14s\n" "$c" "$e" "${m:-0}"
    done
    exit 0
fi

[ -s POTCAR ] || { echo "Build POTCAR first, e.g. cat \$VASP_PP_PATH/potpaw_PBE/Fe_pv/POTCAR > POTCAR"; exit 1; }

for c in FM AFM NM; do
    mkdir -p $c
    cp POSCAR KPOINTS POTCAR job.slurm $c/
    case $c in
        FM)  sed 's/^MAGMOM .*/MAGMOM  = 3.0 3.0/' INCAR > $c/INCAR ;;
        AFM) sed 's/^MAGMOM .*/MAGMOM  = 3.0 -3.0/' INCAR > $c/INCAR ;;
        NM)  sed -e 's/^ISPIN .*/ISPIN   = 1/' -e '/^MAGMOM /d' INCAR > $c/INCAR ;;
    esac
    (cd $c && sbatch job.slurm)
done
