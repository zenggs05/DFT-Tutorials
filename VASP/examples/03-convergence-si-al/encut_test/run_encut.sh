#!/bin/bash
# run_encut.sh - ENCUT convergence scan for bulk Si (one directory per cutoff)
#
# Usage (inside encut_test/, POTCAR must exist here):
#   bash run_encut.sh                       # ideal structure (POSCAR)
#   bash run_encut.sh POSCAR_displaced      # displaced atom -> force convergence
#
# Environment variables:
#   ENCUTS   list of cutoffs in eV   (default "250 275 300 325 350 375 400 450 500")
#   MODE     local | slurm | collect (default local)
#   NP       MPI ranks for MODE=local (default 4)
#   VASP_CMD full command for MODE=local (default "mpirun -np $NP vasp_std")
#
# Result: encut_<POSCAR name>.dat with the columns
#   ENCUT natoms E0 E0_per_atom TOTEN EENTRO pressure_kB maxF_eVA nkpts time_s
# Analyse with: python ../plot_conv.py encut_POSCAR.dat

HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=../common.sh
source "$HERE/../common.sh"
cd "$HERE" || exit 1

POSCAR_FILE=${1:-POSCAR}
ENCUTS=${ENCUTS:-"250 275 300 325 350 375 400 450 500"}
TAG=$(basename "$POSCAR_FILE")
OUT=encut_${TAG}.dat

[ -s "$POSCAR_FILE" ] || { echo "ERROR: $POSCAR_FILE not found"; exit 1; }
check_potcar || exit 1

if [ "$MODE" != "collect" ]; then
    for e in $ENCUTS; do
        d=${TAG}_ENCUT_$e
        mkdir -p "$d"
        sed "s/__ENCUT__/$e/" INCAR.template > "$d/INCAR"
        cp "$POSCAR_FILE" "$d/POSCAR"
        cp KPOINTS POTCAR "$d/"
        run_or_submit "$d" || exit 1
    done
fi

if [ "$MODE" = "slurm" ]; then
    echo "Jobs submitted. When they are done run:  MODE=collect bash run_encut.sh $POSCAR_FILE"
    exit 0
fi

n=$(natoms "$POSCAR_FILE")
echo "# ENCUT natoms E0 E0_per_atom TOTEN EENTRO pressure_kB maxF_eVA nkpts time_s" > "$OUT"
for e in $ENCUTS; do
    d=${TAG}_ENCUT_$e
    read -r e0 f ts p fm nk t <<< "$(outcar_row "$d")"
    if ! grep -q "General timing and accounting" "$d/OUTCAR" 2>/dev/null; then
        echo "WARNING: $d did not finish normally" >&2
    fi
    epa=$(awk -v e="$e0" -v n="$n" 'BEGIN{ if (e == "NaN") print "NaN"; else printf "%.8f", e / n }')
    echo "$e $n $e0 $epa $f $ts $p $fm $nk $t" >> "$OUT"
done
column -t "$OUT" 2>/dev/null || cat "$OUT"
echo "Written $OUT. Next: python ../plot_conv.py $OUT"
