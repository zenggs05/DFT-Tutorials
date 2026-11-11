#!/bin/bash
# run_kpts.sh - k-point convergence scan for bulk Si
#
# Usage (inside kpoint_test/, POTCAR must exist here):
#   bash run_kpts.sh                        # Gamma-centred N x N x N meshes
#   SCAN=kspacing bash run_kpts.sh          # no KPOINTS file, KSPACING in INCAR
#
# Environment variables:
#   SCAN      mesh | kspacing                     (default mesh)
#   MESHES    list of N for N x N x N meshes      (default "2 4 6 8 10 12 14 16")
#   SPACINGS  list of KSPACING values in 1/A      (default "0.5 0.4 0.3 0.25 0.2 0.15")
#   MODE      local | slurm | collect             (default local)
#   NP, VASP_CMD as in ../common.sh
#
# Result: kpts_mesh.dat or kpts_kspacing.dat with the columns
#   X natoms E0 E0_per_atom TOTEN EENTRO pressure_kB maxF_eVA nkpts time_s
# where X = N (mesh) or KSPACING. Analyse with
#   python ../plot_conv.py kpts_mesh.dat
#   python ../plot_conv.py kpts_kspacing.dat --ref min     (smallest spacing = reference)

HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=../common.sh
source "$HERE/../common.sh"
cd "$HERE" || exit 1

SCAN=${SCAN:-mesh}
MESHES=${MESHES:-"2 4 6 8 10 12 14 16"}
SPACINGS=${SPACINGS:-"0.5 0.4 0.3 0.25 0.2 0.15"}
case "$SCAN" in
    mesh)     LIST=$MESHES ;;
    kspacing) LIST=$SPACINGS ;;
    *) echo "ERROR: SCAN must be mesh or kspacing"; exit 1 ;;
esac
OUT=kpts_${SCAN}.dat

check_potcar || exit 1

if [ "$MODE" != "collect" ]; then
    for x in $LIST; do
        d=${SCAN}_$x
        mkdir -p "$d"
        cp INCAR POSCAR POTCAR "$d/"
        if [ "$SCAN" = "mesh" ]; then
            printf 'Gamma-centred %sx%sx%s\n0\nGamma\n  %s %s %s\n  0 0 0\n' \
                "$x" "$x" "$x" "$x" "$x" "$x" > "$d/KPOINTS"
        else
            rm -f "$d/KPOINTS"            # KSPACING is ignored if KPOINTS exists
            printf 'KSPACING = %s\nKGAMMA = .TRUE.\n' "$x" >> "$d/INCAR"
        fi
        run_or_submit "$d" || exit 1
    done
fi

if [ "$MODE" = "slurm" ]; then
    echo "Jobs submitted. When they are done run:  MODE=collect SCAN=$SCAN bash run_kpts.sh"
    exit 0
fi

n=$(natoms POSCAR)
echo "# X natoms E0 E0_per_atom TOTEN EENTRO pressure_kB maxF_eVA nkpts time_s" > "$OUT"
for x in $LIST; do
    d=${SCAN}_$x
    read -r e0 f ts p fm nk t <<< "$(outcar_row "$d")"
    if ! grep -q "General timing and accounting" "$d/OUTCAR" 2>/dev/null; then
        echo "WARNING: $d did not finish normally" >&2
    fi
    epa=$(awk -v e="$e0" -v n="$n" 'BEGIN{ if (e == "NaN") print "NaN"; else printf "%.8f", e / n }')
    echo "$x $n $e0 $epa $f $ts $p $fm $nk $t" >> "$OUT"
done
column -t "$OUT" 2>/dev/null || cat "$OUT"
if [ "$SCAN" = "kspacing" ]; then
    echo "Mesh information written by VASP (format may differ between versions):"
    grep -H -i "generate k-points\|k-mesh" kspacing_*/OUTCAR 2>/dev/null | head -20
    echo "Written $OUT. Next: python ../plot_conv.py $OUT --ref min"
else
    echo "Written $OUT. Next: python ../plot_conv.py $OUT"
fi
