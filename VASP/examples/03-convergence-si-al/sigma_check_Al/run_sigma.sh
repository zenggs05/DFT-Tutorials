#!/bin/bash
# run_sigma.sh - smearing (ISMEAR/SIGMA) and k-mesh test for fcc Al (metal)
#
# Usage (inside sigma_check_Al/, POTCAR for Al must exist here):
#   bash run_sigma.sh
#
# Environment variables:
#   ISMEARS  smearing methods to test        (default "1 0")
#   SIGMAS   smearing widths in eV           (default "0.05 0.1 0.2 0.3")
#   MESHES   Gamma-centred N x N x N meshes  (default "12 16 20")
#   REF_N    mesh of the tetrahedron reference run, ISMEAR = -5 (default 24)
#   MODE     local | slurm | collect         (default local)
#   NP, VASP_CMD as in ../common.sh
#
# Result: sigma_Al.dat with the columns
#   ISMEAR SIGMA N natoms E0 TOTEN EENTRO TS_meV_per_atom dE0_meV_vs_ref nkpts
# A setting is acceptable when |T*S| < 1 meV/atom and E0 is close to the
# tetrahedron (ISMEAR = -5) reference obtained with the densest mesh.

HERE=$(cd "$(dirname "$0")" && pwd)
# shellcheck source=../common.sh
source "$HERE/../common.sh"
cd "$HERE" || exit 1

ISMEARS=${ISMEARS:-"1 0"}
SIGMAS=${SIGMAS:-"0.05 0.1 0.2 0.3"}
MESHES=${MESHES:-"12 16 20"}
REF_N=${REF_N:-24}
OUT=sigma_Al.dat

check_potcar || exit 1

make_dir() {   # make_dir ISMEAR SIGMA N
    local d="ismear$1_sigma$2_k$3"
    mkdir -p "$d"
    sed -e "s/__ISMEAR__/$1/" -e "s/__SIGMA__/$2/" INCAR.template > "$d/INCAR"
    cp POSCAR POTCAR "$d/"
    printf 'Gamma-centred %sx%sx%s\n0\nGamma\n  %s %s %s\n  0 0 0\n' "$3" "$3" "$3" "$3" "$3" "$3" > "$d/KPOINTS"
    echo "$d"
}

DIRS=()
DIRS+=("$(make_dir -5 0.05 "$REF_N")")          # reference: tetrahedron + Bloechl, densest mesh
for is in $ISMEARS; do
    for s in $SIGMAS; do
        for k in $MESHES; do
            DIRS+=("$(make_dir "$is" "$s" "$k")")
        done
    done
done

if [ "$MODE" != "collect" ]; then
    for d in "${DIRS[@]}"; do run_or_submit "$d" || exit 1; done
fi
if [ "$MODE" = "slurm" ]; then
    echo "Jobs submitted. When they are done run:  MODE=collect bash run_sigma.sh"
    exit 0
fi

n=$(natoms POSCAR)
ref=$(outcar_row "${DIRS[0]}" | awk '{print $1}')
echo "# ISMEAR SIGMA N natoms E0 TOTEN EENTRO TS_meV_per_atom dE0_meV_vs_ref nkpts" > "$OUT"
for d in "${DIRS[@]}"; do
    is=$(echo "$d" | sed 's/ismear\(-*[0-9]*\)_.*/\1/')
    s=$(echo "$d" | sed 's/.*_sigma\([0-9.]*\)_k.*/\1/')
    k=${d##*_k}
    read -r e0 f ts p fm nk t <<< "$(outcar_row "$d")"
    awk -v is="$is" -v s="$s" -v k="$k" -v n="$n" -v e0="$e0" -v f="$f" -v ts="$ts" -v ref="$ref" -v nk="$nk" '
        BEGIN {
            tsm = (ts == "NaN") ? "NaN" : sprintf("%.3f", 1000 * ts / n)
            de  = (e0 == "NaN" || ref == "NaN") ? "NaN" : sprintf("%.3f", 1000 * (e0 - ref) / n)
            print is, s, k, n, e0, f, ts, tsm, de, nk
        }' >> "$OUT"
done
column -t "$OUT" 2>/dev/null || cat "$OUT"
echo
echo "Rule of thumb: keep |TS_meV_per_atom| < 1 and |dE0_meV_vs_ref| < 1."
awk '!/^#/ && $1 != -5 && $8 != "NaN" { t = $8 < 0 ? -$8 : $8; d = $9 < 0 ? -$9 : $9
     flag = (t < 1 && d < 1) ? "OK" : "--"
     printf "ISMEAR=%-3s SIGMA=%-5s N=%-3s  |TS|=%7.3f meV/atom  |dE0|=%7.3f meV/atom  %s\n", $1, $2, $3, t, d, flag }' "$OUT"
