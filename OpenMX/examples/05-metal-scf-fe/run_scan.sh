#!/bin/bash
# Electronic-temperature and k-point scans for ferromagnetic bcc Fe (OpenMX 4.0).
#
# Usage:
#   bash run_scan.sh temp  [run|collect]   # scan scf.ElectronicTemperature (k-grid as in the input)
#   bash run_scan.sh kgrid [run|collect]   # scan scf.Kgrid n n n (temperature as in the input)
#
# Environment variables (optional):
#   MPI_RUN (default "mpirun -np 4")  OPENMX (default openmx)  NT (default 1)
#   TEMPS   (default "300 600 1000")             electronic temperatures in K
#   KGRIDS  (default "10 12 14 16 18 20 24")     n for an n x n x n grid
#
# Inputs are generated from Fe_bcc_metal_scf.dat with GNU sed.
set -euo pipefail

MPI_RUN=${MPI_RUN:-"mpirun -np 4"}
OPENMX=${OPENMX:-openmx}
NT=${NT:-1}
TEMPS=${TEMPS:-"300 600 1000"}
KGRIDS=${KGRIDS:-"10 12 14 16 18 20 24"}
BASE=Fe_bcc_metal_scf.dat
HA2MEV=27211.386

usage() { echo "usage: bash $0 temp|kgrid [run|collect]" >&2; exit 1; }

mode=${1:-}
action=${2:-run}
case "$mode" in
    temp)  values=$TEMPS;  prefix=fe_T ;;
    kgrid) values=$KGRIDS; prefix=fe_k ;;
    *) usage ;;
esac
case "$action" in run|collect) ;; *) usage ;; esac

if [ "$action" = run ]; then
    for v in $values; do
        inp=${prefix}${v}.dat
        if [ "$mode" = temp ]; then
            sed -e "s/^System.Name .*/System.Name                 ${prefix}${v}/" \
                -e "s/^scf.ElectronicTemperature .*/scf.ElectronicTemperature   ${v}.0/" \
                "$BASE" > "$inp"
        else
            sed -e "s/^System.Name .*/System.Name                 ${prefix}${v}/" \
                -e "s/^scf.Kgrid .*/scf.Kgrid                   ${v} ${v} ${v}/" \
                "$BASE" > "$inp"
        fi
        echo ">>> $MPI_RUN $OPENMX $inp -nt $NT"
        $MPI_RUN $OPENMX "$inp" -nt "$NT" > "${prefix}${v}.std" 2>&1
    done
fi

printf "# %-6s %18s %20s %16s %10s %9s\n" "$mode" "Utot(Ha/atom)" "dE_vs_last(meV)" "spin(muB/cell)" "time(s)" "scf_conv"
last_e=""
for v in $values; do last_e_file=${prefix}${v}.out; done
[ -f "$last_e_file" ] && last_e=$(awk '/Utot\./{x=$2} END{print x}' "$last_e_file")
for v in $values; do
    out=${prefix}${v}.out
    if [ ! -f "$out" ]; then
        echo "WARNING: $out not found, skipped" >&2
        continue
    fi
    e=$(awk '/Utot\./{x=$2} END{print x}' "$out")
    m=$(awk '/Total spin moment/{x=$NF} END{print x}' "$out")
    t=$(awk '/Elapsed.Time\./{x=$2} END{print x}' "$out")
    c=$(grep -o 'scf_conv=[01]' "${prefix}${v}.md" 2>/dev/null | tail -1 | cut -d= -f2 || true)
    awk -v v="$v" -v e="$e" -v r="${last_e:-$e}" -v m="${m:-NA}" -v t="${t:-NA}" -v c="${c:-NA}" -v k="$HA2MEV" \
        'BEGIN{printf "  %-6s %18.10f %20.3f %16s %10s %9s\n", v, e, (e-r)*k, m, t, c}'
done
echo "# Compare energies only between runs with the SAME temperature and k-grid settings,"
echo "# except for the scanned parameter itself. Use only rows with scf_conv=1."
