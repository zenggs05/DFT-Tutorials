#!/bin/bash
# Cutoff-energy and k-point convergence tests for diamond Si (OpenMX 4.0).
#
# Usage:
#   bash run_conv.sh ecut  [run|collect]   # scf.energycutoff scan in conv_ecut/
#   bash run_conv.sh kgrid [run|collect]   # scf.Kgrid scan in conv_kgrid/
#
# Environment variables (optional):
#   MPI_RUN  (default "mpirun -np 4")   OPENMX (default openmx)   NT (default 1)
#   ECUTS    (default "120 150 180 220 260 300 350 400")   cutoff values in Ry
#   KGRIDS   (default "4 6 8 10 12 14 16")                 n for an n x n x n grid
#   THRESH   (default 1.0)                                 tolerance in meV/atom
#
# The summary compares every run with the LAST value of the list (the most accurate one)
# and is also written to conv_<mode>/<mode>_conv.txt.
set -euo pipefail

MPI_RUN=${MPI_RUN:-"mpirun -np 4"}
OPENMX=${OPENMX:-openmx}
NT=${NT:-1}
ECUTS=${ECUTS:-"120 150 180 220 260 300 350 400"}
KGRIDS=${KGRIDS:-"4 6 8 10 12 14 16"}
THRESH=${THRESH:-1.0}
NATOM=2
HERE=$(cd "$(dirname "$0")" && pwd)

usage() { echo "usage: bash $0 ecut|kgrid [run|collect]" >&2; exit 1; }

mode=${1:-}
action=${2:-run}
case "$mode" in
    ecut)  dir=$HERE/conv_ecut;  tmpl=Si_ecut_template.dat; tag=@ECUT@; values=$ECUTS;  prefix=si_ecut_ ;;
    kgrid) dir=$HERE/conv_kgrid; tmpl=Si_k_template.dat;    tag=@K@;    values=$KGRIDS; prefix=si_k ;;
    *) usage ;;
esac
case "$action" in run|collect) ;; *) usage ;; esac
cd "$dir"

if [ "$action" = run ]; then
    for v in $values; do
        inp=${prefix}${v}.dat
        sed "s/${tag}/${v}/g" "$tmpl" > "$inp"
        echo ">>> $MPI_RUN $OPENMX $inp -nt $NT"
        $MPI_RUN $OPENMX "$inp" -nt "$NT" > "${prefix}${v}.std" 2>&1
    done
fi

raw=$(mktemp)
trap 'rm -f "$raw"' EXIT
for v in $values; do
    out=${prefix}${v}.out
    if [ ! -f "$out" ]; then
        echo "WARNING: $out not found, skipped" >&2
        continue
    fi
    e=$(awk '/Utot\./{x=$2} END{print x}' "$out")
    t=$(awk '/Elapsed.Time\./{x=$2} END{print x}' "$out")
    g=$(awk -F'=' '/Num. of grids of a-, b-, and c-axes/{x=$2} END{gsub(/ /,"",x); print x}' "$out")
    c=$(grep -o 'scf_conv=[01]' "${prefix}${v}.md" 2>/dev/null | tail -1 | cut -d= -f2 || true)
    [ -n "$e" ] || { echo "WARNING: no Utot. in $out, skipped" >&2; continue; }
    echo "$v $e ${t:-NA} ${g:-NA} ${c:-NA}" >> "$raw"
done

awk -v n="$NATOM" -v th="$THRESH" -v mode="$mode" '
    { v[NR]=$1; e[NR]=$2; t[NR]=$3; g[NR]=$4; c[NR]=$5 }
    END {
        if (NR == 0) { print "no finished runs found" > "/dev/stderr"; exit 1 }
        ref = e[NR]
        printf "# %-6s %18s %22s %10s %14s %9s %8s\n", mode, "Utot(Ha/cell)", "dE_vs_last(meV/atom)", "time(s)", "grid", "scf_conv", "ok(<" th ")"
        for (i = 1; i <= NR; i++) {
            d = (e[i] - ref) * 27211.386 / n
            ok = (d < th && d > -th) ? "yes" : "no"
            printf "  %-6s %18.10f %22.3f %10s %14s %9s %8s\n", v[i], e[i], d, t[i], g[i], c[i], ok
        }
    }' "$raw" | tee "${mode}_conv.txt"
