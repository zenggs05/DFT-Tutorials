#!/usr/bin/env bash
# Two-step CP2K grid convergence scan: CUTOFF first, then REL_CUTOFF.
#
# Usage:
#   bash run_scan.sh cutoff  [prepare|run|collect]   # scan CUTOFF at fixed REL_CUTOFF
#   bash run_scan.sh rel     [prepare|run|collect]   # scan REL_CUTOFF at fixed CUTOFF
#
# Optional environment variables:
#   CUTOFFS="200 250 300 350 400 500 600 800"   REL=60       (step 1)
#   RELS="30 40 50 60 70 80 100"                CUT=400      (step 2, use the step-1 result)
#   CP2K="mpirun -np 8 cp2k.psmp"   NATOMS=8
#
# collect writes cutoff.dat / rel.dat with columns:  value  E(hartree)  NATOMS
# and prints the number of Gaussians mapped on each grid ("count for grid") and the
# "Electronic density on regular grids" line (integral and deviation from the exact count).
# Plot:  python ../plot_convergence.py cutoff.dat --unit hartree --xlabel "CUTOFF (Ry)"
set -euo pipefail

MODE=${1:-}
ACTION=${2:-run}
CUTOFFS=${CUTOFFS:-"200 250 300 350 400 500 600 800"}
RELS=${RELS:-"30 40 50 60 70 80 100"}
REL=${REL:-60}
CUT=${CUT:-400}
CP2K=${CP2K:-"mpirun -np 8 cp2k.psmp"}
NATOMS=${NATOMS:-8}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

case "$MODE" in
  cutoff) values=$CUTOFFS; dat=cutoff.dat ;;
  rel)    values=$RELS;    dat=rel.dat ;;
  *) sed -n '2,15p' "$0"; exit 1 ;;
esac

dirname_of() { if [ "$MODE" = cutoff ]; then echo "cut_${1}_rel_${REL}"; else echo "cut_${CUT}_rel_${1}"; fi; }

if [ "$ACTION" = prepare ] || [ "$ACTION" = run ]; then
  for v in $values; do
    d=$(dirname_of "$v")
    mkdir -p "$d"
    if [ "$MODE" = cutoff ]; then c=$v; r=$REL; else c=$CUT; r=$v; fi
    sed -E -e "s/^@SET[[:space:]]+CUTOFF[[:space:]].*/@SET CUTOFF ${c}/" \
           -e "s/^@SET[[:space:]]+REL_CUTOFF[[:space:]].*/@SET REL_CUTOFF ${r}/" \
           template.inp > "$d/input.inp"
    if [ "$ACTION" = run ]; then
      echo ">>> CUTOFF=$c Ry REL_CUTOFF=$r Ry"
      ( cd "$d" && $CP2K -i input.inp -o output.out )
    fi
  done
  if [ "$ACTION" = prepare ]; then echo "inputs prepared"; exit 0; fi
fi

: > "$dat"
for v in $values; do
  o=$(dirname_of "$v")/output.out
  [ -f "$o" ] || { echo "WARNING: $o missing" >&2; continue; }
  grep -q "SCF run NOT converged" "$o" && echo "WARNING: SCF not converged in $o" >&2
  e=$(grep "ENERGY| Total FORCE_EVAL" "$o" | tail -1 | awk '{print $NF}')
  counts=$(grep "count for grid" "$o" | awk '{printf "%s ", $5}')
  # integrated electron density on the grids and its deviation from the exact count
  rho=$(grep "Electronic density on regular grids" "$o" | tail -1 | awk '{print $(NF-1), $NF}')
  [ -n "$e" ] || { echo "WARNING: no energy in $o" >&2; continue; }
  echo "$v $e $NATOMS" >> "$dat"
  printf "%6s  %20s  grids: %s  rho(grid), error: %s\n" "$v" "$e" "$counts" "${rho:-n/a}"
done
echo "written: $dat (energies in hartree)"
