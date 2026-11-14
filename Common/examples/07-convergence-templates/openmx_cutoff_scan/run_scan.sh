#!/usr/bin/env bash
# OpenMX convergence scan of scf.energycutoff or scf.Kgrid.
#
# Usage:
#   bash run_scan.sh ecut  [prepare|run|collect]
#   bash run_scan.sh kgrid [prepare|run|collect]
#
# Optional environment variables:
#   ECUTS="120 150 180 220 260 300 350 400"   KFIX="8 8 8"   (ecut scan, fixed k mesh)
#   KGRIDS="4 6 8 10 12 14 16"                EFIX=220       (kgrid scan, n x n x n, fixed cutoff)
#   MPI_RUN="mpirun -np 4"  OPENMX=openmx  NT=1  NATOMS=2
#
# collect writes ecut.dat / kgrid.dat:  value  Utot(hartree)  NATOMS
# Plot:  python ../plot_convergence.py ecut.dat --unit hartree --xlabel "scf.energycutoff (Ry)"
set -euo pipefail

MODE=${1:-}
ACTION=${2:-run}
ECUTS=${ECUTS:-"120 150 180 220 260 300 350 400"}
KGRIDS=${KGRIDS:-"4 6 8 10 12 14 16"}
KFIX=${KFIX:-"8 8 8"}
EFIX=${EFIX:-220}
MPI_RUN=${MPI_RUN:-"mpirun -np 4"}
OPENMX=${OPENMX:-openmx}
NT=${NT:-1}
NATOMS=${NATOMS:-2}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

case "$MODE" in
  ecut)  values=$ECUTS ;;
  kgrid) values=$KGRIDS ;;
  *) sed -n '2,14p' "$0"; exit 1 ;;
esac
mkdir -p "scan_$MODE"
cd "scan_$MODE"

if [ "$ACTION" = prepare ] || [ "$ACTION" = run ]; then
  for v in $values; do
    name=${MODE}_$v
    if [ "$MODE" = ecut ]; then e=$v; k=$KFIX; else e=$EFIX; k="$v $v $v"; fi
    sed -e "s/@ECUT@/$e/" -e "s/@KGRID@/$k/" -e "s/@NAME@/$name/" ../template.dat > "$name.dat"
    if [ "$ACTION" = run ]; then
      echo ">>> $name (scf.energycutoff=$e, scf.Kgrid=$k)"
      $MPI_RUN $OPENMX "$name.dat" -nt "$NT" > "$name.std" 2>&1
    fi
  done
  if [ "$ACTION" = prepare ]; then echo "inputs prepared in scan_$MODE/"; exit 0; fi
fi

: > "../$MODE.dat"
for v in $values; do
  out=${MODE}_$v.out
  [ -f "$out" ] || { echo "WARNING: $out missing" >&2; continue; }
  u=$(awk '/Utot\./{x=$2} END{print x}' "$out")
  [ -n "$u" ] || { echo "WARNING: no Utot. in $out" >&2; continue; }
  echo "$v $u $NATOMS" >> "../$MODE.dat"
  printf "%6s  Utot = %s Hartree\n" "$v" "$u"
done
echo "written: $MODE.dat"
