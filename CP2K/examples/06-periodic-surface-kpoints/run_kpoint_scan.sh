#!/usr/bin/env bash
# k-point convergence scan: MONKHORST-PACK NK x NK x NK for NK in $NKS.
# Collects: number of k-points, total energy, dE per atom w.r.t. the densest mesh,
# and the pressure (1/3 trace of the stress tensor, bar).
#
# Usage (all settings are optional environment variables):
#   bash run_kpoint_scan.sh
#   NKS="2 4 6 8 10" CP2K="mpirun -np 8 cp2k.psmp" bash run_kpoint_scan.sh
#   DRYRUN=1 bash run_kpoint_scan.sh        # only create kp_*/input.inp
#   COLLECT_ONLY=1 bash run_kpoint_scan.sh  # only re-collect existing outputs
set -euo pipefail

TEMPLATE=${TEMPLATE:-kpoint_convergence_template.inp}
NKS=${NKS:-"2 3 4 5 6 8 10 12"}
NATOMS=${NATOMS:-2}                        # atoms in the cell of the template
CP2K=${CP2K:-cp2k.psmp}
DRYRUN=${DRYRUN:-0}
COLLECT_ONLY=${COLLECT_ONLY:-0}
OUT=kpoint_scan.csv

[ -f "$TEMPLATE" ] || { echo "ERROR: template $TEMPLATE not found" >&2; exit 1; }
grep -Eq '^@SET[[:space:]]+NK[[:space:]]' "$TEMPLATE" || {
  echo "ERROR: $TEMPLATE has no '@SET NK' line" >&2; exit 1; }

if [ "$COLLECT_ONLY" != "1" ]; then
  for n in $NKS; do
    d=kp_${n}
    mkdir -p "$d"
    sed -E "s/^@SET[[:space:]]+NK[[:space:]].*/@SET NK ${n}/" "$TEMPLATE" > "$d/input.inp"
    [ "$DRYRUN" = "1" ] && continue
    echo "running MONKHORST-PACK ${n}x${n}x${n} in $d"
    ( cd "$d" && $CP2K -i input.inp -o output.out )
  done
  if [ "$DRYRUN" = "1" ]; then echo "DRYRUN: inputs written to kp_*/"; exit 0; fi
fi

tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
for n in $NKS; do
  f=kp_${n}/output.out
  [ -f "$f" ] || { echo "skip kp_${n}: no output.out" >&2; continue; }
  e=$(grep "ENERGY| Total FORCE_EVAL" "$f" | tail -n 1 | awk '{print $NF}')
  [ -n "$e" ] || { echo "skip kp_${n}: no final energy" >&2; continue; }
  nkp=$(grep "BRILLOUIN| List of Kpoints" "$f" | tail -n 1 | awk '{print $NF}')
  p=$(grep "STRESS| 1/3 Trace" "$f" | tail -n 1 | awk '{print $NF}')
  echo "$n ${nkp:-n/a} $e ${p:-n/a}" >> "$tmp"
done
[ -s "$tmp" ] || { echo "ERROR: no usable output found" >&2; exit 1; }

awk -v n="$NATOMS" '
  BEGIN { ha2mev = 27211.386 }
  { k[NR] = $1; nk[NR] = $2; e[NR] = $3; p[NR] = $4 }
  END {
    ref = e[NR]
    print "NK,nkpts,energy_Ha,dE_meV_per_atom,pressure_bar"
    for (i = 1; i <= NR; i++)
      printf "%s,%s,%.10f,%.4f,%s\n", k[i], nk[i], e[i], (e[i] - ref) / n * ha2mev, p[i]
  }' "$tmp" > "$OUT"

if command -v column >/dev/null 2>&1; then column -s, -t "$OUT"; else cat "$OUT"; fi
echo "written: $OUT  (natoms = $NATOMS, reference = densest mesh = last row)"
