#!/usr/bin/env bash
# Step 2 of the two-step grid convergence test:
# scan REL_CUTOFF (Ry) at the CUTOFF chosen in step 1.
#
# Usage (all settings are optional environment variables):
#   CUTOFF=350 bash run_relcutoff_scan.sh
#   CUTOFF=350 RELCUTOFFS="30 40 50 60 70 80" bash run_relcutoff_scan.sh
#   CP2K="mpirun -np 4 cp2k.psmp" NPAR=2 CUTOFF=350 bash run_relcutoff_scan.sh
#   DRYRUN=1 CUTOFF=350 bash run_relcutoff_scan.sh
set -euo pipefail

TEMPLATE=${TEMPLATE:-template.inp}
EXTRA_FILES=${EXTRA_FILES:-si_bulk8.xyz}
CUTOFF=${CUTOFF:-350}                      # result of step 1
RELCUTOFFS=${RELCUTOFFS:-"20 30 40 50 60 70 80 90 100"}
CP2K=${CP2K:-cp2k.psmp}
NPAR=${NPAR:-1}
DRYRUN=${DRYRUN:-0}

[ -f "$TEMPLATE" ] || { echo "ERROR: template $TEMPLATE not found" >&2; exit 1; }
grep -Eq '^@SET[[:space:]]+REL_CUTOFF[[:space:]]' "$TEMPLATE" || {
  echo "ERROR: $TEMPLATE has no '@SET REL_CUTOFF' line" >&2; exit 1; }

NATOMS=${NATOMS:-}
if [ -z "$NATOMS" ]; then
  for f in $EXTRA_FILES; do
    case $f in *.xyz) NATOMS=$(head -n 1 "$f" | tr -d '[:space:]'); break ;; esac
  done
fi
NATOMS=${NATOMS:-1}

running=0
for r in $RELCUTOFFS; do
  d=relcutoff_${r}Ry
  mkdir -p "$d"
  sed -E -e "s/^@SET[[:space:]]+CUTOFF[[:space:]].*/@SET CUTOFF ${CUTOFF}/" \
         -e "s/^@SET[[:space:]]+REL_CUTOFF[[:space:]].*/@SET REL_CUTOFF ${r}/" \
         "$TEMPLATE" > "$d/input.inp"
  for f in $EXTRA_FILES; do cp "$f" "$d/"; done
  [ "$DRYRUN" = "1" ] && continue
  echo "running REL_CUTOFF = ${r} Ry (CUTOFF = ${CUTOFF} Ry) in $d"
  ( cd "$d" && $CP2K -i input.inp -o output.out ) &
  running=$((running + 1))
  if [ "$running" -ge "$NPAR" ]; then wait; running=0; fi
done
wait

if [ "$DRYRUN" = "1" ]; then
  echo "DRYRUN: inputs written to relcutoff_*Ry/, nothing was run."
  exit 0
fi
bash collect_scan.sh relcutoff "$NATOMS"
