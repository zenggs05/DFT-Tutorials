#!/usr/bin/env bash
# Step 1 of the two-step grid convergence test:
# scan CUTOFF (Ry) at a fixed REL_CUTOFF, then collect energies and grid counts.
#
# Usage (all settings are optional environment variables):
#   bash run_cutoff_scan.sh
#   CP2K="mpirun -np 4 cp2k.psmp" NPAR=2 bash run_cutoff_scan.sh
#   REL_CUTOFF=60 CUTOFFS="200 300 400 500" bash run_cutoff_scan.sh
#   TEMPLATE=h2o_template.inp EXTRA_FILES=h2o.xyz bash run_cutoff_scan.sh
#   DRYRUN=1 bash run_cutoff_scan.sh        # only create the directories/inputs
#
# The template must contain the lines "@SET CUTOFF <value>" and
# "@SET REL_CUTOFF <value>"; they are rewritten for every scan point.
set -euo pipefail

TEMPLATE=${TEMPLATE:-template.inp}
EXTRA_FILES=${EXTRA_FILES:-si_bulk8.xyz}   # files copied into every directory
REL_CUTOFF=${REL_CUTOFF:-60}
CUTOFFS=${CUTOFFS:-"150 200 250 300 350 400 450 500 600"}
CP2K=${CP2K:-cp2k.psmp}                    # may contain an MPI launcher
NPAR=${NPAR:-1}                            # number of calculations run at the same time
DRYRUN=${DRYRUN:-0}

[ -f "$TEMPLATE" ] || { echo "ERROR: template $TEMPLATE not found" >&2; exit 1; }
grep -Eq '^@SET[[:space:]]+CUTOFF[[:space:]]' "$TEMPLATE" || {
  echo "ERROR: $TEMPLATE has no '@SET CUTOFF' line" >&2; exit 1; }

# number of atoms (for energies per atom): first line of the first .xyz file
NATOMS=${NATOMS:-}
if [ -z "$NATOMS" ]; then
  for f in $EXTRA_FILES; do
    case $f in *.xyz) NATOMS=$(head -n 1 "$f" | tr -d '[:space:]'); break ;; esac
  done
fi
NATOMS=${NATOMS:-1}

running=0
for c in $CUTOFFS; do
  d=cutoff_${c}Ry
  mkdir -p "$d"
  sed -E -e "s/^@SET[[:space:]]+CUTOFF[[:space:]].*/@SET CUTOFF ${c}/" \
         -e "s/^@SET[[:space:]]+REL_CUTOFF[[:space:]].*/@SET REL_CUTOFF ${REL_CUTOFF}/" \
         "$TEMPLATE" > "$d/input.inp"
  for f in $EXTRA_FILES; do cp "$f" "$d/"; done
  [ "$DRYRUN" = "1" ] && continue
  echo "running CUTOFF = ${c} Ry (REL_CUTOFF = ${REL_CUTOFF} Ry) in $d"
  ( cd "$d" && $CP2K -i input.inp -o output.out ) &
  running=$((running + 1))
  if [ "$running" -ge "$NPAR" ]; then wait; running=0; fi
done
wait

if [ "$DRYRUN" = "1" ]; then
  echo "DRYRUN: inputs written to cutoff_*Ry/, nothing was run."
  exit 0
fi
bash collect_scan.sh cutoff "$NATOMS"
