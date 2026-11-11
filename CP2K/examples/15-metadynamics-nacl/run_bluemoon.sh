#!/bin/bash
# run_bluemoon.sh - run the blue-moon windows of chapter 15 one after another.
# Window k starts from the last frame of window k-1 (the first from nacl_eq.xyz);
# CONSTRAINT_INIT T moves the ions onto the new target before the first step.
#
# Usage:
#   bash run_bluemoon.sh                       # DFT level, default targets
#   LEVEL=xtb STEPS=20000 bash run_bluemoon.sh # cheap GFN1-xTB windows
#   DRYRUN=1 bash run_bluemoon.sh              # only create the directories
# Environment:
#   MPIRUN   (default "mpirun -np 64")   CP2K_EXE (default cp2k.psmp)
#   TARGETS  (default "2.60 2.80 ... 5.20", Angstrom)
#   START    (default nacl_eq.xyz)
set -euo pipefail

LEVEL=${LEVEL:-dft}
STEPS=${STEPS:-10000}
MPIRUN=${MPIRUN:-mpirun -np 64}
CP2K_EXE=${CP2K_EXE:-cp2k.psmp}
START=${START:-nacl_eq.xyz}
TARGETS=${TARGETS:-"2.60 2.80 3.00 3.20 3.40 3.60 3.80 4.00 4.20 4.40 4.60 4.80 5.00 5.20"}
here=$(pwd)

[ -f "$START" ] || { echo "start structure $START not found (run nacl_equil.inp first)" >&2; exit 1; }
prev="$here/$START"

for t in $TARGETS; do
  dir="window_${t}"
  mkdir -p "$dir"
  cp "$here/constrained_md_distance.inp" "$dir/"
  ln -sf "$here/nacl_${LEVEL}_force_eval.inc" "$dir/"
  cp "$prev" "$dir/start.xyz"
  echo "=== $dir  (TARGET = $t A, start from $(basename "$(dirname "$prev")")/$(basename "$prev"))"
  if [ "${DRYRUN:-0}" = "1" ]; then
    prev="$here/$dir/start.xyz"
    continue
  fi
  ( cd "$dir" && $MPIRUN "$CP2K_EXE" -E LEVEL="$LEVEL" -E TARGET="$t" -E STEPS="$STEPS" \
        -E COORD=start.xyz -i constrained_md_distance.inp -o bm.out )
  grep -q "PROGRAM ENDED AT" "$dir/bm.out" || { echo "window $t did not finish, see $dir/bm.out" >&2; exit 1; }
  traj="$dir/nacl_bm_${LEVEL}-pos-1.xyz"
  nat=$(head -n 1 "$traj" | tr -d '[:space:]')
  tail -n $((nat + 2)) "$traj" > "$dir/last.xyz"
  prev="$here/$dir/last.xyz"
done

echo "done. Integrate with: python3 bluemoon_integrate.py window_*/ --pmf"
