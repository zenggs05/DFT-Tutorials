#!/bin/bash
# run_bluemoon.sh - run the blue-moon windows of chapter 15.
#
# CHAIN=1 (default): windows run one after another; window k starts from the
#   last frame of window k-1 (the first from $START).
# CHAIN=0: every window starts from its own window_<t>/start.xyz, prepared e.g.
#   by pick_windows.py from an xTB metadynamics trajectory.  The windows are then
#   independent: run all of them here, or one per job with TARGETS="3.20".
# CONSTRAINT_INIT T moves the ions onto the target before the first step.
#
# Usage:
#   bash run_bluemoon.sh                        # DFT level, default targets, chained
#   LEVEL=xtb STEPS=20000 bash run_bluemoon.sh  # cheap GFN1-xTB windows
#   CHAIN=0 TARGETS="3.20 3.40" bash run_bluemoon.sh
#   DRYRUN=1 bash run_bluemoon.sh               # only create the directories
# Environment:
#   MPIRUN   (default "mpirun -np 64")   CP2K_EXE (default cp2k.psmp)
#   TARGETS  (default "2.60 2.80 ... 5.20", Angstrom; same strings as the
#             window_<t> directory names)
#   START    (default nacl_eq.xyz, only used with CHAIN=1)
set -euo pipefail

LEVEL=${LEVEL:-dft}
STEPS=${STEPS:-10000}
CHAIN=${CHAIN:-1}
MPIRUN=${MPIRUN:-mpirun -np 64}
CP2K_EXE=${CP2K_EXE:-cp2k.psmp}
START=${START:-nacl_eq.xyz}
TARGETS=${TARGETS:-"2.60 2.80 3.00 3.20 3.40 3.60 3.80 4.00 4.20 4.40 4.60 4.80 5.00 5.20"}
here=$(pwd)

if [ "$CHAIN" = "1" ]; then
  [ -f "$START" ] || { echo "start structure $START not found (run nacl_equil.inp first)" >&2; exit 1; }
  prev="$here/$START"
fi

for t in $TARGETS; do
  dir="window_${t}"
  mkdir -p "$dir"
  cp "$here/constrained_md_distance.inp" "$dir/"
  ln -sf "$here/nacl_${LEVEL}_force_eval.inc" "$dir/"
  if [ "$CHAIN" = "1" ]; then
    cp "$prev" "$dir/start.xyz"
    echo "=== $dir  (TARGET = $t A, start from $(basename "$(dirname "$prev")")/$(basename "$prev"))"
  else
    [ -f "$dir/start.xyz" ] || { echo "$dir/start.xyz missing (run pick_windows.py first)" >&2; exit 1; }
    echo "=== $dir  (TARGET = $t A, start from its own start.xyz)"
  fi
  if [ "${DRYRUN:-0}" = "1" ]; then
    prev="$here/$dir/start.xyz"
    continue
  fi
  # CP2K appends to existing print-key files: remove old multipliers so that a
  # re-run window does not mix two trajectories
  rm -f "$dir/bluemoon.LagrangeMultLog"
  ( cd "$dir" && $MPIRUN "$CP2K_EXE" -E LEVEL="$LEVEL" -E TARGET="$t" -E STEPS="$STEPS" \
        -E COORD=start.xyz -i constrained_md_distance.inp -o bm.out )
  grep -q "PROGRAM ENDED AT" "$dir/bm.out" || { echo "window $t did not finish, see $dir/bm.out" >&2; exit 1; }
  traj="$dir/nacl_bm_${LEVEL}-pos-1.xyz"
  nat=$(head -n 1 "$traj" | tr -d '[:space:]')
  tail -n $((nat + 2)) "$traj" > "$dir/last.xyz"
  prev="$here/$dir/last.xyz"
done

echo "done. Integrate with: python3 bluemoon_integrate.py window_*/ --pmf"
