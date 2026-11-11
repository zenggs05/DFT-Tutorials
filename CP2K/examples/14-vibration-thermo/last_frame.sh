#!/bin/bash
# last_frame.sh - copy the last frame of a CP2K xyz trajectory to a new xyz file
# Usage: bash last_frame.sh <traj.xyz> <out.xyz>
#   e.g. bash last_frame.sh opt_initial-pos-1.xyz initial_opt.xyz
set -euo pipefail
traj=${1:?usage: last_frame.sh traj.xyz out.xyz}
out=${2:?usage: last_frame.sh traj.xyz out.xyz}
nat=$(head -n 1 "$traj" | tr -d '[:space:]')
nlines=$(wc -l < "$traj")
if (( nlines % (nat + 2) != 0 )); then
  echo "warning: $traj has $nlines lines, not a multiple of $((nat + 2))" >&2
fi
tail -n $((nat + 2)) "$traj" > "$out"
echo "wrote $out ($nat atoms): $(sed -n 2p "$out")"
