#!/bin/bash
# batch_submit_all_gjf.sh - submit every *.gjf in a folder with g16_slurm.sh, skipping
# inputs that already finished (log ends with "Normal termination") or are already queued.
#
# Usage:  bash batch_submit_all_gjf.sh [folder] [max_jobs_in_queue]
#         bash batch_submit_all_gjf.sh conformers/ 40
# Extra sbatch options can be passed through the SBATCH_OPTS variable, e.g.
#         SBATCH_OPTS="--cpus-per-task=8 --mem=16G" bash batch_submit_all_gjf.sh .
set -uo pipefail

DIR=${1:-.}
MAXQ=${2:-50}
TEMPLATE=$(cd "$(dirname "$0")" && pwd)/g16_slurm.sh
[ -f "$TEMPLATE" ] || { echo "ERROR: $TEMPLATE not found"; exit 1; }
cd "$DIR" || exit 1

shopt -s nullglob
inputs=(*.gjf)
[ ${#inputs[@]} -gt 0 ] || { echo "no .gjf files in $PWD"; exit 0; }

nsub=0; nskip=0
for gjf in "${inputs[@]}"; do
    name=${gjf%.gjf}
    if [ -f "$name.log" ] && tail -n 1 "$name.log" | grep -q "Normal termination"; then
        echo "skip  $gjf (finished)"; nskip=$((nskip + 1)); continue
    fi
    if squeue -h -u "$USER" -o "%j" | grep -qx "$name"; then
        echo "skip  $gjf (already queued)"; nskip=$((nskip + 1)); continue
    fi
    # throttle: do not flood the queue (many sites limit jobs per user)
    while [ "$(squeue -h -u "$USER" | wc -l)" -ge "$MAXQ" ]; do sleep 60; done
    # shellcheck disable=SC2086
    jid=$(sbatch --parsable --job-name="$name" ${SBATCH_OPTS:-} "$TEMPLATE" "$gjf") || { echo "sbatch failed for $gjf"; continue; }
    echo "submit $gjf -> job $jid"
    nsub=$((nsub + 1))
done
echo "submitted $nsub, skipped $nskip"
