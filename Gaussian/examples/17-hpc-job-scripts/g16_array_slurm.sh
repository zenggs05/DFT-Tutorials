#!/bin/bash
# g16_array_slurm.sh - run many Gaussian inputs as ONE Slurm job array.
#
# Usage:  ls conf_*.gjf > joblist.txt
#         sbatch --array=1-$(wc -l < joblist.txt)%10 g16_array_slurm.sh
#           (%10 = at most 10 array tasks running at the same time)
# Each task takes line number $SLURM_ARRAY_TASK_ID of joblist.txt and hands it to
# g16_slurm.sh (which must be in the submit folder, or set G16_TEMPLATE).
#
#SBATCH --job-name=g16-array
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%A_%a.out

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
LIST=${LIST:-joblist.txt}
TEMPLATE=${G16_TEMPLATE:-$PWD/g16_slurm.sh}

INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID:?not an array job}p" "$LIST")
[ -n "$INPUT" ] || { echo "no line $SLURM_ARRAY_TASK_ID in $LIST"; exit 1; }
echo "array task $SLURM_ARRAY_TASK_ID -> $INPUT"
exec bash "$TEMPLATE" "$INPUT"
