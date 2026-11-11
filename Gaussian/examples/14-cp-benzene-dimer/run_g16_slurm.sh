#!/bin/bash
#SBATCH --job-name=g16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out
#
# Slurm script for one Gaussian 16 input (may contain several --Link1-- steps).
# Usage:   sbatch run_g16_slurm.sh input.gjf
# The inputs in this directory use %nprocshared=16 and %mem=16GB, i.e. the
# same number of cores and about 80% of the memory requested above.
# A more complete wrapper is Gaussian/scripts/sub_g16.sh (chapter 17).

set -eo pipefail

INP=${1:?usage: sbatch run_g16_slurm.sh input.gjf}
BASE=${INP%.*}

# --- Gaussian environment: edit for your cluster ----------------------------
module load gaussian/16 2>/dev/null || true          # option A: module
# export g16root=/opt/software                        # option B: manual
# source "$g16root/g16/bsd/g16.profile"
command -v g16 >/dev/null || { echo "g16 not found in PATH" >&2; exit 1; }

# --- private scratch directory, removed when the job ends -------------------
export GAUSS_SCRDIR="${SCR_BASE:-/tmp/$USER}/g16_${SLURM_JOB_ID:-$$}"
mkdir -p "$GAUSS_SCRDIR"
trap 'rm -rf "$GAUSS_SCRDIR"' EXIT

# --- warn if %nprocshared differs from the allocation -----------------------
NP=$(grep -i -m1 '^%nprocshared' "$INP" | cut -d= -f2 | tr -d '[:space:]' || true)
if [ -n "$NP" ] && [ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ "$NP" != "$SLURM_CPUS_PER_TASK" ]; then
    echo "WARNING: %nprocshared=$NP but --cpus-per-task=$SLURM_CPUS_PER_TASK" >&2
fi

echo "start $BASE on $(hostname) at $(date)"
g16 < "$INP" > "$BASE.log" || true

# --- every Link1 step prints its own termination line -----------------------
NSTEP=$(grep -c -- '--Link1--' "$INP" || true)
NSTEP=$((NSTEP + 1))
NOK=$(grep -c "Normal termination" "$BASE.log" || true)
echo "$NOK of $NSTEP job step(s) terminated normally"
if [ "$NOK" -ne "$NSTEP" ]; then
    echo "Problem in $BASE.log; last lines:" >&2
    tail -n 8 "$BASE.log" >&2
    exit 1
fi

# --- formatted checkpoint files for GaussView / Multiwfn --------------------
for chk in $(grep -i '^%chk' "$INP" | cut -d= -f2 | tr -d '\r' | sort -u); do
    [ -f "$chk" ] && formchk "$chk" "${chk%.chk}.fchk" >/dev/null 2>&1 || true
done
echo "done at $(date)"
