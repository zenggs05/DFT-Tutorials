#!/bin/bash
#SBATCH --job-name=g16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#
# Minimal Slurm script for one Gaussian 16 job.
# Usage:   sbatch run_g16.sh input.gjf
# The inputs in this directory use %nprocshared=8 and %mem=12GB,
# i.e. the same core count and ~75% of the memory requested above.
# A full-featured wrapper is Gaussian/scripts/sub_g16.sh (see chapter 17).

set -eo pipefail

INP=${1:?usage: sbatch run_g16.sh input.gjf}
BASE=${INP%.*}

# --- Gaussian environment (adapt to your cluster) ---------------------------
# Option A: environment module
module load gaussian/16 2>/dev/null || true
# Option B: manual setup
# export g16root=/opt/software
# source "$g16root/g16/bsd/g16.profile"

# --- per-job scratch directory, removed on exit ------------------------------
SCR_BASE=${SCR_BASE:-/tmp/$USER}
export GAUSS_SCRDIR="$SCR_BASE/g16_${SLURM_JOB_ID:-$$}"
mkdir -p "$GAUSS_SCRDIR"
trap 'rm -rf "$GAUSS_SCRDIR"' EXIT

# --- sanity check: %nprocshared in the input vs. allocated cores -------------
NP=$(grep -i -m1 '^%nprocshared' "$INP" | cut -d= -f2 | tr -d '[:space:]' || true)
if [ -n "${NP:-}" ] && [ -n "${SLURM_CPUS_PER_TASK:-}" ] && [ "$NP" != "$SLURM_CPUS_PER_TASK" ]; then
    echo "WARNING: %nprocshared=$NP but --cpus-per-task=$SLURM_CPUS_PER_TASK" >&2
fi

echo "Job $BASE started on $(hostname) at $(date)"
g16 < "$INP" > "$BASE.log"

# --- check termination and make a formatted checkpoint -----------------------
if tail -n 1 "$BASE.log" | grep -q "Normal termination"; then
    echo "Normal termination: $BASE.log"
    # formchk every %chk named in the input (Link1 steps may share one file)
    for chk in $(grep -i '^%chk' "$INP" | cut -d= -f2 | tr -d '\r' | sort -u); do
        if [ -f "$chk" ]; then
            formchk "$chk" "${chk%.chk}.fchk" > /dev/null 2>&1 || echo "formchk failed for $chk" >&2
        fi
    done
else
    echo "Error termination, last lines of $BASE.log:" >&2
    tail -n 5 "$BASE.log" >&2
    exit 1
fi
