#!/bin/bash
# run_g16_slurm.sh - minimal Slurm script for the small water examples.
# Usage:  sbatch run_g16_slurm.sh h2o_b3lyp_opt.gjf
# Keep %nprocshared in the .gjf equal to --cpus-per-task, and %mem at
# about 75-85 % of --mem. A full-featured template is in chapter 17.
#SBATCH --job-name=g16-h2o
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=%x-%j.out

set -u
INPUT=${1:-h2o_b3lyp_opt.gjf}

# --- load Gaussian (site specific: use 'module load' OR source the profile)
# module load gaussian/16-C.02
export g16root=${g16root:-/opt}
source "$g16root/g16/bsd/g16.profile"

# --- per-job scratch directory on a local disk, removed at exit
SCR_BASE=${SCR_BASE:-/tmp}
export GAUSS_SCRDIR="$SCR_BASE/$USER/g16_${SLURM_JOB_ID:-$$}"
mkdir -p "$GAUSS_SCRDIR"
trap 'rm -rf "$GAUSS_SCRDIR"' EXIT

cd "${SLURM_SUBMIT_DIR:-.}"
echo "host: $(hostname)  input: $INPUT  scratch: $GAUSS_SCRDIR"
g16 < "$INPUT" > "${INPUT%.*}.log"

tail -n 1 "${INPUT%.*}.log"
