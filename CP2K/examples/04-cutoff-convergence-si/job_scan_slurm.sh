#!/bin/bash
# Slurm job: run both steps of the CUTOFF / REL_CUTOFF scan on one node.
# Each calculation uses all MPI ranks of the job; calculations run one after another.
#
#   sbatch job_scan_slurm.sh
#   sbatch --export=ALL,CUTOFF2=400 job_scan_slurm.sh   # CUTOFF used in step 2
#
#SBATCH --job-name=cp2k-cutoff
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail

# ---- environment: adapt to your cluster (see chapter 1 of the CP2K tutorial) ----
# module load cp2k/2026.2
# source /path/to/your/cp2k/environment/script
# export CP2K_DATA_DIR=/path/to/cp2k/data
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

cd "${SLURM_SUBMIT_DIR:-.}"
command -v cp2k.psmp >/dev/null || { echo "cp2k.psmp not found in PATH" >&2; exit 1; }

export CP2K="srun cp2k.psmp"
export NPAR=1

bash run_cutoff_scan.sh
CUTOFF=${CUTOFF2:-350} bash run_relcutoff_scan.sh

echo "done: see cutoff_scan.csv and relcutoff_scan.csv"
