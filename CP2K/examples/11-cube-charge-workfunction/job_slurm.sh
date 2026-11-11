#!/bin/bash
# Slurm job for the chapter-11 examples.
# Usage:  sbatch job_slurm.sh h2o_density_mo_cubes.inp
#         sbatch job_slurm.sh pt111_workfunction.inp
#SBATCH --job-name=cp2k-cube
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.log

set -euo pipefail

INP=${1:-h2o_density_mo_cubes.inp}
OUT=${INP%.inp}.out

# --- environment: adapt ONE of the following lines to your cluster ---------
# module load cp2k/2026.2
# source /path/to/cp2k/tools/toolchain/install/setup    # toolchain build
# source /path/to/cp2k/install/cp2k_env                 # build_cp2k.sh build
# conda activate cp2k_env                               # conda-forge build
# ---------------------------------------------------------------------------

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export CP2K_DATA_DIR=${CP2K_DATA_DIR:-/path/to/cp2k/data}

srun cp2k.psmp -i "$INP" -o "$OUT"

if grep -q "PROGRAM ENDED AT" "$OUT"; then
    echo "finished: $OUT"
    ls -1 ./*.cube ./*.molden 2>/dev/null || true
    grep "Fermi energy:" "$OUT" | tail -n 1 || true
else
    echo "CP2K did not finish normally, check $OUT" >&2
    exit 1
fi
