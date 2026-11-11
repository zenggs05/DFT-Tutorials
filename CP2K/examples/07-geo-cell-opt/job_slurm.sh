#!/bin/bash
#SBATCH --job-name=cp2k-ch07
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#
# Generic Slurm job for the chapter-7 examples.
# Usage:  sbatch --export=ALL,INP=si_cell_opt.inp job_slurm.sh
# Adjust partition/account lines and the environment block for your cluster.

set -euo pipefail

INP=${INP:-h2o_geo_opt_bfgs.inp}
OUT=${INP%.inp}.out
CP2K_EXE=${CP2K_EXE:-cp2k.psmp}

# ---- environment: keep exactly one of these, matching how CP2K was installed ----
# source /path/to/cp2k/install/cp2k_env      # toolchain / build_cp2k.sh build
# module load cp2k/2026.2                    # cluster module
# source ~/miniconda3/bin/activate cp2k_env  # conda-forge build
# export CP2K_DATA_DIR=/path/to/cp2k/data    # only if your build does not set it

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

echo "CP2K binary : $(command -v "$CP2K_EXE")"
echo "Input       : $INP"
echo "MPI x OMP   : ${SLURM_NTASKS:-1} x ${OMP_NUM_THREADS}"

srun "$CP2K_EXE" -i "$INP" -o "$OUT"

if grep -q "PROGRAM ENDED AT" "$OUT"; then
    echo "CP2K finished normally: $OUT"
else
    echo "CP2K did not finish. Continue from the *-1.restart file (see README.md)." >&2
    exit 1
fi
