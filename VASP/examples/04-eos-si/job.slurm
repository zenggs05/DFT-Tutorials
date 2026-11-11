#!/bin/bash
#SBATCH --job-name=si-eos
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

# Submitted from inside each V_* directory by "bash run_eos.sh submit"
module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}

mpirun -np ${SLURM_NTASKS:-8} $VASP_EXE > vasp.out 2>&1
