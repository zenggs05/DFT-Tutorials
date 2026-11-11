#!/bin/bash
#SBATCH --job-name=mos2-strain
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
#SBATCH --partition=normal
# Adjust partition/account/module lines for your cluster.

# module load vasp/6.6.1        # example only
ulimit -s unlimited
export OMP_NUM_THREADS=1
srun ${VASP:-vasp_std} > vasp.out 2>&1   # or: mpirun -np ${SLURM_NTASKS} vasp_std
