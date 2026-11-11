#!/bin/bash
#SBATCH --job-name=co-box
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
# Gamma-only executable: about 2x faster, half the memory, KPOINTS must be Gamma 1 1 1
VASP_EXE=${VASP_EXE:-vasp_gam}

mpirun -np ${SLURM_NTASKS:-8} $VASP_EXE > vasp.out 2>&1

grep -q "reached required accuracy" OUTCAR && echo "relaxation converged" || echo "NOT converged"
grep "F=" OSZICAR | tail -1
