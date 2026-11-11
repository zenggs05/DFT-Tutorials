#!/bin/bash
#SBATCH -J h2o_nve
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH -o slurm-%j.out

# Example Slurm script for OpenMX 4.0 (edit module names / paths for your cluster)
# module load openmx/4.0
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_STACKSIZE=512M

NAME=H2O_NVE
# OpenMX appends to .ene/.md: start from clean files for a fresh run
rm -f ${NAME}.ene ${NAME}.md

mpirun -np ${SLURM_NTASKS} openmx ${NAME}.dat -nt ${OMP_NUM_THREADS} > ${NAME}.std

python3 nve_check.py ${NAME} > nve_check.log 2>&1
