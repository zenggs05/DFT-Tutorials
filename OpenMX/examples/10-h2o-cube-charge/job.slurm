#!/bin/bash
#SBATCH -J omx-h2o-cube
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH -t 00:30:00
#SBATCH -o slurm-%j.out

# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

INPUT=H2O_cube.dat
mpirun -np ${SLURM_NTASKS} openmx ${INPUT} -nt ${OMP_NUM_THREADS} > h2o.std 2>&1

# quick checks
grep "Utot\." h2o.out
grep -A 8 "Mulliken populations" h2o.out | head -12
grep -A 14 "Voronoi charges" h2o.out
