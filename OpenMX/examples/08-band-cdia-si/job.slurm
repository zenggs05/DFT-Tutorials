#!/bin/bash
#SBATCH -J omx_band
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH -o slurm-%j.out

# Band structures of diamond and silicon with OpenMX 4.0 (edit modules/paths).
# bandgnu13 must be compiled once:  gcc bandgnu13.c -lm -o bandgnu13  (in openmx4.0/source)
# module load openmx/4.0
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

for NAME in Cdia_band Si_band; do
  mpirun -np ${SLURM_NTASKS} openmx ${NAME}.dat -nt ${OMP_NUM_THREADS} > ${NAME}.std
done

# System.Name is "cdia" and "si" in the two inputs
for SYS in cdia si; do
  bandgnu13 ${SYS}.Band
  python3 band_parse.py ${SYS}.Band --emin -15 --emax 10 --plot ${SYS}_band.png > ${SYS}_gap.txt
done
