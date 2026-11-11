#!/bin/bash
#SBATCH -J omx_band2
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 04:00:00
#SBATCH -o slurm-%j.out

# Graphene and ferromagnetic bcc Fe band structures with OpenMX 4.0 (edit modules/paths).
# module load openmx/4.0
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mpirun -np ${SLURM_NTASKS} openmx Graphene_band.dat    -nt ${OMP_NUM_THREADS} > Graphene_band.std
mpirun -np ${SLURM_NTASKS} openmx Fe_bcc_spin_band.dat -nt ${OMP_NUM_THREADS} > Fe_bcc_spin_band.std

for SYS in graphene fe; do
  bandgnu13 ${SYS}.Band
  python3 ../08-band-cdia-si/band_parse.py ${SYS}.Band --emin -10 --emax 6 --plot ${SYS}_band.png > ${SYS}_gap.txt
done
