#!/bin/bash
#SBATCH -J omx-gaas-soc
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out
#
# GaAs band structure with and without SOC, then SOC splittings at Gamma and L.
# bandgnu13 is compiled separately: gcc bandgnu13.c -lm -o bandgnu13

# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}

mpirun -np ${NP} openmx GaAs_SOC_band.dat   -nt ${OMP_NUM_THREADS} > gaas_soc.std   2>&1
mpirun -np ${NP} openmx GaAs_noSOC_band.dat -nt ${OMP_NUM_THREADS} > gaas_nosoc.std 2>&1

python3 soc_splitting.py gaas_soc.Band   --at G L
python3 soc_splitting.py gaas_nosoc.Band --at G L

if command -v bandgnu13 > /dev/null; then
  bandgnu13 gaas_soc.Band      # -> gaas_soc.GNUBAND, gaas_soc.BANDDAT1
  bandgnu13 gaas_nosoc.Band
fi
