#!/bin/bash
#SBATCH -J omx-mos2-soc
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out
#
# 1) SCF + bands with SOC  2) splittings at K  3) spin texture along the path with kSpin
# kSpin is built once:  cd openmx4.0/source && make kSpin   (copy it to ../work if needed)

# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}
KSPIN=${KSPIN:-$HOME/openmx4.0/work/kSpin}

mpirun -np ${NP} openmx MoS2_SOC_band.dat -nt ${OMP_NUM_THREADS} > mos2_soc.std 2>&1
python3 ../12-gaas-soc-band/soc_splitting.py mos2_soc.Band --at G K Kp

if [ -x "${KSPIN}" ]; then
  mpirun -np ${NP} ${KSPIN} MoS2_SOC_band.dat > kspin.log 2>&1
  ls mos2_soc_BD.* | head
fi
