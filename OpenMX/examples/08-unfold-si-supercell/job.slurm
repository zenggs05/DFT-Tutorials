#!/bin/bash
#SBATCH -J omx_unfold
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 04:00:00
#SBATCH -o slurm-%j.out

# Band unfolding of Si (perfect 8-atom cell and Si7Ge1) with OpenMX 4.0 (edit modules/paths).
# module load openmx/4.0
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mpirun -np ${SLURM_NTASKS} openmx Si_supercell_unfold.dat -nt ${OMP_NUM_THREADS} > si8_unfold.std
mpirun -np ${SLURM_NTASKS} openmx Si7Ge1_unfold.dat       -nt ${OMP_NUM_THREADS} > si7ge1_unfold.std

BAND=../08-band-cdia-si/si.BANDDAT1
OPT=""
[ -f ${BAND} ] && OPT="--band ${BAND}"
python3 plot_unfold.py si8_unfold    ${OPT} --out si8_unfold.png
python3 plot_unfold.py si7ge1_unfold ${OPT} --out si7ge1_unfold.png
