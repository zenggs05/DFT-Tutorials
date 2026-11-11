#!/bin/bash
#SBATCH -J omx_pdos
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 06:00:00
#SBATCH -o slurm-%j.out

# PDOS of monolayer MoS2 and AFM NiO with OpenMX 4.0 (edit modules/paths).
# DosMain is built with "make DosMain" in openmx4.0/source.
# module load openmx/4.0
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
TOOLS=../09-dos-pdos-cdia

mpirun -np ${SLURM_NTASKS} openmx MoS2_monolayer_pdos.dat -nt ${OMP_NUM_THREADS} > mos2.std
mpirun -np ${SLURM_NTASKS} openmx NiO_AFM_pdos.dat        -nt ${OMP_NUM_THREADS} > nio_afm.std

bash ${TOOLS}/run_dosmain.sh mos2    tetra 0.1
bash ${TOOLS}/run_dosmain.sh nio_afm tetra 0.1
python3 ${TOOLS}/sum_pdos.py mos2    --dat MoS2_monolayer_pdos.dat --resolve-m --plot --emin -8 --emax 5 > mos2_pdos.txt
python3 ${TOOLS}/sum_pdos.py nio_afm --dat NiO_AFM_pdos.dat --by atom --resolve-m --plot --emin -9 --emax 5 > nio_pdos.txt
