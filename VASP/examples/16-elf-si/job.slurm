#!/bin/bash
#SBATCH --job-name=si-elf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1

[ -s POTCAR ] || { echo "POTCAR missing: cat \$VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR" >&2; exit 1; }
mpirun -np ${SLURM_NTASKS:-8} ${VASP_EXE:-vasp_std} > vasp.out 2>&1

grep -q "General timing" OUTCAR && ls -lh ELFCAR CHGCAR
grep -E "ELF|KPAR" vasp.out | head
