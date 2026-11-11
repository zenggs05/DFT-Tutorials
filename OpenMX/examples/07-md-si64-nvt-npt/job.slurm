#!/bin/bash
#SBATCH -J omx_md
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH -o slurm-%j.out

# Run one MD input of this example. Select it with NAME, e.g.
#   sbatch job.slurm                                   # Si64_NVT_VS (default)
#   sbatch --export=ALL,NAME=Si64_NVT_NH job.slurm
#   sbatch --export=ALL,NAME=Si8_NPT_VS_PR job.slurm
# Set FRESH=0 to keep existing .ene/.md (e.g. when continuing from a .dat# file).

# module load openmx/4.0        # edit for your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_STACKSIZE=512M

NAME=${NAME:-Si64_NVT_VS}
INPUT=${INPUT:-${NAME}.dat}
FRESH=${FRESH:-1}

if [ "${FRESH}" = "1" ]; then
  # OpenMX appends to .ene/.md/.npt.ene: remove them for a fresh run
  rm -f ${NAME}.ene ${NAME}.md ${NAME}.npt.ene
fi

mpirun -np ${SLURM_NTASKS} openmx ${INPUT} -nt ${OMP_NUM_THREADS} > ${NAME}.std

python3 md_analysis.py ${NAME} --skip 200 > ${NAME}_analysis.log 2>&1
