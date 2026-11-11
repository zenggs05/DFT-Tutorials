#!/bin/bash
#SBATCH --job-name=zno-hse-lcalceps
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out
#SBATCH --partition=normal
# Adjust partition/account/module lines for your cluster.

# module load vasp/6.6.1        # example only; GPU builds older than 6.6.1 give
#                               # wrong Born charges with LCALCEPS (known issue 93)
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP=${VASP:-vasp_std}
RUN="srun ${VASP}"              # or: mpirun -np ${SLURM_NTASKS} ${VASP}

if [ ! -s POTCAR ]; then
    PP=${VASP_PP_PATH:?please set VASP_PP_PATH}
    cat "$PP/potpaw_PBE/Zn/POTCAR" "$PP/potpaw_PBE/O/POTCAR" > POTCAR
fi

cp INCAR.pbe INCAR
${RUN} > vasp_pbe.out 2>&1
cp OUTCAR OUTCAR.pbe

cp INCAR.hse_lcalceps INCAR
${RUN} > vasp_hse.out 2>&1
grep -A 5 "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects)" OUTCAR
grep -A 14 "BORN EFFECTIVE CHARGES (including local field effects)" OUTCAR
