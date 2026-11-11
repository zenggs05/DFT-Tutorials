#!/bin/bash
#SBATCH --job-name=si-mlff
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out
#SBATCH --partition=normal
# Train -> refit -> run in one job.  Adjust partition/account/module lines.
# The MLFF code is MPI-parallel only: keep OMP_NUM_THREADS=1.
set -e
# module load vasp/6.6.1        # example only; needs a scaLAPACK build
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP=${VASP:-vasp_std}
RUN="srun ${VASP}"              # or: mpirun -np ${SLURM_NTASKS} ${VASP}

if [ ! -s POTCAR ]; then
    cat "${VASP_PP_PATH:?please set VASP_PP_PATH}/potpaw_PBE/Si/POTCAR" > POTCAR
fi
mkdir -p 1_train 2_refit 3_run

# ---- stage 1: on-the-fly training ----
cp INCAR_train 1_train/INCAR
cp POSCAR KPOINTS POTCAR 1_train/
(cd 1_train && ${RUN} > vasp.out 2>&1)
python3 parse_ml_logfile.py -f 1_train/ML_LOGFILE | tail -20

# ---- stage 2: refit for fast prediction ----
cp INCAR_refit 2_refit/INCAR
cp POSCAR KPOINTS POTCAR 2_refit/
cp 1_train/ML_ABN 2_refit/ML_AB
(cd 2_refit && ${RUN} > vasp.out 2>&1)
grep "ERR" 2_refit/ML_LOGFILE | grep -v "#" | tail -1

# ---- stage 3: production with the force field only ----
cp INCAR_run 3_run/INCAR
cp 1_train/CONTCAR 3_run/POSCAR
cp KPOINTS POTCAR 3_run/
cp 2_refit/ML_FFN 3_run/ML_FF
(cd 3_run && ${RUN} > vasp.out 2>&1)
python3 parse_ml_logfile.py -f 3_run/ML_LOGFILE | tail -5
