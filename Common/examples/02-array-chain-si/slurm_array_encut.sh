#!/bin/bash
# =====================================================================
# ENCUT convergence scan of Si as a Slurm job array.
# Task i runs ENCUT = VALUES[i] in its own directory encut_<E>/.
#   sbatch slurm_array_encut.sh
#   bash collect_encut.sh            # after all tasks have finished
# =====================================================================
#SBATCH --job-name=si-encut
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --array=0-6%4               # 7 values, at most 4 running at once
#SBATCH --output=%x-%A_%a.out

VALUES=(250 300 350 400 450 500 550)

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}" || exit 1
E=${VALUES[$SLURM_ARRAY_TASK_ID]}
[ -n "${E}" ] || { echo "no value for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
[ -f POTCAR ] || { echo "POTCAR missing (see README)" >&2; exit 1; }

D=encut_${E}
mkdir -p "${D}"
cp POSCAR KPOINTS POTCAR "${D}/"
sed "s/@ENCUT@/${E}/g" INCAR.encut.template > "${D}/INCAR"

cd "${D}" || exit 1
srun --ntasks="${SLURM_NTASKS}" --cpus-per-task=1 vasp_std > vasp.out 2>&1
grep "free  energy   TOTEN" OUTCAR | tail -n 1
