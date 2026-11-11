#!/bin/bash
# =====================================================================
# Generic Slurm job array: run the same calculation for a list of
# parameter values (ENCUT, k-mesh, lattice constant, CUTOFF, ...).
#
# Layout expected in the submit directory:
#   template/            all input files; one of them contains the
#                        placeholder @VALUE@ (e.g. template/INCAR)
#   slurm_array.sh       this script
# Every array task creates run_<VALUE>/ from template/ and runs there.
#
#   sbatch slurm_array.sh                  # 7 tasks, max 4 at a time
#   sbatch --array=2,5 slurm_array.sh      # rerun only tasks 2 and 5
# =====================================================================
#SBATCH --job-name=scan
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --array=0-6%4               # task ids 0..6, at most 4 running
#SBATCH --output=%x-%A_%a.out       # %A = array master id, %a = task id

# ----------------------------- user settings -------------------------
VALUES=(250 300 350 400 450 500 550)   # one entry per array task id
PLACEHOLDER="@VALUE@"
EXE="vasp_std"
ARGS=""

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}" || exit 1
ID=${SLURM_ARRAY_TASK_ID}
if [ "${ID}" -ge "${#VALUES[@]}" ]; then
    echo "array id ${ID} has no value in VALUES" >&2; exit 1
fi
VAL=${VALUES[$ID]}
RUN="run_${VAL}"

mkdir -p "${RUN}"
cp -r template/. "${RUN}/"
# replace the placeholder in every text file of the run directory
grep -rl "${PLACEHOLDER}" "${RUN}" | xargs -r sed -i "s/${PLACEHOLDER}/${VAL}/g"
[ -f POTCAR ] && [ ! -f "${RUN}/POTCAR" ] && cp POTCAR "${RUN}/"

cd "${RUN}" || exit 1
echo "task ${ID}: value ${VAL} in $(pwd) on ${SLURM_JOB_NODELIST}"
srun --ntasks="${SLURM_NTASKS}" --cpus-per-task=1 ${EXE} ${ARGS} > run.log 2>&1
echo "exit code $?"
