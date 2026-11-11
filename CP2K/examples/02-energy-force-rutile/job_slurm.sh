#!/bin/bash
# job_slurm.sh - Slurm job script for the rutile TiO2 examples (CP2K 2026.2, cp2k.psmp)
#
# Submit (default input: input_file.inp):
#     sbatch job_slurm.sh
# Choose another input without editing the script:
#     sbatch --export=ALL,INP=rutile_kpoints.inp job_slurm.sh
#     sbatch --export=ALL,INP=rutile_supercell_gamma.inp job_slurm.sh
#
# Layout: MPI ranks = nodes * ntasks-per-node, OpenMP threads = cpus-per-task.
#SBATCH --job-name=rutile-ef
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log
##SBATCH --partition=compute       # remove one '#' and set your partition
##SBATCH --account=myproject       # remove one '#' if your site needs an account

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# ---- CP2K environment: keep ONE of the following lines -------------------
# module load cp2k/2026.2
# source $HOME/src/cp2k-2026.2/install/cp2k_env
# eval "$(conda shell.bash hook)" && conda activate cp2k_env
# ---------------------------------------------------------------------------

CP2K_EXE=${CP2K_EXE:-cp2k.psmp}
INP=${INP:-input_file.inp}
OUT=${INP%.inp}.out

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}   # srun >= 22.05 does not inherit it
ulimit -s unlimited

command -v "${CP2K_EXE}" > /dev/null || { echo "ERROR: ${CP2K_EXE} not in PATH" >&2; exit 1; }
[[ -f "${INP}" ]] || { echo "ERROR: input ${INP} not found" >&2; exit 1; }

# cp2k -o appends to an existing file: keep the old output under another name
[[ -f "${OUT}" ]] && mv "${OUT}" "${OUT}.$(date +%Y%m%d%H%M%S)"

echo "Job ${SLURM_JOB_ID:-local} on $(hostname): ${SLURM_NTASKS:-1} MPI x ${OMP_NUM_THREADS} OMP"
echo "Input ${INP} -> output ${OUT}"
srun --cpu-bind=cores "${CP2K_EXE}" -i "${INP}" -o "${OUT}"
rc=$?

if grep -q "PROGRAM ENDED AT" "${OUT}"; then
  echo "CP2K finished normally."
  grep "ENERGY| Total FORCE_EVAL" "${OUT}" | tail -1
else
  echo "CP2K did NOT finish normally (exit code ${rc}); check ${OUT} and the Slurm log." >&2
fi
exit ${rc}
