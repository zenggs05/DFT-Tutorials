#!/bin/bash
# Slurm job for the chapter-5 SCF examples.
#   sbatch job_slurm.sh                      # run all four inputs in the right order
#   sbatch job_slurm.sh al_metal_smear.inp   # run only the given input(s)
#
#SBATCH --job-name=cp2k-scf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail

# ---- environment: adapt to your cluster (see chapter 1 of the CP2K tutorial) ----
# module load cp2k/2026.2
# source /path/to/your/cp2k/environment/script
# export CP2K_DATA_DIR=/path/to/cp2k/data
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

cd "${SLURM_SUBMIT_DIR:-.}"
command -v cp2k.psmp >/dev/null || { echo "cp2k.psmp not found in PATH" >&2; exit 1; }

# h2o_ot.inp must run before restart_from_wfn.inp (it writes h2o_ot-RESTART.wfn)
if [ $# -gt 0 ]; then
  inputs="$*"
else
  inputs="h2o_ot.inp restart_from_wfn.inp si_diag_broyden.inp al_metal_smear.inp"
fi

for inp in $inputs; do
  out=${inp%.inp}.out
  echo "=== $inp -> $out  ($(date))"
  srun cp2k.psmp -i "$inp" -o "$out"
  if grep -q "PROGRAM ENDED AT" "$out"; then
    grep -q "SCF run converged" "$out" && echo "    SCF converged" || echo "    WARNING: no 'SCF run converged' in $out"
  else
    echo "    ERROR: $inp did not finish normally, check $out" >&2
  fi
done
