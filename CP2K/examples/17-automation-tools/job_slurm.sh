#!/bin/bash
# Generic Slurm template for one CP2K input (MPI x OpenMP). Chapter 17 of the CP2K tutorial.
#   sbatch job_slurm.sh my.inp                 # output: my.out
#   sbatch --ntasks-per-node=16 --cpus-per-task=4 job_slurm.sh my.inp
#   sbatch job_slurm.sh --ase ase_h2o_opt.py   # run an ASE driver script instead
#
#SBATCH --job-name=cp2k
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail

# ---- environment: adapt to your cluster ----
# module load cp2k/2026.2
# source /path/to/cp2k/install/cp2k_env
# conda activate cp2k_env
# export CP2K_DATA_DIR=/path/to/cp2k/data
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}   # srun >= 22.05 no longer inherits --cpus-per-task
export OMP_PLACES=cores
export OMP_PROC_BIND=close
ulimit -s unlimited || true

cd "${SLURM_SUBMIT_DIR:-.}"
command -v cp2k.psmp >/dev/null || { echo "cp2k.psmp not found in PATH" >&2; exit 1; }

if [ "${1:-}" = "--ase" ]; then
  # ASE drives CP2K through the shell mode; stdin must reach MPI rank 0
  script=${2:?usage: job_slurm.sh --ase script.py}
  export ASE_CP2K_COMMAND=${ASE_CP2K_COMMAND:-"mpirun -np ${SLURM_NTASKS:-1} cp2k.psmp -s"}
  python "$script"
  exit $?
fi

inp=${1:?usage: sbatch job_slurm.sh input.inp}
out=${inp%.inp}.out
# "cp2k -c" prints "SUCCESS, the input could be parsed correctly." on success
# (output captured first: with pipefail, "cmd | grep -q" can fail through SIGPIPE)
chk=$(cp2k.psmp -c "$inp" 2>&1 || true)
if ! grep -q "SUCCESS, the input could be parsed" <<< "$chk"; then
  echo "syntax check failed for $inp:" >&2
  echo "$chk" | tail -20 >&2
  exit 1
fi
rm -f "$out"                          # -o appends to an existing file

echo "start $(date)  ranks=${SLURM_NTASKS:-1}  threads=$OMP_NUM_THREADS  host=$(hostname)"
srun --cpu-bind=cores cp2k.psmp -i "$inp" -o "$out" || true

if grep -q "PROGRAM ENDED AT" "$out"; then
  echo "finished normally $(date)"
  grep -m1 "DBCSR| MPI: Number of processes" "$out" || true
  grep "ENERGY| Total FORCE_EVAL" "$out" | tail -1 || true
else
  echo "CP2K did not finish: check the tail of $out" >&2
  ls -t ./*-1.restart ./*-RESTART.wfn ./*-RESTART.kp 2>/dev/null | head -3 | \
    sed 's/^/restart candidate: /' >&2 || true
  exit 1
fi
