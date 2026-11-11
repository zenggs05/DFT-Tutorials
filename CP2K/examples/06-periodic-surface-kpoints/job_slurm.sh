#!/bin/bash
# Slurm job for the chapter-6 examples.
#   sbatch job_slurm.sh pt111_slab_dipole.inp      # one input
#   sbatch job_slurm.sh                            # default: the three single-point inputs
#   sbatch --export=ALL,KSCAN=1 job_slurm.sh       # additionally run the k-point scan
#
#SBATCH --job-name=cp2k-periodic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail

# ---- environment: adapt to your cluster (see chapter 1 of the CP2K tutorial) ----
# module load cp2k/2026.2
# source /path/to/your/cp2k/environment/script
# export CP2K_DATA_DIR=/path/to/cp2k/data
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

cd "${SLURM_SUBMIT_DIR:-.}"
command -v cp2k.psmp >/dev/null || { echo "cp2k.psmp not found in PATH" >&2; exit 1; }

if [ $# -gt 0 ]; then
  inputs="$*"
else
  inputs="h2o_isolated_mt.inp si_primitive_kpoints.inp pt111_slab_dipole.inp"
fi

for inp in $inputs; do
  out=${inp%.inp}.out
  echo "=== $inp -> $out  ($(date))"
  srun cp2k.psmp -i "$inp" -o "$out"
  grep -q "PROGRAM ENDED AT" "$out" || echo "    ERROR: $inp did not finish normally" >&2
done

if [ "${KSCAN:-0}" = "1" ]; then
  CP2K="srun cp2k.psmp" bash run_kpoint_scan.sh
fi
