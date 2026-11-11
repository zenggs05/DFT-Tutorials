#!/bin/bash
#SBATCH --job-name=cp2k-ch14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%j.out
#
# Slurm job for the chapter-14 examples (vibrational analysis).
# Usage (one step per job):
#   sbatch --export=ALL,STEP=h2o_opt job_slurm.sh
#   sbatch --export=ALL,STEP=h2o_vib job_slurm.sh
#   sbatch --export=ALL,STEP=co_opt  job_slurm.sh
#   sbatch --export=ALL,STEP=co_vib  job_slurm.sh
#   sbatch --export=ALL,STEP=ts_freq job_slurm.sh   # needs ts_dimer.xyz (chapter 13)
# NPROC_REP is chosen so that, for H2O and CO/Pt, all displaced geometries of
# one sign run at once (ts_freq: 6 replicas, i.e. 5 rounds per sign).

set -euo pipefail

STEP=${STEP:-h2o_opt}
CP2K_EXE=${CP2K_EXE:-cp2k.psmp}

# ---- environment: keep exactly one of these ----
# source /path/to/cp2k/install/cp2k_env
# module load cp2k/2026.2
# source ~/miniconda3/bin/activate cp2k_env
# export CP2K_DATA_DIR=/path/to/cp2k/data

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
cd "${SLURM_SUBMIT_DIR:-.}"
command -v "$CP2K_EXE" >/dev/null || { echo "$CP2K_EXE not found in PATH" >&2; exit 1; }
NTASKS=${SLURM_NTASKS:-48}

# number of displaced coordinates per sign -> ranks per replica
case "$STEP" in
  h2o_opt) ARGS=(-i h2o_opt_tight.inp); OUT=h2o_opt_tight.out ;;
  h2o_vib) NPR=$(( NTASKS >= 9 ? NTASKS / 9 : 1 ))
           ARGS=(-E NPROC_REP=${NPR} -i h2o_vib_ir.inp); OUT=h2o_vib_ir.out ;;
  co_opt)  ARGS=(-i co_pt111_opt.inp); OUT=co_pt111_opt.out ;;
  co_vib)  NPR=$(( NTASKS >= 6 ? NTASKS / 6 : 1 ))
           ARGS=(-E NPROC_REP=${NPR} -i co_on_pt_partial_hessian.inp); OUT=co_on_pt_partial_hessian.out ;;
  ts_freq) NPR=$(( NTASKS >= 6 ? NTASKS / 6 : 1 ))
           ARGS=(-E NPROC_REP=${NPR} -i ts_freq_check.inp); OUT=ts_freq_check.out ;;
  *) echo "unknown STEP=$STEP" >&2; exit 1 ;;
esac

echo "STEP=$STEP  ranks=$NTASKS  args: ${ARGS[*]}"
srun "$CP2K_EXE" "${ARGS[@]}" -o "$OUT"

if grep -q "PROGRAM ENDED AT" "$OUT"; then
  echo "finished: $OUT"
  grep "VIB|Frequency" "$OUT" || true
else
  echo "CP2K did not finish normally, check $OUT" >&2
  exit 1
fi
