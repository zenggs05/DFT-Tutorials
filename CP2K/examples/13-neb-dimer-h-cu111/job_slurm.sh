#!/bin/bash
#SBATCH --job-name=cp2k-ch13
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#
# Slurm job for the chapter-13 example (H on Cu(111): NEB + Dimer).
# Usage (one step per job, in this order):
#   sbatch --export=ALL,STEP=opt_initial job_slurm.sh
#   sbatch --export=ALL,STEP=opt_final   job_slurm.sh
#   sbatch --export=ALL,STEP=neb         job_slurm.sh   # 56 ranks = 7 replicas x 8
#   sbatch --export=ALL,STEP=vib         job_slurm.sh
#   sbatch --export=ALL,STEP=dimer       job_slurm.sh
# Adjust partition/account and the environment block for your cluster.

set -euo pipefail

STEP=${STEP:-opt_initial}
CP2K_EXE=${CP2K_EXE:-cp2k.psmp}

# ---- environment: keep exactly one of these ----
# source /path/to/cp2k/install/cp2k_env      # toolchain / build_cp2k.sh build
# module load cp2k/2026.2                    # cluster module
# source ~/miniconda3/bin/activate cp2k_env  # conda-forge build
# export CP2K_DATA_DIR=/path/to/cp2k/data    # only if your build does not set it

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
cd "${SLURM_SUBMIT_DIR:-.}"
command -v "$CP2K_EXE" >/dev/null || { echo "$CP2K_EXE not found in PATH" >&2; exit 1; }

NTASKS=${SLURM_NTASKS:-56}
case "$STEP" in
  opt_initial) ARGS=(-i geo_opt_endpoint.inp); OUT=opt_initial.out ;;
  opt_final)   ARGS=(-E COORD=final.xyz -E PROJ=opt_final -i geo_opt_endpoint.inp); OUT=opt_final.out ;;
  neb)
    # NPROC_REP = ranks per replica; 7 replicas must share NTASKS evenly
    NPR=$(( NTASKS / 7 ))
    ARGS=(-E NPROC_REP=${NPR} -i h_diffusion_cineb.inp); OUT=h_diffusion_cineb.out ;;
  vib)
    NPR=$(( NTASKS / 6 ))
    ARGS=(-E NPROC_REP=${NPR} -i ts_guess_vib.inp); OUT=ts_guess_vib.out ;;
  dimer)       ARGS=(-i dimer_ts.inp); OUT=dimer_ts.out ;;
  *) echo "unknown STEP=$STEP" >&2; exit 1 ;;
esac

echo "STEP=$STEP  ranks=$NTASKS  args: ${ARGS[*]}"
srun "$CP2K_EXE" "${ARGS[@]}" -o "$OUT"

if grep -q "PROGRAM ENDED AT" "$OUT"; then
  echo "finished: $OUT"
else
  echo "CP2K did not finish normally, check $OUT (and *-BAND*.out for NEB)" >&2
  exit 1
fi
