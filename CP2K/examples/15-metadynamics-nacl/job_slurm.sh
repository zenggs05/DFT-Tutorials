#!/bin/bash
#SBATCH --job-name=cp2k-ch15
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%j.out
#
# Slurm job for the chapter-15 examples (NaCl in water: metadynamics, PLUMED,
# blue moon).  Usage:
#   sbatch --export=ALL,STEP=equil,LEVEL=xtb    job_slurm.sh
#   sbatch --export=ALL,STEP=cvmon,LEVEL=xtb    job_slurm.sh   # unbiased CV fluctuations
#   sbatch --export=ALL,STEP=wtmetad,LEVEL=xtb,STEPS=200000 job_slurm.sh   # 100 ps pre-sampling
#   sbatch --export=ALL,STEP=wtmetad,LEVEL=dft  job_slurm.sh
#   sbatch --export=ALL,STEP=wt2d,LEVEL=xtb     job_slurm.sh   # 2D: d(Na-Cl) + CN(Na-O)
#   sbatch --export=ALL,STEP=plumed,LEVEL=dft   job_slurm.sh   # needs PLUMED build
#   sbatch --export=ALL,STEP=bluemoon,LEVEL=dft job_slurm.sh   # all windows in sequence
#   sbatch --export=ALL,STEP=bluemoon,LEVEL=dft,CHAIN=0,TARGETS=3.20 job_slurm.sh
#          (one independent window whose start.xyz came from pick_windows.py)
# STEPS (optional) overrides the default number of MD steps of the chosen input.
# Long AIMD runs: continue from <PROJECT>-1.restart (see README.md).

set -euo pipefail

STEP=${STEP:-equil}
LEVEL=${LEVEL:-xtb}
CP2K_EXE=${CP2K_EXE:-cp2k.psmp}

# ---- environment: keep exactly one of these ----
# source /path/to/cp2k/install/cp2k_env
# module load cp2k/2026.2
# source ~/miniconda3/bin/activate cp2k_env
# export CP2K_DATA_DIR=/path/to/cp2k/data

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
cd "${SLURM_SUBMIT_DIR:-.}"
command -v "$CP2K_EXE" >/dev/null || { echo "$CP2K_EXE not found in PATH" >&2; exit 1; }

EXTRA=()
if [ -n "${STEPS:-}" ]; then EXTRA=(-E STEPS="$STEPS"); fi

run() {  # run <input> <output> [extra -E options]
  local inp=$1 out=$2; shift 2
  srun "$CP2K_EXE" -E LEVEL="$LEVEL" ${EXTRA[@]+"${EXTRA[@]}"} "$@" -i "$inp" -o "$out"
  grep -q "PROGRAM ENDED AT" "$out" || { echo "CP2K did not finish: $out" >&2; exit 1; }
}

case "$STEP" in
  equil)    run nacl_equil.inp nacl_equil_${LEVEL}.out
            bash last_frame.sh nacl_equil_${LEVEL}-pos-1.xyz nacl_eq.xyz ;;
  cvmon)    run cv_monitor.inp cv_monitor_${LEVEL}.out ;;
  wtmetad)  run nacl_dissociation_wtmetad.inp wtmetad_${LEVEL}.out ;;
  wt2d)     run nacl_wtmetad_2d.inp wt2d_${LEVEL}.out ;;
  plumed)   "$CP2K_EXE" -v | grep -q plumed2 || { echo "this CP2K has no PLUMED support" >&2; exit 1; }
            run plumed_metad.inp plumed_metad_${LEVEL}.out ;;
  bluemoon) MPIRUN="srun" LEVEL="$LEVEL" CP2K_EXE="$CP2K_EXE" bash run_bluemoon.sh ;;
  *) echo "unknown STEP=$STEP" >&2; exit 1 ;;
esac
echo "STEP=$STEP LEVEL=$LEVEL finished"
