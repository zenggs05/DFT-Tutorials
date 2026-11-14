#!/bin/bash
# GPU job for MLIP training on Slurm. Choose the code with TASK:
#   sbatch --export=ALL,TASK=nep   job_mlip_gpu.sh     # GPUMD nep: nep.in + train.xyz [+ test.xyz]
#   sbatch --export=ALL,TASK=mace  job_mlip_gpu.sh     # mace_run_train --config mace_config.yaml
#   sbatch --export=ALL,TASK=dp    job_mlip_gpu.sh     # DeePMD-kit (PyTorch backend)
#SBATCH --job-name=mlip-train
#SBATCH --partition=gpu            # EDIT: GPU partition of your cluster
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out

set -euo pipefail
TASK=${TASK:-nep}
cd "$SLURM_SUBMIT_DIR"
# module load cuda                 # EDIT: cluster specific
# source ~/miniforge3/etc/profile.d/conda.sh && conda activate mlip
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
nvidia-smi || true

case "$TASK" in
  nep)
    # export PATH=/path/to/GPUMD/src:$PATH   # EDIT: directory with the nep executable
    nep > nep.log 2>&1 ;;
  mace)
    mace_run_train --config mace_config.yaml > mace.log 2>&1 ;;
  dp)
    dp --pt train deepmd_input.json > dp_train.log 2>&1
    dp --pt freeze -o frozen_model.pth >> dp_train.log 2>&1 ;;
  *)
    echo "unknown TASK=$TASK" >&2; exit 1 ;;
esac
