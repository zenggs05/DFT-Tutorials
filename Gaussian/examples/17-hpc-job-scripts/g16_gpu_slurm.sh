#!/bin/bash
# g16_gpu_slurm.sh - Slurm template for a Gaussian 16 GPU job (Linux only; HF/DFT energies,
# gradients and frequencies of LARGER molecules; small jobs and MP2/CCSD gain little).
# Supported GPUs: K40, K80, P100 (Rev. B.01+), V100 (Rev. C.01+), A100 (Rev. C.02 only).
#
# Usage:   sbatch g16_gpu_slurm.sh job.gjf
#
# Gaussian rules (see chapter 17):
#   * GAUSS_CDEF (= %CPU) lists ALL cores used, INCLUDING the GPU control cores;
#   * GAUSS_GDEF (= %GPUCPU) is "gpu-list=control-cpu-list"; one control core per GPU,
#     preferably on the CPU socket the GPU is attached to (check: nvidia-smi topo -m);
#   * memory is shared equally among threads, so %mem must be large:
#     threads x (memory one GPU needs), e.g. ~12 GB per thread for 16 GB GPUs.
#
#SBATCH --job-name=g16-gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:2
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
##SBATCH --partition=gpu

set -uo pipefail

INPUT=${1:?usage: sbatch g16_gpu_slurm.sh job.gjf}
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
LOG=${INPUT%.*}.log

# module load gaussian/16-C.02
export g16root=${g16root:-/opt}
# shellcheck disable=SC1091
source "$g16root/g16/bsd/g16.profile" || exit 1

# ---------- cores ----------
expand_cpus() {
    local IFS=, r
    for r in $1; do
        if [[ $r == *-* ]]; then seq "${r%-*}" "${r#*-}"; else echo "$r"; fi
    done
}
CPULIST=$(awk '/^Cpus_allowed_list/ {print $2}' /proc/self/status)
mapfile -t CPUS < <(expand_cpus "$CPULIST")
export GAUSS_CDEF=$CPULIST

# ---------- GPUs and their control cores ----------
# Inside the job Slurm usually exposes only the allocated GPUs (CUDA_VISIBLE_DEVICES),
# numbered from 0. Verify with the nvidia-smi output printed below.
NGPU=$(nvidia-smi -L 2>/dev/null | wc -l)
[ "$NGPU" -ge 1 ] || { echo "ERROR: no GPU visible"; exit 1; }
GPULIST=$(seq -s, 0 $(( NGPU - 1 )))
# Naive choice: the first NGPU allocated cores control the GPUs. Better: pick cores whose
# "CPU Affinity" in 'nvidia-smi topo -m' matches each GPU, and set CTRL by hand.
CTRL=$(IFS=,; echo "${CPUS[*]:0:$NGPU}")
export GAUSS_GDEF="${GPULIST}=${CTRL}"

# ---------- memory (~80% of --mem) ----------
MEM_MB=${SLURM_MEM_PER_NODE:-245760}
export GAUSS_MDEF="$(( MEM_MB * 80 / 100 ))MB"

# ---------- scratch ----------
export GAUSS_SCRDIR="${SCR_BASE:-/tmp}/$USER/g16_${SLURM_JOB_ID:-$$}"
mkdir -p "$GAUSS_SCRDIR" || exit 1
trap 'rm -rf -- "${GAUSS_SCRDIR:?}"' EXIT

echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
nvidia-smi -L
nvidia-smi topo -m 2>/dev/null | head -n $(( NGPU + 2 ))
echo "GAUSS_CDEF=$GAUSS_CDEF  GAUSS_GDEF=$GAUSS_GDEF  GAUSS_MDEF=$GAUSS_MDEF"

g16 < "$INPUT" > "$LOG"
tail -n 1 "$LOG" | grep -q "Normal termination" && exit 0
grep -B 6 "Error termination" "$LOG" | tail -n 14
exit 1
