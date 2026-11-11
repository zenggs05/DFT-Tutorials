#!/bin/bash
# g16_linda_slurm.sh - Slurm template for a multi-node Gaussian 16 job with Linda.
# Requires the Linda version of Gaussian (separate licence; Rev. C.01 needs Linda 9.2),
# the same $g16root path on every node, and password-less ssh between the nodes.
# Worth it for large HF/DFT energies, gradients and frequencies; CCSD/EOM-CCSD are SMP only.
#
# Usage:   sbatch g16_linda_slurm.sh job.gjf
# One Linda worker per node; each worker uses all cores of its node (shared memory).
# %CPU and %mem apply PER NODE.
#
#SBATCH --job-name=g16-linda
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --exclusive
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out

set -uo pipefail

INPUT=${1:?usage: sbatch g16_linda_slurm.sh job.gjf}
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
LOG=${INPUT%.*}.log

# module load gaussian/16-C.02-linda
export g16root=${g16root:-/opt}
# shellcheck disable=SC1091
source "$g16root/g16/bsd/g16.profile" || exit 1

# ---------- worker list: every node of the allocation ----------
NODES=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | paste -sd, -)
export GAUSS_WDEF=$NODES          # same as %LindaWorkers=node1,node2
export GAUSS_SDEF=ssh             # start workers with ssh (same as %UseSSH)
# For numerical-derivative jobs (e.g. Freq=Anharm) Gaussian wants one extra worker on the
# master node:  export GAUSS_WDEF="$(hostname -s):2,<other nodes>"

# ---------- per-node cores and memory ----------
NCPU=${SLURM_CPUS_PER_TASK:-32}
export GAUSS_CDEF="0-$(( NCPU - 1 ))"             # --exclusive: whole node, so 0..N-1 is safe
MEM_MB=${SLURM_MEM_PER_NODE:-122880}
export GAUSS_MDEF="$(( MEM_MB * 80 / 100 ))MB"

# ---------- scratch: must exist on EVERY node (local disk) ----------
export GAUSS_SCRDIR="${SCR_BASE:-/tmp}/$USER/g16_${SLURM_JOB_ID:-$$}"
srun --ntasks-per-node=1 mkdir -p "$GAUSS_SCRDIR"
trap 'srun --ntasks-per-node=1 rm -rf -- "${GAUSS_SCRDIR:?}"' EXIT
# Linda debugging output if workers do not start:
# export GAUSS_LFLAGS="-vv"

echo "workers: $GAUSS_WDEF | per node: GAUSS_CDEF=$GAUSS_CDEF GAUSS_MDEF=$GAUSS_MDEF"
g16 < "$INPUT" > "$LOG"
tail -n 1 "$LOG" | grep -q "Normal termination" && exit 0
grep -B 6 "Error termination" "$LOG" | tail -n 14
exit 1
