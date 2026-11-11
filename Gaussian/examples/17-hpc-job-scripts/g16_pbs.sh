#!/bin/bash
# g16_pbs.sh - PBS/Torque template: one Gaussian 16 job on one node (shared-memory parallel).
#
# Usage:   qsub -v INPUT=job.gjf g16_pbs.sh
#          qsub -v INPUT=big.gjf -l nodes=1:ppn=32,mem=120gb g16_pbs.sh
# PBS Pro uses the select syntax instead of nodes/ppn, e.g.
#          #PBS -l select=1:ncpus=16:mem=64gb
# Keep %nprocshared/%mem OUT of the .gjf; this script sets GAUSS_PDEF and GAUSS_MDEF.
#
#PBS -N g16
#PBS -l nodes=1:ppn=16
#PBS -l mem=64gb
#PBS -l walltime=24:00:00
#PBS -j oe

set -uo pipefail

INPUT=${INPUT:?submit with: qsub -v INPUT=job.gjf g16_pbs.sh}
cd "${PBS_O_WORKDIR:-$PWD}" || exit 1
[ -f "$INPUT" ] || { echo "ERROR: $INPUT not found in $PWD"; exit 1; }
NAME=${INPUT%.*}
LOG=$NAME.log
JOBID=${PBS_JOBID%%.*}
SCR_BASE=${SCR_BASE:-/tmp}

# ---------- Gaussian environment (site specific) ----------
# module load gaussian/16-C.02
export g16root=${g16root:-/opt}
# shellcheck disable=SC1091
source "$g16root/g16/bsd/g16.profile" || { echo "ERROR: cannot source g16.profile"; exit 1; }

# ---------- cores and memory ----------
# Torque sets PBS_NUM_PPN; PBS Pro sets NCPUS. Edit MEM_GB to match the -l mem request.
NCPU=${PBS_NUM_PPN:-${NCPUS:-16}}
MEM_GB=${MEM_GB:-64}
export GAUSS_PDEF=$NCPU                                 # thread count (no pinning)
export GAUSS_MDEF="$(( MEM_GB * 1024 * 80 / 100 ))MB"  # ~80% of the request
# If the batch system confines the job to a cpuset, pinning is safer:
# CPULIST=$(awk '/^Cpus_allowed_list/ {print $2}' /proc/self/status); export GAUSS_CDEF=$CPULIST; unset GAUSS_PDEF

if grep -qiE '^[[:space:]]*%(nprocshared|nproc|cpu|mem)[[:space:]]*=' "$INPUT"; then
    echo "WARNING: Link 0 resource lines in $INPUT override GAUSS_PDEF/GAUSS_MDEF:"
    grep -inE '^[[:space:]]*%(nprocshared|nproc|cpu|mem)[[:space:]]*=' "$INPUT"
fi

# ---------- private scratch ----------
export GAUSS_SCRDIR="$SCR_BASE/$USER/g16_${JOBID:-$$}"
mkdir -p "$GAUSS_SCRDIR" || { echo "ERROR: cannot create $GAUSS_SCRDIR"; exit 1; }
trap 'rm -rf -- "${GAUSS_SCRDIR:?}"' EXIT

# ---------- run ----------
echo "job ${PBS_JOBID:-local} on $(hostname) | $INPUT | cores=$GAUSS_PDEF mem=$GAUSS_MDEF"
g16 < "$INPUT" > "$LOG"

if tail -n 1 "$LOG" | grep -q "Normal termination"; then
    echo "OK: Normal termination"
    exit 0
else
    grep -B 6 "Error termination" "$LOG" | tail -n 14
    exit 1
fi
