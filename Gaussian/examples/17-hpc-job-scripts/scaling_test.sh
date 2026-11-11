#!/bin/bash
# scaling_test.sh - measure how a Gaussian job scales with the number of cores, inside ONE
# allocation, so that every run sees the same hardware.
#
# Usage:  sbatch scaling_test.sh benzene_scaling_freq.gjf "1 2 4 8 16 32"
# Output: scaling_<name>.dat with cores, wall seconds, speedup and parallel efficiency.
# The input must NOT contain %nprocshared/%cpu (they would override GAUSS_CDEF).
#
#SBATCH --job-name=g16-scaling
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --exclusive
#SBATCH --time=12:00:00
#SBATCH --output=%x-%j.out

set -uo pipefail
INPUT=${1:?usage: sbatch scaling_test.sh input.gjf "1 2 4 8 16 32"}
NLIST=${2:-1 2 4 8 16 32}
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
NAME=${INPUT%.*}

# module load gaussian/16-C.02
export g16root=${g16root:-/opt}
# shellcheck disable=SC1091
source "$g16root/g16/bsd/g16.profile" || exit 1

expand_cpus() {
    local IFS=, r
    for r in $1; do
        if [[ $r == *-* ]]; then seq "${r%-*}" "${r#*-}"; else echo "$r"; fi
    done
}
mapfile -t CPUS < <(expand_cpus "$(awk '/^Cpus_allowed_list/ {print $2}' /proc/self/status)")
MEM_MB=${SLURM_MEM_PER_NODE:-65536}
export GAUSS_MDEF="$(( MEM_MB * 80 / 100 ))MB"       # same memory for every run
export GAUSS_SCRDIR="${SCR_BASE:-/tmp}/$USER/g16_${SLURM_JOB_ID:-$$}"
mkdir -p "$GAUSS_SCRDIR" || exit 1
trap 'rm -rf -- "${GAUSS_SCRDIR:?}"' EXIT

OUT=scaling_$NAME.dat
printf "%6s %10s %8s %10s\n" cores wall_s speedup efficiency > "$OUT"
T1=""
for n in $NLIST; do
    [ "$n" -le "${#CPUS[@]}" ] || { echo "skip $n cores (only ${#CPUS[@]} allocated)"; continue; }
    GAUSS_CDEF=$(IFS=,; echo "${CPUS[*]:0:$n}")
    export GAUSS_CDEF
    log=${NAME}_n${n}.log
    t0=$(date +%s.%N)
    g16 < "$INPUT" > "$log"
    t1=$(date +%s.%N)
    tail -n 1 "$log" | grep -q "Normal termination" || echo "WARNING: $log did not terminate normally"
    wall=$(awk -v a="$t0" -v b="$t1" 'BEGIN {printf "%.1f", b - a}')
    [ -n "$T1" ] || T1=$(awk -v w="$wall" -v n="$n" 'BEGIN {print w * n}')   # serial-time estimate from the first run
    awk -v n="$n" -v w="$wall" -v t1="$T1" \
        'BEGIN {s = t1 / w; printf "%6d %10.1f %8.2f %10.2f\n", n, w, s, s / n}' >> "$OUT"
done
cat "$OUT"
