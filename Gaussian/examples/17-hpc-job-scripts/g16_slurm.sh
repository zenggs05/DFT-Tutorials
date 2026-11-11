#!/bin/bash
# g16_slurm.sh - Slurm template: one Gaussian 16 job, one node, shared-memory parallel.
#
# Usage:   sbatch g16_slurm.sh job.gjf
#          sbatch --cpus-per-task=32 --mem=120G --time=48:00:00 g16_slurm.sh big.gjf
#
# Policy:  keep %nprocshared/%cpu/%mem OUT of the .gjf. This script sets
#          GAUSS_CDEF (the core list Slurm gave us) and GAUSS_MDEF (~80% of the job memory).
#          Link 0 lines in the input have HIGHER priority than these variables,
#          so the script warns if it finds them.
# Exit code: 0 only if the log ends with "Normal termination" (needed for --dependency=afterok).
#
#SBATCH --job-name=g16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
##SBATCH --partition=cpu
##SBATCH --hint=nomultithread

set -uo pipefail

INPUT=${1:?usage: sbatch g16_slurm.sh job.gjf}
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1
[ -f "$INPUT" ] || { echo "ERROR: $INPUT not found in $PWD"; exit 1; }
NAME=${INPUT%.*}
LOG=$NAME.log
MEM_FRACTION=${MEM_FRACTION:-80}      # percent of the Slurm memory given to Gaussian (%mem)
SCR_BASE=${SCR_BASE:-/tmp}            # local disk; many sites offer /scratch/$USER or $TMPDIR
KEEP_SCR_ON_FAIL=${KEEP_SCR_ON_FAIL:-0}

# ---------- 1. Gaussian environment (site specific) ----------
# module load gaussian/16-C.02
export g16root=${g16root:-/opt}
# shellcheck disable=SC1091
source "$g16root/g16/bsd/g16.profile" || { echo "ERROR: cannot source g16.profile"; exit 1; }

# ---------- 2. cores: pin threads to exactly the CPUs Slurm allocated ----------
expand_cpus() {                        # "0-3,8,10-11" -> one CPU number per line
    local IFS=, r
    for r in $1; do
        if [[ $r == *-* ]]; then seq "${r%-*}" "${r#*-}"; else echo "$r"; fi
    done
}
NCPU=${SLURM_CPUS_PER_TASK:-1}
CPULIST=$(awk '/^Cpus_allowed_list/ {print $2}' /proc/self/status)
NLIST=$(expand_cpus "$CPULIST" | wc -l)
if [ -n "$CPULIST" ] && [ "$NLIST" -eq "$NCPU" ]; then
    export GAUSS_CDEF=$CPULIST
    unset GAUSS_PDEF
else
    # e.g. hyperthreads visible in the cgroup: fall back to a plain thread count
    echo "NOTE: allowed CPUs '$CPULIST' ($NLIST) != cpus-per-task ($NCPU); using GAUSS_PDEF=$NCPU"
    unset GAUSS_CDEF
    export GAUSS_PDEF=$NCPU
fi

# ---------- 3. memory: ~80% of what Slurm granted ----------
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
    MEM_MB=$SLURM_MEM_PER_NODE
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
    MEM_MB=$(( SLURM_MEM_PER_CPU * NCPU ))
else
    MEM_MB=4096
fi
export GAUSS_MDEF="$(( MEM_MB * MEM_FRACTION / 100 ))MB"

# ---------- 4. sanity checks on the input ----------
if grep -q $'\r' "$INPUT"; then
    echo "NOTE: CRLF line endings found, converting (backup: $INPUT.crlf.bak)"
    sed -i.crlf.bak 's/\r$//' "$INPUT"
fi
if [ "$(tail -c 2 "$INPUT" | od -An -c | tr -d ' ')" != '\n\n' ]; then
    echo "NOTE: input did not end with a blank line, appending one"
    echo >> "$INPUT"
fi
if grep -qiE '^[[:space:]]*%(nprocshared|nproc|cpu|mem)[[:space:]]*=' "$INPUT"; then
    echo "WARNING: Link 0 resource lines in $INPUT override GAUSS_CDEF/GAUSS_MDEF:"
    grep -inE '^[[:space:]]*%(nprocshared|nproc|cpu|mem)[[:space:]]*=' "$INPUT"
    echo "         make sure they match --cpus-per-task=$NCPU and --mem (${MEM_MB} MB)"
fi

# ---------- 5. private scratch directory, removed at exit ----------
export GAUSS_SCRDIR="$SCR_BASE/$USER/g16_${SLURM_JOB_ID:-$$}"
mkdir -p "$GAUSS_SCRDIR" || { echo "ERROR: cannot create $GAUSS_SCRDIR"; exit 1; }
cleanup() {
    if tail -n 1 "$LOG" 2>/dev/null | grep -q "Normal termination" || [ "$KEEP_SCR_ON_FAIL" != 1 ]; then
        rm -rf -- "${GAUSS_SCRDIR:?}"
    else
        echo "scratch kept for inspection/restart: $(hostname):$GAUSS_SCRDIR"
    fi
}
trap cleanup EXIT
trap 'echo "caught SIGTERM/SIGINT (time limit or scancel)"; kill "${G16PID:-0}" 2>/dev/null; exit 143' TERM INT

# ---------- 6. run (in background + wait, so the traps fire immediately) ----------
echo "job ${SLURM_JOB_ID:-local} on $(hostname) | input $INPUT"
echo "GAUSS_CDEF=${GAUSS_CDEF:-unset} GAUSS_PDEF=${GAUSS_PDEF:-unset} GAUSS_MDEF=$GAUSS_MDEF"
echo "GAUSS_SCRDIR=$GAUSS_SCRDIR ($(df -h "$GAUSS_SCRDIR" | awk 'NR==2 {print $4}') free)"
START=$(date +%s)
g16 < "$INPUT" > "$LOG" &
G16PID=$!
wait "$G16PID"
echo "wall time: $(( $(date +%s) - START )) s"

# ---------- 7. check the result, formchk every named chk ----------
if tail -n 1 "$LOG" | grep -q "Normal termination"; then
    echo "OK: $LOG ended with Normal termination"
    grep -ioE '^[[:space:]]*%chk=[^[:space:]]+' "$INPUT" | cut -d= -f2 | sort -u | while read -r chk; do
        [[ $chk == *.chk ]] || chk=$chk.chk
        [ -f "$chk" ] && formchk "$chk" "${chk%.chk}.fchk" > /dev/null 2>&1 && echo "formchk: ${chk%.chk}.fchk"
    done
    exit 0
else
    echo "FAILED: last lines of $LOG"
    grep -B 6 "Error termination" "$LOG" | tail -n 14
    exit 1
fi
