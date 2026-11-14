#!/usr/bin/env bash
# =============================================================================================
# sub_g16.sh - submit Gaussian 16 jobs to Slurm or PBS with consistent cores / memory / scratch.
#
# For every input (.gjf/.com) the script writes a job script <stem>.sub_g16.sh next to the input
# and submits it (sbatch or qsub).  The generated job script
#   * loads Gaussian (module and/or g16.profile), sets ulimit -s unlimited,
#   * creates a private GAUSS_SCRDIR per job and removes it on exit, scancel/qdel or timeout
#     (keep it for failed jobs with --keep-scr-on-fail),
#   * runs g16 in the background and waits, so signals are handled immediately,
#   * checks "Normal termination", runs formchk on every %chk of the input,
#   * on failure prints the error termination with the lines above it and exits non-zero
#     (so that afterok dependencies stop).
#
# Resources:
#   default    : cores / memory are READ from %nprocshared (%nproc) and %mem of the input; the
#                scheduler request is cores = nprocshared, memory = %mem / 0.8 (rounded up).
#                Missing Link 0 lines are supplied through GAUSS_PDEF / GAUSS_MDEF (80 % of --mem).
#   --rewrite  : REWRITE %nprocshared / %mem in every job step of the input to match -c / -m
#                (%mem = 80 % of -m); the original is kept as <input>.bak.
#
# Usage examples
#   bash sub_g16.sh water.gjf                          # resources from the input
#   bash sub_g16.sh -c 16 -m 32G -t 48:00:00 ts.gjf     # explicit request (warns on mismatch)
#   bash sub_g16.sh -c 32 -m 120G --rewrite big.gjf    # make %nprocshared/%mem match the request
#   bash sub_g16.sh -a -c 8 -m 16G                     # every *.gjf / *.com in this directory
#   bash sub_g16.sh -s pbs -q normal job.gjf           # PBS/Torque (nodes=1:ppn=N)
#   bash sub_g16.sh -s pbspro job.gjf                  # PBS Pro (select=1:ncpus=N:mem=...)
#   bash sub_g16.sh -d afterok:123456 freq.gjf         # dependency on another job
#   bash sub_g16.sh -n job.gjf                         # dry run: write the job script only
#
# Site settings: edit the USER SETTINGS block below or use the environment variables
#   SUB_G16_MODULE, g16root, SCR_BASE, SUB_G16_PARTITION, SUB_G16_ACCOUNT.
#
# Dependencies: bash, Slurm (sbatch) or PBS (qsub), Gaussian 16 (g16, formchk) on the compute nodes.
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
# =============================================================================================
set -euo pipefail

# ------------------------------ USER SETTINGS (edit for your cluster) ------------------------
G16_MODULE="${SUB_G16_MODULE:-}"          # e.g. "gaussian/16-C.02"; empty = no module load
G16ROOT="${g16root:-/opt}"                # directory that contains g16/ (g16.profile lives in g16/bsd)
SCR_BASE_DEFAULT="${SCR_BASE:-/tmp}"      # node-local disk for scratch; per-job dir is created below it
DEF_CPUS=8                                # used when neither -c nor %nprocshared is given
DEF_MEM_GB=16                             # used when neither -m nor %mem is given
DEF_TIME="24:00:00"
DEF_PARTITION="${SUB_G16_PARTITION:-}"
DEF_ACCOUNT="${SUB_G16_ACCOUNT:-}"
MEM_FRACTION=80                           # %mem (or GAUSS_MDEF) = this percentage of the request
# ---------------------------------------------------------------------------------------------

usage() {
    cat <<'EOF'
Usage: bash sub_g16.sh [options] input.gjf [more.gjf ...]
       bash sub_g16.sh [options] -a [DIR ...]

Resources
  -c, --cpus N          cores (Slurm --cpus-per-task / PBS ppn); default: %nprocshared of the input
  -m, --mem SIZE        memory for the scheduler, e.g. 32G, 64000M; default: %mem / 0.8
  -t, --time HH:MM:SS   wall time (default 24:00:00)
  -q, --partition NAME  Slurm partition / PBS queue
  -A, --account NAME    account / project
  -J, --name NAME       job name (default: input stem; one input only)
      --rewrite         rewrite %nprocshared and %mem of the input to match -c/-m (backup .bak)
      --pin             Slurm: pin threads to the allocated CPUs via GAUSS_CDEF (Cpus_allowed_list)
Scheduler
  -s, --scheduler S     slurm | pbs (Torque nodes/ppn) | pbspro (select) ; default: auto
  -d, --dependency DEP  e.g. afterok:123456 (passed as --dependency / -W depend=)
  -x, --extra "ARGS"    extra arguments for sbatch/qsub, e.g. "--exclusive"
Files
  -a, --all             submit every *.gjf and *.com in the given directories (default: .)
      --scr-base DIR    parent of the per-job scratch directory (default: $SCR_BASE or /tmp)
      --keep-scr-on-fail  keep the scratch directory when the job fails
      --no-formchk      do not run formchk after a successful job
      --module NAME     module to load (default: $SUB_G16_MODULE)
      --g16root DIR     Gaussian root (default: $g16root or /opt)
  -n, --dry-run         write <stem>.sub_g16.sh but do not submit
  -h, --help            this help

The generated job script is <input stem>.sub_g16.sh in the directory of the input; the Gaussian
output is <input stem>.log.  Link 0 lines in the input override GAUSS_PDEF/GAUSS_MDEF.
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }

# ------------------------------ argument parsing --------------------------------------------
CPUS="" MEM="" TIME="$DEF_TIME" PART="$DEF_PARTITION" ACCOUNT="$DEF_ACCOUNT" NAME=""
REWRITE=0 PIN=0 SCHED="" DEP="" EXTRA="" ALL=0 SCRBASE="$SCR_BASE_DEFAULT" KEEPFAIL=0
FORMCHK=1 DRYRUN=0
ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -c|--cpus) CPUS="${2:?missing value for $1}"; shift 2 ;;
        -m|--mem) MEM="${2:?missing value for $1}"; shift 2 ;;
        -t|--time) TIME="${2:?missing value for $1}"; shift 2 ;;
        -q|--partition|--queue) PART="${2:?missing value for $1}"; shift 2 ;;
        -A|--account) ACCOUNT="${2:?missing value for $1}"; shift 2 ;;
        -J|--name) NAME="${2:?missing value for $1}"; shift 2 ;;
        --rewrite) REWRITE=1; shift ;;
        --pin) PIN=1; shift ;;
        -s|--scheduler) SCHED="${2:?missing value for $1}"; shift 2 ;;
        -d|--dependency) DEP="${2:?missing value for $1}"; shift 2 ;;
        -x|--extra) EXTRA="${2:?missing value for $1}"; shift 2 ;;
        -a|--all) ALL=1; shift ;;
        --scr-base) SCRBASE="${2:?missing value for $1}"; shift 2 ;;
        --keep-scr-on-fail) KEEPFAIL=1; shift ;;
        --no-formchk) FORMCHK=0; shift ;;
        --module) G16_MODULE="${2:?missing value for $1}"; shift 2 ;;
        --g16root) G16ROOT="${2:?missing value for $1}"; shift 2 ;;
        -n|--dry-run) DRYRUN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        --) shift; ARGS+=("$@"); break ;;
        -*) usage >&2; die "unknown option $1" ;;
        *) ARGS+=("$1"); shift ;;
    esac
done

[ -n "$CPUS" ] && { [[ "$CPUS" =~ ^[1-9][0-9]*$ ]] || die "-c needs a positive integer"; }
[[ "$TIME" =~ ^([0-9]+-)?[0-9]+(:[0-9]{2}){0,2}$ ]] || die "bad --time '$TIME' (use HH:MM:SS)"

if [ -z "$SCHED" ]; then
    if command -v sbatch >/dev/null 2>&1; then SCHED=slurm
    elif command -v qsub >/dev/null 2>&1; then SCHED=pbs
    elif [ "$DRYRUN" -eq 1 ]; then SCHED=slurm
    else die "neither sbatch nor qsub found; use -s slurm|pbs|pbspro (and -n for a dry run)"
    fi
fi
case "$SCHED" in slurm|pbs|pbspro) ;; *) die "unknown scheduler '$SCHED'" ;; esac

# ------------------------------ helpers -----------------------------------------------------
# memory string -> MB.  Scheduler style (32G, 64000M, 64gb; bare number = MB) when $2=sched,
# Gaussian %mem style (8GB, 800MB, 1GW; bare number = 8-byte words) when $2=g16.
to_mb() {
    local v u
    v=$(echo "$1" | tr 'a-z' 'A-Z' | tr -d ' ')
    [[ "$v" =~ ^([0-9]+(\.[0-9]+)?)([KMGT]?)([BW]?)$ ]] || return 1
    u="${BASH_REMATCH[3]}${BASH_REMATCH[4]}"
    if [ -z "$u" ]; then
        if [ "$2" = g16 ]; then u=W; else u=M; fi
    fi
    awk -v n="${BASH_REMATCH[1]}" -v u="$u" 'BEGIN{
        f = 1.0
        if (u ~ /^K/) f = 1.0/1024; else if (u ~ /^G/) f = 1024; else if (u ~ /^T/) f = 1048576
        else if (u ~ /^[BW]$/) f = 1.0/1048576
        if (u ~ /W$/) f *= 8
        printf "%d\n", n * f + 0.5 }'
}

# all Link 0 values of one keyword (case-insensitive), one per line
link0_values() {
    tr -d '\r' < "$2" | awk -v key="$1" '
        { line = tolower($0); sub(/^[ \t]+/, "", line) }
        index(line, "%" key "=") == 1 { split($0, a, "="); v = a[2]; gsub(/[ \t]/, "", v); print v }'
}

count_cpus() {   # "0-3,8,10-11" -> 7
    echo "$1" | tr ',' '\n' | awk -F- 'NF==2{n+=$2-$1+1; next} NF==1 && $1!=""{n++} END{print n+0}'
}

sanitize_name() {
    local s
    s=$(echo "$1" | tr -c 'A-Za-z0-9_.-' '_' | sed 's/_*$//')
    [[ "$s" =~ ^[A-Za-z] ]] || s="g$s"
    if [ "$SCHED" != slurm ]; then s="${s:0:15}"; fi
    echo "$s"
}

# rewrite %nprocshared/%nproc/%cpu and %mem in every job step (Link0 block before each route)
rewrite_link0() {
    local f="$1" np="$2" memmb="$3"
    cp -p "$f" "$f.bak"
    awk -v np="$np" -v mem="${memmb}MB" '
        BEGIN { in0 = 1; seen_np = 0; seen_mem = 0 }
        function flush_missing() {
            if (!seen_np) print "%nprocshared=" np
            if (!seen_mem) print "%mem=" mem
        }
        {
            l = $0; sub(/\r$/, "", l); t = tolower(l); sub(/^[ \t]+/, "", t)
            if (t ~ /^--link1--/) { print l; in0 = 1; seen_np = 0; seen_mem = 0; next }
            if (in0 && t ~ /^%(nprocshared|nproc|cpu)=/) {
                if (!seen_np) print "%nprocshared=" np; seen_np = 1; next }
            if (in0 && t ~ /^%mem=/) { if (!seen_mem) print "%mem=" mem; seen_mem = 1; next }
            if (in0 && t ~ /^#/) { flush_missing(); in0 = 0 }
            print l
        }' "$f.bak" > "$f"
    echo "  rewrote Link 0 of $f: %nprocshared=$np %mem=${memmb}MB (backup $f.bak)"
}

# ------------------------------ collect inputs ----------------------------------------------
INPUTS=()
if [ "$ALL" -eq 1 ]; then
    [ ${#ARGS[@]} -eq 0 ] && ARGS=(.)
    for d in "${ARGS[@]}"; do
        [ -d "$d" ] || die "-a expects directories, '$d' is not one"
        while IFS= read -r f; do INPUTS+=("$f"); done < <(find "$d" -maxdepth 1 -type f \( -name '*.gjf' -o -name '*.com' \) | sort)
    done
else
    INPUTS=("${ARGS[@]+"${ARGS[@]}"}")
fi
[ ${#INPUTS[@]} -gt 0 ] || { usage >&2; die "no input files"; }
if [ -n "$NAME" ] && [ ${#INPUTS[@]} -gt 1 ]; then die "-J/--name needs a single input"; fi
if [ "$REWRITE" -eq 1 ] && { [ -z "$CPUS" ] || [ -z "$MEM" ]; }; then
    die "--rewrite needs both -c and -m"
fi

# ------------------------------ per input ---------------------------------------------------
NSUB=0
for INP in "${INPUTS[@]}"; do
    [ -f "$INP" ] || { warn "$INP not found - skipped"; continue; }
    DIR=$(cd "$(dirname "$INP")" && pwd)
    BASE=$(basename "$INP")
    STEM="${BASE%.*}"
    INPUT="$DIR/$BASE"
    echo "== $INPUT"

    # input hygiene: CRLF and the final blank line
    if grep -q $'\r' "$INPUT"; then
        cp -p "$INPUT" "$INPUT.crlf.bak"
        tr -d '\r' < "$INPUT.crlf.bak" > "$INPUT"
        echo "  converted CRLF line endings (backup $BASE.crlf.bak)"
    fi
    if [ -n "$(tail -n 1 "$INPUT" | tr -d ' \t')" ]; then
        printf '\n\n' >> "$INPUT"; echo "  appended the missing final blank line"
    fi

    # Link 0 resources of the input
    NP_IN=$(link0_values nprocshared "$INPUT"; link0_values nproc "$INPUT")
    NP_IN=$(echo "$NP_IN" | awk 'NF && $1+0>m{m=$1+0} END{if(m)print m}')
    CPU_IN=$(link0_values cpu "$INPUT" | head -n 1)
    if [ -n "$CPU_IN" ]; then
        NC=$(count_cpus "$CPU_IN")
        [ -z "$NP_IN" ] || [ "$NC" -gt "$NP_IN" ] && NP_IN=$NC
    fi
    MEM_IN_MB=""
    while IFS= read -r m; do
        [ -n "$m" ] || continue
        mb=$(to_mb "$m" g16) || { warn "cannot parse %mem=$m"; continue; }
        if [ -z "$MEM_IN_MB" ] || [ "$mb" -gt "$MEM_IN_MB" ]; then MEM_IN_MB=$mb; fi
    done < <(link0_values mem "$INPUT")

    # scheduler request
    NCPU=${CPUS:-${NP_IN:-$DEF_CPUS}}
    if [ -n "$MEM" ]; then
        REQ_MB=$(to_mb "$MEM" sched) || die "bad --mem '$MEM'"
    elif [ -n "$MEM_IN_MB" ] && [ "$REWRITE" -eq 0 ]; then
        REQ_MB=$(( (MEM_IN_MB * 100 + MEM_FRACTION - 1) / MEM_FRACTION ))
        REQ_MB=$(( (REQ_MB + 1023) / 1024 * 1024 ))
    else
        REQ_MB=$(( DEF_MEM_GB * 1024 ))
    fi
    G16_MB=$(( REQ_MB * MEM_FRACTION / 100 ))

    if [ "$REWRITE" -eq 1 ]; then
        [ -n "$CPU_IN" ] && warn "%cpu=$CPU_IN replaced by %nprocshared=$NCPU"
        rewrite_link0 "$INPUT" "$NCPU" "$G16_MB"
        NP_IN=$NCPU; MEM_IN_MB=$G16_MB
    else
        if [ -n "$NP_IN" ] && [ "$NP_IN" -ne "$NCPU" ]; then
            warn "input asks for $NP_IN cores but $NCPU are requested (use --rewrite or change -c)"
        fi
        if [ -n "$MEM_IN_MB" ] && [ "$MEM_IN_MB" -gt $(( REQ_MB * 90 / 100 )) ]; then
            warn "%mem (${MEM_IN_MB} MB) is above 90% of the request (${REQ_MB} MB): risk of OOM kill"
        fi
    fi
    [ -n "$CPU_IN" ] && [ "$REWRITE" -eq 0 ] && \
        warn "input pins cores with %cpu=$CPU_IN; on shared nodes prefer %nprocshared (or --pin)"

    # %chk files for formchk
    CHKS=$(link0_values chk "$INPUT" | awk 'NF{ if ($0 !~ /\.chk$/) $0 = $0 ".chk"; print }' | sort -u | tr '\n' ' ')

    JNAME=$(sanitize_name "${NAME:-$STEM}")
    JOB="$DIR/$STEM.sub_g16.sh"
    LOG="$STEM.log"
    REQ_GB=$(( (REQ_MB + 1023) / 1024 ))

    {
        echo "#!/usr/bin/env bash"
        echo "# generated by sub_g16.sh (DFT-Tutorials) on $(date '+%Y-%m-%d %H:%M') for $BASE"
        case "$SCHED" in
            slurm)
                echo "#SBATCH --job-name=$JNAME"
                echo "#SBATCH --nodes=1"
                echo "#SBATCH --ntasks=1"
                echo "#SBATCH --cpus-per-task=$NCPU"
                echo "#SBATCH --mem=${REQ_MB}M"
                echo "#SBATCH --time=$TIME"
                echo "#SBATCH --output=$DIR/%x-%j.out"
                [ -n "$PART" ] && echo "#SBATCH --partition=$PART"
                [ -n "$ACCOUNT" ] && echo "#SBATCH --account=$ACCOUNT"
                ;;
            pbs)
                echo "#PBS -N $JNAME"
                echo "#PBS -l nodes=1:ppn=$NCPU"
                echo "#PBS -l mem=${REQ_GB}gb"
                echo "#PBS -l walltime=$TIME"
                echo "#PBS -j oe"
                [ -n "$PART" ] && echo "#PBS -q $PART"
                [ -n "$ACCOUNT" ] && echo "#PBS -A $ACCOUNT"
                ;;
            pbspro)
                echo "#PBS -N $JNAME"
                echo "#PBS -l select=1:ncpus=$NCPU:mem=${REQ_GB}gb"
                echo "#PBS -l walltime=$TIME"
                echo "#PBS -j oe"
                [ -n "$PART" ] && echo "#PBS -q $PART"
                [ -n "$ACCOUNT" ] && echo "#PBS -A $ACCOUNT"
                ;;
        esac
        cat <<EOF

set -uo pipefail
# ---- values filled in by sub_g16.sh -------------------------------------------------------
WORKDIR="$DIR"
INPUT="$BASE"
LOG="$LOG"
NCPU=$NCPU
G16_MB=$G16_MB
CHKS="$CHKS"
G16_MODULE="$G16_MODULE"
G16ROOT="$G16ROOT"
SCR_BASE="$SCRBASE"
KEEP_SCR_ON_FAIL=$KEEPFAIL
RUN_FORMCHK=$FORMCHK
PIN=$PIN
EOF
        cat <<'EOF'
# --------------------------------------------------------------------------------------------
cd "$WORKDIR" || exit 1
if [ -n "$G16_MODULE" ]; then module load "$G16_MODULE" || exit 1; fi
if ! command -v g16 >/dev/null 2>&1; then
    export g16root="$G16ROOT"
    # shellcheck disable=SC1091
    source "$g16root/g16/bsd/g16.profile" || { echo "ERROR: cannot source $g16root/g16/bsd/g16.profile" >&2; exit 1; }
fi
command -v g16 >/dev/null 2>&1 || { echo "ERROR: g16 not found" >&2; exit 1; }
ulimit -s unlimited 2>/dev/null || true

# threads and memory used when the input has no %nprocshared / %mem (Link 0 wins otherwise)
export GAUSS_PDEF=$NCPU
export GAUSS_MDEF="${G16_MB}MB"
if [ "$PIN" -eq 1 ] && [ -r /proc/self/status ]; then
    CPULIST=$(awk '/^Cpus_allowed_list/ {print $2}' /proc/self/status)
    NLIST=$(echo "$CPULIST" | tr ',' '\n' | awk -F- 'NF==2{n+=$2-$1+1; next} NF==1 && $1!=""{n++} END{print n+0}')
    if [ -n "$CPULIST" ] && [ "$NLIST" -eq "$NCPU" ]; then
        unset GAUSS_PDEF; export GAUSS_CDEF="$CPULIST"
    fi
fi

JOBID="${SLURM_JOB_ID:-${PBS_JOBID:-$$}}"; JOBID="${JOBID%%.*}"
export GAUSS_SCRDIR="$SCR_BASE/${USER:-$(id -un)}/g16_${JOBID}"
mkdir -p "$GAUSS_SCRDIR" || { echo "ERROR: cannot create $GAUSS_SCRDIR" >&2; exit 1; }
STATUS=1
G16PID=""
cleanup() {
    if [ -n "$G16PID" ] && kill -0 "$G16PID" 2>/dev/null; then kill -TERM "$G16PID" 2>/dev/null; wait "$G16PID" 2>/dev/null; fi
    if [ "$STATUS" -ne 0 ] && [ "$KEEP_SCR_ON_FAIL" -eq 1 ]; then
        echo "scratch kept for inspection: $GAUSS_SCRDIR"
    else
        rm -rf -- "${GAUSS_SCRDIR:?}"
    fi
}
trap cleanup EXIT
trap 'echo "signal received - stopping g16" >&2; exit 143' TERM INT

echo "host $(hostname)  job $JOBID  start $(date)"
echo "input $INPUT  cores $NCPU  GAUSS_MDEF ${G16_MB}MB  scratch $GAUSS_SCRDIR"
df -h "$GAUSS_SCRDIR" 2>/dev/null | tail -n 1
g16 < "$INPUT" > "$LOG" 2>&1 &
G16PID=$!
wait "$G16PID"
RC=$?
G16PID=""

if [ "$RC" -eq 0 ] && tail -n 3 "$LOG" | grep -q "Normal termination" && ! grep -q "Error termination" "$LOG"; then
    STATUS=0
    echo "Normal termination ($(grep -c 'Normal termination' "$LOG") job step(s)) at $(date)"
    if [ "$RUN_FORMCHK" -eq 1 ]; then
        for c in $CHKS; do
            if [ -f "$c" ]; then formchk "$c" "${c%.chk}.fchk" >/dev/null && echo "formchk $c -> ${c%.chk}.fchk"
            else echo "WARNING: $c not found (no formchk)"; fi
        done
    fi
else
    echo "FAILED (g16 exit code $RC) - see $LOG" >&2
    LN=$(grep -n "Error termination" "$LOG" | head -n 1 | cut -d: -f1)
    if [ -n "$LN" ]; then
        START=$(( LN > 25 ? LN - 25 : 1 ))
        sed -n "${START},$((LN + 2))p" "$LOG" >&2
    else
        echo "no 'Error termination' line: killed (time/memory limit?) or g16 did not start" >&2
        tail -n 20 "$LOG" >&2
    fi
    echo "diagnose with: python g16_error_doctor.py $LOG" >&2
fi
exit "$STATUS"
EOF
    } > "$JOB"
    chmod +x "$JOB"
    echo "  request: $NCPU cores, ${REQ_MB} MB, $TIME ($SCHED); g16 memory ${G16_MB} MB; job script $JOB"

    # submit
    CMD=()
    if [ "$SCHED" = slurm ]; then
        CMD=(sbatch --chdir="$DIR")
        [ -n "$DEP" ] && CMD+=(--dependency="$DEP")
    else
        CMD=(qsub -d "$DIR")
        [ "$SCHED" = pbspro ] && CMD=(qsub)
        [ -n "$DEP" ] && CMD+=(-W "depend=$DEP")
    fi
    if [ -n "$EXTRA" ]; then
        read -r -a XARR <<< "$EXTRA"
        CMD+=("${XARR[@]}")
    fi
    CMD+=("$JOB")
    if [ "$DRYRUN" -eq 1 ]; then
        echo "  dry run: ${CMD[*]}"
    else
        (cd "$DIR" && "${CMD[@]}") || { warn "submission of $BASE failed"; continue; }
        NSUB=$((NSUB + 1))
    fi
done
[ "$DRYRUN" -eq 1 ] || echo "submitted $NSUB job(s)"
