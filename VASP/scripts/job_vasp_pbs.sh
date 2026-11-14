#!/usr/bin/env bash
#PBS -N vasp
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#PBS -l walltime=24:00:00
#PBS -j oe
##PBS -q workq                       # <- remove one '#' and edit for your cluster
##PBS -A myproject
## Torque instead of PBS Pro: replace the 'select' line by
## #PBS -l nodes=1:ppn=64
#
# job_vasp_pbs.sh - PBS Pro / Torque job template for VASP 6.x (vasp_std / vasp_gam / vasp_ncl)
#
# Usage:
#   qsub job_vasp_pbs.sh                                  # executable chosen automatically
#   qsub -v VARIANT=ncl job_vasp_pbs.sh                   # std | gam | ncl | auto
#   qsub -v NTHREADS=2 -l select=2:ncpus=64:mpiprocs=32:ompthreads=2 job_vasp_pbs.sh
#   qsub -v RESTART=1 job_vasp_pbs.sh                     # continue a relaxation: CONTCAR -> POSCAR
#   bash job_vasp_pbs.sh -h                               # this help
#   (qsub does not pass command-line arguments to the script: use -v VAR=value)
#
# Executable selection (VARIANT=auto):
#   vasp_ncl  if INCAR sets LSORBIT or LNONCOLLINEAR = .TRUE.
#   vasp_gam  if KPOINTS is a Gamma-only 1x1x1 Gamma-centred mesh
#   vasp_std  otherwise
#
# MPI x OpenMP: number of MPI ranks = lines of $PBS_NODEFILE (PBS Pro: 'mpiprocs' per chunk,
#   Torque: 'ppn'); threads per rank = NTHREADS (default 1).  PBS Pro derives OMP_NUM_THREADS
#   from 'ompthreads' (which defaults to ncpus), so this script sets it explicitly from NTHREADS:
#   keep mpiprocs * NTHREADS = ncpus.  With more than one thread VASP uses NCORE = 1 internally
#   (VASP wiki "Combining MPI and OpenMP"); OMP_STACKSIZE=512m, OMP_PLACES=cores,
#   OMP_PROC_BIND=close follow the same page.
#
# What the script does: loads the environment (edit block below), "ulimit -s unlimited",
#   pre-checks the inputs (missing files, CR/TAB in INCAR, POSCAR/POTCAR element order,
#   NCORE/KPAR vs. number of ranks), runs VASP with stdout in vasp.out, then checks
#   "General timing and accounting", "reached required accuracy", SCF steps that hit NELM
#   and known error strings.  If CHECK_VASP points to check_vasp.py it is called as well.
#   Exit code 0 = finished normally, 1 = problem.
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

usage() {
    sed -n '/^# job_vasp_pbs.sh - /,/^# Part of DFT-Tutorials/p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}
case "${1:-}" in -h|--help) usage 0 ;; esac

set -euo pipefail

# ======================= user-editable settings ===========================
VARIANT=${1:-${VARIANT:-auto}}           # auto | std | gam | ncl
NTHREADS=${NTHREADS:-1}                  # OpenMP threads per MPI rank (builds with -D_OPENMP)
VASP_BIN_DIR=${VASP_BIN_DIR:-}           # e.g. $HOME/software/vasp.6.5.1/bin (empty = use PATH)
LAUNCHER=${LAUNCHER:-mpirun}             # mpirun of the MPI library VASP was built with
USE_HOSTFILE=${USE_HOSTFILE:-1}          # 1: pass -machinefile $PBS_NODEFILE (needed without tm support)
MPI_EXTRA=${MPI_EXTRA:-}                 # extra launcher options, e.g. "--report-bindings"
STDOUT_FILE=${STDOUT_FILE:-vasp.out}     # VASP standard output
RESTART=${RESTART:-0}                    # 1: copy a non-empty CONTCAR to POSCAR before the run
CHECK_VASP=${CHECK_VASP:-}               # path to check_vasp.py (optional post-run report)
PYTHON=${PYTHON:-python3}

# ---- environment: keep exactly ONE block that matches your installation ----
# module purge
# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl     # Intel oneAPI build
# module load nvhpc hdf5                                                   # NVIDIA HPC SDK (GPU) build
# export PATH=$HOME/software/vasp.6.5.1/bin:$PATH
# ==========================================================================

die()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }

cd "${PBS_O_WORKDIR:-$PWD}"

# ---- INCAR helpers ---------------------------------------------------------
# value of an INCAR tag (last occurrence, comments removed, 'A = 1; B = 2' lines split)
incar_get() {
    [ -f INCAR ] || return 0
    tr -d '\r' < INCAR | sed 's/[#!].*//' | tr ';' '\n' \
        | awk -v t="$1" '
              { n = index($0, "="); if (n == 0) next
                k = substr($0, 1, n - 1); gsub(/[ \t]/, "", k)
                if (toupper(k) == toupper(t)) { v = substr($0, n + 1); gsub(/^[ \t]+|[ \t]+$/, "", v); last = v } }
              END { if (last != "") print last }'
}
is_true() { case "${1:-}" in .[Tt]*|[Tt]*) return 0 ;; *) return 1 ;; esac; }

# ---- choose the executable --------------------------------------------------
if [ "$VARIANT" = "auto" ]; then
    VARIANT=std
    if is_true "$(incar_get LSORBIT)" || is_true "$(incar_get LNONCOLLINEAR)"; then
        VARIANT=ncl
    elif [ -f KPOINTS ] && [ -z "$(incar_get KSPACING)" ]; then
        if tr -d '\r' < KPOINTS | awk 'NR==2{a=$1} NR==3{b=toupper(substr($1,1,1))} NR==4{m=$1" "$2" "$3}
                                         END{exit !(a=="0" && b=="G" && m=="1 1 1")}'; then
            VARIANT=gam
        fi
    fi
fi
case "$VARIANT" in
    std|gam|ncl) VASP_EXE="vasp_$VARIANT" ;;
    *) die "unknown VARIANT '$VARIANT' (use std, gam, ncl or auto)" ;;
esac
[ -n "$VASP_BIN_DIR" ] && VASP_EXE="$VASP_BIN_DIR/$VASP_EXE"
command -v "$VASP_EXE" >/dev/null 2>&1 || die "$VASP_EXE not found - edit the environment block of $0"

# ---- parallel environment --------------------------------------------------
if [ -n "${PBS_NODEFILE:-}" ] && [ -f "$PBS_NODEFILE" ]; then
    NRANKS=$(wc -l < "$PBS_NODEFILE")
    NNODES=$(sort -u "$PBS_NODEFILE" | wc -l)
else
    NRANKS=${NRANKS:-1}; NNODES=1
    warn "PBS_NODEFILE not set (not running under PBS?) - using NRANKS=$NRANKS"
fi
[[ "$NTHREADS" =~ ^[0-9]+$ ]] && [ "$NTHREADS" -ge 1 ] || die "NTHREADS must be a positive integer"
export OMP_NUM_THREADS=$NTHREADS
export OMP_STACKSIZE=${OMP_STACKSIZE:-512m}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
ulimit -s unlimited 2>/dev/null || warn "could not set 'ulimit -s unlimited' (risk of segmentation faults)"

# ---- input checks ------------------------------------------------------------
if [ "$RESTART" = "1" ] && [ -s CONTCAR ]; then
    n=1; while [ -e "POSCAR.$n" ]; do n=$((n + 1)); done
    cp POSCAR "POSCAR.$n"; cp CONTCAR POSCAR
    echo "RESTART=1: POSCAR saved as POSCAR.$n, CONTCAR copied to POSCAR"
fi
for f in INCAR POSCAR POTCAR; do
    [ -s "$f" ] || die "$f is missing or empty in $(pwd)"
done
if [ ! -s KPOINTS ] && [ -z "$(incar_get KSPACING)" ]; then
    warn "no KPOINTS file and no KSPACING in INCAR (VASP then uses KSPACING = 0.5)"
fi
if grep -q $'\t' INCAR || grep -q $'\r' INCAR; then
    warn "INCAR contains TAB or CR (Windows) characters - run: dos2unix INCAR"
fi
pos_el=$(tr -d '\r' < POSCAR | sed -n '6p' | tr -s ' \t' '\n' | sed '/^$/d' | sed 's#[/_].*##' | tr '\n' ' ')
pot_el=$(grep -E '^[[:space:]]*VRHFIN[[:space:]]*=' POTCAR | sed -E 's/.*=[[:space:]]*([A-Za-z]+).*/\1/' | tr '\n' ' ')
if echo "$pos_el" | grep -qE '^[A-Z][a-z]? '; then
    [ "$pos_el" = "$pot_el" ] || warn "POSCAR species ($pos_el) and POTCAR order ($pot_el) differ"
fi
NCORE=$(incar_get NCORE); KPAR=$(incar_get KPAR)
if [ "$NTHREADS" -gt 1 ] && [ -n "$NCORE" ] && [ "$NCORE" != "1" ]; then
    echo "NOTE: $NTHREADS OpenMP threads per rank -> VASP uses NCORE = 1 (INCAR NCORE = $NCORE ignored)"
fi
if [ "$NTHREADS" -eq 1 ] && [[ "$NCORE" =~ ^[0-9]+$ ]] && [ "$NCORE" -gt 0 ] && [ $((NRANKS % NCORE)) -ne 0 ]; then
    warn "NCORE = $NCORE does not divide the number of MPI ranks ($NRANKS)"
fi
if [[ "$KPAR" =~ ^[0-9]+$ ]] && [ "$KPAR" -gt 0 ] && [ $((NRANKS % KPAR)) -ne 0 ]; then
    warn "KPAR = $KPAR does not divide the number of MPI ranks ($NRANKS)"
fi

# ---- launcher command ------------------------------------------------------
RUN=("$LAUNCHER" -np "$NRANKS")
if [ "$USE_HOSTFILE" = "1" ] && [ -n "${PBS_NODEFILE:-}" ] && [ -f "$PBS_NODEFILE" ]; then
    RUN+=(-machinefile "$PBS_NODEFILE")
fi
if "$LAUNCHER" --version 2>&1 | grep -qi "open mpi"; then
    RUN+=(-x OMP_NUM_THREADS -x OMP_STACKSIZE -x OMP_PLACES -x OMP_PROC_BIND)
    [ "$NTHREADS" -gt 1 ] && RUN+=(--map-by "slot:PE=$NTHREADS" --bind-to core)
elif [ "$NTHREADS" -gt 1 ]; then
    export I_MPI_PIN_DOMAIN=omp            # Intel MPI placement
    export I_MPI_PIN=yes
fi
# shellcheck disable=SC2206
[ -n "$MPI_EXTRA" ] && RUN+=($MPI_EXTRA)

echo "========================================================"
echo "Job ID        : ${PBS_JOBID:-none}   Nodes: $NNODES ($(sort -u "${PBS_NODEFILE:-/dev/null}" 2>/dev/null | tr '\n' ' '))"
echo "Directory     : $(pwd)"
echo "Executable    : $(command -v "$VASP_EXE")"
echo "MPI x OpenMP  : $NRANKS x $NTHREADS   (NCORE=${NCORE:-default} KPAR=${KPAR:-default})"
echo "Command       : ${RUN[*]} $VASP_EXE > $STDOUT_FILE"
echo "Start         : $(date)"
echo "========================================================"

set +e
"${RUN[@]}" "$VASP_EXE" > "$STDOUT_FILE" 2>&1
rc=$?
set -e
echo "End           : $(date)   exit code $rc"

# ---------------------------- post-run check ------------------------------
status=0
if grep -q "General timing and accounting" OUTCAR 2>/dev/null; then
    echo "OK: VASP finished normally"
    grep "Elapsed time (sec)" OUTCAR | tail -n 1 || true
else
    status=1
    echo "PROBLEM: OUTCAR has no 'General timing and accounting' (crash, error or wall-time limit)" >&2
fi
grep "energy(sigma->0)" OUTCAR 2>/dev/null | tail -n 1 || true
NSW=$(incar_get NSW); IBRION=$(incar_get IBRION); NELM=$(incar_get NELM)
NELM=${NELM:-60}
if [ -f OSZICAR ]; then
    nion=$(grep -c " F= " OSZICAR || true)
    # SCF steps per ionic step (lines like 'DAV:', 'RMM:', 'CGA:'); count steps that hit NELM
    nhit=$(awk -v nelm="$NELM" '/^[A-Za-z]{2,4}:/{n++} / F= /{if (n >= nelm) h++; n = 0} END{print h + 0}' OSZICAR)
    echo "ionic steps   : $nion   (SCF reaching NELM = $NELM in $nhit of them)"
    [ "$nhit" -gt 0 ] && warn "$nhit ionic step(s) did not converge electronically (raise NELM, change ALGO/AMIX)"
    if [[ "${NSW:-0}" =~ ^[0-9]+$ ]] && [ "${NSW:-0}" -gt 0 ] && [[ "${IBRION:-0}" =~ ^[123]$ ]]; then
        if grep -q "reached required accuracy" OUTCAR 2>/dev/null; then
            echo "OK: relaxation reached required accuracy"
        else
            status=1
            echo "NOTE: relaxation NOT converged ($nion / NSW = $NSW steps) - resubmit with qsub -v RESTART=1" >&2
        fi
    fi
fi
pat='VERY BAD NEWS|ZBRENT|BRMIX|EDDDAV|EDDRMM|ZPOTRF|Sub-Space-Matrix|internal error|BZINTS|Tetrahedron method fails|Error reading item|BAD TERMINATION|forrtl: severe|Segmentation fault|oom-kill|Out Of Memory'
if grep -E -q "$pat" "$STDOUT_FILE" 2>/dev/null; then
    echo "Messages found in $STDOUT_FILE:" >&2
    grep -E -h "$pat" "$STDOUT_FILE" 2>/dev/null | sort | uniq -c | head -n 10 >&2 || true
fi
if [ -n "$CHECK_VASP" ] && [ -f "$CHECK_VASP" ]; then
    "$PYTHON" "$CHECK_VASP" . -o "$STDOUT_FILE" || true
fi
[ "$rc" -ne 0 ] && status=1
exit "$status"
