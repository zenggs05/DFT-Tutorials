#!/usr/bin/env bash
#SBATCH --job-name=vasp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64          # MPI ranks per node
#SBATCH --cpus-per-task=1             # OpenMP threads per rank (only for builds with -D_OPENMP)
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=cpu              # <- remove one '#' and edit for your cluster
##SBATCH --account=myproject
##SBATCH --mem=0                      # request all memory of the node
##SBATCH --gres=gpu:4                 # GPU (OpenACC) build: 1 MPI rank per GPU, e.g. ntasks-per-node=4
#
# job_vasp_slurm.sh - Slurm job template for VASP 6.x (vasp_std / vasp_gam / vasp_ncl)
#
# Usage:
#   sbatch job_vasp_slurm.sh                  # executable chosen automatically (VARIANT=auto)
#   sbatch job_vasp_slurm.sh ncl              # force vasp_ncl (std | gam | ncl | auto)
#   sbatch --export=ALL,VARIANT=gam,LAUNCHER=srun job_vasp_slurm.sh
#   sbatch -N 2 --ntasks-per-node=32 --cpus-per-task=2 job_vasp_slurm.sh    # MPI x OpenMP
#   sbatch --export=ALL,RESTART=1 job_vasp_slurm.sh    # continue a relaxation: CONTCAR -> POSCAR
#   bash job_vasp_slurm.sh -h                 # this help
#   (conv_test.sh / phonopy_workflow.sh copy this file into every directory as job.sh)
#
# Executable selection (VARIANT=auto):
#   vasp_ncl  if INCAR sets LSORBIT or LNONCOLLINEAR = .TRUE.
#   vasp_gam  if KPOINTS is a Gamma-only 1x1x1 Gamma-centred mesh (molecules, large cells)
#   vasp_std  otherwise
#
# MPI x OpenMP: total cores = nodes * ntasks-per-node * cpus-per-task (physical cores, no
#   oversubscription).  With more than one OpenMP thread VASP sets NCORE = 1 internally
#   (VASP wiki "Combining MPI and OpenMP"); the threads take over the role of NCORE.
#   OMP_STACKSIZE=512m, OMP_PLACES=cores, OMP_PROC_BIND=close follow the same wiki page.
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
    sed -n '/^# job_vasp_slurm.sh - /,/^# Part of DFT-Tutorials/p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}
case "${1:-}" in -h|--help) usage 0 ;; esac

set -euo pipefail

# ======================= user-editable settings ===========================
VARIANT=${1:-${VARIANT:-auto}}           # auto | std | gam | ncl
VASP_BIN_DIR=${VASP_BIN_DIR:-}           # e.g. $HOME/software/vasp.6.5.1/bin (empty = use PATH)
LAUNCHER=${LAUNCHER:-mpirun}             # mpirun (Intel MPI / Open MPI) or srun
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

cd "${SLURM_SUBMIT_DIR:-$PWD}"

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
        # Gamma-only: line 2 = 0 (automatic), line 3 starts with G, mesh 1 1 1
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
NRANKS=${SLURM_NTASKS:-1}
NTHREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$NTHREADS
export OMP_STACKSIZE=${OMP_STACKSIZE:-512m}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export SRUN_CPUS_PER_TASK=$NTHREADS            # Slurm >= 22.05: srun does not inherit --cpus-per-task
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
# element order POSCAR (line 6, VASP 5 format) vs POTCAR (VRHFIN = El: ...)
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
if [ "$LAUNCHER" = "srun" ]; then
    RUN=(srun --cpu-bind=cores)
else
    RUN=("$LAUNCHER" -np "$NRANKS")
    if [ "$NTHREADS" -gt 1 ]; then
        if "$LAUNCHER" --version 2>&1 | grep -qi "open mpi"; then
            RUN+=(--map-by "slot:PE=$NTHREADS" --bind-to core)     # Open MPI placement
        else
            export I_MPI_PIN_DOMAIN=omp                             # Intel MPI placement
            export I_MPI_PIN=yes
        fi
    fi
fi
# shellcheck disable=SC2206
[ -n "$MPI_EXTRA" ] && RUN+=($MPI_EXTRA)

echo "========================================================"
echo "Job ID        : ${SLURM_JOB_ID:-none}   Nodes: ${SLURM_JOB_NODELIST:-$(hostname)}"
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
            echo "NOTE: relaxation NOT converged ($nion / NSW = $NSW steps) - continue with RESTART=1" >&2
        fi
    fi
fi
pat='VERY BAD NEWS|ZBRENT|BRMIX|EDDDAV|EDDRMM|ZPOTRF|Sub-Space-Matrix|internal error|BZINTS|Tetrahedron method fails|Error reading item|BAD TERMINATION|forrtl: severe|Segmentation fault|oom-kill|Out Of Memory|DUE TO TIME LIMIT'
if grep -E -q "$pat" "$STDOUT_FILE" "slurm-${SLURM_JOB_ID:-x}.out" 2>/dev/null; then
    echo "Messages found in the output:" >&2
    grep -E -h "$pat" "$STDOUT_FILE" "slurm-${SLURM_JOB_ID:-x}.out" 2>/dev/null | sort | uniq -c | head -n 10 >&2 || true
fi
if [ -n "$CHECK_VASP" ] && [ -f "$CHECK_VASP" ]; then
    "$PYTHON" "$CHECK_VASP" . -o "$STDOUT_FILE" || true
fi
[ "$rc" -ne 0 ] && status=1
exit "$status"
