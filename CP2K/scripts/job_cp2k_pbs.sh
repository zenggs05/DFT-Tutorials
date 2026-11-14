#!/usr/bin/env bash
#PBS -N cp2k
#PBS -l select=1:ncpus=32:mpiprocs=32:ompthreads=1
#PBS -l walltime=24:00:00
#PBS -j oe
##PBS -q normal                        # <- remove one '#' and edit for your cluster
##PBS -A myproject
## Torque / OpenPBS without 'select': use instead
##PBS -l nodes=1:ppn=32
#
# job_cp2k_pbs.sh - PBS Pro / OpenPBS / Torque job template for CP2K (cp2k.psmp, MPI x OpenMP)
#
# Usage:
#   qsub -v INP=my.inp job_cp2k_pbs.sh                   # PBS does not pass script arguments
#   qsub -v INP=my.inp,OUT=my.out,NTHREADS=2 job_cp2k_pbs.sh
#   bash job_cp2k_pbs.sh -h                              # this help
#
# MPI x OpenMP: with PBS Pro "select=N:ncpus=C:mpiprocs=M:ompthreads=T" use C = M * T;
# the script takes the number of MPI ranks from $PBS_NODEFILE (one line per rank) and the
# threads from NTHREADS (default: $OMP_NUM_THREADS set by ompthreads, else 1).
# Torque: "-l nodes=1:ppn=32" gives 32 lines (one per core) in $PBS_NODEFILE; submit with
# -v NODEFILE_LINES=cores,NTHREADS=T to start 32/T ranks.
# NEB: total ranks must be NUMBER_OF_REPLICA * NPROC_REP.
#
# What the script does: loads the CP2K environment (edit block below), sets OpenMP variables
# and "ulimit -s unlimited", keeps an existing output as <out>.<n>.bak, runs
# "mpirun -np <ranks> cp2k.psmp -i INP -o OUT", then checks for "PROGRAM ENDED AT" / ABORT /
# "SCF run NOT converged" and prints how to restart from <PROJECT>-1.restart.
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

usage() {
    sed -n '/^# job_cp2k_pbs.sh - /,/^# Part of DFT-Tutorials/p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}
case "${1:-}" in -h|--help) usage 0 ;; esac

set -euo pipefail

# ======================= user-editable settings ===========================
INP=${1:-${INP:-input.inp}}              # CP2K input file
OUT=${OUT:-${INP%.*}.out}                # CP2K output file
CP2K_EXE=${CP2K_EXE:-cp2k.psmp}          # cp2k.psmp (MPI+OpenMP); cp2k.ssmp = OpenMP only
NTHREADS=${NTHREADS:-${OMP_NUM_THREADS:-1}}   # OpenMP threads per MPI rank
MPIRUN=${MPIRUN:-mpirun}                 # MPI launcher (mpirun / mpiexec)
MPI_EXTRA=${MPI_EXTRA:-}                 # extra launcher options, e.g. "--bind-to core --map-by slot:PE=2"
NODEFILE_LINES=${NODEFILE_LINES:-ranks}  # "ranks" (PBS Pro select:mpiprocs) or "cores" (Torque nodes:ppn)

# ---- environment: keep exactly ONE block that matches your installation ----
# source "$HOME/cp2k-2026.2/install/cp2k_env"     # toolchain build (install_cp2k_toolchain.sh)
# module load cp2k/2026.2                         # cluster module
# source "$HOME/miniconda3/bin/activate" cp2k     # conda-forge package
# export CP2K_DATA_DIR=/path/to/cp2k/data         # only if your build does not set it
# ==========================================================================

cd "${PBS_O_WORKDIR:-$PWD}"

[ -f "$INP" ] || { echo "ERROR: input file '$INP' not found in $(pwd)" >&2; exit 1; }
command -v "$CP2K_EXE" >/dev/null 2>&1 || {
    echo "ERROR: $CP2K_EXE not in PATH - edit the environment block of $0" >&2; exit 1; }

if [ -n "${PBS_NODEFILE:-}" ] && [ -f "$PBS_NODEFILE" ]; then
    NSLOTS=$(wc -l < "$PBS_NODEFILE")
    MACHINEFILE_OPT=(-machinefile "$PBS_NODEFILE")
else
    NSLOTS=${NSLOTS:-1}
    MACHINEFILE_OPT=()
fi
# PBS Pro "mpiprocs=M": the node file has one line per MPI rank   -> NODEFILE_LINES=ranks
# Torque "nodes=N:ppn=P": the node file has one line per core      -> NODEFILE_LINES=cores
if [ "$NODEFILE_LINES" = "cores" ]; then
    NP=$(( NSLOTS / NTHREADS ))
else
    NP=$NSLOTS
fi
[ "$NP" -ge 1 ] || NP=1

export OMP_NUM_THREADS=$NTHREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=close
ulimit -s unlimited 2>/dev/null || echo "WARNING: could not set 'ulimit -s unlimited'" >&2

if [ -f "$OUT" ]; then
    n=1; while [ -e "$OUT.$n.bak" ]; do n=$((n + 1)); done
    mv "$OUT" "$OUT.$n.bak"
    echo "previous output moved to $OUT.$n.bak"
fi

PROJECT=$(grep -iE '^[[:space:]]*PROJECT(_NAME)?[[:space:]]' "$INP" | head -n 1 | awk '{print $2}' || true)
echo "========================================================"
echo "Job ID        : ${PBS_JOBID:-none}   Host: $(hostname)"
echo "CP2K binary   : $(command -v "$CP2K_EXE")"
echo "Input/Output  : $INP -> $OUT   (PROJECT ${PROJECT:-?})"
echo "MPI x OpenMP  : $NP x $OMP_NUM_THREADS   (node-file lines: $NSLOTS)"
echo "Start         : $(date)"
echo "========================================================"

set +e
# shellcheck disable=SC2086
"$MPIRUN" -np "$NP" ${MACHINEFILE_OPT[@]+"${MACHINEFILE_OPT[@]}"} $MPI_EXTRA "$CP2K_EXE" -i "$INP" -o "$OUT"
rc=$?
set -e
echo "End           : $(date)   exit code $rc"

# ---------------------------- post-run check ------------------------------
status=0
if [ -f "$OUT" ] && grep -q "PROGRAM ENDED AT" "$OUT"; then
    echo "OK: CP2K finished normally ($OUT)"
    grep "ENERGY| Total FORCE_EVAL" "$OUT" | tail -n 1 || true
    grep -q "GEOMETRY OPTIMIZATION COMPLETED" "$OUT" && echo "Optimization converged."
    grep -q "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED" "$OUT" && \
        echo "NOTE: optimization stopped at MAX_ITER - not converged."
    nwarn=$(grep -c "\*\*\* WARNING in" "$OUT" || true)
    [ "${nwarn:-0}" -gt 0 ] && echo "NOTE: $nwarn warning(s); list them with: python cp2k_outparse.py $OUT"
else
    status=1
    echo "PROBLEM: no 'PROGRAM ENDED AT' in $OUT (crash, ABORT or walltime limit)." >&2
    if [ -f "$OUT" ]; then
        grep -A8 "\[ABORT\]" "$OUT" | head -n 10 >&2 || true
        grep -q "SCF run NOT converged" "$OUT" && \
            echo "-> SCF not converged: see the SCF chapter (OT/OUTER_SCF, mixing, smearing, restart .wfn)." >&2
    fi
    nstart=$(grep -c "PROGRAM STARTED AT" "$OUT" 2>/dev/null || true)
    [ "${nstart:-0}" -gt 1 ] && echo "-> 'PROGRAM STARTED AT' appears $nstart times: MPI mismatch (several serial runs)." >&2
    if [ -n "$PROJECT" ] && [ -f "${PROJECT}-1.restart" ]; then
        echo "-> To continue, add to a copy of $INP:" >&2
        echo "     &EXT_RESTART" >&2
        echo "       RESTART_FILE_NAME ${PROJECT}-1.restart" >&2
        echo "     &END EXT_RESTART" >&2
        echo "   or write the last structure into a new input: python restart_update.py ${PROJECT}-1.restart -i $INP" >&2
        ls "${PROJECT}"-RESTART.wfn "${PROJECT}"-RESTART.kp 2>/dev/null | sed 's/^/   wavefunction available: /' >&2 || true
    fi
fi
[ "$rc" -ne 0 ] && status=1
exit "$status"
