#!/usr/bin/env bash
#SBATCH --job-name=cp2k
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32          # MPI ranks per node
#SBATCH --cpus-per-task=1             # OpenMP threads per MPI rank (cp2k.psmp)
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal           # <- remove one '#' and edit for your cluster
##SBATCH --account=myproject
##SBATCH --mem=0                      # all memory of the node (hybrid functionals!)
#
# job_cp2k_slurm.sh - Slurm job template for CP2K (cp2k.psmp, MPI x OpenMP)
#
# Usage:
#   sbatch job_cp2k_slurm.sh my.inp                      # input as argument
#   sbatch --export=ALL,INP=my.inp job_cp2k_slurm.sh     # or as variable (OUT=... optional)
#   sbatch -N 2 --ntasks-per-node=16 --cpus-per-task=4 job_cp2k_slurm.sh my.inp
#   bash job_cp2k_slurm.sh -h                            # this help
#
# MPI x OpenMP: total cores = nodes * ntasks-per-node * cpus-per-task.  Start with
# 1-2 threads per rank; hybrid functionals / large memory often run better with fewer ranks
# and more threads.  NEB: the total number of ranks must be NUMBER_OF_REPLICA * NPROC_REP.
#
# What the script does: loads the CP2K environment (edit block below), sets OpenMP variables
# and "ulimit -s unlimited", keeps an existing output as <out>.<n>.bak, runs
# "srun cp2k.psmp -i INP -o OUT", then checks for "PROGRAM ENDED AT" / ABORT /
# "SCF run NOT converged" and prints how to restart from <PROJECT>-1.restart.
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License

usage() {
    sed -n '/^# job_cp2k_slurm.sh - /,/^# Part of DFT-Tutorials/p' "$0" | sed 's/^# \{0,1\}//'
    exit "${1:-0}"
}
case "${1:-}" in -h|--help) usage 0 ;; esac

set -euo pipefail

# ======================= user-editable settings ===========================
INP=${1:-${INP:-input.inp}}              # CP2K input file
OUT=${OUT:-${INP%.*}.out}                # CP2K output file
CP2K_EXE=${CP2K_EXE:-cp2k.psmp}          # cp2k.psmp (MPI+OpenMP); cp2k.ssmp = OpenMP only
LAUNCHER=${LAUNCHER:-srun}               # srun (recommended) or "mpirun"

# ---- environment: keep exactly ONE block that matches your installation ----
# source "$HOME/cp2k-2026.2/install/cp2k_env"     # toolchain build (install_cp2k_toolchain.sh)
# module load cp2k/2026.2                         # cluster module
# source "$HOME/miniconda3/bin/activate" cp2k     # conda-forge package (conda install -c conda-forge cp2k)
# export CP2K_DATA_DIR=/path/to/cp2k/data         # only if your build does not set it
# ==========================================================================

cd "${SLURM_SUBMIT_DIR:-$PWD}"

[ -f "$INP" ] || { echo "ERROR: input file '$INP' not found in $(pwd)" >&2; exit 1; }
command -v "$CP2K_EXE" >/dev/null 2>&1 || {
    echo "ERROR: $CP2K_EXE not in PATH - edit the environment block of $0" >&2; exit 1; }

NTHREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=$NTHREADS
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export SRUN_CPUS_PER_TASK=$NTHREADS     # Slurm >= 22.05: srun does not inherit --cpus-per-task
ulimit -s unlimited 2>/dev/null || echo "WARNING: could not set 'ulimit -s unlimited'" >&2

# keep an existing output instead of overwriting it
if [ -f "$OUT" ]; then
    n=1; while [ -e "$OUT.$n.bak" ]; do n=$((n + 1)); done
    mv "$OUT" "$OUT.$n.bak"
    echo "previous output moved to $OUT.$n.bak"
fi

PROJECT=$(grep -iE '^[[:space:]]*PROJECT(_NAME)?[[:space:]]' "$INP" | head -n 1 | awk '{print $2}' || true)
echo "========================================================"
echo "Job ID        : ${SLURM_JOB_ID:-none}   Nodes: ${SLURM_JOB_NODELIST:-$(hostname)}"
echo "CP2K binary   : $(command -v "$CP2K_EXE")"
echo "Input/Output  : $INP -> $OUT   (PROJECT ${PROJECT:-?})"
echo "MPI x OpenMP  : ${SLURM_NTASKS:-1} x $OMP_NUM_THREADS"
echo "Start         : $(date)"
echo "========================================================"

set +e
if [ "$LAUNCHER" = "srun" ]; then
    srun --cpu-bind=cores "$CP2K_EXE" -i "$INP" -o "$OUT"
else
    $LAUNCHER -np "${SLURM_NTASKS:-1}" "$CP2K_EXE" -i "$INP" -o "$OUT"
fi
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
    echo "PROBLEM: no 'PROGRAM ENDED AT' in $OUT (crash, ABORT or wall-time limit)." >&2
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
