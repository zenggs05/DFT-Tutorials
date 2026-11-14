#!/usr/bin/env bash
#SBATCH --job-name=openmx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16          # MPI processes per node
#SBATCH --cpus-per-task=2             # OpenMP threads per MPI process (= openmx -nt)
#SBATCH --time=24:00:00
#SBATCH --partition=cpu               # EDIT: site specific
#SBATCH --output=%x-%j.log
##SBATCH --account=myproject          # EDIT: uncomment if your site needs it
##SBATCH --mem=0                      # all memory of the node
#
# job_omx_slurm.sh -- Slurm job template for OpenMX 4.0 (MPI x OpenMP hybrid).
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
#
# Usage
#   sbatch job_omx_slurm.sh input.dat
#   sbatch --nodes=2 --ntasks-per-node=32 --cpus-per-task=2 job_omx_slurm.sh si64.dat
#   INPUT=fe.dat OMX_DATA=/opt/openmx4.0/DFT_DATA19 sbatch job_omx_slurm.sh
#   NP=4 NT=1 bash job_omx_slurm.sh Methane.dat      # also works without Slurm (workstation)
#   bash job_omx_slurm.sh -h
#
# Rules followed (OpenMX manual Sec. 28)
#   * total cores = MPI processes x threads; OMP_NUM_THREADS is kept equal to -nt
#   * ulimit -s unlimited and OMP_STACKSIZE avoid segmentation faults of large runs
#   * success = the line "The calculation was normally finished." in System.Name.std
#     (the exit code of mpirun is not reliable); in addition the SCF step count of
#     the last block is compared with scf.maxIter
#   * System.Name_rst/ (restart files) is kept unless CLEAN_RST=1
#
# Dependencies: bash, awk, an OpenMX executable and its MPI runtime.
#
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: sbatch job_omx_slurm.sh [input.dat]        (or: NP=4 NT=1 bash job_omx_slurm.sh input.dat)

Slurm resources: edit the #SBATCH lines (or override them on the sbatch command line);
MPI processes = SLURM_NTASKS, threads = SLURM_CPUS_PER_TASK (passed as -nt).

Environment variables (all optional):
  INPUT       input file if no argument is given           (default input.dat)
  OPENMX      OpenMX executable                             (default openmx)
  OMX_DATA    rewrite DATA.PATH of the input to this dir   (default: keep the input)
  LAUNCH      mpirun | srun                                 (default mpirun)
  MPI_EXTRA   extra launcher options                        (default none)
  CLEAN_RST   1 = remove System.Name_rst/ after success     (default 0 = keep)
  USE_DEB     1 = add /opt/openmx4.0/runtime/lib (deb package) to LD_LIBRARY_PATH
  NP, NT      MPI processes / threads outside Slurm         (default 4 / 1)

Exit status: 0 = finished and last SCF converged, 1 = run failed, 2 = SCF hit scf.maxIter.
EOF
}
case "${1:-}" in -h|--help) usage; exit 0 ;; esac

# ============================ user-editable variables ============================
INPUT=${1:-${INPUT:-input.dat}}        # OpenMX input file
OPENMX=${OPENMX:-openmx}               # executable (full path if not in PATH)
OMX_DATA=${OMX_DATA:-}                 # if set: DATA.PATH of the input is rewritten to this directory
LAUNCH=${LAUNCH:-mpirun}               # mpirun | srun (srun needs MPI built with Slurm PMI support)
MPI_EXTRA=${MPI_EXTRA:-}               # extra launcher options, e.g. "--bind-to core"
CLEAN_RST=${CLEAN_RST:-0}              # 1 = delete System.Name_rst/ after a successful run
USE_DEB=${USE_DEB:-0}                  # 1 = runtime libraries of the official .deb in /opt/openmx4.0
# module purge                                                             # EDIT: environment
# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl     # Intel toolchain
# module load gcc openmpi openblas scalapack fftw                          # GNU toolchain
# source /opt/intel/oneapi/setvars.sh
# ==================================================================================

if [ "$USE_DEB" = 1 ]; then
  # fixes "libmpifort.so.12: cannot open shared object file" of the .deb tools
  export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:${LD_LIBRARY_PATH:-}
fi
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then cd "$SLURM_SUBMIT_DIR"; fi

NP=${SLURM_NTASKS:-${NP:-4}}           # MPI processes
NT=${SLURM_CPUS_PER_TASK:-${NT:-1}}    # OpenMP threads per process
export OMP_NUM_THREADS=$NT             # must equal -nt
export OMP_STACKSIZE=${OMP_STACKSIZE:-512M}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MKL_NUM_THREADS=$NT
export OPENBLAS_NUM_THREADS=1          # avoid nested threading inside OpenBLAS
export I_MPI_PIN_DOMAIN=omp            # Intel MPI: one pinning domain per process
ulimit -s unlimited 2>/dev/null || echo "warning: could not set 'ulimit -s unlimited'" >&2

# ---------------------------------------------------------------- checks
[ -f "$INPUT" ] || { echo "ERROR: input '$INPUT' not found (usage: sbatch $0 input.dat)" >&2; exit 1; }
if ! command -v "$OPENMX" >/dev/null 2>&1 && [ ! -x "$OPENMX" ]; then
  echo "ERROR: OpenMX executable '$OPENMX' not found; set OPENMX=/path/to/openmx or load the module" >&2
  exit 1
fi
if grep -q $'\r' "$INPUT"; then
  echo "ERROR: $INPUT has Windows (CRLF) line endings; run: dos2unix $INPUT" >&2
  exit 1
fi
if [ -n "$OMX_DATA" ]; then
  sed -i "s#^[[:space:]]*[Dd][Aa][Tt][Aa]\.[Pp][Aa][Tt][Hh][[:space:]].*#DATA.PATH   ${OMX_DATA}#" "$INPUT"
fi
NAME=$(awk 'tolower($1)=="system.name"{print $2; exit}' "$INPUT")
[ -n "$NAME" ] || { echo "ERROR: System.Name not found in $INPUT" >&2; exit 1; }
DPATH=$(awk 'tolower($1)=="data.path"{print $2; exit}' "$INPUT")
DPATH=${DPATH:-../DFT_DATA19}
if [ ! -d "$DPATH/PAO" ] || [ ! -d "$DPATH/VPS" ]; then
  echo "ERROR: DATA.PATH '$DPATH' does not contain PAO/ and VPS/ (set it to .../DFT_DATA19 or use OMX_DATA=...)" >&2
  exit 1
fi
MAXITER=$(awk 'tolower($1)=="scf.maxiter"{print $2; exit}' "$INPUT")
MAXITER=${MAXITER:-40}
STD=${NAME}.std
OUT=${NAME}.out

echo "host $(hostname) | input $INPUT | System.Name $NAME | NP=$NP x NT=$NT | start $(date '+%F %T')"

# ---------------------------------------------------------------- run
set +e
if [ "$LAUNCH" = srun ] && [ -n "${SLURM_JOB_ID:-}" ]; then
  export SRUN_CPUS_PER_TASK=$NT        # Slurm >= 22.05 does not pass --cpus-per-task to srun
  # shellcheck disable=SC2086
  srun --cpus-per-task="$NT" --cpu-bind=cores $MPI_EXTRA "$OPENMX" "$INPUT" -nt "$NT" 2>&1 | tee "$STD"
else
  # shellcheck disable=SC2086
  mpirun -np "$NP" $MPI_EXTRA "$OPENMX" "$INPUT" -nt "$NT" 2>&1 | tee "$STD"
fi
RC=${PIPESTATUS[0]}
set -e

# ---------------------------------------------------------------- post-run check
if ! grep -q "The calculation was normally finished." "$STD"; then
  echo "ERROR: OpenMX did NOT finish normally (launcher exit code $RC); last lines of $STD:" >&2
  tail -n 30 "$STD" >&2
  exit 1
fi
STATUS=0
if [ -f "$OUT" ]; then
  NSCF=$(awk '/SCF history at MD=/{n=0} /^ *SCF= *[0-9]+ +NormRD=/{n++} END{print n+0}' "$OUT")
  UTOT=$(awk '/^ *Utot\./{v=$2} END{print v}' "$OUT")
  TTOT=$(awk '/Elapsed\.Time\./{v=$2} END{print v}' "$OUT")
  echo "finished normally $(date '+%F %T'): Utot = ${UTOT:-n/a} Ha, last SCF block $NSCF steps (scf.maxIter $MAXITER), elapsed ${TTOT:-n/a} s"
  if [ "$NSCF" -ge "$MAXITER" ]; then
    echo "WARNING: the last SCF block used scf.maxIter steps -> probably NOT converged" >&2
    STATUS=2
  fi
else
  echo "WARNING: $OUT not found" >&2
fi
if [ "$CLEAN_RST" = 1 ] && [ "$STATUS" = 0 ]; then
  rm -rf "${NAME}_rst"
fi
exit "$STATUS"
