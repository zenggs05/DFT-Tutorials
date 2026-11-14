#!/usr/bin/env bash
#PBS -N openmx
#PBS -l select=1:ncpus=32:mpiprocs=16:ompthreads=2
#PBS -l walltime=24:00:00
#PBS -q workq
#PBS -j oe
##PBS -A myproject
#
# job_omx_pbs.sh -- PBS Pro / OpenPBS / Torque job template for OpenMX 4.0 (MPI x OpenMP).
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
#
# Usage
#   qsub -v INPUT=input.dat job_omx_pbs.sh
#   qsub -v INPUT=si64.dat,OMX_DATA=/opt/openmx4.0/DFT_DATA19 job_omx_pbs.sh
#   qsub -l select=2:ncpus=64:mpiprocs=32:ompthreads=2 -v INPUT=si64.dat job_omx_pbs.sh
#   bash job_omx_pbs.sh -h
#
# Resource line
#   PBS Pro/OpenPBS: select=<nodes>:ncpus=<cores per node>:mpiprocs=<MPI per node>:ompthreads=<threads>
#                    $PBS_NODEFILE has one line per MPI process, OMP_NUM_THREADS is set from ompthreads
#   Torque:          use "#PBS -l nodes=1:ppn=32" instead, set PBS_FLAVOR=torque and NT below;
#                    $PBS_NODEFILE has one line per core, the script builds a machine file with
#                    (cores per node / NT) slots per node
#
# Success check: "The calculation was normally finished." in System.Name.std plus the SCF step
# count of the last block compared with scf.maxIter.  System.Name_rst/ is kept unless CLEAN_RST=1.
#
# Dependencies: bash, awk, an OpenMX executable and its MPI runtime.
#
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: qsub -v INPUT=input.dat job_omx_pbs.sh     (or: INPUT=x.dat NP=4 NT=1 bash job_omx_pbs.sh)

Environment variables (pass with qsub -v VAR=value,...; all optional):
  INPUT       OpenMX input file                              (default input.dat)
  OPENMX      OpenMX executable                              (default openmx)
  OMX_DATA    rewrite DATA.PATH of the input to this dir    (default: keep the input)
  PBS_FLAVOR  pro | torque  (layout of $PBS_NODEFILE)        (default pro)
  NT          OpenMP threads per MPI process                 (default: OMP_NUM_THREADS or 1)
  NP          MPI processes (default: from $PBS_NODEFILE; 4 outside PBS)
  MPI_EXTRA   extra mpirun options                           (default none)
  CLEAN_RST   1 = remove System.Name_rst/ after success      (default 0 = keep)
  USE_DEB     1 = add /opt/openmx4.0/runtime/lib (deb package) to LD_LIBRARY_PATH

Exit status: 0 = finished and last SCF converged, 1 = run failed, 2 = SCF hit scf.maxIter.
EOF
}
case "${1:-}" in -h|--help) usage; exit 0 ;; esac

# ============================ user-editable variables ============================
INPUT=${1:-${INPUT:-input.dat}}
OPENMX=${OPENMX:-openmx}
OMX_DATA=${OMX_DATA:-}
PBS_FLAVOR=${PBS_FLAVOR:-pro}          # pro (PBS Pro/OpenPBS select=...) | torque (nodes=:ppn=)
NT=${NT:-${OMP_NUM_THREADS:-1}}        # threads per MPI process (= openmx -nt)
MPI_EXTRA=${MPI_EXTRA:-}
CLEAN_RST=${CLEAN_RST:-0}
USE_DEB=${USE_DEB:-0}
# module purge                                                             # EDIT: environment
# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl
# module load gcc openmpi openblas scalapack fftw
# ==================================================================================

if [ "$USE_DEB" = 1 ]; then
  export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:${LD_LIBRARY_PATH:-}
fi
if [ -n "${PBS_O_WORKDIR:-}" ]; then cd "$PBS_O_WORKDIR"; fi

# ---------------------------------------------------------------- MPI layout
MACHINEFILE=""
if [ -n "${PBS_NODEFILE:-}" ] && [ -f "$PBS_NODEFILE" ]; then
  MACHINEFILE=$(mktemp "${TMPDIR:-/tmp}/omx_hosts.XXXXXX")
  trap 'rm -f "$MACHINEFILE"' EXIT
  if [ "$PBS_FLAVOR" = torque ]; then
    # one line per core -> (cores / NT) MPI slots per host
    awk -v nt="$NT" '{c[$1]++} END{for (h in c) for (i = 0; i < int(c[h] / nt); i++) print h}' \
        "$PBS_NODEFILE" | sort > "$MACHINEFILE"
  else
    cp "$PBS_NODEFILE" "$MACHINEFILE"   # PBS Pro: one line per MPI process (mpiprocs)
  fi
  NP=${NP:-$(wc -l < "$MACHINEFILE" | tr -d " ")}
else
  NP=${NP:-4}
fi
[ "$NP" -ge 1 ] || { echo "ERROR: no MPI slots (check select/mpiprocs or NT)" >&2; exit 1; }

export OMP_NUM_THREADS=$NT
export OMP_STACKSIZE=${OMP_STACKSIZE:-512M}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export MKL_NUM_THREADS=$NT
export OPENBLAS_NUM_THREADS=1
export I_MPI_PIN_DOMAIN=omp
ulimit -s unlimited 2>/dev/null || echo "warning: could not set 'ulimit -s unlimited'" >&2

# ---------------------------------------------------------------- checks
[ -f "$INPUT" ] || { echo "ERROR: input '$INPUT' not found (qsub -v INPUT=file.dat ...)" >&2; exit 1; }
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
  echo "ERROR: DATA.PATH '$DPATH' does not contain PAO/ and VPS/" >&2
  exit 1
fi
MAXITER=$(awk 'tolower($1)=="scf.maxiter"{print $2; exit}' "$INPUT")
MAXITER=${MAXITER:-40}
STD=${NAME}.std
OUT=${NAME}.out

echo "host $(hostname) | input $INPUT | System.Name $NAME | NP=$NP x NT=$NT | start $(date '+%F %T')"

# ---------------------------------------------------------------- run
set +e
if [ -n "$MACHINEFILE" ]; then
  # shellcheck disable=SC2086
  mpirun -np "$NP" -machinefile "$MACHINEFILE" $MPI_EXTRA "$OPENMX" "$INPUT" -nt "$NT" 2>&1 | tee "$STD"
else
  # shellcheck disable=SC2086
  mpirun -np "$NP" $MPI_EXTRA "$OPENMX" "$INPUT" -nt "$NT" 2>&1 | tee "$STD"
fi
RC=${PIPESTATUS[0]}
set -e

# ---------------------------------------------------------------- post-run check
if ! grep -q "The calculation was normally finished." "$STD"; then
  echo "ERROR: OpenMX did NOT finish normally (mpirun exit code $RC); last lines of $STD:" >&2
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
