#!/usr/bin/env bash
# check_install.sh - sanity check of a CP2K installation.
#
# Usage:
#     bash check_install.sh [cp2k_executable]
#
# Environment variables (all optional):
#     CP2K_EXE         executable to test (default: first of cp2k.psmp cp2k.ssmp cp2k in PATH)
#     NP               number of MPI ranks for the smoke test (default 1)
#     MPI_LAUNCHER     launcher used when NP > 1 (default: mpirun; e.g. "srun" or "mpiexec")
#     OMP_NUM_THREADS  OpenMP threads per rank (default 1)
#     CP2K_DATA_DIR    data directory; checked if set
#
# What is checked:
#     1. the executable exists and prints its version (cp2k -v)
#     2. important cp2kflags (omp, parallel, libint, libxc, fftw3, spglib)
#     3. key data files in CP2K_DATA_DIR (if set)
#     4. syntax check of h2o_smoke_test.inp (cp2k -c)
#     5. smoke test run: SCF convergence, total energy, normal termination
#     6. MPI sanity: DBCSR must report exactly NP processes, printed once

set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INP="${HERE}/h2o_smoke_test.inp"
NP=${NP:-1}
MPI_LAUNCHER=${MPI_LAUNCHER:-mpirun}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

nfail=0
ok()   { echo "[ OK ] $*"; }
info() { echo "[INFO] $*"; }
warn() { echo "[WARN] $*"; }
fail() { echo "[FAIL] $*"; nfail=$((nfail + 1)); }

# ---------------------------------------------------------------- 1. binary
EXE=${1:-${CP2K_EXE:-}}
if [[ -z "${EXE}" ]]; then
  for cand in cp2k.psmp cp2k.ssmp cp2k; do
    if command -v "${cand}" > /dev/null 2>&1; then EXE=${cand}; break; fi
  done
fi
if [[ -z "${EXE}" ]] || ! command -v "${EXE}" > /dev/null 2>&1; then
  fail "no CP2K executable found (set CP2K_EXE or source env_cp2k.sh)"
  exit 1
fi
ok "executable : $(command -v "${EXE}")"

VERSION_TXT=$("${EXE}" -v 2>&1)
VLINE=$(echo "${VERSION_TXT}" | grep -m1 "CP2K version")
FLAGS=$(echo "${VERSION_TXT}" | grep -m1 "cp2kflags" | sed 's/.*cp2kflags: *//')
if [[ -n "${VLINE}" ]]; then
  ok "version    : $(echo "${VLINE}" | sed 's/^ *//')"
else
  fail "'${EXE} -v' did not print a version (output below)"
  echo "${VERSION_TXT}" | head -20
fi
info "cp2kflags  : ${FLAGS}"

# ---------------------------------------------------------------- 2. flags
for f in omp libint libxc fftw3 spglib; do
  if [[ " ${FLAGS} " == *" ${f} "* ]]; then
    ok "feature ${f}"
  else
    warn "feature ${f} not compiled in"
  fi
done
if [[ " ${FLAGS} " == *" parallel "* ]]; then
  ok "MPI build (parallel)"
else
  info "serial build (no MPI): NP is forced to 1"
  NP=1
fi

# ---------------------------------------------------------------- 3. data dir
if [[ -n "${CP2K_DATA_DIR:-}" ]]; then
  if [[ -d "${CP2K_DATA_DIR}" ]]; then
    ok "CP2K_DATA_DIR = ${CP2K_DATA_DIR}"
    for f in BASIS_MOLOPT BASIS_MOLOPT_UZH POTENTIAL_UZH GTH_POTENTIALS BASIS_ADMM_UZH dftd3.dat t_c_g.dat; do
      if [[ -e "${CP2K_DATA_DIR}/${f}" ]]; then
        ok "  data file ${f}"
      else
        warn "  data file ${f} missing"
      fi
    done
    for f in BASIS_MOLOPT_UZH_2026.2 POTENTIAL_UZH_2026.2; do
      [[ -e "${CP2K_DATA_DIR}/${f}" ]] && info "  versioned snapshot ${f} present (CP2K >= 2026.2 data)"
    done
  else
    fail "CP2K_DATA_DIR=${CP2K_DATA_DIR} is not a directory"
  fi
else
  info "CP2K_DATA_DIR not set: the compiled-in data path will be used (see output header)"
fi

# ---------------------------------------------------------------- 4. syntax
WORK=$(mktemp -d "${TMPDIR:-/tmp}/cp2k_check.XXXXXX")
cp "${INP}" "${WORK}/"
cd "${WORK}" || exit 1
if "${EXE}" -c h2o_smoke_test.inp 2>&1 | grep -q "SUCCESS, the input could be parsed correctly"; then
  ok "syntax check (-c) passed"
else
  fail "syntax check (-c) failed:"
  "${EXE}" -c h2o_smoke_test.inp 2>&1 | tail -20
fi

# ---------------------------------------------------------------- 5. run
OUT=h2o_smoke_test.out
rm -f "${OUT}"          # cp2k -o appends to existing files
t0=$(date +%s)
if [[ "${NP}" -gt 1 ]]; then
  case "${MPI_LAUNCHER}" in
    srun*) ${MPI_LAUNCHER} -n "${NP}" "${EXE}" -i h2o_smoke_test.inp -o "${OUT}" 2> run.err ;;
    *)     ${MPI_LAUNCHER} -np "${NP}" "${EXE}" -i h2o_smoke_test.inp -o "${OUT}" 2> run.err ;;
  esac
else
  "${EXE}" -i h2o_smoke_test.inp -o "${OUT}" 2> run.err
fi
rc=$?
t1=$(date +%s)
info "run finished with exit code ${rc} in $((t1 - t0)) s (work dir: ${WORK})"

if [[ ! -s "${OUT}" ]]; then
  fail "no output written; stderr:"
  tail -20 run.err
  exit 1
fi

DDIR=$(grep -m1 "CP2K| Data directory path" "${OUT}" | awk '{print $NF}')
[[ -n "${DDIR}" ]] && info "data directory used by CP2K: ${DDIR}"

SCF=$(grep -m1 "SCF run converged in" "${OUT}")
if [[ -n "${SCF}" ]]; then ok "SCF converged      : $(echo "${SCF}" | sed 's/^ *//')"; else fail "SCF did not converge"; fi

EN=$(grep -m1 "ENERGY| Total FORCE_EVAL" "${OUT}")
if [[ -n "${EN}" ]]; then ok "total energy       : $(echo "${EN}" | sed 's/^ *//')"; else fail "no total energy line"; fi

if grep -q "PROGRAM ENDED AT" "${OUT}"; then ok "program ended normally"; else fail "PROGRAM ENDED AT not found (crash or abort)"; fi

# ---------------------------------------------------------------- 6. MPI sanity
nline=$(grep -c "DBCSR| MPI: Number of processes" "${OUT}")
nproc=$(grep -m1 "DBCSR| MPI: Number of processes" "${OUT}" | awk '{print $NF}')
if [[ "${nline}" -gt 1 ]]; then
  fail "DBCSR MPI line printed ${nline} times: independent serial copies are running (MPI mismatch)"
  info "compare: which ${MPI_LAUNCHER}  vs  ldd \$(command -v ${EXE}) | grep -i mpi"
elif [[ -n "${nproc}" && "${nproc}" != "${NP}" ]]; then
  fail "MPI ranks reported by DBCSR: ${nproc} (expected ${NP})"
elif [[ -n "${nproc}" ]]; then
  ok "MPI ranks reported by DBCSR: ${nproc} (expected ${NP})"
fi

nwarn=$(grep -m1 "The number of warnings for this run is" "${OUT}" | awk '{print $NF}')
[[ -n "${nwarn}" ]] && info "number of warnings: ${nwarn}"

echo "----------------------------------------------------------------"
if [[ "${nfail}" -eq 0 ]]; then
  echo "RESULT: PASS  (output kept in ${WORK}/${OUT})"
else
  echo "RESULT: ${nfail} check(s) FAILED  (see ${WORK})"
  exit 1
fi
