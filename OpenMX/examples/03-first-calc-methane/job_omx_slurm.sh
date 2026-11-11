#!/bin/bash
#SBATCH --job-name=omx-ch03
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --output=%x-%j.log
#
# Run the chapter-3 molecules (CH4, H2O, CO) one after another.
# Works as a Slurm batch script (sbatch job_omx_slurm.sh) and also as a plain
# bash script on a workstation (bash job_omx_slurm.sh).
#
# Environment variables you may set:
#   OPENMX    path of the openmx executable           (default: openmx)
#   OMX_DATA  absolute path of DFT_DATA19; if set, the DATA.PATH line of every
#             input is rewritten before the run      (default: unset)
#   INPUTS    list of input files                     (default: the three molecules)
#   NP        number of MPI processes outside Slurm   (default: 4)
#   MPIRUN    MPI launcher                            (default: mpirun)
#
set -euo pipefail

# --- load your environment here (examples, adapt to your cluster) ---
# module load intel-oneapi-compilers intel-oneapi-mkl intel-oneapi-mpi
# source /opt/intel/oneapi/setvars.sh
# export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:${LD_LIBRARY_PATH:-}   # deb package tools

OPENMX=${OPENMX:-openmx}
MPIRUN=${MPIRUN:-mpirun}
INPUTS=${INPUTS:-"Methane_manual_check.dat Methane.dat H2O.dat CO.dat"}
NP=${SLURM_NTASKS:-${NP:-4}}
NT=${SLURM_CPUS_PER_TASK:-1}

export OMP_NUM_THREADS=${NT}
export OMP_STACKSIZE=512M
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=${NT}
ulimit -s unlimited 2>/dev/null || echo "warning: could not set unlimited stack size"

if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then cd "${SLURM_SUBMIT_DIR}"; fi

if ! command -v "${OPENMX}" >/dev/null 2>&1 && [ ! -x "${OPENMX}" ]; then
  echo "error: openmx executable '${OPENMX}' not found (set OPENMX=/full/path/openmx)" >&2
  exit 1
fi

for inp in ${INPUTS}; do
  if [ ! -f "${inp}" ]; then echo "skip: ${inp} not found"; continue; fi

  if [ -n "${OMX_DATA:-}" ]; then
    sed -i "s#^[[:space:]]*DATA\.PATH.*#DATA.PATH                     ${OMX_DATA}#" "${inp}"
  fi
  datapath=$(awk 'tolower($1)=="data.path"{print $2; exit}' "${inp}")
  if [ ! -d "${datapath}/PAO" ] || [ ! -d "${datapath}/VPS" ]; then
    echo "error: DATA.PATH '${datapath}' in ${inp} does not contain PAO/ and VPS/" >&2
    echo "       edit the file or run with OMX_DATA=/abs/path/to/DFT_DATA19" >&2
    exit 1
  fi

  name=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${inp}")
  echo "=== $(date '+%F %T')  running ${inp}  (System.Name=${name}, np=${NP}, nt=${NT})"
  ${MPIRUN} -np "${NP}" "${OPENMX}" "${inp}" -nt "${NT}" > "${name}.std" 2>&1 || {
    echo "openmx returned an error for ${inp}; see ${name}.std" >&2; }
  tail -n 1 "${name}.std"
done

echo "=== all done; summary:"
bash ./check_results.sh || true
