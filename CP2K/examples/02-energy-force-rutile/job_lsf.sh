#!/bin/bash
# job_lsf.sh - LSF job script for the rutile TiO2 examples (CP2K, cp2k.psmp or cp2k.popt)
#
# This is the corrected version of the original cp2k.slurm, which was an LSF
# script despite its name and ran "input.inp" instead of "input_file.inp".
#
# Submit (LSF reads the #BSUB lines from standard input):
#     bsub < job_lsf.sh
# Another input:
#     INP=rutile_kpoints.inp bsub < job_lsf.sh      (if your site forwards the environment)
#     or edit the INP line below.
#BSUB -J rutile-e-f
#BSUB -q 52-medium
#BSUB -W 12:00
#BSUB -n 28
#BSUB -R "span[ptile=28]"
#BSUB -o cp2k.%J.stdout
#BSUB -e cp2k.%J.stderr

set -uo pipefail
cd "${LS_SUBCWD:-$PWD}"

# ---- CP2K environment (edit for your cluster) -----------------------------
# The original script loaded:  intel/17.5.239  mpi/intel/2017.5.239  gcc/5.5.0  cp2k/7.1
# CP2K 7.1 is very old; prefer a current module, e.g.:
# module load cp2k/2026.2
# or a toolchain build:
# source $HOME/src/cp2k-2026.2/install/cp2k_env
# ---------------------------------------------------------------------------

CP2K_EXE=${CP2K_EXE:-cp2k.psmp}
INP=${INP:-input_file.inp}
OUT=${INP%.inp}_${LSB_JOBID:-local}.out
NP=${LSB_DJOB_NUMPROC:-28}

export OMP_NUM_THREADS=1          # pure MPI, same as the original cp2k.popt run
ulimit -s unlimited

[[ -f "${INP}" ]] || { echo "ERROR: input ${INP} not found" >&2; exit 1; }

# Intel MPI: mpiexec.hydra / mpirun; OpenMPI: mpirun. Use the launcher that
# belongs to the MPI library CP2K was built with.
mpirun -np "${NP}" "${CP2K_EXE}" -i "${INP}" -o "${OUT}"
rc=$?

if grep -q "PROGRAM ENDED AT" "${OUT}"; then
  echo "CP2K finished normally."
  grep "ENERGY| Total FORCE_EVAL" "${OUT}" | tail -1
else
  echo "CP2K did NOT finish normally (exit code ${rc}); check ${OUT}." >&2
fi
exit ${rc}
