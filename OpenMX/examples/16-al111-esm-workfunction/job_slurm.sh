#!/bin/bash
#SBATCH --job-name=al111-wf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log
# Usage: sbatch job_slurm.sh Al111_ESM_on1_2E.dat
# Runs OpenMX (MPI x OpenMP), then the work-function analysis if gcube2oned is available.

set -u
INPUT=${1:-Al111_ESM_on1_2E.dat}
OPENMX=${OPENMX:-openmx}                 # e.g. /home/user/openmx4.0/work/openmx
GCUBE2ONED=${GCUBE2ONED:-./gcube2oned}   # compile: gcc gcube2oned.c -lm -o gcube2oned

# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl   # site specific

NP=${SLURM_NTASKS:-4}
NT=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${NT}
export OMP_STACKSIZE=512M
ulimit -s unlimited

NAME=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${INPUT}")
STD=${NAME}.std

mpirun -np "${NP}" "${OPENMX}" "${INPUT}" -nt "${NT}" > "${STD}" 2>&1

# OpenMX returns exit code 0 even after input errors: check the log text instead.
if ! grep -q "The calculation was normally finished." "${STD}"; then
    echo "OpenMX did not finish normally, see ${STD}" >&2
    tail -n 20 "${STD}" >&2
    exit 1
fi
grep -i "Chemical potential (Hartree)" "${NAME}.out" | head -n 1

if [ -x "${GCUBE2ONED}" ] && [ -f "${NAME}.v0.cube" ]; then
    "${GCUBE2ONED}" "${NAME}.v0.cube" 1 > "${NAME}_1d_pot.txt"
    python3 workfunc.py "${NAME}_1d_pot.txt" "${NAME}.out" --plot "${NAME}_wf.png"
fi
