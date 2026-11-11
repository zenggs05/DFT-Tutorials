#!/bin/bash
#PBS -N openmx
#PBS -l select=2:ncpus=64:mpiprocs=32:ompthreads=2
#PBS -l walltime=24:00:00
#PBS -q workq
#PBS -j oe
#
# OpenMX 4.0 hybrid MPI x OpenMP job (PBS Pro / OpenPBS).
#   qsub -v INPUT=input.dat job_pbs.sh
# select=<nodes>:ncpus=<cores per node>:mpiprocs=<MPI per node>:ompthreads=<threads per MPI>
# For Torque use instead:  #PBS -l nodes=2:ppn=64  and set NT by hand.
set -u
cd "${PBS_O_WORKDIR}"
INPUT=${INPUT:-input.dat}
[ -f "${INPUT}" ] || { echo "input ${INPUT} not found"; exit 1; }

# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl     # site specific
if [ "${USE_DEB:-0}" = "1" ]; then
    export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:${LD_LIBRARY_PATH:-}
fi
OPENMX=${OPENMX:-openmx}

NP=$(wc -l < "${PBS_NODEFILE}")       # one line per MPI process with mpiprocs=
NT=${OMP_NUM_THREADS:-${NCPUS:-1}}    # PBS Pro sets OMP_NUM_THREADS from ompthreads
export OMP_NUM_THREADS=${NT}
export OMP_STACKSIZE=512M
ulimit -s unlimited
export I_MPI_PIN_DOMAIN=omp

NAME=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${INPUT}")
STD=${NAME}.std
echo "NP=${NP} NT=${NT} input=${INPUT} start $(date)"

mpirun -np "${NP}" -machinefile "${PBS_NODEFILE}" "${OPENMX}" "${INPUT}" -nt "${NT}" > "${STD}" 2>&1

if grep -q "The calculation was normally finished." "${STD}"; then
    echo "finished normally $(date)"
    grep -h "Utot\." "${NAME}.out" | tail -n 1
else
    echo "OpenMX did NOT finish normally; see ${STD}" >&2
    tail -n 30 "${STD}" >&2
    exit 1
fi
