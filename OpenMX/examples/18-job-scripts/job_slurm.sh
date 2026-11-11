#!/bin/bash
#SBATCH --job-name=openmx
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32          # MPI processes per node
#SBATCH --cpus-per-task=2             # OpenMP threads per MPI process (= -nt)
#SBATCH --time=24:00:00
#SBATCH --partition=cpu               # site specific
#SBATCH --output=%x-%j.log
#
# OpenMX 4.0 hybrid MPI x OpenMP job (Slurm).
#   sbatch job_slurm.sh input.dat
# Total cores = nodes x ntasks-per-node x cpus-per-task (here 2 x 32 x 2 = 128).
#
# Environment switches (export before sbatch, or edit here):
#   OPENMX=/path/to/openmx       executable (default: openmx in PATH)
#   USE_DEB=1                    official Debian package in /opt/openmx4.0
#   LAUNCH=srun                  use srun instead of mpirun (MPI built with Slurm PMI support)
#   KEEP_SCRATCH=0               remove *_rst/ after a successful run
set -u
INPUT=${1:?usage: sbatch job_slurm.sh input.dat}
[ -f "${INPUT}" ] || { echo "input ${INPUT} not found"; exit 1; }

# ---------- environment ----------
# module purge
# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl      # Intel toolchain
# module load gcc openmpi openblas scalapack fftw                           # GNU toolchain
if [ "${USE_DEB:-0}" = "1" ]; then
    # runtime libraries of the official .deb (fixes 'libmpifort.so.12: cannot open shared object file')
    export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:${LD_LIBRARY_PATH:-}
fi
OPENMX=${OPENMX:-openmx}

NP=${SLURM_NTASKS:-1}
NT=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${NT}          # keep identical to -nt (MKL/OpenBLAS threads follow OpenMP)
export OMP_STACKSIZE=512M             # large private stacks for OpenMP threads
export OMP_PLACES=cores
export OMP_PROC_BIND=close
ulimit -s unlimited                   # avoid segmentation faults from deep stacks
# Intel MPI: one pinning domain per MPI process, sized by OMP_NUM_THREADS
export I_MPI_PIN_DOMAIN=omp

NAME=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${INPUT}")
STD=${NAME}.std
echo "host $(hostname)  NP=${NP}  NT=${NT}  input=${INPUT}  name=${NAME}  start $(date)"

# ---------- run ----------
if [ "${LAUNCH:-mpirun}" = "srun" ]; then
    srun --cpu-bind=cores "${OPENMX}" "${INPUT}" -nt "${NT}" > "${STD}" 2>&1
else
    # Open MPI binding (ignored by Intel MPI, which uses I_MPI_PIN_DOMAIN):
    #   mpirun --map-by ppr:${SLURM_NTASKS_PER_NODE}:node:PE=${NT} --bind-to core ...
    mpirun -np "${NP}" "${OPENMX}" "${INPUT}" -nt "${NT}" > "${STD}" 2>&1
fi

# ---------- check (OpenMX returns exit code 0 even after input errors) ----------
if grep -q "The calculation was normally finished." "${STD}"; then
    echo "finished normally $(date)"
    grep -h "Utot\." "${NAME}.out" | tail -n 1
    grep -h "Total Computational Time" "${NAME}.out" | tail -n 1
    if [ "${KEEP_SCRATCH:-1}" = "0" ]; then rm -rf "${NAME}_rst"; fi
else
    echo "OpenMX did NOT finish normally; last lines of ${STD}:" >&2
    tail -n 30 "${STD}" >&2
    exit 1
fi
