#!/bin/bash
# =====================================================================
# Generic Slurm template: hybrid MPI x OpenMP job.
# Typical users: CP2K (cp2k.psmp), OpenMX (openmx -nt N), VASP built
# with OpenMP (and every VASP GPU run: 1 rank per GPU + threads).
#
# cores used per node = ntasks-per-node x cpus-per-task
# Example below: 64-core node = 16 MPI ranks x 4 OpenMP threads.
# Choose threads per rank so that one rank stays inside one NUMA domain
# (lscpu | grep -i numa) and benchmark 1, 2, 4, 8 threads per rank.
# =====================================================================
#SBATCH --job-name=dft-hybrid
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16         # MPI ranks per node
#SBATCH --cpus-per-task=4            # OpenMP threads per rank
#SBATCH --time=24:00:00
#SBATCH --mem=0
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
##SBATCH --account=my_project
##SBATCH --exclusive

# ----------------------------- user settings -------------------------
CODE="cp2k"               # cp2k | openmx | vasp
INPUT="input.inp"         # cp2k: input.inp, openmx: Si.dat, vasp: (none)
LOG="run.log"

# ----------------------------- environment ---------------------------
module purge
# module load gcc openmpi openblas scalapack fftw        # site specific
# module load cp2k/2026.2      or: source /path/to/cp2k/tools/toolchain/install/setup

ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_PLACES=cores           # one thread per physical core
export OMP_PROC_BIND=close        # keep a rank's threads next to each other
export OMP_STACKSIZE=512m         # private stack of every OpenMP thread

cd "${SLURM_SUBMIT_DIR}" || exit 1
echo "Job ${SLURM_JOB_ID}: ${SLURM_NTASKS} ranks x ${OMP_NUM_THREADS} threads on ${SLURM_JOB_NODELIST}"

case "${CODE}" in
    cp2k)   CMD="cp2k.psmp -i ${INPUT}" ;;
    openmx) CMD="openmx ${INPUT} -nt ${OMP_NUM_THREADS}" ;;   # -nt = OMP threads
    vasp)   CMD="vasp_std" ;;   # VASP sets NCORE=1 itself when threads > 1
    *) echo "unknown CODE=${CODE}" >&2; exit 1 ;;
esac

# Option A (recommended): Slurm launcher; binding handled by Slurm.
srun --ntasks="${SLURM_NTASKS}" --cpus-per-task="${OMP_NUM_THREADS}" \
     --cpu-bind=cores ${CMD} > "${LOG}" 2>&1
status=$?

# Option B: Open MPI's mpirun (4.x syntax; 5.x writes the modifier as pe=N)
# mpirun -np ${SLURM_NTASKS} --map-by ppr:${SLURM_NTASKS_PER_NODE}:node:PE=${OMP_NUM_THREADS} \
#        --bind-to core -x OMP_NUM_THREADS -x OMP_PLACES -x OMP_PROC_BIND \
#        -x OMP_STACKSIZE --report-bindings ${CMD} > "${LOG}" 2>&1
# Option C: Intel MPI's mpirun (one pinning domain per rank sized by OMP)
# export I_MPI_PIN_DOMAIN=omp
# mpirun -np ${SLURM_NTASKS} ${CMD} > "${LOG}" 2>&1

echo "exit code ${status}"
exit ${status}
