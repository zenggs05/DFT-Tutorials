#!/bin/bash
# =====================================================================
# Generic Slurm template: pure MPI job (VASP / CP2K / OpenMX / any MPI
# code). One MPI rank per physical core, 1 OpenMP thread per rank.
#
#   sbatch slurm_mpi.sh                 # submit
#   sbatch --nodes=2 slurm_mpi.sh       # command-line options override
#                                       # the #SBATCH lines below
# Lines starting with "##SBATCH" are switched off (comment).
# =====================================================================
#SBATCH --job-name=dft-mpi
#SBATCH --partition=cpu              # sinfo -s lists the partitions
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64         # = physical cores per node
#SBATCH --cpus-per-task=1            # pure MPI: 1 core per rank
#SBATCH --time=24:00:00              # wall time limit (D-HH:MM:SS also ok)
#SBATCH --mem=0                      # 0 = all memory of the node (site
                                     # dependent; else e.g. --mem=180G)
#SBATCH --output=%x-%j.out           # %x = job name, %j = job id
#SBATCH --error=%x-%j.err
##SBATCH --account=my_project        # needed on many clusters
##SBATCH --qos=normal
##SBATCH --exclusive                 # do not share the node
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=you@example.com
##SBATCH --signal=B:USR1@300         # send USR1 to this script 5 min
                                     # before the time limit (see trap)

# ----------------------------- user settings -------------------------
EXE="vasp_std"            # vasp_std | vasp_gam | vasp_ncl | cp2k.psmp | openmx
ARGS=""                   # CP2K: "-i input.inp"   OpenMX: "Si.dat"
LOG="run.log"             # stdout+stderr of the program
LAUNCHER="srun"           # srun (Slurm native) or mpirun (MPI library)

# ----------------------------- environment ---------------------------
# Load EXACTLY the compiler/MPI/library modules the program was built
# with (ask your admins or run "module spider vasp").
module purge
# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl
# module load vasp/6.5.1
# source ~/miniforge3/bin/activate dft      # only if Python is needed

ulimit -s unlimited       # large stack: avoids SIGSEGV in Fortran codes
export OMP_NUM_THREADS=1  # pure MPI: never let MKL/OpenMP oversubscribe

cd "${SLURM_SUBMIT_DIR}" || exit 1

# Refuse to run with Windows (CRLF) line endings in the input files
# (-U only matters for Git-Bash/Windows grep; it is a no-op on Linux)
crlf=$(grep -lU $'\r' INCAR POSCAR KPOINTS *.inp *.dat *.gjf *.com 2>/dev/null)
if [ -n "${crlf}" ]; then
    echo "ERROR: CRLF line endings in: ${crlf}  -> run dos2unix" >&2
    exit 1
fi

# Optional: copy results back / write a marker when the time limit is near
trap 'echo "USR1 received: time limit is near" >> "${LOG}"' USR1

echo "Job ${SLURM_JOB_ID} on ${SLURM_JOB_NODELIST}"
echo "Nodes=${SLURM_JOB_NUM_NODES} ntasks=${SLURM_NTASKS} cpus/task=${SLURM_CPUS_PER_TASK:-1}"
echo "Start: $(date)"

# ----------------------------- run -----------------------------------
if [ "${LAUNCHER}" = "srun" ]; then
    # Passing --cpus-per-task explicitly works with every Slurm version
    # (22.05-23.02 did not inherit it from sbatch).
    srun --ntasks="${SLURM_NTASKS}" --cpus-per-task="${SLURM_CPUS_PER_TASK:-1}" \
         ${EXE} ${ARGS} > "${LOG}" 2>&1 &
else
    mpirun -np "${SLURM_NTASKS}" ${EXE} ${ARGS} > "${LOG}" 2>&1 &
fi
PID=$!
# "wait" returns early when a trapped signal arrives; keep waiting until
# the program itself has finished.
while true; do
    wait "${PID}"; status=$?
    kill -0 "${PID}" 2>/dev/null || break
done

echo "End: $(date), exit code ${status}"
tail -n 5 "${LOG}"
exit ${status}
