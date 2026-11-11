#!/bin/bash
# =====================================================================
# Generic PBS template (PBS Pro / OpenPBS syntax), pure MPI job.
# Torque users: see the commented "Torque" lines.
#
#   qsub pbs_mpi.sh
#   qstat -u $USER          qdel <jobid>          qstat -f <jobid>
# =====================================================================
#PBS -N dft-mpi
#PBS -q workq                                  # queue name (qstat -Q)
#PBS -l select=2:ncpus=64:mpiprocs=64:ompthreads=1:mem=200gb
#                  ^ chunks (nodes)  ^ MPI ranks per chunk
#PBS -l walltime=24:00:00
#PBS -j oe                                     # merge stdout and stderr
#PBS -o dft-mpi.log
##PBS -A my_project                            # account / project
##PBS -m ae -M you@example.com                 # mail on abort / end
##PBS -l place=scatter:excl                    # one chunk per node, exclusive
# --- Torque equivalents ---------------------------------------------
##PBS -l nodes=2:ppn=64
##PBS -l walltime=24:00:00

EXE="vasp_std"
ARGS=""
LOG="run.log"

module purge
# module load intel-oneapi-mpi intel-oneapi-mkl vasp/6.5.1

ulimit -s unlimited
export OMP_NUM_THREADS=1

# PBS starts the script in $HOME: go to the submission directory
cd "${PBS_O_WORKDIR}" || exit 1

NP=$(wc -l < "${PBS_NODEFILE}")      # one line per MPI slot
echo "Job ${PBS_JOBID}: ${NP} MPI ranks on $(sort -u "${PBS_NODEFILE}" | tr '\n' ' ')"

# Most MPI libraries detect PBS automatically; the hostfile is explicit
# here for clarity (Intel MPI: -machinefile, Open MPI: --hostfile).
mpirun -np "${NP}" -machinefile "${PBS_NODEFILE}" ${EXE} ${ARGS} > "${LOG}" 2>&1
status=$?
echo "exit code ${status}"
exit ${status}
