#!/bin/bash
# =====================================================================
# Generic IBM Spectrum LSF template, pure MPI job.
#
#   bsub < lsf_mpi.sh        (note the "<": LSF reads #BSUB lines from stdin)
#   bjobs                    bkill <jobid>          bjobs -l <jobid>
#   bpeek <jobid>            (look at stdout of a running job)
# =====================================================================
#BSUB -J dft-mpi
#BSUB -q normal                      # queue (bqueues lists them)
#BSUB -n 128                         # TOTAL number of MPI slots
#BSUB -R "span[ptile=64]"            # 64 slots per host -> 2 hosts
#BSUB -W 24:00                       # wall time HH:MM
#BSUB -o %J.out                      # %J = job id
#BSUB -e %J.err
##BSUB -R "rusage[mem=3000]"         # memory per slot; unit is site
                                     # dependent (LSF_UNIT_FOR_LIMITS)
##BSUB -x                            # exclusive hosts
##BSUB -P my_project

EXE="vasp_std"
ARGS=""
LOG="run.log"

module purge
# module load intel-oneapi-mpi intel-oneapi-mkl vasp/6.5.1

ulimit -s unlimited
export OMP_NUM_THREADS=1

# LSF starts the job in the submission directory (LS_SUBCWD)
cd "${LS_SUBCWD:-$PWD}" || exit 1

NP=${LSB_DJOB_NUMPROC:-128}
echo "Job ${LSB_JOBID}: ${NP} slots on: ${LSB_HOSTS}"

# Intel MPI and Open MPI both read the LSF allocation; some sites
# provide a wrapper (mpirun.lsf, blaunch) - check the local guide.
mpirun -np "${NP}" ${EXE} ${ARGS} > "${LOG}" 2>&1
status=$?
echo "exit code ${status}"
exit ${status}
