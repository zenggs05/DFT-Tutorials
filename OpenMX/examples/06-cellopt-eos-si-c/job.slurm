#!/bin/bash
#SBATCH --job-name=omx-cell
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal            # edit for your cluster

# Usage:
#   sbatch job.slurm Si_RFC5_cellopt.dat
#   sbatch job.slurm Cdia_EvsLC.dat
module purge
# module load openmx/4.0.1                     # edit for your cluster
# Debian package users: export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_STACKSIZE=512M
export OPENMX=${OPENMX:-openmx}
export NT=$OMP_NUM_THREADS
export MPI_RUN="mpirun -np ${SLURM_NTASKS:-8}"

target=${1:?usage: sbatch job.slurm <input.dat | driver.sh [args]>}
shift
case "$target" in
    *.dat) $MPI_RUN $OPENMX "$target" -nt $NT > "${target%.dat}.std" 2>&1 ;;
    *.sh)  bash "$target" "$@" ;;
    *)     echo "unknown target: $target" >&2; exit 1 ;;
esac
