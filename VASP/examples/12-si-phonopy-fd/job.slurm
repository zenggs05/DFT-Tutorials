#!/bin/bash
#SBATCH --job-name=si-phon
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster
#
# Submit from a calculation directory, e.g.  cd disp-001 && sbatch ../job.slurm
module purge
# module load vasp/6.5.1          # edit for your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
mpirun -np ${SLURM_NTASKS:-32} $VASP_EXE > vasp.out 2>&1
grep -q "</modeling>" vasprun.xml && echo "vasprun.xml complete" || echo "vasprun.xml incomplete!"
