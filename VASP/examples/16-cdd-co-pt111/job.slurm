#!/bin/bash
#SBATCH --job-name=cdd-copt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=06:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster
##SBATCH --account=myproject     # edit for your cluster

# ---- environment (edit for your cluster) ----
module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
export VASP_EXE=${VASP_EXE:-vasp_std}
export NP=${SLURM_NTASKS:-32}

# three single points (AB, slab, CO) followed by the difference
bash run_cdd.sh prepare
bash run_cdd.sh run
bash run_cdd.sh analyze
