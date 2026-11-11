#!/bin/bash
#SBATCH --job-name=bader-copt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster
##SBATCH --account=myproject     # edit for your cluster

# ---- environment (edit for your cluster) ----
module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}

# POTCAR must exist and follow the element order of POSCAR (Pt C O)
if [ ! -s POTCAR ]; then
    echo "POTCAR missing: cat \$VASP_PP_PATH/potpaw_PBE/{Pt,C,O}/POTCAR > POTCAR" >&2
    exit 1
fi

mpirun -np ${SLURM_NTASKS:-32} $VASP_EXE > vasp.out 2>&1

# Bader post-processing (needs chgsum.pl and bader in PATH)
bash run_bader.sh
