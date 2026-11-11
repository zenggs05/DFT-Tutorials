#!/bin/bash
#SBATCH --job-name=mos2-2l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --partition=normal
#SBATCH --output=slurm-%j.out
#
# Minimal Slurm job script for this example.
# Adjust partition, core count, module names and VASP_EXE to your cluster.
# A more complete template lives in ../../scripts/job_vasp_slurm.sh

module purge
module load vasp/6.6.1            # site specific: change to your module name

export OMP_NUM_THREADS=1
ulimit -s unlimited

VASP_EXE=${VASP_EXE:-vasp_std} # vasp_std | vasp_gam | vasp_ncl

# Basic sanity checks before wasting queue time
for f in INCAR POSCAR KPOINTS POTCAR; do
    [ -s "$f" ] || { echo "missing $f in $(pwd)"; exit 1; }
done

mpirun -np ${SLURM_NTASKS} ${VASP_EXE} > vasp.out 2>&1

grep -q "reached required accuracy" OUTCAR && echo "ionic relaxation converged"
grep "F=" OSZICAR | tail -1
