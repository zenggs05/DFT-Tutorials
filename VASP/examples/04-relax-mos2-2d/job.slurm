#!/bin/bash
#SBATCH --job-name=mos2-2d
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.5.1        # LATTICE_CONSTRAINTS needs >= 6.4.3 with IBRION=1/2
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}

mpirun -np ${SLURM_NTASKS:-16} $VASP_EXE > vasp.out 2>&1

# quick report
grep -q "reached required accuracy" OUTCAR && echo "relaxation converged" || echo "NOT converged"
echo "final lattice vectors (a3 must still be 0 0 20):"
sed -n 3,5p CONTCAR
grep "in kB" OUTCAR | tail -1
