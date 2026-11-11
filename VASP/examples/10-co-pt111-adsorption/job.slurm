#!/bin/bash
#SBATCH --job-name=co-pt111
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster
##SBATCH --account=myproject     # edit for your cluster
#
# Submit from inside one calculation directory, e.g.
#   cd 1_slab && sbatch ../job.slurm
# Gas molecules use the Gamma-only executable:
#   cd 2_co_gas && VASP_EXE=vasp_gam sbatch --export=ALL ../job.slurm

module purge
# module load vasp/6.5.1          # edit for your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}

for f in INCAR POSCAR KPOINTS POTCAR; do
    [ -f "$f" ] || { echo "missing $f in $(pwd)"; exit 1; }
done
mpirun -np ${SLURM_NTASKS:-32} $VASP_EXE > vasp.out 2>&1

grep -q "reached required accuracy" OUTCAR && echo "relaxation converged" || \
    echo "no 'reached required accuracy' (normal for NSW<=1 runs; otherwise check OSZICAR)"
grep "energy(sigma->0)" OUTCAR | tail -1
