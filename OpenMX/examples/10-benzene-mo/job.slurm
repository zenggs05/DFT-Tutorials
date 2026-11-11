#!/bin/bash
#SBATCH -J omx-benzene-mo
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH -t 00:30:00
#SBATCH -o slurm-%j.out

# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mpirun -np ${SLURM_NTASKS} openmx Benzene_MO.dat -nt ${OMP_NUM_THREADS} > benzene.std 2>&1

ls -l benzene.homo* benzene.lumo*
grep "HOMO =" benzene.out
grep -A 22 "Eigenvalues (Hartree)" benzene.out | head -26
for f in benzene.homo0_0.cube benzene.homo0_1.cube benzene.lumo0_0.cube benzene.lumo0_1.cube; do
  echo "$f : $(head -n 1 $f)"
done
