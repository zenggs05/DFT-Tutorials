#!/bin/bash
#SBATCH -J omx-gra-pz
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out
#
# Graphene: SCF + 2 pz MLWFs, then hopping table and band comparison.

# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-8}

mpirun -np ${NP} openmx Graphene_pz_MLWF.dat -nt ${OMP_NUM_THREADS} > graphene_pz.std 2>&1

grep "Selected bands within the Outer and Inner Windows" -A 6 graphene_pz.std
grep "CONVERGENCE ACHIEVED" graphene_pz.std | tail -n 2
python3 hwr_hopping.py graphene_pz.HWR --top 12
python3 ../14-si-mlwf-sp3/band_compare.py graphene_pz.Band graphene_pz.Wannier_Band \
        --emin -2.0 --emax 1.5 --plot graphene_pz_bands.png --ylim -10 10
