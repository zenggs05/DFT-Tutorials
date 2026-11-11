#!/bin/bash
#SBATCH -J omx_cohp_cddf
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH -t 06:00:00
#SBATCH -o slurm-%j.out

# COHP/COOP and optical properties (CDDF) of Si with OpenMX 4.0 (edit modules/paths).
# The post-processing code is built with "make cohp" in openmx4.0/source.
# module load openmx/4.0
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# ---- COHP / COOP ----
mpirun -np ${SLURM_NTASKS} openmx Si_COHP.dat -nt ${OMP_NUM_THREADS} > si_cohp.std
for PAIR in 1 2 3; do
  # answers: tetrahedron(1), number of energy meshes, pair index
  printf '1\n3000\n%d\n' ${PAIR} | cohp si_cohp.cohp > cohp_pair${PAIR}.log
done
python3 plot_cohp_cddf.py cohp si_cohp_i1_s0_tetra.cohp si_cohp_i2_s0_tetra.cohp \
        si_cohp_i3_s0_tetra.cohp --out si_cohp.png > cohp_summary.txt

# ---- optical conductivity / dielectric function ----
mpirun -np ${SLURM_NTASKS} openmx Si_CDDF.dat -nt ${OMP_NUM_THREADS} > si_cddf.std
python3 plot_cohp_cddf.py optics si_cddf --out si_optics.png > optics_summary.txt
