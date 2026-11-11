#!/bin/bash
#SBATCH --job-name=water-npt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out
#SBATCH --partition=normal
# Adjust partition/account/module lines for your cluster.

# module load vasp/6.6.1        # example only
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP=${VASP:-vasp_gam}
RUN="srun ${VASP}"              # or: mpirun -np ${SLURM_NTASKS} ${VASP}

if [ ! -s POTCAR ]; then
    PP=${VASP_PP_PATH:?please set VASP_PP_PATH}
    cat "$PP/potpaw_PBE/O/POTCAR" "$PP/potpaw_PBE/H/POTCAR" > POTCAR
fi
[ -s POSCAR ] || python3 build_water_box.py

${RUN} > vasp.out 2>&1

# density along the run: rho = M(amu) x 1.66054 / V(A^3); M = 32 x 18.015 amu
grep "volume of cell" OUTCAR | awk 'NR>1{print NR-1, 957.27/$5}' > density.dat
tail -3 density.dat
