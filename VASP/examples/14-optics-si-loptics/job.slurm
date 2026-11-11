#!/bin/bash
#SBATCH --job-name=si-loptics
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#SBATCH --partition=normal
# Adjust partition/account/module lines for your cluster.

# module load vasp/6.6.1        # example only
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP=${VASP:-vasp_std}
RUN="srun ${VASP}"              # or: mpirun -np ${SLURM_NTASKS} ${VASP}

if [ ! -s POTCAR ]; then
    cat "${VASP_PP_PATH:?please set VASP_PP_PATH}/potpaw_PBE/Si/POTCAR" > POTCAR
fi

# step 1: ground state
cp INCAR.scf INCAR
${RUN} > vasp_scf.out 2>&1
cp OUTCAR OUTCAR.scf

# step 2: empty bands + LOPTICS (reads WAVECAR)
cp INCAR.optics INCAR
${RUN} > vasp_optics.out 2>&1
grep "NBANDS=" OUTCAR | head -1
python3 plot_optics.py --tauc direct
