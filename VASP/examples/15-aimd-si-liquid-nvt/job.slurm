#!/bin/bash
#SBATCH --job-name=si-liquid-aimd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --partition=normal
# Adjust partition/account/module lines for your cluster.

# module load vasp/6.6.1        # example only
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP=${VASP:-vasp_gam}          # Gamma-only build for a Gamma-point supercell
RUN="srun ${VASP}"              # or: mpirun -np ${SLURM_NTASKS} ${VASP}

if [ ! -s POTCAR ]; then
    cat "${VASP_PP_PATH:?please set VASP_PP_PATH}/potpaw_PBE/Si/POTCAR" > POTCAR
fi
cp POSCAR POSCAR.orig

# stage 1: melt
cp INCAR.melt INCAR
${RUN} > vasp_melt.out 2>&1
for f in OUTCAR OSZICAR XDATCAR REPORT; do cp $f $f.melt 2>/dev/null; done
cp CONTCAR POSCAR

# stage 2: production
cp INCAR.prod INCAR
${RUN} > vasp_prod.out 2>&1

python3 analyze_md.py --potim 2.0 --skip 500 --pair Si-Si
