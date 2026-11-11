#!/bin/bash
#SBATCH --job-name=gaas-lepsilon
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
    PP=${VASP_PP_PATH:?please set VASP_PP_PATH}
    cat "$PP/potpaw_PBE/Ga_d/POTCAR" "$PP/potpaw_PBE/As/POTCAR" > POTCAR
fi

# Part A: electronic dielectric tensor and Born charges
mkdir -p A_lepsilon && cp INCAR POSCAR KPOINTS POTCAR A_lepsilon/
(cd A_lepsilon && ${RUN} > vasp.out 2>&1)

# Part B: + ionic contribution (IBRION = 8)
mkdir -p B_ionic && cp POSCAR KPOINTS POTCAR B_ionic/ && cp INCAR.ionic B_ionic/INCAR
(cd B_ionic && ${RUN} > vasp.out 2>&1)

grep -A 5 "MACROSCOPIC STATIC DIELECTRIC TENSOR" A_lepsilon/OUTCAR
grep -A 5 "MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION" B_ionic/OUTCAR
