#!/bin/bash
#SBATCH --job-name=al-elastic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
#SBATCH --partition=normal
# Adjust partition/account/module lines for your cluster.

# module load vasp/6.6.1        # example only; use your site module
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP=${VASP:-vasp_std}
RUN="srun ${VASP}"              # or: mpirun -np ${SLURM_NTASKS} ${VASP}

# POTCAR is never stored in the repository: build it from your own library
if [ ! -s POTCAR ]; then
    cat "${VASP_PP_PATH:?please set VASP_PP_PATH}/potpaw_PBE/Al/POTCAR" > POTCAR
fi
cp POSCAR POSCAR.orig

# ---- step 1: ISIF=3 relaxation, two rounds against Pulay stress ----
cp INCAR.relax INCAR
for i in 1 2; do
    ${RUN} > vasp_relax${i}.out 2>&1
    cp OUTCAR OUTCAR.relax${i}
    cp CONTCAR POSCAR
    grep "volume of cell" OUTCAR | tail -1
done

# ---- step 2: elastic tensor (IBRION=6, ISIF=3) ----
cp INCAR.elastic INCAR
${RUN} > vasp_elastic.out 2>&1
grep -A 10 "TOTAL ELASTIC MODULI" OUTCAR
python3 voigt_vrh.py -f OUTCAR --poscar POSCAR > elastic_summary.txt
cat elastic_summary.txt
