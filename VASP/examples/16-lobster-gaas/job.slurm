#!/bin/bash
#SBATCH --job-name=gaas-lobster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.5.1
ulimit -s unlimited
VASP_EXE=${VASP_EXE:-vasp_std}           # never vasp_gam for LOBSTER
LOBSTER_EXE=${LOBSTER_EXE:-lobster}      # e.g. lobster-5.1.1 or lobster-6.0.0, as downloaded

[ -s POTCAR ] || { echo "POTCAR missing: cat \$VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR" >&2; exit 1; }

# 1) VASP single point
export OMP_NUM_THREADS=1
mpirun -np ${SLURM_NTASKS:-16} $VASP_EXE > vasp.out 2>&1
grep -q "General timing" OUTCAR || { echo "VASP did not finish" >&2; exit 1; }
if grep -q "number of bands has been changed" OUTCAR vasp.out; then
    echo "WARNING: VASP changed NBANDS - adjust NBANDS or the parallel setup" >&2
fi

# 2) LOBSTER (OpenMP parallel)
export OMP_NUM_THREADS=${SLURM_NTASKS:-16}
$LOBSTER_EXE > lobster.out 2>&1
grep -i "spilling" lobsterout
head -5 ICOHPLIST.lobster
