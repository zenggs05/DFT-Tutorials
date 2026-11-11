#!/bin/bash
#SBATCH --job-name=zno-r2scan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.6.1        # avoid 6.6.0 for spin-polarized meta-GGA forces (known issue 94)
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
MPI_RUN="mpirun -np ${SLURM_NTASKS:-16}"

# The POTCAR must contain kinetic-energy-density information for meta-GGA
if [ "$(grep -c 'kinetic energy-density' POTCAR)" -eq 0 ]; then
    echo "POTCAR has no kinetic energy-density data: not suitable for r2SCAN"; exit 1
fi

[ -f POSCAR.orig ] || cp POSCAR POSCAR.orig
# Two rounds: the second one removes the Pulay-stress bias of the first cell relaxation
for round in 1 2; do
    $MPI_RUN $VASP_EXE > vasp.out.$round 2>&1
    mkdir -p round$round && cp OSZICAR OUTCAR CONTCAR round$round/
    echo "round $round: $(grep -c ' F= ' OSZICAR) ionic steps, $(grep 'volume of cell' OUTCAR | tail -1)"
    cp CONTCAR POSCAR
done
