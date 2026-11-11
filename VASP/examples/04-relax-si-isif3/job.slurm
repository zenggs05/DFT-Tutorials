#!/bin/bash
#SBATCH --job-name=si-isif3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster
##SBATCH --account=myproject     # edit for your cluster

# ---- environment (edit for your cluster) ----
module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
MPI_RUN="mpirun -np ${SLURM_NTASKS:-16}"

# Keep the original starting structure
[ -f POSCAR.orig ] || cp POSCAR POSCAR.orig

# Relax repeatedly (CONTCAR -> POSCAR) until one round needs <= 3 ionic steps.
# Every new run rebuilds the plane-wave basis for the new cell, which removes
# the Pulay-stress bias of the previous round.
for round in 1 2 3 4; do
    $MPI_RUN $VASP_EXE > vasp.out.$round 2>&1
    nstep=$(grep -c " F= " OSZICAR)
    vol=$(grep "volume of cell" OUTCAR | tail -1 | awk '{print $NF}')
    pres=$(grep "external pressure" OUTCAR | tail -1 | awk '{print $4}')
    conv="no"
    grep -q "reached required accuracy" OUTCAR && conv="yes"
    echo "round $round: ionic steps = $nstep, volume = $vol A^3, pressure = $pres kB, converged = $conv" | tee -a relax_rounds.log
    mkdir -p round$round
    cp OSZICAR OUTCAR CONTCAR INCAR round$round/
    cp CONTCAR POSCAR
    if [ "$conv" = "yes" ] && [ "$nstep" -le 3 ]; then
        echo "volume converged after round $round" | tee -a relax_rounds.log
        break
    fi
done
