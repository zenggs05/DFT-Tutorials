#!/bin/bash
#SBATCH --job-name=si-mbj
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.6.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}          # KPOINTS_OPT needs vasp_std
MPI_RUN="mpirun -np ${SLURM_NTASKS:-16}"

# step 1: PBE -> WAVECAR (KPOINTS_OPT ignored via LKPOINTS_OPT = .FALSE.)
cp INCAR.pbe INCAR
$MPI_RUN $VASP_EXE > vasp_pbe.out 2>&1
mkdir -p pbe && cp OSZICAR OUTCAR vasprun.xml pbe/

# step 2: mBJ from the PBE WAVECAR, band path from KPOINTS_OPT
cp INCAR.mbj INCAR
$MPI_RUN $VASP_EXE > vasp_mbj.out 2>&1

echo "self-consistent c of mBJ:"
grep CMBJ OUTCAR | tail -1
grep -A12 "Band structure" OUTCAR
