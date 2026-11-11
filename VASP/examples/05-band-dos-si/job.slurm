#!/bin/bash
#SBATCH --job-name=si-band-dos
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

# Runs all four steps in sequence. Submit from VASP/examples/05-band-dos-si
# after building ./POTCAR (see README.md).
module purge
# module load vasp/6.6.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}          # KPOINTS_OPT does not work with vasp_gam
MPI_RUN="mpirun -np ${SLURM_NTASKS:-8}"
TOP=$(pwd)

[ -f POTCAR ] || { echo "POTCAR missing in $TOP"; exit 1; }

# step 1: SCF
cd "$TOP/1-scf"
cp "$TOP/POTCAR" .
$MPI_RUN $VASP_EXE > vasp.out 2>&1
grep "E-fermi" OUTCAR | tail -1
[ -s CHGCAR ] || { echo "SCF did not write CHGCAR"; exit 1; }

# step 2: band structure, ICHARG=11
cd "$TOP/2-band-icharg11"
cp "$TOP/POTCAR" "$TOP/1-scf/POSCAR" "$TOP/1-scf/CHGCAR" .
$MPI_RUN $VASP_EXE > vasp.out 2>&1

# step 3: DOS on a dense mesh, ICHARG=11
cd "$TOP/3-dos"
cp "$TOP/POTCAR" "$TOP/1-scf/POSCAR" "$TOP/1-scf/CHGCAR" .
$MPI_RUN $VASP_EXE > vasp.out 2>&1

# step 4: SCF + band path in one run (KPOINTS_OPT, VASP >= 6.3)
cd "$TOP/4-band-kpoints-opt"
cp "$TOP/POTCAR" "$TOP/1-scf/POSCAR" .
$MPI_RUN $VASP_EXE > vasp.out 2>&1

cd "$TOP"
echo "done. Fermi energy of the SCF run:"
grep "E-fermi" 1-scf/OUTCAR | tail -1
