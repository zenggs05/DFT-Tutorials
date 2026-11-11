#!/bin/bash
#SBATCH --job-name=neb-pt100
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster
#
# Run in the NEB root directory (INCAR, KPOINTS, POTCAR and 00 .. NN/POSCAR).
# The total number of MPI ranks must be a multiple of IMAGES: here
# 3 nodes x 32 ranks = 96 ranks = 3 images x 32 ranks per image.
# The executable must be compiled with the VTST code (IOPT, LCLIMB, ICHAIN).

module purge
# module load vasp/6.5.1-vtst      # edit for your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
NP=${SLURM_NTASKS:-96}

IMAGES=$(grep -E "^[[:space:]]*IMAGES" INCAR | head -1 | sed 's/.*=//; s/[#!].*//' | tr -d ' ')
[ -n "$IMAGES" ] || { echo "IMAGES not set in INCAR"; exit 1; }
if [ $((NP % IMAGES)) -ne 0 ]; then
    echo "ERROR: $NP MPI ranks cannot be split into $IMAGES images"; exit 1
fi
for f in INCAR KPOINTS POTCAR; do [ -f $f ] || { echo "missing $f"; exit 1; }; done

mpirun -np $NP $VASP_EXE > vasp.out 2>&1

# quick summary (VTST scripts must be in PATH)
command -v nebef.pl > /dev/null && nebef.pl
