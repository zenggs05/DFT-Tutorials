#!/bin/bash
# PBS/Torque job script for a CPU build of VASP (pure MPI).
# Torque resource syntax is active; the PBS Pro equivalent is given below it.
#PBS -N vasp_cpu
#PBS -q workq
#PBS -l nodes=1:ppn=32
## PBS Pro: replace the line above by
## #PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -l walltime=24:00:00
#PBS -j oe

module purge
module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl
export PATH=$HOME/software/vasp.6.6.1/bin:$PATH
VASP_EXE=${VASP_EXE:-vasp_std}

cd "$PBS_O_WORKDIR" || exit 1
ulimit -s unlimited
export OMP_NUM_THREADS=1

for f in INCAR POSCAR POTCAR; do
    [ -s "$f" ] || { echo "ERROR: $f is missing or empty"; exit 1; }
done

NP=$(wc -l < "$PBS_NODEFILE")
echo "Running $VASP_EXE on $NP cores, start $(date)"
mpirun -np "$NP" -machinefile "$PBS_NODEFILE" "$VASP_EXE" > vasp.out 2>&1
echo "End $(date)"

grep -q "General timing and accounting" OUTCAR 2>/dev/null \
    && echo "VASP finished normally" \
    || echo "VASP did NOT finish normally, check vasp.out and OUTCAR"
