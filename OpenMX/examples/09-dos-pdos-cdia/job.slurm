#!/bin/bash
#SBATCH -J omx_dos
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH -o slurm-%j.out

# DOS/PDOS of diamond with OpenMX 4.0 (edit modules/paths).
# DosMain is built with "make DosMain" in openmx4.0/source.
# module load openmx/4.0
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mpirun -np ${SLURM_NTASKS} openmx Cdia_dos.dat -nt ${OMP_NUM_THREADS} > cdia_dos.std

# interactive answers stored in files ...
DosMain cdia_dos.Dos.val cdia_dos.Dos.vec < dosmain_answers_dos.txt > dosmain_dos.log
DosMain cdia_dos.Dos.val cdia_dos.Dos.vec < dosmain_answers.txt     > dosmain_pdos.log
# ... or the same with the generic driver:  bash run_dosmain.sh cdia_dos tetra

python3 sum_pdos.py cdia_dos --dat Cdia_dos.dat --plot --emin -28 --emax 15 > sum_pdos.txt
