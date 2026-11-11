#!/bin/bash
#SBATCH -J omx-fe-jx
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH -o slurm-%j.out
#
# 1) collinear SCF with HS.fileout on  2) jx (MPI over k points)  3) mean-field Tc
# jx must be compiled once:  cd openmx4.0/source && make jx   (copied to ../work)

# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}
JX=${JX:-$HOME/openmx4.0/work/jx}     # path to the jx executable

mpirun -np ${NP} openmx Fe_bcc_jx.dat -nt ${OMP_NUM_THREADS} > fe_jx.std 2>&1
grep "Total spin moment" fe_jx.out
ls -l fe_jx.scfout

mpirun -np ${NP} ${JX} fe_jx.scfout Fe_bcc_jx.config > jx.log 2>&1
python3 jx_tc_mf.py jx.log Fe_bcc_jx.dat
