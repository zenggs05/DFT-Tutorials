#!/bin/bash
#SBATCH -J omx-fe-fm
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 01:00:00
#SBATCH -o slurm-%j.out
#
# Runs the ferromagnetic calculation and a non-magnetic reference (fe_nm)
# with identical basis, cutoff, k-mesh and smearing, then prints both energies.

# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}

mpirun -np ${NP} openmx Fe_bcc_FM.dat -nt ${OMP_NUM_THREADS} > fe_fm.std 2>&1

# non-magnetic reference: spin polarization off, 8.0 + 8.0 electrons
sed -e 's/^System.Name .*/System.Name                 fe_nm/' \
    -e 's/^scf.SpinPolarization .*/scf.SpinPolarization        off/' \
    -e 's/10\.0  6\.0$/8.0  8.0/' Fe_bcc_FM.dat > Fe_bcc_NM.dat
mpirun -np ${NP} openmx Fe_bcc_NM.dat -nt ${OMP_NUM_THREADS} > fe_nm.std 2>&1

echo "== total energies (Hartree)"
grep "Utot\." fe_fm.out fe_nm.out
echo "== spin moments (muB)"
grep "Total spin moment" fe_fm.out
grep "by Voronoi charges" fe_fm.out
grep -A 6 "Mulliken populations" fe_fm.out
