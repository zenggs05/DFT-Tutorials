#!/bin/bash
#SBATCH -J omx-bi2se3
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 04:00:00
#SBATCH -o slurm-%j.out
#
# 1) bands with SOC  2) bands without SOC (same j-dependent VPS, SOC off)
# 3) SCF with HS.fileout on  4) Z2 invariant with Z2FH
# Z2FH is built once:  cd openmx4.0/source && make Z2FH   (copied to ../work)

# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}
Z2FH=${Z2FH:-$HOME/openmx4.0/work/Z2FH}

mpirun -np ${NP} openmx Bi2Se3_SOC_band.dat -nt ${OMP_NUM_THREADS} > bi2se3_soc.std 2>&1

sed -e 's/^System.Name .*/System.Name                 bi2se3_nosoc/' \
    -e 's/^scf.SpinOrbit.Coupling .*/scf.SpinOrbit.Coupling      off/' \
    Bi2Se3_SOC_band.dat > Bi2Se3_noSOC_band.dat
mpirun -np ${NP} openmx Bi2Se3_noSOC_band.dat -nt ${OMP_NUM_THREADS} > bi2se3_nosoc.std 2>&1

python3 ../12-gaas-soc-band/soc_splitting.py bi2se3_soc.Band   --at G Z F L
python3 ../12-gaas-soc-band/soc_splitting.py bi2se3_nosoc.Band --at G Z F L

mpirun -np ${NP} openmx Bi2Se3_Z2.dat -nt ${OMP_NUM_THREADS} > bi2se3_z2.std 2>&1
if [ -x "${Z2FH}" ]; then
  mpirun -np ${NP} ${Z2FH} bi2se3_z2.scfout < Z2FH.in > Z2FH.out 2>&1
  cat Z2.dat
fi
