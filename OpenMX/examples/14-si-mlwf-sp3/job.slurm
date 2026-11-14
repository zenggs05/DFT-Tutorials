#!/bin/bash
#SBATCH -J omx-si-mlwf
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH -o slurm-%j.out
#
# Si MLWFs: 8 sp3-type functions (disentangled) and 4 bond-centred valence
# functions, then the Wannier90 interface run. Each input does its own SCF.

# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-8}

mpirun -np ${NP} openmx Si_MLWF_sp3.dat   -nt ${OMP_NUM_THREADS} > si_mlwf_sp3.std   2>&1
mpirun -np ${NP} openmx Si_MLWF_bond4.dat -nt ${OMP_NUM_THREADS} > si_mlwf_bond4.std 2>&1
mpirun -np ${NP} openmx Si_MLWF_W90.dat   -nt ${OMP_NUM_THREADS} > si_w90.std        2>&1

# quick checks (python3 with numpy is needed for band_compare.py)
python3 wf_centers.py si_mlwf_sp3.std
python3 wf_centers.py si_mlwf_bond4.std
python3 band_compare.py si_mlwf_sp3.Band   si_mlwf_sp3.Wannier_Band   --emin -14 --emax 1
python3 band_compare.py si_mlwf_bond4.Band si_mlwf_bond4.Wannier_Band --emin -14 --emax 0

# optional: Wannier90 (needs wannier90.x in PATH)
if command -v wannier90.x > /dev/null && [ -f si_w90.win ]; then
  cat w90_extra.win >> si_w90.win
  wannier90.x si_w90
fi
