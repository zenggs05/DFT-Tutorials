#!/bin/bash
#SBATCH -J omx-gaas-cwf
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH -o slurm-%j.out
#
# GaAs closest Wannier functions: AO- and HO-guided runs, analysis, optional
# window scan that re-uses the SCF restart files (scf.restart pp).

# module load openmx/4.0.1          # adapt to your cluster (CWF needs OpenMX 4.0)
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-8}

mpirun -np ${NP} openmx GaAs_CWF_AO.dat -nt ${OMP_NUM_THREADS} > gaas_cwf_ao.std 2>&1
mpirun -np ${NP} openmx GaAs_CWF_HO.dat -nt ${OMP_NUM_THREADS} > gaas_cwf_ho.std 2>&1

python3 cwf_analyze.py gaas_cwf_ao
python3 cwf_analyze.py gaas_cwf_ho
python3 ../14-si-mlwf-sp3/band_compare.py gaas_cwf_ao.Band gaas_cwf_ao.CWF.Band --emin -13.8 --emax 1.0
python3 ../14-si-mlwf-sp3/band_compare.py gaas_cwf_ho.Band gaas_cwf_ho.CWF.Band --emin -13.8 --emax 1.0

# optional window scan (about one diagonalization per window)
# NP=${NP} NT=${OMP_NUM_THREADS} bash scan_cwf_window.sh
