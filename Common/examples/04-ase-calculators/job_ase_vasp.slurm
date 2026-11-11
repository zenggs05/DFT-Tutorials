#!/bin/bash
#SBATCH --job-name=ase-vasp-si
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --partition=cpu
#SBATCH --output=%x-%j.out

# Run an ASE driver script inside a Slurm allocation.
# The Python process itself is serial; ASE launches VASP with the MPI
# command stored in ASE_VASP_COMMAND for every energy/force evaluation.

set -euo pipefail
ulimit -s unlimited
export OMP_NUM_THREADS=1

# --- edit to match your cluster -------------------------------------------
module purge
module load intel-oneapi vasp/6.5.1          # example module names
source ~/miniconda3/etc/profile.d/conda.sh
conda activate dft                             # env with ase >= 3.23
export VASP_PP_PATH=$HOME/pp                   # contains potpaw_PBE/ (or potpaw_PBE.64/)
# ---------------------------------------------------------------------------

export ASE_VASP_COMMAND="srun --mpi=pmix -n ${SLURM_NTASKS} vasp_std"
# On clusters without srun integration use instead:
# export ASE_VASP_COMMAND="mpirun -np ${SLURM_NTASKS} vasp_std"

python ase_vasp_relax.py --mode vasp

# Quick success check (the script also prints it)
grep -H "reached required accuracy" si_relax_vasp/OUTCAR || echo "WARNING: relaxation not converged"
