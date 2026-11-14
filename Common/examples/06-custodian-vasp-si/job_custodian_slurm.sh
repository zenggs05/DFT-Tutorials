#!/bin/bash
#SBATCH --job-name=si-custodian
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --partition=cpu
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

# custodian runs on the first node of the allocation and launches VASP itself.
# Adjust module names and the conda environment to your cluster.

set -eo pipefail     # no -u: module and conda activation scripts use unset variables
ulimit -s unlimited
export OMP_NUM_THREADS=1

module purge
module load intel-oneapi vasp/6.5.1          # example module names

source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate wf                            # env with custodian (Python >= 3.11)

cd "$SLURM_SUBMIT_DIR"

# POTCAR is never stored in the repository: build it from your licensed library
if [ ! -f POTCAR ]; then
    : "${VASP_PP_PATH:?POTCAR missing and VASP_PP_PATH is not set}"
    cat "$VASP_PP_PATH/potpaw_PBE/Si/POTCAR" > POTCAR
fi

# Remove Windows line endings that would make VASP ignore INCAR tags
command -v dos2unix >/dev/null && dos2unix -q INCAR KPOINTS POSCAR

export VASP_CMD="srun vasp_std"              # or: mpirun -np $SLURM_NTASKS vasp_std
# Wall time in seconds, keep it equal to --time above. Slurm does not export the
# limit to the job, so custodian's WalltimeHandler needs it explicitly.
export CUSTODIAN_WALLTIME=7200
python custodian_vasp_run.py

# Alternative: YAML driven run
#   cstdn run custodian_spec.yaml
