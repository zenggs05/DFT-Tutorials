#!/bin/bash
# Slurm job for the chapter-12 AIMD examples (one input per job).
# Usage:
#   jid=$(sbatch --parsable job_slurm.sh water32_nvt_csvr.inp)
#   sbatch --dependency=afterok:$jid job_slurm.sh md_restart.inp
#   sbatch --dependency=afterok:$jid job_slurm.sh water32_nve.inp
#   sbatch --dependency=afterok:$jid job_slurm.sh water32_npt_i.inp
#   (optional) water32_nvt_nose.inp, reftraj_eval.inp, reftraj_npt.inp
#SBATCH --job-name=cp2k-aimd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.log

set -euo pipefail

INP=${1:-water32_nvt_csvr.inp}
OUT=${INP%.inp}.out

# --- environment: adapt ONE of the following lines to your cluster ---------
# module load cp2k/2026.2
# source /path/to/cp2k/tools/toolchain/install/setup    # toolchain build
# source /path/to/cp2k/install/cp2k_env                 # build_cp2k.sh build
# conda activate cp2k_env                               # conda-forge build
# ---------------------------------------------------------------------------

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export CP2K_DATA_DIR=${CP2K_DATA_DIR:-/path/to/cp2k/data}

# never overwrite the output of a previous attempt
if [ -e "$OUT" ]; then
    mv "$OUT" "${OUT}.$(date +%Y%m%d%H%M%S)"
fi

srun cp2k.psmp -i "$INP" -o "$OUT"

if grep -q "PROGRAM ENDED AT" "$OUT"; then
    echo "finished: $OUT"
else
    echo "CP2K stopped early (walltime or SCF failure?): restart from the"
    echo "latest <PROJECT>-1.restart with &EXT_RESTART, see md_restart.inp" >&2
    exit 1
fi
