#!/bin/bash
# Slurm job for the chapter-10 examples.
# Usage:  sbatch job_slurm.sh h2o_pdos_2026.inp
#         sbatch job_slurm.sh si_band_structure.inp
#SBATCH --job-name=cp2k-dos
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.log

set -euo pipefail

INP=${1:-h2o_pdos_2026.inp}
OUT=${INP%.inp}.out

# --- environment: adapt ONE of the following lines to your cluster ---------
# module load cp2k/2026.2
# source /path/to/cp2k/tools/toolchain/install/setup    # toolchain build
# source /path/to/cp2k/install/cp2k_env                 # build_cp2k.sh build
# conda activate cp2k_env                               # conda-forge build
# ---------------------------------------------------------------------------

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export CP2K_DATA_DIR=${CP2K_DATA_DIR:-/path/to/cp2k/data}

# the band-structure file is opened in APPEND mode: start from a clean file
if grep -q "FILE_NAME" "$INP"; then
    BS=$(awk 'toupper($1)=="FILE_NAME"{print $2}' "$INP" | head -n 1)
    [ -n "$BS" ] && rm -f "$BS"
fi

srun cp2k.psmp -i "$INP" -o "$OUT"

if grep -q "PROGRAM ENDED AT" "$OUT"; then
    echo "finished: $OUT"
    ls -1 ./*.dos ./*.pdos ./*.bs 2>/dev/null || true
else
    echo "CP2K did not finish normally, check $OUT" >&2
    exit 1
fi
