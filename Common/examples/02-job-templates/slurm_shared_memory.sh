#!/bin/bash
# =====================================================================
# Generic Slurm template: shared-memory (single node, many threads),
# with a per-job scratch directory on the node-local disk.
# Typical users: Gaussian 16 (g16), ORCA-like serial+threads codes,
# Multiwfn, phonopy/Python post-processing.
#
# Rule of thumb for Gaussian:
#   %NProcShared (or %CPU) = --cpus-per-task
#   %Mem                   = about 75-85 % of --mem
# A full Gaussian wrapper is in Gaussian/scripts/sub_g16.sh.
#
#   sbatch slurm_shared_memory.sh
# =====================================================================
#SBATCH --job-name=g16
#SBATCH --partition=cpu
#SBATCH --nodes=1                    # shared memory: never more than 1 node
#SBATCH --ntasks=1                   # one process ...
#SBATCH --cpus-per-task=16           # ... with 16 threads/cores
#SBATCH --mem=40G                    # memory for the whole job
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out

# ----------------------------- user settings -------------------------
INPUT="benzene.gjf"                  # Gaussian input in the submit dir
BASENAME="${INPUT%.*}"

# ----------------------------- environment ---------------------------
module purge
# module load gaussian/g16c02          # or set g16root as below
# export g16root=/opt/gaussian
# source "$g16root/g16/bsd/g16.profile"

ulimit -s unlimited
cd "${SLURM_SUBMIT_DIR}" || exit 1

# Per-job scratch on the local disk (fast, not shared). Use the site's
# variable if it exists ($TMPDIR, $SLURM_TMPDIR, /scratch/$USER ...).
SCR_BASE=${TMPDIR:-/tmp}
export GAUSS_SCRDIR="${SCR_BASE}/${USER}/g16_${SLURM_JOB_ID}"
mkdir -p "${GAUSS_SCRDIR}" || { echo "cannot create ${GAUSS_SCRDIR}" >&2; exit 1; }
# remove the scratch directory whatever happens (normal end, error, scancel)
trap 'rm -rf "${GAUSS_SCRDIR}"' EXIT

# Consistency check between the input and the allocation
nproc_in=$(grep -i -m1 -E '^%nprocshared' "${INPUT}" | cut -d= -f2 | tr -d '[:space:]')
if [ -n "${nproc_in}" ] && [ "${nproc_in}" != "${SLURM_CPUS_PER_TASK}" ]; then
    echo "WARNING: %NProcShared=${nproc_in} but --cpus-per-task=${SLURM_CPUS_PER_TASK}" >&2
fi

echo "Job ${SLURM_JOB_ID} on $(hostname), ${SLURM_CPUS_PER_TASK} cores, scratch ${GAUSS_SCRDIR}"
df -h "${SCR_BASE}" | tail -n 1

g16 < "${INPUT}" > "${BASENAME}.log" 2>&1
status=$?

if grep -q "Normal termination" "${BASENAME}.log"; then
    echo "Normal termination"
else
    echo "Gaussian did NOT terminate normally (exit ${status}):" >&2
    grep -B2 -A2 "Error termination" "${BASENAME}.log" | tail -n 8 >&2
fi
exit ${status}
