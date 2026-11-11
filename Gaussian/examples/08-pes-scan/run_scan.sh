#!/bin/bash
#SBATCH --job-name=g16_scan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --output=%x-%j.out
#
# Chapter 8: run the PES scan examples one after another.
# Usage:  sbatch run_scan.sh                    (all five inputs)
#         sbatch run_scan.sh butane_dihedral_relaxed_scan   (only one input, no .gjf suffix)
# %nprocshared/%mem in the inputs: 8 cores / 8GB (4 / 4GB for the small H2O2 job).

set -u

# ---- site specific: load Gaussian 16 (edit for your cluster) ----
# module load gaussian/16-C.02
export g16root=${g16root:-/opt/g16root}
source "$g16root/g16/bsd/g16.profile"

export GAUSS_SCRDIR=${SCRATCH_BASE:-/tmp}/$USER/g16_${SLURM_JOB_ID:-$$}
mkdir -p "$GAUSS_SCRDIR"
trap 'rm -rf "$GAUSS_SCRDIR"' EXIT

cd "${SLURM_SUBMIT_DIR:-.}"

if [ $# -gt 0 ]; then
    jobs="$*"
else
    jobs="h2o2_rigid_zmatrix_scan butane_dihedral_relaxed_scan sn2_bond_scan_ts_guess two_dimensional_scan gic_bond_difference_scan"
fi

for name in $jobs; do
    echo "[$(date '+%F %T')] start $name"
    g16 < "$name.gjf" > "$name.log"
    if tail -n 1 "$name.log" | grep -q "Normal termination"; then
        echo "[$(date '+%F %T')] done  $name"
    else
        echo "WARNING: $name.log did not end with Normal termination" >&2
        grep "Error termination" "$name.log" >&2
    fi
    # relaxed scans: number of converged scan points
    if grep -q "ModRedundant\|AddGIC" "$name.gjf"; then
        echo "  converged scan points: $(grep -c 'Stationary point found' "$name.log")"
    fi
done
