#!/bin/bash
#SBATCH --job-name=g16_solv
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#
# Chapter 9: implicit solvation examples (SMD / IEFPCM).
# Usage:  sbatch run_solvation.sh                       (all inputs)
#         sbatch run_solvation.sh phenol_smd_water_opt_freq   (selected inputs, no .gjf suffix)

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
    jobs="phenol_smd_water_opt_freq phenoxide_smd_water_opt_freq solvation_free_energy_gas_and_smd custom_solvent_generic_read redox_potential_cycle_neutral_and_cation"
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
done

echo "---- SCF energies (Hartree) ----"
grep "SCF Done" ./*.log | awk '{print $1, $6}'
echo "---- Free energies (Hartree, 1 atm thermal corrections) ----"
grep "Sum of electronic and thermal Free Energies" ./*.log
