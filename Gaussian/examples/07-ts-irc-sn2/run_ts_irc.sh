#!/bin/bash
#SBATCH --job-name=sn2_ts_irc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#
# Chapter 7 workflow: TS (Opt+Freq) -> check 1 imaginary frequency -> IRC
#                     -> IRC endpoint Opt+Freq -> separated reactants.
# %nprocshared=8 and %mem=8GB in the .gjf files match the request above
# (Gaussian %mem is kept at ~80% of the Slurm --mem).

set -u

# ---- site specific: load Gaussian 16 (edit for your cluster) ----
# module load gaussian/16-C.02
export g16root=${g16root:-/opt/g16root}
source "$g16root/g16/bsd/g16.profile"

# ---- private scratch directory, removed when the job ends ----
export GAUSS_SCRDIR=${SCRATCH_BASE:-/tmp}/$USER/g16_${SLURM_JOB_ID:-$$}
mkdir -p "$GAUSS_SCRDIR"
trap 'rm -rf "$GAUSS_SCRDIR"' EXIT

cd "${SLURM_SUBMIT_DIR:-.}"

run_g16() {
    local name=$1
    echo "[$(date '+%F %T')] start $name"
    g16 < "$name.gjf" > "$name.log"
    if ! tail -n 1 "$name.log" | grep -q "Normal termination"; then
        echo "ERROR: $name.log did not end with Normal termination" >&2
        grep "Error termination" "$name.log" >&2
        exit 1
    fi
    echo "[$(date '+%F %T')] done  $name"
}

# 1) TS optimization + frequencies
run_g16 sn2_cl_ch3br_ts_calcfc

nimag=$(grep "imaginary frequencies" sn2_cl_ch3br_ts_calcfc.log | tail -n 1 | awk '{print $2}')
nimag=${nimag:-0}
echo "Number of imaginary frequencies at the TS: $nimag"
if [ "$nimag" != "1" ]; then
    echo "ERROR: expected exactly 1 imaginary frequency, found $nimag. Stop before IRC." >&2
    exit 1
fi

# 2) IRC in both directions (reads force constants from sn2_ts.chk)
run_g16 sn2_irc_both_directions

# 3) Optimize the IRC endpoints (replace coordinates with your own IRC endpoints first!)
run_g16 irc_endpoint_opt_freq

# 4) Separated reactants (reference state for the barrier)
run_g16 sn2_separated_reactants

echo "Free energies (Hartree):"
grep "Sum of electronic and thermal Free Energies" ./*.log
