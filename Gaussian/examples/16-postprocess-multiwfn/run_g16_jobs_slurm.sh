#!/bin/bash
# run_g16_jobs_slurm.sh - run the Gaussian jobs of this example in order on one node.
# Usage:  sbatch run_g16_jobs_slurm.sh
# %nprocshared=8 / %mem=8GB in the inputs match --cpus-per-task=8 / --mem=10G below.
#SBATCH --job-name=g16-ch16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out

set -u
# module load gaussian/16-C.02          # site specific: use module OR source the profile
export g16root=${g16root:-/opt}
source "$g16root/g16/bsd/g16.profile"

export GAUSS_SCRDIR=${SCR_BASE:-/tmp}/$USER/g16_${SLURM_JOB_ID:-$$}
mkdir -p "$GAUSS_SCRDIR"
trap 'rm -rf "$GAUSS_SCRDIR"' EXIT
cd "${SLURM_SUBMIT_DIR:-.}"

run() {
    g16 < "$1.gjf" > "$1.log"
    if tail -n 1 "$1.log" | grep -q "Normal termination"; then
        echo "$1: Normal termination"
    else
        echo "$1: FAILED, see $1.log" >&2
        exit 1
    fi
}
run phenol_opt_sp          # must finish before phenol_td (it reads phenol.chk)
run phenol_td
run water_dimer
run wfx_output_example

for c in phenol phenol_td water_dimer; do formchk "$c.chk" "$c.fchk"; done
grep -c "imaginary frequencies" phenol_opt_sp.log water_dimer.log || true
