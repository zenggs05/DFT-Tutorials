#!/bin/bash
# run_workflow_slurm.sh - run the whole Diels-Alder workflow on one node, stopping at the first
# failed check:  minima opt+freq -> TS opt+freq (exactly 1 imaginary frequency) -> IRC
#                -> high-level single points -> GoodVibes.
# The .gjf files carry no %nprocshared/%mem: cores and memory come from GAUSS_CDEF/GAUSS_MDEF.
# CREST (crest_conformers.sh) is a separate, optional first step.
#
# Usage:  sbatch run_workflow_slurm.sh
#
#SBATCH --job-name=da-workflow
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}" || exit 1

# module load gaussian/16-C.02
export g16root=${g16root:-/opt}
# shellcheck disable=SC1091
source "$g16root/g16/bsd/g16.profile" || exit 1

CPULIST=$(awk '/^Cpus_allowed_list/ {print $2}' /proc/self/status)
if [ -n "$CPULIST" ]; then export GAUSS_CDEF=$CPULIST; else export GAUSS_PDEF=${SLURM_CPUS_PER_TASK:-1}; fi
export GAUSS_MDEF="$(( ${SLURM_MEM_PER_NODE:-32768} * 80 / 100 ))MB"
export GAUSS_SCRDIR="${SCR_BASE:-/tmp}/$USER/g16_${SLURM_JOB_ID:-$$}"
mkdir -p "$GAUSS_SCRDIR" || exit 1
trap 'rm -rf -- "${GAUSS_SCRDIR:?}"' EXIT
echo "GAUSS_CDEF=${GAUSS_CDEF:-} GAUSS_PDEF=${GAUSS_PDEF:-} GAUSS_MDEF=$GAUSS_MDEF"

run() {                                   # run NAME.gjf, stop the workflow if it fails
    local name=$1
    if [ -f "$name.log" ] && tail -n 1 "$name.log" | grep -q "Normal termination"; then
        echo "skip $name (already finished)"; return 0
    fi
    echo "run  $name  $(date '+%F %T')"
    g16 < "$name.gjf" > "$name.log"
    if ! tail -n 1 "$name.log" | grep -q "Normal termination"; then
        echo "FAILED: $name"; grep -B 6 "Error termination" "$name.log" | tail -n 12; exit 1
    fi
}
nimag() { grep "imaginary frequencies" "$1" | tail -n 1 | awk '{print $2}' | grep -E '^[0-9]+$' || echo 0; }

# 1) minima
for m in reactant_opt reactant2_opt product_opt; do
    run "$m"
    [ "$(nimag "$m.log")" -eq 0 ] || { echo "STOP: $m.log is not a minimum (imaginary frequencies)"; exit 1; }
done

# 2) transition state
run ts_opt
n=$(nimag ts_opt.log)
[ "$n" -eq 1 ] || { echo "STOP: ts_opt.log has $n imaginary frequencies, expected 1"; exit 1; }
grep -m 1 "Frequencies --" ts_opt.log      # the first value should be negative (imaginary)

# 3) IRC (reads ts_opt.chk); check both ends afterwards in GaussView or with log2xyz.py
run irc
grep -A 4 "Summary of reaction path following" irc.log | head -n 6

# 4) high-level single points <root>_TZ.gjf, generated from sp_highlevel.gjf
bash make_sp_inputs.sh
for r in reactant_opt reactant2_opt ts_opt product_opt; do run "${r}_TZ"; done

# 5) thermochemistry (qRRHO, 1 M) and the profile, if GoodVibes is installed
if command -v goodvibes >/dev/null 2>&1; then
    bash goodvibes_run.sh
else
    echo "GoodVibes not found: run 'bash goodvibes_run.sh' where it is installed"
fi
