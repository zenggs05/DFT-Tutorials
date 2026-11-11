#!/bin/bash
# run_error_cases.sh - run every deliberately broken input (and its fixed version) and print
# the Link that failed together with the lines just above "Error termination".
# All jobs are tiny (seconds); run them on an interactive/debug node, not on a login node.
#
# Usage:  bash run_error_cases.sh            (Gaussian environment must be loaded)
set -uo pipefail
command -v g16 >/dev/null 2>&1 || { echo "ERROR: g16 not found (source g16.profile first)"; exit 1; }
export GAUSS_SCRDIR=${GAUSS_SCRDIR:-/tmp/$USER/g16_errors_$$}
mkdir -p "$GAUSS_SCRDIR"

# the L101 example is stored as .txt so that editors/linters do not "repair" the missing blank line
cp err_l101_missing_blank_line.gjf.txt err_l101_missing_blank_line.gjf

for gjf in err_*.gjf; do
    log=${gjf%.gjf}.log
    g16 < "$gjf" > "$log" 2>&1
    echo "=================== $gjf"
    if tail -n 1 "$log" | grep -q "Normal termination"; then
        echo "Normal termination"
    else
        # the failing link, then the real message a few lines above it
        grep -m 1 "Error termination" "$log" || echo "(no 'Error termination' line: killed or crashed?)"
        grep -B 8 -m 1 "Error termination" "$log" | head -n 8
    fi
done
echo
echo "Explanations: tutorials/18-errors-capstone.md, or try ../../scripts/g16_error_doctor.py (-h for usage)."
