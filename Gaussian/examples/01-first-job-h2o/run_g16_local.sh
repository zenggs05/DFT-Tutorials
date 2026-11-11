#!/usr/bin/env bash
# run_g16_local.sh - run Gaussian 16 inputs one after another on a local
# Linux/WSL machine and report whether each job terminated normally.
#
# Usage:
#   bash run_g16_local.sh                 # run every *.gjf in this directory
#   bash run_g16_local.sh a.gjf b.gjf     # run only the given inputs
#
# Requirements: g16 in PATH (source $g16root/g16/bsd/g16.profile first).
# Do NOT run long jobs on a shared login node; use run_g16_slurm.sh there.

set -u

# ---------------------------------------------------------------- checks
if ! command -v g16 >/dev/null 2>&1; then
    echo "ERROR: g16 not found in PATH." >&2
    echo "       export g16root=...; source \$g16root/g16/bsd/g16.profile" >&2
    exit 1
fi

# scratch directory: keep the user's setting, otherwise use /tmp
if [ -z "${GAUSS_SCRDIR:-}" ]; then
    export GAUSS_SCRDIR="/tmp/${USER:-g16user}/g16scr"
fi
if ! mkdir -p "$GAUSS_SCRDIR"; then
    echo "ERROR: cannot create GAUSS_SCRDIR=$GAUSS_SCRDIR" >&2
    exit 1
fi
echo "GAUSS_SCRDIR = $GAUSS_SCRDIR"

# ---------------------------------------------------------------- inputs
if [ "$#" -gt 0 ]; then
    inputs=("$@")
else
    shopt -s nullglob
    inputs=(*.gjf *.com)
    shopt -u nullglob
fi
if [ "${#inputs[@]}" -eq 0 ]; then
    echo "No .gjf/.com files found." >&2
    exit 1
fi

n_ok=0
n_bad=0
for inp in "${inputs[@]}"; do
    if [ ! -f "$inp" ]; then
        echo "skip: $inp does not exist"
        continue
    fi

    # Windows CR characters break blank-line detection in Gaussian input
    if grep -q $'\r' "$inp"; then
        sed -i 's/\r$//' "$inp"
        echo "note: converted CRLF line endings to LF in $inp"
    fi
    # Gaussian input must end with a blank line
    if [ -n "$(tail -c 1 "$inp")" ] || [ -n "$(tail -n 1 "$inp")" ]; then
        echo "warning: $inp does not end with a blank line (expect 'End of file in ZSymb')"
    fi

    base="${inp%.*}"
    log="$base.log"
    echo "=== $(date '+%F %T')  running $inp -> $log"
    t0=$(date +%s)
    g16 < "$inp" > "$log" 2>&1
    t1=$(date +%s)

    n_norm=$(grep -c "Normal termination" "$log")
    if grep -q "Error termination" "$log"; then
        status="FAILED"
        n_bad=$((n_bad + 1))
    elif tail -n 3 "$log" | grep -q "Normal termination"; then
        status="OK"
        n_ok=$((n_ok + 1))
    else
        status="INCOMPLETE"
        n_bad=$((n_bad + 1))
    fi

    echo "    status       : $status  (normal-termination lines: $n_norm, wall: $((t1 - t0)) s)"
    last_e=$(awk '/SCF Done/{e=$5} END{print e}' "$log")
    [ -n "$last_e" ] && echo "    last SCF E   : $last_e Hartree"
    if grep -q "Stationary point found" "$log"; then
        echo "    optimization : stationary point found after $(grep -c 'SCF Done' "$log") SCF calculations"
    fi
    if [ "$status" = "FAILED" ]; then
        echo "    failing link : $(grep 'Error termination' "$log" | tail -n 1 | sed 's/.*\/\(l[0-9]*\.exe\).*/\1/')"
        echo "    lines above the error message:"
        grep -B 6 "Error termination" "$log" | head -n 6 | sed 's/^/      | /'
    fi
done

echo "=== finished: $n_ok OK, $n_bad failed/incomplete"
[ "$n_bad" -eq 0 ]
