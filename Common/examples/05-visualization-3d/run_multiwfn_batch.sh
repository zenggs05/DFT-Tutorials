#!/bin/bash
# run_multiwfn_batch.sh - apply one Multiwfn command stream to many wavefunction files
#
# Usage:
#   bash run_multiwfn_batch.sh <cmds.txt> <file1> [file2 ...]
#   bash run_multiwfn_batch.sh ../../../Gaussian/examples/16-postprocess-multiwfn/multiwfn_esp_surface_cmds.txt mols/*.fchk
#   NTHREADS=16 bash run_multiwfn_batch.sh my_cmds.txt *.molden
#
# Inputs: anything Multiwfn can load (.fchk/.fch, .wfx/.wfn, .molden, .mwfn, .cube/.cub,
# CP2K molden files, ...). Gaussian .chk files are converted with formchk first.
# Each input gets its own folder <name>_mwfn/ containing the screen output and all
# files Multiwfn wrote (Multiwfn uses fixed output names such as surfanalysis.txt or
# *.cub, so runs must not share a directory).
#
# Command streams: one answer per line, exactly the keystrokes of an interactive
# session. Write them by doing the analysis once interactively. Ready-made, tested
# streams for Gaussian are in Gaussian/examples/16-postprocess-multiwfn/.
set -uo pipefail

CMDS=${1:-}
shift || true
[ -n "$CMDS" ] && [ -f "$CMDS" ] && [ $# -gt 0 ] || { sed -n '2,17p' "$0"; exit 1; }
command -v Multiwfn >/dev/null 2>&1 || { echo "ERROR: Multiwfn not in PATH (use the noGUI build on clusters)"; exit 1; }
CMDS=$(cd "$(dirname "$CMDS")" && pwd)/$(basename "$CMDS")

# Stack settings: without them large systems end in "Segmentation fault"
ulimit -s unlimited 2>/dev/null || echo "WARNING: could not set ulimit -s unlimited"
export OMP_STACKSIZE=${OMP_STACKSIZE:-200M}
NTHREADS=${NTHREADS:-${SLURM_CPUS_PER_TASK:-4}}

# Private settings.ini with the requested number of threads (Multiwfn reads
# settings.ini from the current directory first; check the thread count it prints)
SETTINGS_SRC=${Multiwfnpath:-$(dirname "$(command -v Multiwfn)")}/settings.ini

for f in "$@"; do
    [ -f "$f" ] || { echo "skip (not found): $f"; continue; }
    name=$(basename "$f"); stem=${name%.*}
    absf=$(cd "$(dirname "$f")" && pwd)/$name
    if [ "${name##*.}" = "chk" ]; then
        command -v formchk >/dev/null 2>&1 || { echo "skip $f: formchk not found"; continue; }
        formchk "$absf" "${absf%.chk}.fchk" > /dev/null || { echo "formchk failed: $f"; continue; }
        absf=${absf%.chk}.fchk
    fi
    out="${stem}_mwfn"
    mkdir -p "$out"
    (
        cd "$out" || exit 1
        if [ -f "$SETTINGS_SRC" ]; then
            sed "s/nthreads= *[0-9]*/nthreads= $NTHREADS/" "$SETTINGS_SRC" > settings.ini
        fi
        Multiwfn "$absf" -silent < "$CMDS" > "${stem}.out" 2>&1
    )
    if grep -qi "segmentation fault\|error" "$out/${stem}.out"; then
        status="CHECK"
    else
        status="OK"
    fi
    printf "%-30s %-6s %s\n" "$name" "$status" "$(ls "$out" | grep -v -e '^settings.ini$' | tr '\n' ' ')"
done
echo "Note: a final 'end of file during read' type message only means the command stream ended."
