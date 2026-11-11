#!/bin/bash
# Build a POTCAR in each given directory from the element line (line 6) of its POSCAR.
# The repository never ships POTCAR files (VASP licence); they come from your own
# PAW library pointed to by VASP_PP_PATH.
#
# Usage:
#   export VASP_PP_PATH=/path/to/your/paw/library
#   bash make_potcar.sh 1_slab 2_co_gas 3_co_top
# Options (environment variables):
#   POTDIR=potpaw_PBE.64   sub-directory of the PBE set (default: potpaw_PBE)
#   POT_Li=Li_sv           choose a variant for an element (POT_<Element>=<dir>)
set -e
: "${VASP_PP_PATH:?please export VASP_PP_PATH=/path/to/your/PAW/library}"
PP="$VASP_PP_PATH/${POTDIR:-potpaw_PBE}"
[ $# -gt 0 ] || { echo "usage: bash make_potcar.sh dir1 [dir2 ...]"; exit 1; }
for d in "$@"; do
    [ -f "$d/POSCAR" ] || { echo "skip $d: no POSCAR"; continue; }
    els=$(sed -n 6p "$d/POSCAR")
    : > "$d/POTCAR"
    for e in $els; do
        var="POT_$e"
        sub="${!var:-$e}"
        if [ ! -f "$PP/$sub/POTCAR" ]; then
            echo "ERROR: $PP/$sub/POTCAR not found"; exit 1
        fi
        cat "$PP/$sub/POTCAR" >> "$d/POTCAR"
    done
    echo "== $d: POSCAR elements: $els"
    grep -E "TITEL" "$d/POTCAR" | sed 's/^/   /'
done
