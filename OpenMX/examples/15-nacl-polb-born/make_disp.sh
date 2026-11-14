#!/bin/bash
# Generate displaced NaCl inputs from NaCl_polB.dat.
# Na (atom 1) is moved along Cartesian x by each value in DISPS (Ang).
# Usage:  bash make_disp.sh                 # default: -0.05 0.00 0.05
#         DISPS="-0.02 0.02" bash make_disp.sh
set -euo pipefail
DISPS=${DISPS:-"-0.05 0.00 0.05"}
TEMPLATE=${TEMPLATE:-NaCl_polB.dat}

for u in $DISPS; do
    # tag: -0.05 -> m0.05, 0.05 -> p0.05, 0.00 -> p0.00
    tag=$(awk -v u="$u" 'BEGIN{ if (u<0) printf "m%.2f", -u; else printf "p%.2f", u }')
    out="NaCl_dx_${tag}.dat"
    awk -v u="$u" -v name="nacl_dx_${tag}" '
        tolower($1)=="system.name" { printf "System.Name                 %s\n", name; next }
        $1=="1" && $2=="Na" && NF>=7 {
            printf "  1  Na  %9.5f  %9.5f  %9.5f   %s  %s\n", $3+u, $4, $5, $6, $7; next }
        { print }
    ' "$TEMPLATE" > "$out"
    echo "wrote $out (Na x shifted by $u Ang)"
done
