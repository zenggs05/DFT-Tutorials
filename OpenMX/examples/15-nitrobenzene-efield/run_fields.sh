#!/bin/bash
# Generate (and optionally run) nitrobenzene inputs for several field strengths.
# Field values are in GV/m along the a axis; the zero-field run uses the template.
#
# Usage:
#   bash run_fields.sh                  # only write the inputs
#   RUN=1 NP=8 bash run_fields.sh       # write and run them one after another
#   FIELDS="-1 1" bash run_fields.sh    # custom list
set -euo pipefail

FIELDS=${FIELDS:-"-10 -2 2 10"}
TEMPLATE=${TEMPLATE:-Nitro_Benzene_EF.dat}
OPENMX=${OPENMX:-openmx}
NP=${NP:-4}
NT=${NT:-1}
RUN=${RUN:-0}

tag_of () {   # -10 -> m10, 2 -> p2, 0.5 -> p0.5
    awk -v f="$1" 'BEGIN{ if (f<0) printf "m%g", -f; else printf "p%g", f }'
}

inputs=("$TEMPLATE")
for f in $FIELDS; do
    tag=$(tag_of "$f")
    out="Nitro_Benzene_EF_${tag}.dat"
    awk -v f="$f" -v name="nben_${tag}" '
        tolower($1)=="system.name"        { printf "System.Name                 %s\n", name; next }
        tolower($1)=="scf.electric.field" { printf "scf.Electric.Field          %s 0.0 0.0\n", f; next }
        { print }
    ' "$TEMPLATE" > "$out"
    echo "wrote $out (field $f GV/m along a)"
    inputs+=("$out")
done

if [ "$RUN" = "1" ]; then
    for inp in "${inputs[@]}"; do
        base=${inp%.dat}
        echo "== $(date '+%F %T') running $inp"
        mpirun -np "$NP" "$OPENMX" "$inp" -nt "$NT" > "${base}.std" 2>&1
        grep -q "The calculation was normally finished." "${base}.std" \
            || echo "WARNING: $inp did not finish normally"
    done
fi
