#!/bin/bash
# Counterpoise (Boys-Bernardi) correction for the water dimer with OpenMX ghost atoms.
#
# Usage:
#   bash run_bsse.sh          # run the five inputs, then print the analysis
#   bash run_bsse.sh collect  # only print the analysis from existing .out files
#
# Environment variables (optional):
#   MPI_RUN (default "mpirun -np 4"), OPENMX (default openmx), NT (default 1)
set -euo pipefail

MPI_RUN=${MPI_RUN:-"mpirun -np 4"}
OPENMX=${OPENMX:-openmx}
NT=${NT:-1}

# input file : System.Name
JOBS="H2O_dimer.dat:h2o_dimer
H2O_dimer_BSSE_ghost.dat:h2o_A_ghostB
H2O_dimer_BSSE_ghost2.dat:h2o_B_ghostA
H2O_monoA.dat:h2o_monoA
H2O_monoB.dat:h2o_monoB"

last_utot() {
    if [ ! -f "$1" ]; then echo "ERROR: $1 not found" >&2; exit 1; fi
    awk '/Utot\./{e=$2} END{if (e=="") exit 1; print e}' "$1" \
        || { echo "ERROR: no Utot. line in $1 (did the run finish?)" >&2; exit 1; }
}

if [ "${1:-run}" != "collect" ]; then
    for j in $JOBS; do
        inp=${j%%:*}; name=${j##*:}
        echo ">>> $MPI_RUN $OPENMX $inp -nt $NT"
        $MPI_RUN $OPENMX "$inp" -nt "$NT" > "$name.std" 2>&1
    done
fi

E_AB=$(last_utot h2o_dimer.out)
E_A_AB=$(last_utot h2o_A_ghostB.out)
E_B_AB=$(last_utot h2o_B_ghostA.out)
E_A=$(last_utot h2o_monoA.out)
E_B=$(last_utot h2o_monoB.out)

awk -v ab="$E_AB" -v aab="$E_A_AB" -v bab="$E_B_AB" -v a="$E_A" -v b="$E_B" 'BEGIN{
    c = 627.5095                       # kcal/mol per Hartree
    raw = (ab - a - b) * c
    cp  = (ab - aab - bab) * c
    printf "E_AB(AB)  = %.10f Ha\n", ab
    printf "E_A(AB)   = %.10f Ha   E_A(A) = %.10f Ha   BSSE(A) = %.3f kcal/mol\n", aab, a, (a - aab) * c
    printf "E_B(AB)   = %.10f Ha   E_B(B) = %.10f Ha   BSSE(B) = %.3f kcal/mol\n", bab, b, (b - bab) * c
    printf "dE_raw    = %8.3f kcal/mol  (no counterpoise)\n", raw
    printf "dE_CP     = %8.3f kcal/mol  (counterpoise corrected)\n", cp
    printf "BSSE      = %8.3f kcal/mol  (dE_CP - dE_raw)\n", cp - raw
}'
echo "# SCF convergence flags (expect scf_conv=1):"
grep -H -o "scf_conv=[01]" h2o_*.md || true
