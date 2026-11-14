#!/bin/bash
# Accuracy test of the O(N) solvers against the conventional band calculation (Si64).
#
#   NP=16 NT=1 bash scan_hopping.sh                 # DC, DC-LNO and Krylov, R = 7 8 9 10 11 A
#   NP=16 NT=1 RLIST="7 9 11" SOLVERS="DC" bash scan_hopping.sh
#
# Steps:
#   1. runs Si64_Band_ref.dat once (skipped if si64_band.out already finished),
#   2. for every solver and every orderN.HoppingRanges R writes <name>_R<R>.dat and runs it,
#   3. writes hopping_scan.csv:
#      solver, R_ang, Utot_Hartree, dE_meV_per_atom (vs. band), wall_s, FNAN_SNAN_line
#
# Environment: OPENMX (default openmx), MPIRUN (default mpirun), NP (default 8), NT (default 1)
set -u
OPENMX=${OPENMX:-openmx}
MPIRUN=${MPIRUN:-mpirun}
NP=${NP:-8}
NT=${NT:-1}
RLIST=${RLIST:-"7 8 9 10 11"}
SOLVERS=${SOLVERS:-"DC DC-LNO Krylov"}
NATOM=64
HARTREE_MEV=27211.386

export OMP_NUM_THREADS=${NT}
export OMP_STACKSIZE=${OMP_STACKSIZE:-512M}
ulimit -s unlimited 2>/dev/null || true

finished() { grep -q "The calculation was normally finished." "$1" 2>/dev/null; }
utot()     { awk '/^ *Utot\./{v=$2} END{print v}' "$1" 2>/dev/null; }
walltime() { awk '/Total Computational Time/{v=$NF} END{print v}' "$1" 2>/dev/null; }

# ---- 1. reference ----
if ! finished si64_band.std; then
    echo ">>> reference band calculation"
    ${MPIRUN} -np "${NP}" "${OPENMX}" Si64_Band_ref.dat -nt "${NT}" > si64_band.std 2>&1
    finished si64_band.std || { echo "reference run failed, see si64_band.std"; exit 1; }
fi
EREF=$(utot si64_band.out)
echo "reference Utot = ${EREF} Hartree"

# ---- 2. scan ----
CSV=hopping_scan.csv
echo "solver,R_ang,Utot_Hartree,dE_meV_per_atom,wall_s,FNAN_SNAN_line" > "${CSV}"
for solver in ${SOLVERS}; do
    src=Si64_${solver}.dat
    [ -f "${src}" ] || { echo "missing ${src}"; continue; }
    base=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${src}")
    for R in ${RLIST}; do
        tag=${base}_R${R}
        awk -v n="${tag}" -v r="${R}" '
            tolower($1)=="system.name"          {print "System.Name   " n; next}
            tolower($1)=="ordern.hoppingranges" {print "orderN.HoppingRanges   " r; next}
            {print}' "${src}" > "${tag}.dat"
        echo ">>> ${solver}  R = ${R} A"
        ${MPIRUN} -np "${NP}" "${OPENMX}" "${tag}.dat" -nt "${NT}" > "${tag}.std" 2>&1
        if finished "${tag}.std"; then
            e=$(utot "${tag}.out")
            de=$(awk -v e="${e}" -v r="${EREF}" -v n="${NATOM}" -v h="${HARTREE_MEV}" \
                 'BEGIN{printf "%.3f", (e-r)/n*h}')
        else
            e=FAILED; de=
        fi
        # FNAN/SNAN are reported in the standard output; keep the first line that mentions them
        fs=$(grep -m1 -i "SNAN" "${tag}.std" | tr -s ' ' | tr ',' ' ')
        echo "${solver},${R},${e},${de},$(walltime "${tag}.out"),${fs}" >> "${CSV}"
    done
done
echo
column -s, -t "${CSV}" 2>/dev/null || cat "${CSV}"
echo
echo "Choose the smallest R whose |dE| is below your tolerance (e.g. 1-5 meV/atom),"
echo "then use it in Si512_DC.dat / Si512_DC-LNO.dat."
