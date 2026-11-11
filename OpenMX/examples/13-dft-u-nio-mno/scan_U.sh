#!/bin/bash
# Scan U_eff (Dudarev form) on the Ni 1d orbital of NiO AFM-II.
#
# Usage:
#   bash scan_U.sh                     # default list "0 2 4 5 6 7 8"
#   bash scan_U.sh "0 3 6"             # custom list (eV)
#   RESTART=1 bash scan_U.sh           # reuse the converged state of the previous U
#
# Environment variables:
#   NP      number of MPI processes      (default 4)
#   NT      OpenMP threads per process   (default 1)
#   OPENMX  openmx executable            (default: openmx in PATH)
#   MPIRUN  MPI launcher                 (default: mpirun)
#
# Template: NiO_AFM_U_Dudarev.dat (only the "1d <value>" on the Ni line,
# System.Name and scf.restart are modified).
# Result:   nio_u<U>.out for every U and scan_U_summary.csv (via gap_from_out.py)

set -euo pipefail

ULIST=${1:-"0 2 4 5 6 7 8"}
NP=${NP:-4}
NT=${NT:-1}
OPENMX=${OPENMX:-openmx}
MPIRUN=${MPIRUN:-mpirun}
RESTART=${RESTART:-0}
TEMPLATE=NiO_AFM_U_Dudarev.dat

export OMP_NUM_THREADS=$NT
ulimit -s unlimited 2>/dev/null || true

outs=()
first=1
for U in $ULIST; do
    if [ "$RESTART" = "1" ]; then
        name=nio_scan
        if [ $first -eq 1 ]; then rst=off; else rst=on; fi
    else
        name=nio_u${U}
        rst=off
    fi
    inp=scan_u${U}.dat
    sed -e "s/^System\.Name .*/System.Name                 ${name}/" \
        -e "/^ *Ni  *1s /s/1d  *[0-9.]*/1d ${U}/" \
        -e "s/^scf\.restart .*/scf.restart                 ${rst}/" \
        "$TEMPLATE" > "$inp"

    echo "== U = ${U} eV  (System.Name=${name}, scf.restart=${rst})"
    grep -E "^ *Ni  *1s" "$inp"
    $MPIRUN -np "$NP" "$OPENMX" "$inp" -nt "$NT" > "scan_u${U}.std" 2>&1

    if [ "$RESTART" = "1" ]; then
        cp "${name}.out" "nio_u${U}.out"
    fi
    outs+=("nio_u${U}.out")
    first=0
done

python3 gap_from_out.py --csv scan_U_summary.csv "${outs[@]}"
