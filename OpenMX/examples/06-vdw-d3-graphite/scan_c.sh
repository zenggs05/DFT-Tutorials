#!/bin/bash
# Single-point E(c) scan of AB graphite, independent of MD.Type EvsLC.
# Useful (1) as a fallback if your OpenMX build rejects EvsLC together with DFT-D,
#        (2) to compare PBE with and without DFT-D3 on exactly the same c values.
#
# Usage:
#   DFTD=on  bash scan_c.sh [run|collect]     # PBE-D3(BJ)  (default)
#   DFTD=off bash scan_c.sh [run|collect]     # plain PBE
#
# Environment variables (optional):
#   CVALUES (default "6.2 6.4 6.6 6.7 6.8 7.0 7.2 7.6 8.0")   c in Angstrom
#   MPI_RUN (default "mpirun -np 4")  OPENMX (default openmx)  NT (default 1)
#
# Output: graphite_<on|off>_scan.EvsLC with the same 10-column layout as an EvsLC file,
# so it can be fitted with ../06-cellopt-eos-si-c/fit_evslc.py --mode c
set -euo pipefail

MPI_RUN=${MPI_RUN:-"mpirun -np 4"}
OPENMX=${OPENMX:-openmx}
NT=${NT:-1}
DFTD=${DFTD:-on}
CVALUES=${CVALUES:-"6.2 6.4 6.6 6.7 6.8 7.0 7.2 7.6 8.0"}
BASE=Graphite_D3_opt.dat
TAG=graphite_${DFTD}

case "$DFTD" in on|off) ;; *) echo "DFTD must be on or off" >&2; exit 1 ;; esac
action=${1:-run}

if [ "$action" = run ]; then
    for c in $CVALUES; do
        name=${TAG}_c${c}
        sed -e "s/^System.Name .*/System.Name                 ${name}/" \
            -e "s/^   0.000000   0.000000   6.200000/   0.000000   0.000000   ${c}/" \
            -e "s/^MD.Type .*/MD.Type                     Nomd/" \
            -e "s/^scf.dftD .*/scf.dftD                    ${DFTD}/" \
            "$BASE" > "${name}.dat"
        echo ">>> $MPI_RUN $OPENMX ${name}.dat -nt $NT"
        $MPI_RUN $OPENMX "${name}.dat" -nt "$NT" > "${name}.std" 2>&1
    done
fi

out=${TAG}_scan.EvsLC
echo "# ax ay az bx by bz cx cy cz Etot(Hartree)   from scan_c.sh, DFTD=${DFTD}" > "$out"
for c in $CVALUES; do
    f=${TAG}_c${c}.out
    if [ ! -f "$f" ]; then echo "WARNING: $f not found, skipped" >&2; continue; fi
    e=$(awk '/Utot\./{x=$2} END{print x}' "$f")
    conv=$(grep -o 'scf_conv=[01]' "${TAG}_c${c}.md" 2>/dev/null | tail -1 || true)
    [ "$conv" = "scf_conv=1" ] || echo "WARNING: ${TAG}_c${c} not converged (${conv:-no .md})" >&2
    echo "2.464 0.0 0.0 -1.232 2.133886 0.0 0.0 0.0 $c $e" >> "$out"
done
echo "written: $out"
echo "fit with: python ../06-cellopt-eos-si-c/fit_evslc.py $out --mode c --natoms 4"
