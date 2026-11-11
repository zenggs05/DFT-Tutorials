#!/bin/bash
# Basis-set convergence test for diamond Si with OpenMX 4.0.
#
# Usage:
#   bash run_basis_test.sh scf       # SCF + band for quick / standard / precise, then summary
#   bash run_basis_test.sh eos       # EvsLC lattice scan for each basis (scf.Ngrid fixed)
#   bash run_basis_test.sh collect   # only print the summary of the scf runs
#
# Environment variables (optional):
#   MPI_RUN  MPI launcher            (default: "mpirun -np 4")
#   OPENMX   OpenMX executable       (default: openmx)
#   NT       OpenMP threads/process  (default: 1)
#
# Requires GNU sed (Linux). On macOS use gsed.
set -euo pipefail

MPI_RUN=${MPI_RUN:-"mpirun -np 4"}
OPENMX=${OPENMX:-openmx}
NT=${NT:-1}
NATOM=2
BASES="quick standard precise"
HA2MEV=27211.386

run_one() {   # $1 = input file, $2 = prefix for the captured standard output
    echo ">>> $MPI_RUN $OPENMX $1 -nt $NT"
    $MPI_RUN $OPENMX "$1" -nt "$NT" > "$2.std" 2>&1
}

last_utot() {  # last "Utot." line of an OpenMX .out file (Hartree)
    awk '/Utot\./{e=$2} END{print e}' "$1"
}

collect() {
    local eref=""
    [ -f si_precise.out ] && eref=$(last_utot si_precise.out)
    printf "%-10s %12s %20s %26s\n" "# basis" "nbasis/atom" "Utot(Ha/cell)" "dE_vs_precise(meV/atom)"
    for b in $BASES; do
        f=si_${b}.out
        if [ ! -f "$f" ]; then
            echo "WARNING: $f not found, skipped" >&2
            continue
        fi
        case $b in quick) nb=10 ;; standard) nb=13 ;; precise) nb=22 ;; esac
        e=$(last_utot "$f")
        if [ -n "$eref" ]; then
            awk -v b="$b" -v nb="$nb" -v e="$e" -v r="$eref" -v n="$NATOM" -v c="$HA2MEV" \
                'BEGIN{printf "%-10s %12d %20.10f %26.3f\n", b, nb, e, (e-r)*c/n}'
        else
            printf "%-10s %12d %20.10f %26s\n" "$b" "$nb" "$e" "n/a"
        fi
    done
    echo "# grid check (must be identical for all runs):"
    grep -H "Num. of grids" si_*.out | grep -v _eos || true
}

scf() {
    for b in $BASES; do
        run_one "Si_basis_${b}.dat" "si_${b}"
    done
    collect
}

eos() {
    # Start from a = 5.24 Ang (lattice vectors 2.6200), expand by 0.7 % per step, 12 steps.
    # scf.Ngrid is fixed so that the number of grid points does not jump between steps.
    for b in $BASES; do
        in=Si_basis_${b}.dat
        out=Si_basis_${b}_eos.dat
        sed -e "s/^System.Name .*/System.Name                 si_${b}_eos/" \
            -e "s/2\.7150/2.6200/g" \
            -e "s/^MD.Type .*/MD.Type                     EvsLC/" \
            -e "/^MD.Type/a MD.EvsLC.Step               0.7" \
            -e "/^MD.Type/a MD.maxIter                  12" \
            -e "/^scf.energycutoff/a scf.Ngrid                   32 32 32" \
            -e "s/^Band.dispersion .*/Band.dispersion             off/" \
            "$in" > "$out"
        run_one "$out" "si_${b}_eos"
    done
    echo "Fit each curve with, e.g.:"
    for b in $BASES; do
        echo "  python ../06-cellopt-eos-si-c/fit_evslc.py si_${b}_eos.EvsLC --lattice fcc --natoms 2"
    done
}

case "${1:-}" in
    scf)     scf ;;
    eos)     eos ;;
    collect) collect ;;
    *) echo "usage: bash $0 scf|eos|collect" >&2; exit 1 ;;
esac
