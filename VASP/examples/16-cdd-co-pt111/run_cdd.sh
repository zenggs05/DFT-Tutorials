#!/bin/bash
# Charge density difference workflow for CO/Pt(111):  d_rho = rho(CO/Pt) - rho(Pt slab) - rho(CO)
#
#   bash run_cdd.sh prepare    # make AB/ slab/ co/ with identical cell, INCAR, KPOINTS
#   bash run_cdd.sh run        # run AB, copy its NGXF/NGYF/NGZF into slab/ and co/, run both
#   bash run_cdd.sh analyze    # check grids, compute CHGCAR_cdd and the plane average
#
# Environment: VASP_PP_PATH (for POTCAR), VASP_EXE (default vasp_std), NP (MPI ranks)
set -e
STEP=${1:-all}
VASP_EXE=${VASP_EXE:-vasp_std}
NP=${NP:-${SLURM_NTASKS:-32}}
PP=${VASP_PP_PATH:-}/potpaw_PBE

prepare() {
    mkdir -p AB slab co
    cp POSCAR AB/POSCAR
    python3 split_fragments.py POSCAR --keep Pt  --out slab/POSCAR
    python3 split_fragments.py POSCAR --keep C O --out co/POSCAR
    for d in AB slab co; do cp INCAR KPOINTS $d/; done
    if [ -n "$VASP_PP_PATH" ]; then
        cat $PP/Pt/POTCAR $PP/C/POTCAR $PP/O/POTCAR > AB/POTCAR
        cat $PP/Pt/POTCAR                         > slab/POTCAR
        cat $PP/C/POTCAR $PP/O/POTCAR             > co/POTCAR
        grep -H TITEL AB/POTCAR slab/POTCAR co/POTCAR
    else
        echo "VASP_PP_PATH not set: create AB/POTCAR (Pt C O), slab/POTCAR (Pt), co/POTCAR (C O) yourself"
    fi
}

run_one() {
    (cd "$1" && mpirun -np "$NP" "$VASP_EXE" > vasp.out 2>&1)
    grep -q "General timing" "$1/OUTCAR" || { echo "run in $1 did not finish" >&2; exit 1; }
}

run_all() {
    run_one AB
    # copy the fine grid of AB into the fragments so that all three grids are identical
    grid=$(grep "dimension x,y,z NGXF" AB/OUTCAR | tail -1 | grep -oE "NG[XYZ]F= *[0-9]+" | tr -d ' ')
    echo "fine grid of AB: $grid"
    for d in slab co; do
        sed -i '/^ *NG[XYZ]F *=/d' $d/INCAR
        for kv in $grid; do echo "${kv/=/ = }" >> $d/INCAR; done
    done
    run_one slab
    run_one co
}

analyze() {
    echo "fine FFT grids (must be identical):"
    grep -H "dimension x,y,z NGXF" AB/OUTCAR slab/OUTCAR co/OUTCAR | tail -3
    python3 cdd_planar.py AB/CHGCAR slab/CHGCAR co/CHGCAR --axis z --out CHGCAR_cdd
    # optional cross-check with the VTST perl scripts
    if command -v chgdiff.pl >/dev/null 2>&1; then
        chgdiff.pl AB/CHGCAR slab/CHGCAR >/dev/null && mv CHGCAR_diff CHGCAR_AB-slab
        chgdiff.pl CHGCAR_AB-slab co/CHGCAR >/dev/null && mv CHGCAR_diff CHGCAR_cdd_vtst
        rm -f CHGCAR_AB-slab
        echo "VTST result written to CHGCAR_cdd_vtst"
    fi
}

case $STEP in
    prepare) prepare ;;
    run)     run_all ;;
    analyze) analyze ;;
    all)     prepare; run_all; analyze ;;
    *) echo "usage: bash run_cdd.sh [prepare|run|analyze|all]"; exit 1 ;;
esac
