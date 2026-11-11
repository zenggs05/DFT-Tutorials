#!/bin/bash
# Magnetocrystalline anisotropy of bcc Fe with the force theorem:
#   1) collinear SCF (vasp_std)              -> col/CHGCAR
#   2) non-self-consistent SOC runs (vasp_ncl, ICHARG=11) for SAXIS = 001, 100, 111
#   3) MAE(dir) = E(dir) - E(001)
# Usage:  bash run_mae.sh col | soc | collect
# NP sets the number of MPI ranks; VASP_STD / VASP_NCL the executables.
set -e
NP=${NP:-8}
STD=${VASP_STD:-vasp_std}
NCL=${VASP_NCL:-vasp_ncl}
DIRS="001 100 111"

case "$1" in
col)
    [ -s POTCAR ] || { echo "Build POTCAR: cat \$VASP_PP_PATH/potpaw_PBE/Fe_pv/POTCAR > POTCAR"; exit 1; }
    mkdir -p col
    cp INCAR_collinear col/INCAR; cp POSCAR KPOINTS POTCAR col/
    (cd col && mpirun -np $NP $STD > vasp.out 2>&1)
    grep "F=" col/OSZICAR | tail -1
    ;;
soc)
    [ -s col/CHGCAR ] || { echo "run 'bash run_mae.sh col' first"; exit 1; }
    nb=$(grep -m1 "NBANDS=" col/OUTCAR | sed 's/.*NBANDS= *\([0-9]*\).*/\1/')
    echo "collinear NBANDS = $nb  ->  SOC NBANDS = $((2*nb))"
    for d in $DIRS; do
        mkdir -p soc_$d
        cp POSCAR KPOINTS POTCAR col/CHGCAR soc_$d/
        sed "s/^NBANDS .*/NBANDS     = $((2*nb))/" INCAR_$d > soc_$d/INCAR
        (cd soc_$d && mpirun -np $NP $NCL > vasp.out 2>&1)
    done
    ;;
collect)
    ref=""
    printf "%-6s %18s %14s\n" "SAXIS" "TOTEN (eV)" "MAE (ueV/atom)"
    nat=$(sed -n '7p' POSCAR | awk '{s=0; for(i=1;i<=NF;i++) s+=$i; print s}')
    for d in $DIRS; do
        e=$(grep "free  energy   TOTEN" soc_$d/OUTCAR | tail -1 | awk '{print $(NF-1)}')
        [ -z "$ref" ] && ref=$e
        mae=$(awk -v e=$e -v r=$ref -v n=$nat 'BEGIN{printf "%.2f", (e-r)*1e6/n}')
        printf "%-6s %18s %14s\n" "$d" "$e" "$mae"
    done
    echo "100 vs 001 are equivalent in cubic Fe: their difference measures the numerical noise."
    ;;
*)
    echo "usage: bash run_mae.sh col | soc | collect"; exit 1 ;;
esac
