#!/bin/bash
# Vacuum-thickness convergence for a 2D material.
#   bash vacuum_test.sh prepare   -> c_15 c_20 c_25 c_30 with static INCAR_scf, submitted
#   bash vacuum_test.sh collect   -> total energy per atom and E-fermi vs c
# The geometry is taken from relax/CONTCAR if it exists, otherwise from POSCAR.
set -e
CLIST="15 20 25 30"
SRC=POSCAR; [ -s relax/CONTCAR ] && SRC=relax/CONTCAR

case "$1" in
prepare)
    [ -s POTCAR ] || { echo "Build POTCAR: cat \$VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR"; exit 1; }
    for c in $CLIST; do
        d=c_$c; mkdir -p $d
        python set_vacuum.py $SRC --c $c > $d/POSCAR
        cp INCAR_scf $d/INCAR; cp KPOINTS POTCAR job.slurm $d/
        sed -i 's/^LCHARG .*/LCHARG  = .FALSE./' $d/INCAR
        (cd $d && sbatch job.slurm)
    done
    ;;
collect)
    nat=$(sed -n '7p' $SRC | awk '{s=0; for(i=1;i<=NF;i++) s+=$i; print s}')
    printf "%6s %18s %12s\n" "c(A)" "E0/atom (eV)" "E_F (eV)"
    for c in $CLIST; do
        o=c_$c/OUTCAR; [ -s $o ] || continue
        e=$(grep "energy(sigma->0)" $o | tail -1 | awk -v n=$nat '{printf "%.6f", $NF/n}')
        ef=$(grep "E-fermi" $o | tail -1 | awk '{print $3}')
        printf "%6s %18s %12s\n" "$c" "$e" "$ef"
    done
    echo "Converged when E0/atom changes by < 1 meV between successive c."
    echo "E_F alone is NOT comparable between cells: reference it to the vacuum level"
    echo "(python ../09-pt111-workfunction-dipole/plot_locpot_planar.py c_20/LOCPOT --outcar c_20/OUTCAR)."
    ;;
*)
    echo "usage: bash vacuum_test.sh prepare | collect"; exit 1 ;;
esac
