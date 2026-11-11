#!/bin/bash
# Interlayer binding-energy curve of 2H bilayer MoS2 (rigid-layer scan).
#
#   bash run.sh prepare   -> mono/ and d_<value>/ directories, submitted with sbatch
#   bash run.sh collect   -> binding energy table + binding_curve.png
#
# E_b(d) = [E_bilayer(d) - 2 E_monolayer] / A      (negative = bound)
# All calculations use the SAME cell (a, c), ENCUT, k-mesh and IVDW.
set -e
A=${A:-3.19}            # in-plane lattice constant from your monolayer relaxation
T=${T:-3.13}            # S-S height of one layer (from the relaxed monolayer)
C=${C:-30.0}            # cell length along z
DLIST="5.6 5.8 6.0 6.1 6.2 6.3 6.4 6.6 7.0 7.5 8.5 10.0"
JOB=${JOB:-job.slurm}   # any Slurm script that runs vasp_std in the current directory

case "$1" in
prepare)
    [ -s POTCAR ] || { echo "Build POTCAR: cat \$VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR"; exit 1; }
    [ -s $JOB ] || { echo "no $JOB in this directory"; exit 1; }
    mkdir -p mono
    python make_bilayer.py --mono --a $A --t $T --c $C > mono/POSCAR
    cp INCAR KPOINTS POTCAR $JOB mono/
    (cd mono && sbatch $JOB)
    for d in $DLIST; do
        dir=d_$d; mkdir -p $dir
        python make_bilayer.py --d $d --a $A --t $T --c $C > $dir/POSCAR
        cp INCAR KPOINTS POTCAR $JOB $dir/
        (cd $dir && sbatch $JOB)
    done
    ;;
collect)
    python analyze_binding.py --mono mono/OUTCAR --dirs d_* --a $A
    ;;
*)
    echo "usage: bash run.sh prepare | collect"; exit 1 ;;
esac
