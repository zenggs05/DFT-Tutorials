#!/bin/bash
# phonopy finite-displacement workflow helper (works with phonopy 2.x/3.x and >= 4.0)
#
#   bash run_disp.sh prepare   # supercells with displacements -> disp-XXX/{POSCAR,INCAR,KPOINTS,POTCAR}
#   bash run_disp.sh submit    # sbatch ../job.slurm inside every disp-XXX
#   bash run_disp.sh check     # are all vasprun.xml files complete?
#   bash run_disp.sh collect   # FORCE_SETS from disp-*/vasprun.xml
#
# POSCAR must be the RELAXED unit cell; INCAR/KPOINTS/POTCAR are those for the supercell.
# Supercell size: DIM="3 3 3" bash run_disp.sh prepare
set -e
DIM=${DIM:-"2 2 2"}

# phonopy >= 4.0 moved the setup options (-d, -f, --fc) to "phonopy-init"
if command -v phonopy-init > /dev/null; then
    PHINIT=phonopy-init
else
    PHINIT=phonopy
fi

case "$1" in
prepare)
    for f in POSCAR INCAR KPOINTS POTCAR; do [ -f $f ] || { echo "missing $f"; exit 1; }; done
    $PHINIT -d --dim $DIM -c POSCAR
    for p in POSCAR-[0-9]*; do
        n=${p#POSCAR-}
        mkdir -p disp-$n
        cp "$p" disp-$n/POSCAR
        cp INCAR KPOINTS POTCAR disp-$n/
    done
    echo "created $(ls -d disp-* | wc -l) displacement directories (supercell: SPOSCAR)"
    ;;
submit)
    for d in disp-*/; do (cd "$d" && sbatch ../job.slurm); done
    ;;
check)
    bad=0
    for d in disp-*/; do
        if ! grep -q "</modeling>" "$d/vasprun.xml" 2>/dev/null; then
            echo "not finished: $d"; bad=1
        fi
    done
    [ $bad -eq 0 ] && echo "all vasprun.xml files are complete"
    ;;
collect)
    $PHINIT -f disp-*/vasprun.xml
    if [ "$PHINIT" = "phonopy-init" ]; then
        echo "next: phonopy --config band.conf -p -s"
    else
        echo "next: phonopy -p -s band.conf"
    fi
    ;;
*)
    echo "usage: bash run_disp.sh prepare|submit|check|collect"; exit 1
    ;;
esac
