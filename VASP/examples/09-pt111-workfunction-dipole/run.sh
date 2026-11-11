#!/bin/bash
# Three static runs that illustrate the dipole correction:
#   clean/      relaxed clean Pt(111) (CONTCAR of ../09-pt111-slab), with LDIPOL
#   O_nodip/    O/Pt(111) (O on the top side only) WITHOUT dipole correction
#   O_dip/      O/Pt(111) WITH dipole correction
# Usage: bash run.sh prepare | collect
set -e
case "$1" in
prepare)
    [ -s POTCAR_Pt ] || { echo "Build POTCAR_Pt: cat \$VASP_PP_PATH/potpaw_PBE/Pt/POTCAR > POTCAR_Pt"; exit 1; }
    [ -s POTCAR_PtO ] || { echo "Build POTCAR_PtO: cat \$VASP_PP_PATH/potpaw_PBE/{Pt,O}/POTCAR > POTCAR_PtO"; exit 1; }
    SLAB=../09-pt111-slab/CONTCAR
    [ -s $SLAB ] || { echo "relax the slab first (../09-pt111-slab), or set SLAB=../09-pt111-slab/POSCAR"; exit 1; }
    mkdir -p clean O_nodip O_dip
    cp $SLAB clean/POSCAR;  cp POTCAR_Pt  clean/POTCAR
    cp POSCAR_O_fcc O_nodip/POSCAR; cp POTCAR_PtO O_nodip/POTCAR
    cp POSCAR_O_fcc O_dip/POSCAR;   cp POTCAR_PtO O_dip/POTCAR
    cp INCAR clean/; cp INCAR O_dip/
    sed -e 's/^LDIPOL .*/LDIPOL    = .FALSE./' -e 's/^IDIPOL .*/#IDIPOL   = 3/' \
        -e 's/^DIPOL .*/#DIPOL    = 0.5 0.5 0.5/' -e 's/^LVACPOTAV .*/#LVACPOTAV = .TRUE./' \
        INCAR > O_nodip/INCAR
    for d in clean O_nodip O_dip; do
        cp KPOINTS job.slurm $d/
        (cd $d && sbatch job.slurm)
    done
    ;;
collect)
    for d in clean O_nodip O_dip; do
        echo "=================== $d"
        [ -s $d/LOCPOT ] || { echo "no LOCPOT"; continue; }
        (cd $d && python ../plot_locpot_planar.py LOCPOT --outcar OUTCAR --no-plot | sed -n '4,11p')
        grep "upper" $d/OUTCAR 2>/dev/null | tail -1 || true
    done
    ;;
*)
    echo "usage: bash run.sh prepare | collect"; exit 1 ;;
esac
