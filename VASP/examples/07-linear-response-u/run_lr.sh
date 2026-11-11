#!/bin/bash
# Linear-response U (Cococcioni & de Gironcoli) workflow for VASP, LDAUTYPE = 3.
#
#   bash run_lr.sh gs          # 1) plain-DFT ground state in gs/
#   bash run_lr.sh response    # 2) for each V: nscf_V (ICHARG=11) and scf_V
#   bash run_lr.sh fit         # 3) linear fits -> chi0, chi, U = 1/chi - 1/chi0
#
# Runs are executed directly with mpirun (interactive/allocated node).
# To use the queue instead, replace the RUN line by a call to your job script.
set -e
RUN="mpirun -np ${NP:-32} ${VASP_EXE:-vasp_std}"
VLIST="-0.20 -0.10 -0.05 0.05 0.10 0.20"

[ -s POTCAR ] || { echo "Build POTCAR: cat \$VASP_PP_PATH/potpaw_PBE/{Ni_pv,Ni_pv,O}/POTCAR > POTCAR"; exit 1; }

case "$1" in
gs)
    mkdir -p gs; cp INCAR.gs gs/INCAR; cp POSCAR KPOINTS POTCAR gs/
    (cd gs && $RUN > vasp.out 2>&1)
    ;;
response)
    [ -s gs/CHGCAR ] || { echo "run 'bash run_lr.sh gs' first"; exit 1; }
    for V in $VLIST; do
        for mode in nscf scf; do
            d=${mode}_${V}
            mkdir -p $d
            cp POSCAR KPOINTS POTCAR $d/
            cp gs/WAVECAR $d/
            sed "s/VPERT/${V}/g" INCAR.lr > $d/INCAR
            if [ $mode == nscf ]; then
                cp gs/CHGCAR $d/
                sed -i 's/^#ICHARG  = 11/ICHARG   = 11/' $d/INCAR
            fi
            echo "running $d"
            (cd $d && $RUN > vasp.out 2>&1)
            rm -f $d/WAVECAR
        done
    done
    ;;
fit)
    python fit_u.py --gs gs/OUTCAR --dirs nscf_* scf_*
    ;;
*)
    echo "usage: bash run_lr.sh gs | response | fit"; exit 1 ;;
esac
