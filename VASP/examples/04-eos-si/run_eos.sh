#!/bin/bash
# E-V scan driver for VASP.
#   bash run_eos.sh prep      create V_<scale> directories (scaled POSCAR + INCAR/KPOINTS/POTCAR)
#   bash run_eos.sh submit    sbatch ../job.slurm inside every V_* directory
#   bash run_eos.sh run       run serially in the current allocation (uses $MPI_RUN and $VASP_EXE)
#   bash run_eos.sh collect   write eos.dat: volume(A^3/cell)  energy(sigma->0, eV/cell)
# Volume scale factors V/V_ref; edit as needed (keep the minimum near the middle).
SCALES=${SCALES:-"0.94 0.96 0.98 1.00 1.02 1.04 1.06"}
MODE=${1:-help}
set -e

case "$MODE" in
prep)
    for f in INCAR KPOINTS POSCAR POTCAR; do
        [ -f "$f" ] || { echo "missing $f (POTCAR must be built by you, see README)"; exit 1; }
    done
    for s in $SCALES; do
        d="V_$s"
        mkdir -p "$d"
        cp INCAR KPOINTS POTCAR "$d"/
        # Line 2 of POSCAR is the scaling factor. A volume scale s means a length
        # scale s^(1/3). A negative scaling factor is a target volume: multiply by s.
        awk -v s="$s" 'NR==2 { if ($1 < 0) printf "%.10f\n", $1*s; else printf "%.10f\n", $1*s^(1.0/3.0); next } { print }' \
            POSCAR > "$d"/POSCAR
        echo "prepared $d (scale factor line: $(sed -n 2p "$d"/POSCAR))"
    done
    ;;
submit)
    for d in V_*; do
        (cd "$d" && sbatch ../job.slurm)
    done
    ;;
run)
    MPI_RUN=${MPI_RUN:-"mpirun -np 4"}
    VASP_EXE=${VASP_EXE:-vasp_std}
    for d in V_*; do
        echo "running $d"
        (cd "$d" && $MPI_RUN $VASP_EXE > vasp.out 2>&1)
    done
    ;;
collect)
    rm -f eos.dat
    for d in V_*; do
        if [ ! -f "$d/OUTCAR" ] || ! grep -q "Total CPU time used" "$d/OUTCAR"; then
            echo "warning: $d not finished, skipped" >&2
            continue
        fi
        V=$(grep "volume of cell" "$d/OUTCAR" | tail -1 | awk '{print $NF}')
        E=$(grep "energy(sigma->0)" "$d/OUTCAR" | tail -1 | awk '{print $NF}')
        echo "$V $E" >> eos.dat
    done
    sort -g eos.dat -o eos.dat
    echo "# volume(A^3)  energy(eV)"
    cat eos.dat
    echo "next: python fit_bm.py eos.dat --lattice fcc --plot"
    ;;
*)
    sed -n 2,7p "$0"
    ;;
esac
