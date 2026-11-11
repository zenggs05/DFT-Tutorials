#!/bin/bash
# Cu(100) surface energy: layer convergence, relaxed and unrelaxed.
#   bash run.sh bulk      -> bulk/  : ISIF=3 lattice constant (fcc primitive)
#   bash run.sh slabs     -> obulk/, slab_N/, slab_N_unrel/ built from bulk/CONTCAR, submitted
#   bash run.sh collect   -> surface_energy.py
# Every run uses the same ENCUT, POTCAR, smearing and in-plane k mesh.
set -e
LAYERS="3 5 7 9 11"
JOB=job.slurm

case "$1" in
bulk)
    [ -s POTCAR ] || { echo "Build POTCAR: cat \$VASP_PP_PATH/potpaw_PBE/Cu_pv/POTCAR > POTCAR"; exit 1; }
    mkdir -p bulk
    cp POSCAR_bulk bulk/POSCAR; cp KPOINTS_bulk bulk/KPOINTS; cp INCAR_bulk bulk/INCAR
    cp POTCAR $JOB bulk/
    (cd bulk && sbatch $JOB)
    echo "after it finishes, re-run it once from CONTCAR (Pulay stress), then: bash run.sh slabs"
    ;;
slabs)
    [ -s bulk/CONTCAR ] || { echo "bulk/CONTCAR missing"; exit 1; }
    python make_slabs.py bulk/CONTCAR --layers $LAYERS --vacuum 15 --kin 16 --kz 12
    cp INCAR_static obulk/INCAR; cp POTCAR $JOB obulk/; (cd obulk && sbatch $JOB)
    for n in $LAYERS; do
        cp INCAR_slab slab_$n/INCAR;          cp POTCAR $JOB slab_$n/
        cp INCAR_static slab_${n}_unrel/INCAR; cp POTCAR $JOB slab_${n}_unrel/
        (cd slab_$n && sbatch $JOB)
        (cd slab_${n}_unrel && sbatch $JOB)
    done
    ;;
collect)
    python surface_energy.py --nfit 3
    ;;
*)
    echo "usage: bash run.sh bulk | slabs | collect"; exit 1 ;;
esac
