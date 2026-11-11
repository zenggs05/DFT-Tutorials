#!/bin/bash
# Compare plain PBE, PBE-D3(BJ) and optB86b-vdW for graphite:
#   bash run_compare.sh relax    -> PBE/relax, D3BJ/relax, optB86b/relax  (ISIF = 3)
#   bash run_compare.sh static   -> <method>/graphite (static on CONTCAR)
#                                   <method>/graphene (monolayer with the same a)
#   bash run_compare.sh collect  -> lattice constants and interlayer binding energy
# Jobs are submitted with job.slurm (sbatch); wait for each stage to finish.
set -e
METHODS="PBE D3BJ optB86b"

static_incar () {   # turn a relaxation INCAR into a static one
    sed -e 's/^IBRION .*/IBRION  = -1/' -e 's/^ISIF .*/ISIF    = 2/' \
        -e 's/^NSW .*/NSW     = 0/' -e '/^EDIFFG/d' "$1"
}

case "$1" in
relax)
    [ -s POTCAR ] || { echo "Build POTCAR: cat \$VASP_PP_PATH/potpaw_PBE/C/POTCAR > POTCAR"; exit 1; }
    for m in $METHODS; do
        d=$m/relax; mkdir -p $d
        cp POSCAR KPOINTS POTCAR job.slurm $d/
        cp INCAR_$m $d/INCAR
        (cd $d && sbatch job.slurm)
    done
    ;;
static)
    for m in $METHODS; do
        [ -s $m/relax/CONTCAR ] || { echo "$m/relax/CONTCAR missing"; exit 1; }
        # graphite at the relaxed geometry
        d=$m/graphite; mkdir -p $d
        cp $m/relax/CONTCAR $d/POSCAR; cp KPOINTS POTCAR job.slurm $d/
        static_incar INCAR_$m > $d/INCAR
        # isolated graphene layer with the same in-plane lattice
        d=$m/graphene; mkdir -p $d
        python make_graphene.py $m/relax/CONTCAR --c 20.0 > $d/POSCAR
        printf "graphene: Gamma-centred, 1 k-point along z\n0\nGamma\n15 15 1\n0 0 0\n" > $d/KPOINTS
        cp POTCAR job.slurm $d/
        static_incar INCAR_$m > $d/INCAR
        (cd $m/graphite && sbatch job.slurm)
        (cd $m/graphene && sbatch job.slurm)
    done
    ;;
collect)
    python analyze_graphite.py $METHODS
    ;;
*)
    echo "usage: bash run_compare.sh relax | static | collect"; exit 1 ;;
esac
