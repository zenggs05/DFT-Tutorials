#!/bin/bash
# Set up the four collinear configurations needed to extract J1 and J2 of NiO:
#   FM_R   : ferromagnetic, rhombohedral cell (POSCAR,      KPOINTS)
#   AFM2   : AFM-II,        rhombohedral cell (POSCAR,      KPOINTS)
#   AFM1   : AFM-I,         tetragonal cell   (POSCAR_AFM1, KPOINTS_AFM1)
#   FM_T   : ferromagnetic, tetragonal cell   (consistency check: E/f.u. == FM_R)
# All runs share INCAR (ENCUT, U, smearing ...); only MAGMOM differs.
# Usage: bash run_exchange.sh          -> prepare and submit
#        bash run_exchange.sh collect  -> run exchange_J.py on the results
set -e
if [ "$1" == "collect" ]; then
    python exchange_J.py --fm FM_R/OUTCAR --afm2 AFM2/OUTCAR --afm1 AFM1/OUTCAR --fm-check FM_T/OUTCAR
    exit 0
fi
[ -s POTCAR ] || { echo "Build POTCAR first: cat \$VASP_PP_PATH/potpaw_PBE/{Ni_pv,O}/POTCAR > POTCAR"; exit 1; }

setup () {   # $1 dir, $2 POSCAR, $3 KPOINTS, $4 MAGMOM string
    mkdir -p $1
    cp $2 $1/POSCAR; cp $3 $1/KPOINTS; cp POTCAR job.slurm $1/
    sed "s/^MAGMOM .*/MAGMOM    = $4/" INCAR > $1/INCAR
    (cd $1 && sbatch job.slurm)
}
setup FM_R POSCAR       KPOINTS       "2.0 2.0 2*0.0"
setup AFM2 POSCAR       KPOINTS       "2.0 -2.0 2*0.0"
setup AFM1 POSCAR_AFM1  KPOINTS_AFM1  "2.0 -2.0 2*0.0"
setup FM_T POSCAR_AFM1  KPOINTS_AFM1  "2.0 2.0 2*0.0"
