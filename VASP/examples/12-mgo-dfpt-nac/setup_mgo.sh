#!/bin/bash
# Prepare the two VASP runs of the MgO example.
#   bash setup_mgo.sh          -> dfpt/ (supercell, IBRION=8) and born/ (unit cell, LEPSILON)
#   bash setup_mgo.sh post     -> FORCE_CONSTANTS, BORN and the band structure
# POSCAR-unitcell must be the relaxed conventional cell; POTCAR must exist here.
set -e
if command -v phonopy-init > /dev/null; then PHINIT=phonopy-init; NEW=1; else PHINIT=phonopy; NEW=0; fi
PA="0 1/2 1/2 1/2 0 1/2 1/2 1/2 0"

if [ "$1" != "post" ]; then
    for f in POSCAR-unitcell POTCAR INCAR_dfpt KPOINTS_dfpt INCAR_born KPOINTS_born; do
        [ -f $f ] || { echo "missing $f"; exit 1; }
    done
    $PHINIT -d --dim 2 2 2 -c POSCAR-unitcell     # writes SPOSCAR and phonopy_disp.yaml
    mkdir -p dfpt born
    cp SPOSCAR dfpt/POSCAR;          cp INCAR_dfpt dfpt/INCAR; cp KPOINTS_dfpt dfpt/KPOINTS; cp POTCAR dfpt/
    cp POSCAR-unitcell born/POSCAR;  cp INCAR_born born/INCAR; cp KPOINTS_born born/KPOINTS; cp POTCAR born/
    echo "now run VASP in dfpt/ and in born/ (e.g. cd dfpt && sbatch ../../12-si-phonopy-fd/job.slurm)"
    exit 0
fi

grep -q 'name="hessian"' dfpt/vasprun.xml || { echo "dfpt/vasprun.xml has no hessian"; exit 1; }
$PHINIT --fc dfpt/vasprun.xml                   # -> FORCE_CONSTANTS
(cd born && phonopy-vasp-born --pa "$PA") > BORN # -> BORN (eps_inf and Z*)
cat BORN
if [ $NEW -eq 1 ]; then
    phonopy --config band.conf -p -s
else
    phonopy -c POSCAR-unitcell -p -s band.conf
fi
