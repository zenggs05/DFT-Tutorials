#!/bin/bash
# =====================================================================
# Same ENCUT scan as slurm_array_encut.sh, for IBM Spectrum LSF.
#   bsub < lsf_array_encut.sh
#   bjobs -A                    # one summary line per array
# LSF array indices start at 1 here: [1-7]; %4 = at most 4 running.
# %J = job id, %I = array index in the output file names.
# =====================================================================
#BSUB -J "si-encut[1-7]%4"
#BSUB -q normal
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -W 00:30
#BSUB -o encut-%J_%I.out
#BSUB -e encut-%J_%I.err

VALUES=(250 300 350 400 450 500 550)

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1

cd "${LS_SUBCWD:-$PWD}" || exit 1
E=${VALUES[$((LSB_JOBINDEX - 1))]}
[ -n "${E}" ] || { echo "no value for index ${LSB_JOBINDEX}" >&2; exit 1; }
[ -f POTCAR ] || { echo "POTCAR missing (see README)" >&2; exit 1; }

D=encut_${E}
mkdir -p "${D}"
cp POSCAR KPOINTS POTCAR "${D}/"
sed "s/@ENCUT@/${E}/g" INCAR.encut.template > "${D}/INCAR"
cd "${D}" || exit 1
mpirun -np "${LSB_DJOB_NUMPROC:-8}" vasp_std > vasp.out 2>&1
grep "free  energy   TOTEN" OUTCAR | tail -n 1
