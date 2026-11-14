#!/bin/bash
# =====================================================================
# Same ENCUT scan as slurm_array_encut.sh, for PBS Pro / OpenPBS.
#   qsub pbs_array_encut.sh
#   qstat -t -u $USER           # -t lists the sub-jobs of the array
# Unlike the Slurm/LSF versions, no "at most 4 at once" limit is set
# here; ask your admins whether your PBS version accepts "-J 0-6%4".
# Torque uses a different syntax: "#PBS -t 0-6%4" and $PBS_ARRAYID
# (see the commented lines).
# =====================================================================
#PBS -N si-encut
#PBS -q workq
#PBS -l select=1:ncpus=8:mpiprocs=8:ompthreads=1
#PBS -l walltime=00:30:00
#PBS -J 0-6
#PBS -j oe
# --- Torque ----------------------------------------------------------
##PBS -l nodes=1:ppn=8
##PBS -t 0-6%4

VALUES=(250 300 350 400 450 500 550)

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1

cd "${PBS_O_WORKDIR}" || exit 1
ID=${PBS_ARRAY_INDEX:-${PBS_ARRAYID}}      # PBS Pro : Torque
E=${VALUES[$ID]}
[ -n "${E}" ] || { echo "no value for index ${ID}" >&2; exit 1; }
[ -f POTCAR ] || { echo "POTCAR missing (see README)" >&2; exit 1; }

D=encut_${E}
mkdir -p "${D}"
cp POSCAR KPOINTS POTCAR "${D}/"
sed "s/@ENCUT@/${E}/g" INCAR.encut.template > "${D}/INCAR"
cd "${D}" || exit 1
NP=$(wc -l < "${PBS_NODEFILE}")
mpirun -np "${NP}" vasp_std > vasp.out 2>&1
grep "free  energy   TOTEN" OUTCAR | tail -n 1
