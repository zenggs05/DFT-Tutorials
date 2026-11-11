#!/bin/bash
# Slurm job array: one VASP calculation per sub-directory listed in dirs.txt
#   ls -d strain_* > dirs.txt          (one directory per line)
#   N=$(wc -l < dirs.txt); sbatch --array=1-${N}%10 array_job.slurm
# '%10' limits the number of tasks running at the same time.
#SBATCH --job-name=vasp-array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%A_%a.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}

DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" dirs.txt)
[ -n "$DIR" ] && [ -d "$DIR" ] || { echo "no directory for task ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
cd "$DIR" || exit 1

# skip directories that already finished successfully (makes resubmission safe)
if [ -f OUTCAR ] && grep -q "General timing" OUTCAR; then
    echo "$DIR already finished, skipping"
    exit 0
fi

for f in INCAR POSCAR KPOINTS POTCAR; do
    [ -s "$f" ] || { echo "$DIR: missing $f" >&2; exit 1; }
done

mpirun -np ${SLURM_NTASKS:-32} $VASP_EXE > vasp.out 2>&1

# short health check (the repository script VASP/scripts/check_vasp.py gives a full report)
if grep -q "General timing" OUTCAR; then
    echo "$DIR: finished, $(grep 'F=' OSZICAR | tail -1)"
else
    echo "$DIR: NOT finished - inspect vasp.out" >&2
fi
