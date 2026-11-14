#!/bin/bash
# Run one scan point per Slurm array task.
#
#   ls -d encut_scan/encut_* > dirs.txt          # one directory per line
#   sbatch --array=1-$(wc -l < dirs.txt)%4 slurm_array_dirs.sh
#
# "%4" limits the number of tasks running at the same time.
# Set RUN to the program call (it is executed inside each directory), e.g.
#   sbatch --export=ALL,RUN="mpirun vasp_std" ...       (VASP)
#   sbatch --export=ALL,RUN="mpirun cp2k.psmp -i input.inp -o output.out" ...
#SBATCH --job-name=conv-scan
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%A_%a.out

set -euo pipefail
ulimit -s unlimited
export OMP_NUM_THREADS=1
# module load vasp/6.5.1          # EDIT: same modules as used for compiling

RUN=${RUN:-"mpirun vasp_std"}
DIRLIST=${DIRLIST:-dirs.txt}
cd "$SLURM_SUBMIT_DIR"
d=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$DIRLIST")
[ -n "$d" ] && [ -d "$d" ] || { echo "no directory for task $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
echo "task $SLURM_ARRAY_TASK_ID -> $d on $(hostname)"
cd "$d"
$RUN > run.log 2>&1
