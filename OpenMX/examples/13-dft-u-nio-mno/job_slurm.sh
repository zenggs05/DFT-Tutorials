#!/bin/bash
#SBATCH --job-name=omx-ch13-dftu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.log

# Slurm job for chapter 13 (DFT+U, NiO/MnO).
# Runs every input listed in INPUTS one after another, then summarizes.
# Override on submission, e.g.:  sbatch --export=ALL,INPUTS="MnO_AFM_U.dat" job_slurm.sh

# module load openmx/4.0            # EDIT: load your OpenMX / MPI environment
# export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH   # deb package only

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_STACKSIZE=512M

OPENMX=${OPENMX:-openmx}
NP=${SLURM_NTASKS:-4}
INPUTS=${INPUTS:-"NiO_AFM_GGA.dat NiO_AFM_U_Dudarev.dat NiO_AFM_UJ_cFLL.dat MnO_AFM_U.dat"}

for inp in $INPUTS; do
    base=${inp%.dat}
    echo "== $(date '+%F %T') running $inp with $NP MPI x $OMP_NUM_THREADS threads"
    mpirun -np "$NP" "$OPENMX" "$inp" -nt "$OMP_NUM_THREADS" > "${base}.std" 2>&1
done

# System.Name differs from the input file name; summarize every .out present
python3 gap_from_out.py *.out
