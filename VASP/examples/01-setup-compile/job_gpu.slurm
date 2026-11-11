#!/bin/bash
# Slurm job script for the OpenACC GPU build of VASP (NVIDIA GPUs).
# Rules: one MPI rank per GPU, OpenMP threads for the remaining CPU cores,
# INCAR: KPAR = number of GPUs, NSIM = 16-32 (NCORE is reset to 1 anyway).
#SBATCH --job-name=vasp_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4          # = number of GPUs per node
#SBATCH --cpus-per-task=8            # OpenMP threads per MPI rank
#SBATCH --gres=gpu:4
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# ---------------- environment (site specific) ----------------
module purge
module load nvhpc/25.1
export PATH=$HOME/software/vasp.6.6.1-gpu/bin:$PATH
# Root of the HPC SDK (the variable name differs between module systems)
NVROOT=${NVHPC_ROOT:-$(dirname "$(dirname "$(dirname "$(which nvfortran)")")")}
export LD_LIBRARY_PATH=$NVROOT/compilers/extras/qd/lib:$LD_LIBRARY_PATH
VASP_EXE=${VASP_EXE:-vasp_std}

ulimit -s unlimited
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
export MKL_THREADING_LAYER=INTEL     # only relevant for builds linked against MKL

# ---------------- checks ----------------
cd "$SLURM_SUBMIT_DIR" || exit 1
for f in INCAR POSCAR POTCAR; do
    [ -s "$f" ] || { echo "ERROR: $f is missing or empty"; exit 1; }
done
grep -E "^ *(KPAR|NSIM|NCORE)" INCAR || echo "Hint: set KPAR = number of GPUs and NSIM = 16-32"
nvidia-smi --query-gpu=index,name,memory.total --format=csv

echo "Job $SLURM_JOB_ID: $SLURM_NTASKS ranks x $OMP_NUM_THREADS threads, start $(date)"
mpirun -np "$SLURM_NTASKS" \
       --map-by ppr:"$SLURM_NTASKS_PER_NODE":node:PE="$SLURM_CPUS_PER_TASK" \
       --bind-to core --report-bindings \
       "$VASP_EXE" > vasp.out 2>&1
echo "End $(date)"

grep -q "General timing and accounting" OUTCAR 2>/dev/null \
    && echo "VASP finished normally" \
    || echo "VASP did NOT finish normally, check vasp.out and OUTCAR"
