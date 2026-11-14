#!/bin/bash
#SBATCH --job-name=omx-api
#SBATCH --nodes=1
#SBATCH --ntasks=16                   # total MPI ranks, split evenly among the inputs
#SBATCH --cpus-per-task=1             # OpenMP threads per rank (= --nt)
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log
#
# Run all *.dat files of this directory as one MPI job through the OpenMX Python API.
#   sbatch job_api_slurm.sh                      # every *.dat in the directory
#   sbatch job_api_slurm.sh a.dat b.dat          # selected inputs
# The number of inputs must not exceed --ntasks.
set -u
# module load ...                      # the SAME MPI/compiler/MKL environment used to build OpenMX
# source ~/venv-openmx/bin/activate    # Python with mpi4py built against that MPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_STACKSIZE=512M
ulimit -s unlimited

shopt -s nullglob
if [ $# -gt 0 ]; then inputs=("$@"); else inputs=(*.dat); fi
[ ${#inputs[@]} -gt 0 ] || { echo "no input files"; exit 1; }
echo "inputs: ${inputs[*]}  ranks: ${SLURM_NTASKS:-?}  threads: ${OMP_NUM_THREADS}"

# pyopenmx_api*.so must be importable (same directory or PYTHONPATH)
export PYTHONPATH=${PWD}:${PYTHONPATH:-}
mpirun -np "${SLURM_NTASKS:-1}" python3 api_batch.py --nt "${OMP_NUM_THREADS}" "${inputs[@]}"

for f in "${inputs[@]}"; do
    name=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${f}")
    if [ -f "${name}.out" ]; then
        printf "%-30s %s\n" "${f}" "$(awk '/^ *Utot\./{v=$2} END{print "Utot = " v}' "${name}.out")"
    else
        printf "%-30s %s\n" "${f}" "no ${name}.out"
    fi
done
