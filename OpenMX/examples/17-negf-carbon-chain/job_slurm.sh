#!/bin/bash
#SBATCH --job-name=negf-chain
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.log
# Usage: sbatch job_slurm.sh [input1.dat input2.dat ...]
# Runs the given inputs one after another and stops at the first failure.
# Default sequence: Step 1 (lead) -> Step 2+3 at 0 V -> Step 3 only -> Step 2+3 at 0.5 V.

set -u
INPUTS=${*:-"Lead-Chain.dat NEGF-Chain.dat NEGF-Chain_step3.dat NEGF-Chain_bias0.5.dat"}
OPENMX=${OPENMX:-openmx}
# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl   # site specific

NP=${SLURM_NTASKS:-8}
NT=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${NT}
export OMP_STACKSIZE=512M
ulimit -s unlimited

for INPUT in ${INPUTS}; do
    NAME=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${INPUT}")
    LOG=${INPUT%.dat}.std
    echo "=== ${INPUT} (System.Name ${NAME}) $(date)"
    mpirun -np "${NP}" "${OPENMX}" "${INPUT}" -nt "${NT}" > "${LOG}" 2>&1
    # OpenMX exits with status 0 even after input errors: check the text.
    if ! grep -q "The calculation was normally finished." "${LOG}"; then
        echo "${INPUT} did not finish normally:" >&2
        tail -n 20 "${LOG}" >&2
        exit 1
    fi
    grep -h "average current for up and down spins" "${NAME}.current" 2>/dev/null
    grep -h "average conductance for up and down spins" "${NAME}.conductance" 2>/dev/null
done
