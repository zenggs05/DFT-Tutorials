#!/bin/bash
#SBATCH --job-name=omx-ch15-polB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log

# Chapter 15: Born effective charge of Na in NaCl from Berry-phase polarization.
# 1) generate displaced inputs  2) SCF with HS.fileout on  3) polB on each scfout
# 4) finite-difference Born charge.
# Override on submission, e.g.:  sbatch --export=ALL,DISPS="-0.02 0.00 0.02" job_slurm.sh

# module load openmx/4.0            # EDIT: load your OpenMX / MPI environment

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

OPENMX=${OPENMX:-openmx}
POLB=${POLB:-polB}                  # EDIT: full path if polB is not in PATH
NP=${SLURM_NTASKS:-4}
export DISPS=${DISPS:-"-0.05 0.00 0.05"}

bash make_disp.sh

for inp in NaCl_dx_*.dat; do
    base=${inp%.dat}
    name=$(awk 'tolower($1)=="system.name"{print $2}' "$inp")
    echo "== $(date '+%F %T') SCF $inp"
    mpirun -np "$NP" "$OPENMX" "$inp" -nt "$OMP_NUM_THREADS" > "${base}.std" 2>&1
    if ! grep -q "The calculation was normally finished." "${base}.std"; then
        echo "SCF of $inp did not finish normally, skipping polB"; continue
    fi
    echo "== $(date '+%F %T') polB ${name}.scfout"
    mpirun -np "$NP" "$POLB" "${name}.scfout" < polB.in > "polB_${base}.out" 2>&1
done

python3 born_charge.py --atom Na --dir x \
    $(for u in $DISPS; do
          tag=$(awk -v u="$u" 'BEGIN{ if (u<0) printf "m%.2f", -u; else printf "p%.2f", u }')
          echo "${u}:polB_NaCl_dx_${tag}.out"
      done)
