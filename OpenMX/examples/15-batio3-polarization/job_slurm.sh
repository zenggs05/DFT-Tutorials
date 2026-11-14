#!/bin/bash
#SBATCH --job-name=omx-ch15-bto
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.log

# Chapter 15: spontaneous polarization of tetragonal BaTiO3 by the Berry phase.
#   STAGE=relax : relax the polar structure (BaTiO3_tetra_relax.dat)
#   STAGE=path  : build the lambda path, SCF + polB on every point, analyse
# Submit:  sbatch --export=ALL,STAGE=relax job_slurm.sh
#          sbatch --export=ALL,STAGE=path  job_slurm.sh

# module load openmx/4.0            # EDIT: load your OpenMX / MPI environment

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

OPENMX=${OPENMX:-openmx}
POLB=${POLB:-polB}                  # EDIT: full path if polB is not in PATH
NP=${SLURM_NTASKS:-4}
STAGE=${STAGE:-path}
LAMBDAS=${LAMBDAS:-"0 0.25 0.5 0.75 1.0"}

run_omx () {
    local inp=$1
    local base=${inp%.dat}
    echo "== $(date '+%F %T') OpenMX $inp"
    mpirun -np "$NP" "$OPENMX" "$inp" -nt "$OMP_NUM_THREADS" > "${base}.std" 2>&1
    grep -q "The calculation was normally finished." "${base}.std"
}

if [ "$STAGE" = "relax" ]; then
    run_omx BaTiO3_tetra_relax.dat || { echo "relaxation failed"; exit 1; }
    grep -A8 "Fractional coordinates of the final structure" bto_relax.out
    exit 0
fi

# path stage: use the relaxed end point if the relaxation output exists
FROM=()
if [ -f bto_relax.out ]; then
    FROM=(--from-out bto_relax.out)
fi
# shellcheck disable=SC2086
python3 make_path.py BaTiO3_cubic_ref.dat BaTiO3_tetra_polB.dat --lambdas $LAMBDAS ${FROM[@]+"${FROM[@]}"}

args=()
for inp in BaTiO3_lam*.dat; do
    name=$(awk 'tolower($1)=="system.name"{print $2}' "$inp")
    lam=${name#bto_lam}
    if ! run_omx "$inp"; then
        echo "SCF of $inp failed, skipped"
        continue
    fi
    echo "== $(date '+%F %T') polB ${name}.scfout"
    mpirun -np "$NP" "$POLB" "${name}.scfout" < polB.in > "polB_${name}.out" 2>&1
    args+=("${lam}:polB_${name}.out")
done

python3 polar_path.py "${args[@]}"
