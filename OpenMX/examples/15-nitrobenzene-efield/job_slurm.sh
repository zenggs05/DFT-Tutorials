#!/bin/bash
#SBATCH --job-name=omx-ch15-efield
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log

# Chapter 15: nitrobenzene in a uniform external field (scf.Electric.Field).
# Runs F = 0, +-2, +-10 GV/m, then dipoles/polarizability and difference densities.
# Override on submission, e.g.:  sbatch --export=ALL,FIELDS="-1 1" job_slurm.sh

# module load openmx/4.0            # EDIT: load your OpenMX / MPI environment

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

ulimit -s unlimited 2>/dev/null || true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

export OPENMX=${OPENMX:-openmx}
export NP=${SLURM_NTASKS:-4}
export NT=$OMP_NUM_THREADS
export FIELDS=${FIELDS:-"-10 -2 2 10"}
DIFF_GCUBE=${DIFF_GCUBE:-diff_gcube}   # EDIT: path to the compiled diff_gcube

RUN=1 bash run_fields.sh

# dipoles and polarizability (small fields only enter the linear fit)
args=("0:nben_f0.out")
for f in $FIELDS; do
    tag=$(awk -v f="$f" 'BEGIN{ if (f<0) printf "m%g", -f; else printf "p%g", f }')
    [ -f "nben_${tag}.out" ] && args+=("${f}:nben_${tag}.out")
done
python3 dipole_alpha.py "${args[@]}"

# field-induced density change: rho(F) - rho(0)
if command -v "$DIFF_GCUBE" > /dev/null 2>&1; then
    for tag in p10 m10; do
        [ -f "nben_${tag}.tden.cube" ] || continue
        "$DIFF_GCUBE" "nben_${tag}.tden.cube" nben_f0.tden.cube "drho_${tag}.cube"
        echo "wrote drho_${tag}.cube"
    done
else
    echo "diff_gcube not found; compile it with: gcc diff_gcube.c -lm -o diff_gcube"
fi
