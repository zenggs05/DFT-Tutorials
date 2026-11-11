#!/bin/bash
# =====================================================================
# One step of the relax -> static -> band chain for Si (VASP).
# The step is chosen with the environment variable STEP:
#   sbatch --export=ALL,STEP=relax  job_step.sh
#   sbatch --export=ALL,STEP=static job_step.sh
#   sbatch --export=ALL,STEP=band   job_step.sh
# chain_relax_static_band.sh submits all three with dependencies.
# The script exits with a non-zero code when the step did not finish
# properly, so that "afterok" dependencies stop the chain.
# =====================================================================
#SBATCH --job-name=si-step
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out

module purge
# module load vasp/6.5.1          # site specific
ulimit -s unlimited
export OMP_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}" || exit 1
STEP=${STEP:?set STEP=relax|static|band with sbatch --export}
[ -f POTCAR ] || { echo "POTCAR missing (see README)" >&2; exit 1; }

run_vasp() {
    srun --ntasks="${SLURM_NTASKS}" --cpus-per-task=1 vasp_std > vasp.out 2>&1
}

finished_ok() {   # VASP writes this block only at a normal end
    grep -q "General timing and accounting" OUTCAR 2>/dev/null
}

case "${STEP}" in
relax)
    mkdir -p 1-relax && cd 1-relax || exit 1
    cp ../POSCAR ../KPOINTS ../POTCAR .
    cp ../INCAR.relax INCAR
    run_vasp
    finished_ok || exit 2
    grep -q "reached required accuracy" OUTCAR || {
        echo "relaxation NOT converged (NSW reached?)" >&2; exit 3; }
    ;;
static)
    mkdir -p 2-static && cd 2-static || exit 1
    cp ../1-relax/CONTCAR POSCAR
    cp ../KPOINTS ../POTCAR .
    cp ../INCAR.static INCAR
    run_vasp
    finished_ok || exit 2
    [ -s CHGCAR ] || { echo "CHGCAR not written" >&2; exit 3; }
    ;;
band)
    mkdir -p 3-band && cd 3-band || exit 1
    cp ../2-static/POSCAR ../2-static/CHGCAR ../POTCAR .
    cp ../KPOINTS.band KPOINTS
    cp ../INCAR.band INCAR
    run_vasp
    finished_ok || exit 2
    ;;
*)
    echo "unknown STEP=${STEP}" >&2; exit 1 ;;
esac

grep "free  energy   TOTEN" OUTCAR | tail -n 1
echo "step ${STEP} finished OK in $(pwd)"
