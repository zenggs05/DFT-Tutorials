#!/bin/bash
# =====================================================================
# Slurm template for VASP 6.4-6.6 (CPU build, pure MPI).
#
#   sbatch slurm_vasp.sh                    # vasp_std
#   sbatch --export=ALL,VARIANT=gam slurm_vasp.sh   # Gamma-only build
#   sbatch --nodes=2 slurm_vasp.sh          # command line overrides #SBATCH
#
# Before the run the script checks that
#   * all four input files exist and have LF line endings,
#   * NCORE divides the number of MPI ranks per node,
#   * KPAR divides the total number of MPI ranks.
# The full-featured script of this repository is
# VASP/scripts/job_vasp_slurm.sh (GPU, OpenMP and restart options).
# =====================================================================
#SBATCH --job-name=vasp
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64         # = physical cores per node
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mem=0                      # all memory of the node (site dependent)
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
##SBATCH --account=my_project
##SBATCH --exclusive

# ----------------------------- user settings -------------------------
VARIANT=${VARIANT:-std}   # std | gam | ncl
LOG="vasp.out"

# ----------------------------- environment ---------------------------
module purge
# module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl
# module load vasp/6.5.1

ulimit -s unlimited
export OMP_NUM_THREADS=1

cd "${SLURM_SUBMIT_DIR}" || exit 1

case "${VARIANT}" in
    std|gam|ncl) EXE="vasp_${VARIANT}" ;;
    *) echo "VARIANT must be std, gam or ncl" >&2; exit 1 ;;
esac
command -v "${EXE}" >/dev/null || { echo "${EXE} not in PATH (module load?)" >&2; exit 1; }

# ----------------------------- sanity checks -------------------------
for f in INCAR POSCAR KPOINTS POTCAR; do
    [ -s "$f" ] || { echo "ERROR: $f missing or empty" >&2; exit 1; }
done
if grep -qU $'\r' INCAR POSCAR KPOINTS; then
    echo "ERROR: CRLF line endings in INCAR/POSCAR/KPOINTS -> dos2unix" >&2; exit 1
fi

# read an integer tag from INCAR (ignores comments; empty if not set)
incar_tag() {
    sed 's/[#!].*//' INCAR | tr ';' '\n' | \
        awk -F= -v k="$1" 'toupper($1) ~ "^[ \t]*"k"[ \t]*$" {gsub(/[ \t]/,"",$2); print $2}' | tail -n 1
}
NTOT=${SLURM_NTASKS}
# SLURM_NTASKS_PER_NODE is only set when --ntasks-per-node was given;
# otherwise use the average number of ranks per node
NPN=${SLURM_NTASKS_PER_NODE:-$(( NTOT / ${SLURM_JOB_NUM_NODES:-1} ))}
NCORE=$(incar_tag NCORE); KPAR=$(incar_tag KPAR)
echo "ranks total=${NTOT} per node=${NPN} NCORE=${NCORE:-1} KPAR=${KPAR:-1}"
if [ -n "${NCORE}" ] && [ "${NCORE}" -gt 0 ] && [ $(( NPN % NCORE )) -ne 0 ]; then
    echo "WARNING: NCORE=${NCORE} does not divide ${NPN} ranks per node" >&2
fi
if [ -n "${KPAR}" ] && [ "${KPAR}" -gt 0 ] && [ $(( NTOT % KPAR )) -ne 0 ]; then
    echo "WARNING: KPAR=${KPAR} does not divide ${NTOT} ranks" >&2
fi

# ----------------------------- run -----------------------------------
echo "Job ${SLURM_JOB_ID} ${EXE} on ${SLURM_JOB_NODELIST}, start $(date)"
srun --ntasks="${NTOT}" --cpus-per-task=1 --cpu-bind=cores "${EXE}" > "${LOG}" 2>&1
# mpirun alternative (needs the MPI that VASP was built with):
# mpirun -np "${NTOT}" "${EXE}" > "${LOG}" 2>&1
status=$?

echo "End $(date), exit code ${status}"
if grep -q "General timing and accounting" OUTCAR 2>/dev/null; then
    grep "free  energy   TOTEN" OUTCAR | tail -n 1
    grep -q "reached required accuracy" OUTCAR && echo "ionic relaxation converged"
else
    echo "OUTCAR has no timing block: VASP did not finish normally" >&2
    tail -n 20 "${LOG}" >&2
    [ ${status} -eq 0 ] && status=2
fi
exit ${status}
