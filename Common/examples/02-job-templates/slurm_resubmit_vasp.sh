#!/bin/bash
# =====================================================================
# Long VASP relaxation / MD that does not fit into one wall-time slot:
# stop VASP cleanly before the limit and resubmit the same script.
#
# How it works
#   --signal=B:USR1@900  Slurm sends USR1 to this batch shell 15 min
#                        before the time limit.
#   trap ... USR1        we write STOPCAR with LSTOP = .TRUE.; VASP then
#                        finishes the current ionic step, writes CONTCAR,
#                        WAVECAR, CHGCAR and stops normally.
#   next round           CONTCAR -> POSCAR, the job resubmits itself
#                        (at most MAX_ROUND times) until the relaxation
#                        is converged.
#
#   sbatch slurm_resubmit_vasp.sh
# LSTOP acts only at the end of the current ionic step: the lead time
# (@900 s) must be longer than one ionic step, otherwise Slurm kills
# VASP first. Slurm may send the signal up to 60 s earlier than asked.
# Keep LWAVE = .TRUE. (default): with the default ISTART VASP restarts
# from the WAVECAR of the previous round when it exists.
# Some sites forbid sbatch on compute nodes: then resubmit by hand.
# For MD (IBRION=0) the appended XDATCAR pieces are kept as XDATCAR.N.
# =====================================================================
#SBATCH --job-name=vasp-long
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --signal=B:USR1@900

MAX_ROUND=5
EXE="vasp_std"
SCRIPT="slurm_resubmit_vasp.sh"   # name of THIS file in the submit dir

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR}" || exit 1

ROUND=$(cat .round 2>/dev/null || echo 0)
ROUND=$((ROUND + 1)); echo "${ROUND}" > .round
echo "round ${ROUND}/${MAX_ROUND}, job ${SLURM_JOB_ID}"

# start from the last geometry of the previous round
if [ "${ROUND}" -gt 1 ] && [ -s CONTCAR ]; then
    cp POSCAR "POSCAR.$((ROUND - 1))"
    cp CONTCAR POSCAR
    for f in OUTCAR OSZICAR XDATCAR vasprun.xml; do
        [ -f "$f" ] && mv "$f" "$f.$((ROUND - 1))"
    done
fi
rm -f STOPCAR

stop_vasp() {
    echo "USR1: time limit is near, writing STOPCAR (LSTOP)"
    echo "LSTOP = .TRUE." > STOPCAR
}
trap stop_vasp USR1

srun --ntasks="${SLURM_NTASKS}" --cpus-per-task=1 "${EXE}" > "vasp.out.${ROUND}" 2>&1 &
PID=$!
# "wait" returns when a trapped signal arrives: keep waiting for VASP
while true; do
    wait "${PID}"; status=$?
    kill -0 "${PID}" 2>/dev/null || break
done
echo "VASP exit code ${status}"

if grep -q "reached required accuracy" OUTCAR 2>/dev/null; then
    echo "converged in round ${ROUND}"; rm -f .round; exit 0
fi
if [ -f STOPCAR ] && [ "${ROUND}" -lt "${MAX_ROUND}" ]; then
    echo "not converged yet: resubmitting"
    # $0 is a spool copy on the compute node: resubmit the original file
    sbatch "${SLURM_SUBMIT_DIR}/${SCRIPT}"
    exit 0
fi
echo "stopped without convergence (round ${ROUND}); check OUTCAR" >&2
exit 3
