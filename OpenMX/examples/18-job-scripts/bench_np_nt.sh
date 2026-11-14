#!/bin/bash
# Scan MPI processes x OpenMP threads for one OpenMX input and tabulate the timings.
#
#   bash bench_np_nt.sh input.dat "4x1 2x2 1x4 8x1 4x2 2x4"
#   sbatch --nodes=1 --ntasks=8 --cpus-per-task=1 --wrap "bash bench_np_nt.sh input.dat '8x1 4x2 2x4'"
#
# Each entry "NPxNT" runs   mpirun -np NP openmx <copy of input> -nt NT
# with System.Name changed to <name>_npNP_ntNT, so the runs do not overwrite each other.
# Results go to bench_<name>.csv:
#   np, nt, cores, total_time_s, diag_time_s, scf_steps, Utot_Hartree, status
# Use a SHORT test (e.g. scf.maxIter 5-10, level.of.fileout 0): only the relative
# timings matter, and Utot must be identical for all runs that finished.
#
# Environment: OPENMX (default: openmx), MPIRUN (default: mpirun), MPI_EXTRA (extra mpirun options)
set -u
INPUT=${1:?usage: bash bench_np_nt.sh input.dat "NPxNT NPxNT ..."}
LIST=${2:-"1x1 2x1 4x1 2x2"}
OPENMX=${OPENMX:-openmx}
MPIRUN=${MPIRUN:-mpirun}
[ -f "${INPUT}" ] || { echo "input ${INPUT} not found"; exit 1; }

NAME=$(awk 'tolower($1)=="system.name"{print $2; exit}' "${INPUT}")
[ -n "${NAME}" ] || { echo "System.Name not found in ${INPUT}"; exit 1; }
CSV=bench_${NAME}.csv
echo "np,nt,cores,total_time_s,diag_time_s,scf_steps,Utot_Hartree,status" > "${CSV}"

ulimit -s unlimited 2>/dev/null || echo "warning: could not set ulimit -s unlimited"
export OMP_STACKSIZE=${OMP_STACKSIZE:-512M}

for item in ${LIST}; do
    np=${item%x*}
    nt=${item#*x}
    tag=${NAME}_np${np}_nt${nt}
    dat=${tag}.dat
    # same input, new System.Name (case-insensitive keyword match)
    awk -v n="${tag}" 'tolower($1)=="system.name"{print "System.Name   " n; next} {print}' "${INPUT}" > "${dat}"
    export OMP_NUM_THREADS=${nt}
    echo ">>> ${tag}: ${MPIRUN} -np ${np} ${OPENMX} ${dat} -nt ${nt}"
    # MPI_EXTRA is intentionally unquoted so that several options can be given
    # shellcheck disable=SC2086
    ${MPIRUN} -np "${np}" ${MPI_EXTRA:-} "${OPENMX}" "${dat}" -nt "${nt}" > "${tag}.std" 2>&1
    out=${tag}.out
    if grep -q "The calculation was normally finished." "${tag}.std" && [ -f "${out}" ]; then
        status=ok
    else
        status=FAILED
    fi
    ttot=$(awk '/Total Computational Time/{v=$NF} END{print v}' "${out}" 2>/dev/null)
    # first "Diagonalization" line of the timing table (max over ranks = last column)
    tdiag=$(awk '/^ *Diagonalization/{v=$NF} END{print v}' "${out}" 2>/dev/null)
    nscf=$(grep -c "^ *SCF= *[0-9]" "${out}" 2>/dev/null)
    utot=$(awk '/^ *Utot\./{v=$2} END{print v}' "${out}" 2>/dev/null)
    echo "${np},${nt},$((np*nt)),${ttot:-},${tdiag:-},${nscf:-0},${utot:-},${status}" >> "${CSV}"
done

echo
column -s, -t "${CSV}" 2>/dev/null || cat "${CSV}"
echo
echo "Check: all 'ok' rows must have the same Utot (to ~1e-8 Hartree); pick the fastest"
echo "combination per core count and remember that memory per node also matters."
