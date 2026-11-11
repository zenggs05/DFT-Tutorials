#!/bin/bash
# Benchmark several MPI x OpenMP layouts of cp2k.psmp on one full node (Slurm).
# Chapter 18 of the CP2K tutorial.
#
#   sbatch job_mpi_omp_hybrid.sh                          # default input: broken_scf_fixed.inp
#   sbatch job_mpi_omp_hybrid.sh my.inp "64x1 32x2 16x4 8x8"
#
# For every layout RANKSxTHREADS the input is run in its own directory bench_<layout>/
# and the wall time of the run (last column of the "CP2K" line of the timing report)
# is collected in bench_summary.dat. Use a short but representative input (a few SCF
# steps / MD steps of the production system), not a toy molecule.
#
#SBATCH --job-name=cp2k-bench
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log

set -uo pipefail

INP=${1:-broken_scf_fixed.inp}
NCORES=${SLURM_CPUS_ON_NODE:-$(nproc)}
LAYOUTS=${2:-"${NCORES}x1 $((NCORES/2))x2 $((NCORES/4))x4 $((NCORES/8))x8"}

# ---- environment: adapt to your cluster ----
# module load cp2k/2026.2
# source /path/to/cp2k/install/cp2k_env
# export CP2K_DATA_DIR=/path/to/cp2k/data
export OMP_PLACES=cores
export OMP_PROC_BIND=close
ulimit -s unlimited || true

cd "${SLURM_SUBMIT_DIR:-.}"
command -v cp2k.psmp >/dev/null || { echo "cp2k.psmp not found in PATH" >&2; exit 1; }
cp2k.psmp -c "$INP" > /dev/null || { echo "syntax check of $INP failed" >&2; exit 1; }

printf "%-10s %6s %8s %14s %s\n" layout ranks threads walltime_s status > bench_summary.dat
for lay in $LAYOUTS; do
  np=${lay%x*}; nt=${lay#*x}
  [ "$np" -ge 1 ] || continue
  d=bench_${lay}
  rm -rf "$d"; mkdir -p "$d"; cp "$INP" "$d/"
  # copy coordinate files etc. referenced by the input if needed:
  # cp *.xyz "$d/" 2>/dev/null || true
  (
    cd "$d"
    export OMP_NUM_THREADS=$nt
    export SRUN_CPUS_PER_TASK=$nt
    srun --nodes=1 --ntasks="$np" --cpus-per-task="$nt" --cpu-bind=cores \
         cp2k.psmp -i "$INP" -o run.out > srun.log 2>&1
  )
  out=$d/run.out
  if grep -q "PROGRAM ENDED AT" "$out" 2>/dev/null; then
    t=$(awk '$1 == "CP2K" && NF >= 7 {t = $NF} END {print t}' "$out")
    st=ok
  else
    t=NA; st=failed
  fi
  # sanity check: the run must report the intended number of MPI ranks
  nr=$(grep -m1 "DBCSR| MPI: Number of processes" "$out" 2>/dev/null | awk '{print $NF}')
  [ "${nr:-0}" = "$np" ] || st="$st(MPI_ranks_reported=${nr:-none})"
  printf "%-10s %6s %8s %14s %s\n" "$lay" "$np" "$nt" "$t" "$st" >> bench_summary.dat
done
cat bench_summary.dat
