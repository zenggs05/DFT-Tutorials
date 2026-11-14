#!/bin/bash
# Generic PBS (PBS Pro / OpenPBS) template for one CP2K input. Chapter 17 of the CP2K tutorial.
#   qsub -v INP=my.inp job_pbs.sh                 # output: my.out
#   qsub -v INP=my.inp,NTHREADS=4 job_pbs.sh      # MPI x OpenMP (adapt select line too)
#   qsub -v ASE_SCRIPT=ase_h2o_opt.py job_pbs.sh  # run an ASE driver script instead
#
# Resource line: mpiprocs x ompthreads must equal ncpus of the chunk.
# Torque users: replace the select line by  #PBS -l nodes=1:ppn=32
#
#PBS -N cp2k
#PBS -l select=1:ncpus=32:mpiprocs=32:ompthreads=1
#PBS -l walltime=04:00:00
#PBS -j oe

set -euo pipefail

# ---- environment: adapt to your cluster ----
# module load cp2k/2026.2
# source /path/to/cp2k/install/cp2k_env
# conda activate cp2k_env
# export CP2K_DATA_DIR=/path/to/cp2k/data
export OMP_NUM_THREADS=${NTHREADS:-1}
export OMP_PLACES=cores
export OMP_PROC_BIND=close
ulimit -s unlimited || true

cd "${PBS_O_WORKDIR:-.}"
command -v cp2k.psmp >/dev/null || { echo "cp2k.psmp not found in PATH" >&2; exit 1; }

# number of MPI ranks = lines of the node file (one line per mpiproc in PBS Pro)
if [ -n "${PBS_NODEFILE:-}" ] && [ -f "$PBS_NODEFILE" ]; then
  NP=$(wc -l < "$PBS_NODEFILE")
else
  NP=1
fi

if [ -n "${ASE_SCRIPT:-}" ]; then
  # ASE drives CP2K through the shell mode; stdin must reach MPI rank 0
  export ASE_CP2K_COMMAND=${ASE_CP2K_COMMAND:-"mpirun -np $NP cp2k.psmp -s"}
  python "$ASE_SCRIPT"
  exit $?
fi

inp=${INP:?usage: qsub -v INP=input.inp job_pbs.sh}
out=${inp%.inp}.out
# "cp2k -c" prints "SUCCESS, the input could be parsed correctly." on success
# (output captured first: with pipefail, "cmd | grep -q" can fail through SIGPIPE)
chk=$(cp2k.psmp -c "$inp" 2>&1 || true)
if ! grep -q "SUCCESS, the input could be parsed" <<< "$chk"; then
  echo "syntax check failed for $inp:" >&2
  echo "$chk" | tail -20 >&2
  exit 1
fi
rm -f "$out"                          # -o appends to an existing file

echo "start $(date)  ranks=$NP  threads=$OMP_NUM_THREADS  host=$(hostname)"
mpirun -np "$NP" cp2k.psmp -i "$inp" -o "$out" || true

if grep -q "PROGRAM ENDED AT" "$out"; then
  echo "finished normally $(date)"
  grep -m1 "DBCSR| MPI: Number of processes" "$out" || true
  grep "ENERGY| Total FORCE_EVAL" "$out" | tail -1 || true
else
  echo "CP2K did not finish: check the tail of $out" >&2
  ls -t ./*-1.restart ./*-RESTART.wfn ./*-RESTART.kp 2>/dev/null | head -3 | \
    sed 's/^/restart candidate: /' >&2 || true
  exit 1
fi
