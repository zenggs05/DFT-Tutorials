#!/bin/bash
# Slurm job for the chapter-16 examples (xTB, DFTB, QM/MM).
#   sbatch job_slurm.sh                        # run the default list below
#   sbatch job_slurm.sh dftb_md.inp            # run only the given input(s)
#   TB_FLAGS="-E NEWTB=0" sbatch job_slurm.sh  # CP2K 2025.2 / 2026.1 syntax for tblite
#
#SBATCH --job-name=cp2k-ch16
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.log

set -euo pipefail

# ---- environment: adapt to your cluster (see chapter 1 of the CP2K tutorial) ----
# module load cp2k/2026.2
# source /path/to/cp2k/install/cp2k_env
# export CP2K_DATA_DIR=/path/to/cp2k/data
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-1}
ulimit -s unlimited || true

cd "${SLURM_SUBMIT_DIR:-.}"
command -v cp2k.psmp >/dev/null || { echo "cp2k.psmp not found in PATH" >&2; exit 1; }

if [ $# -gt 0 ]; then
  inputs="$*"
else
  inputs="ethanol_gfn1_native.inp ice_gfn2_tblite_kpoints.inp dftb_md.inp water_qmmm_tip3p.inp"
  # the AMBER variant needs the tleap output first
  [ -f water125.prmtop ] && inputs="$inputs water_qmmm_amber.inp"
fi

# tblite is optional at compile time: skip the GFN2 example if it is missing
has_tblite=0
cp2k.psmp -v 2>/dev/null | grep -qi tblite && has_tblite=1

# extra preprocessor flags for the tblite example: "-E NEWTB=0" on CP2K 2025.2 / 2026.1
TB_FLAGS=${TB_FLAGS:-}

for inp in $inputs; do
  extra=""
  if [ "$inp" = "ice_gfn2_tblite_kpoints.inp" ]; then
    if [ "$has_tblite" -eq 0 ]; then
      echo "=== skip $inp: this cp2k.psmp was built without tblite"
      continue
    fi
    extra=$TB_FLAGS
  fi
  out=${inp%.inp}.out
  rm -f "$out"                     # -o appends to an existing file
  echo "=== $inp -> $out  ($(date))"
  # shellcheck disable=SC2086
  srun cp2k.psmp $extra -i "$inp" -o "$out" || true
  if grep -q "PROGRAM ENDED AT" "$out"; then
    echo "    finished normally"
    grep "ENERGY| Total FORCE_EVAL" "$out" | tail -1 || true
  else
    echo "    ERROR: $inp did not finish, check the end of $out" >&2
  fi
done
