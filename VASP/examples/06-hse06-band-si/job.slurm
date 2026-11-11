#!/bin/bash
#SBATCH --job-name=si-hse06
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

# Runs the HSE06 workflow for Si. Submit from VASP/examples/06-hse06-band-si
# after building ./POTCAR. Select steps with e.g.:  sbatch --export=ALL,STEPS="0 3" job.slurm
module purge
# module load vasp/6.6.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
MPI_RUN="mpirun -np ${SLURM_NTASKS:-16}"
STEPS=${STEPS:-"0 1 2 3"}
TOP=$(pwd)
[ -f POTCAR ] || { echo "POTCAR missing in $TOP"; exit 1; }

for step in $STEPS; do
  case $step in
  0)  # PBE SCF on the uniform mesh -> WAVECAR, IBZKPT
      cd "$TOP/0-pbe-scf" && cp "$TOP/POSCAR" "$TOP/POTCAR" .
      $MPI_RUN $VASP_EXE > vasp.out 2>&1 ;;
  1)  # HSE06 SCF from the PBE WAVECAR
      cd "$TOP/1-hse-scf" && cp "$TOP/POSCAR" "$TOP/POTCAR" "$TOP/0-pbe-scf/WAVECAR" .
      $MPI_RUN $VASP_EXE > vasp.out 2>&1
      grep -A12 "Band structure" OUTCAR ;;
  2)  # zero-weight method: build KPOINTS, PBE with the same list, then HSE06
      cd "$TOP/2-hse-band-zero-weight" && cp "$TOP/POSCAR" "$TOP/POTCAR" .
      python3 make_zero_weight_kpoints.py --ibzkpt "$TOP/0-pbe-scf/IBZKPT" --kpath KPATH --out KPOINTS
      cp INCAR.pbe INCAR && $MPI_RUN $VASP_EXE > vasp_pbe.out 2>&1
      cp INCAR.hse INCAR && $MPI_RUN $VASP_EXE > vasp_hse.out 2>&1 ;;
  3)  # KPOINTS_OPT method (VASP >= 6.3)
      cd "$TOP/3-hse-band-kpoints-opt" && cp "$TOP/POSCAR" "$TOP/POTCAR" "$TOP/0-pbe-scf/WAVECAR" .
      $MPI_RUN $VASP_EXE > vasp.out 2>&1
      grep -A12 "Band structure" OUTCAR ;;
  esac
done
cd "$TOP"
