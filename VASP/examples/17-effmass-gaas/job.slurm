#!/bin/bash
#SBATCH --job-name=gaas-effmass
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
NP=${SLURM_NTASKS:-16}

[ -s POTCAR ] || { echo "POTCAR missing: cat \$VASP_PP_PATH/potpaw_PBE/{Ga_d,As}/POTCAR > POTCAR" >&2; exit 1; }

# 1) SCF
mkdir -p scf && cp INCAR_scf scf/INCAR && cp KPOINTS_scf scf/KPOINTS && cp POSCAR POTCAR scf/
(cd scf && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)

# 2) non-SCF on the finite-difference k-points
mkdir -p fd && cp INCAR_nscf fd/INCAR && cp KPOINTS_fd fd/KPOINTS && cp POSCAR POTCAR fd/
cp scf/CHGCAR fd/
(cd fd && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)

# 3) effective masses of CBM and VBM at Gamma
(cd fd && python3 ../fit_effmass.py --eigenval EIGENVAL --poscar POSCAR | tee effmass.txt)
