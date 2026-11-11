#!/bin/bash
#SBATCH --job-name=si-unfold
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=08:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
NP=${SLURM_NTASKS:-32}

[ -s POTCAR ] || { echo "POTCAR missing: cat \$VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR" >&2; exit 1; }

# 0) supercell with one vacancy (ASE) and the unfolding k-points (easyunfold)
[ -s POSCAR_sc ] || python3 make_vacancy.py
[ -s KPOINTS_easyunfold ] || easyunfold generate POSCAR_prim POSCAR_sc KPOINTS_band --matrix "3 3 3"

# 1) supercell SCF
mkdir -p scf && cp INCAR_scf scf/INCAR && cp KPOINTS_scf scf/KPOINTS && cp POSCAR_sc scf/POSCAR && cp POTCAR scf/
(cd scf && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)

# 2) non-SCF on the k-points required for unfolding, keep WAVECAR
mkdir -p band && cp INCAR_band band/INCAR && cp KPOINTS_easyunfold band/KPOINTS
cp POSCAR_sc band/POSCAR && cp POTCAR scf/CHGCAR band/
(cd band && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)
EF=$(grep E-fermi band/OUTCAR | tail -1 | awk '{print $3}')

# 3) spectral weights and effective band structure
easyunfold unfold calculate band/WAVECAR
easyunfold unfold plot --eref "$EF" --emin -4 --emax 4 --out-file unfold.png
