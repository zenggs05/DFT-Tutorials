#!/bin/bash
#SBATCH --job-name=cu-fermi
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.5.1
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
NP=${SLURM_NTASKS:-16}

[ -s POTCAR ] || { echo "POTCAR missing: cat \$VASP_PP_PATH/potpaw_PBE/Cu/POTCAR > POTCAR" >&2; exit 1; }

mkdir -p scf && cp INCAR_scf scf/INCAR && cp KPOINTS_scf scf/KPOINTS && cp POSCAR POTCAR scf/
(cd scf && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)

mkdir -p fs && cp INCAR_fs fs/INCAR && cp KPOINTS_fs fs/KPOINTS && cp POSCAR POTCAR scf/CHGCAR fs/
(cd fs && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)
grep E-fermi fs/OUTCAR | tail -1
grep "irreducible k-points" fs/OUTCAR | head -1

# Fermi surface with IFermi (pip install ifermi); interactive HTML needs no extra package
cd fs && ifermi info && ifermi plot --output cu_fermi_surface.html
