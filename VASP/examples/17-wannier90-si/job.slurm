#!/bin/bash
#SBATCH --job-name=si-w90
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal      # edit for your cluster

module purge
# module load vasp/6.5.1   (must be built with -DVASP2WANNIER90)
# module load wannier90/3.1.0
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
W90_EXE=${W90_EXE:-wannier90.x}
NP=${SLURM_NTASKS:-16}

[ -s POTCAR ] || { echo "POTCAR missing: cat \$VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR" >&2; exit 1; }

# 1) SCF with symmetry (cheap)
mkdir -p scf && cp INCAR_scf scf/INCAR && cp KPOINTS_scf scf/KPOINTS && cp POSCAR POTCAR scf/
(cd scf && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)
EF=$(grep E-fermi scf/OUTCAR | tail -1 | awk '{print $3}')
echo "E-fermi (SCF) = $EF eV"

# 2) non-SCF on the full 6x6x6 mesh, VASP writes wannier90.amn/.mmn/.eig/.win
mkdir -p wannier && cp INCAR_wannier wannier/INCAR && cp KPOINTS_wan wannier/KPOINTS
cp POSCAR POTCAR scf/CHGCAR wannier/
cp wannier90.win wannier/wannier90.win
# energy windows relative to the Fermi energy of this calculation
FROZ=$(awk -v e="$EF" 'BEGIN{printf "%.3f", e+2.0}')
WIN=$(awk -v e="$EF" 'BEGIN{printf "%.3f", e+11.0}')
sed -i "s/^dis_froz_max .*/dis_froz_max    = $FROZ/; s/^dis_win_max .*/dis_win_max     = $WIN/" wannier/wannier90.win
grep -E "^dis_(froz|win)_max" wannier/wannier90.win
(cd wannier && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)
ls -l wannier/wannier90.{amn,mmn,eig,win}

# 3) standalone Wannier90 (serial is enough for Si)
(cd wannier && $W90_EXE wannier90)
grep -A 12 "Final State" wannier/wannier90.wout
ls wannier/wannier90_band.dat wannier/wannier90_hr.dat
