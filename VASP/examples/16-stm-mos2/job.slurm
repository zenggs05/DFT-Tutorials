#!/bin/bash
#SBATCH --job-name=mos2-stm
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

[ -s POTCAR ] || { echo "POTCAR missing: cat \$VASP_PP_PATH/potpaw_PBE/{Mo_pv,S}/POTCAR > POTCAR" >&2; exit 1; }

# 1) self-consistent run that writes WAVECAR
mkdir -p scf && cp INCAR_scf scf/INCAR && cp POSCAR KPOINTS POTCAR scf/
(cd scf && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)
grep "E-fermi" scf/OUTCAR | tail -1
grep -E "val. band max|cond. band min|fundamental gap" scf/OUTCAR | tail -3

# 2) partial charges of occupied and empty states (no SCF, very fast)
for w in occ empty; do
    mkdir -p parchg_$w
    cp INCAR_parchg_$w parchg_$w/INCAR
    cp POSCAR KPOINTS POTCAR parchg_$w/
    ln -sf ../scf/WAVECAR parchg_$w/WAVECAR
    (cd parchg_$w && mpirun -np $NP $VASP_EXE > vasp.out 2>&1)
    ls -l parchg_$w/PARCHG* 2>/dev/null
done

# 3) STM images from PARCHG (numpy + matplotlib)
for w in occ empty; do
    (cd parchg_$w && \
     python3 ../stm_from_parchg.py PARCHG --mode height  --height 3.0 --repeat 6 -o stm_height.png && \
     python3 ../stm_from_parchg.py PARCHG --mode current --iso 1e-4   --repeat 6 -o stm_current.png)
done
