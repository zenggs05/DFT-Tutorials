#!/bin/bash
#SBATCH --job-name=gaas-band-dos
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
##SBATCH --partition=normal          # edit for your cluster

# SCF -> band (ICHARG=11, line mode) -> DOS (ICHARG=11, ISMEAR=-5) for GaAs.
# Submit from Common/examples/05-band-dos-plotting after creating ./POTCAR
# (Ga_d + As, PBE). Set RUN_HSE=1 to add the HSE06 zero-weight band run.
module purge
# module load vasp/6.6.1                  # edit: the module of your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=1
VASP_EXE=${VASP_EXE:-vasp_std}
MPI_RUN="mpirun -np ${SLURM_NTASKS:-16}"
RUN_HSE=${RUN_HSE:-0}
TOP=$(pwd)

[ -f POTCAR ] || { echo "POTCAR missing in $TOP"; exit 1; }

prep () {   # prep <dir> <INCAR> <KPOINTS>
    mkdir -p "$TOP/$1"
    cp "$TOP/$2" "$TOP/$1/INCAR"
    cp "$TOP/$3" "$TOP/$1/KPOINTS"
    cp "$TOP/POSCAR" "$TOP/POTCAR" "$TOP/$1/"
}

# 1) SCF
prep scf INCAR_scf KPOINTS_scf
cd "$TOP/scf" && $MPI_RUN $VASP_EXE > vasp.out 2>&1
grep -q "reached required accuracy\|aborting loop because EDIFF is reached" OUTCAR \
    || echo "WARNING: check SCF convergence in scf/OUTCAR"
[ -s CHGCAR ] || { echo "SCF did not write CHGCAR"; exit 1; }
grep "E-fermi" OUTCAR | tail -1

# 2) band structure along the line-mode path
prep band INCAR_band KPOINTS_band_linemode
cp "$TOP/scf/CHGCAR" "$TOP/band/"
cd "$TOP/band" && $MPI_RUN $VASP_EXE > vasp.out 2>&1

# 3) DOS / PDOS on a dense mesh
prep dos INCAR_dos KPOINTS_dos
cp "$TOP/scf/CHGCAR" "$TOP/dos/"
cd "$TOP/dos" && $MPI_RUN $VASP_EXE > vasp.out 2>&1

# 4) optional: HSE06 band structure (self-consistent, zero-weight path points)
if [ "$RUN_HSE" = "1" ]; then
    prep hse_band INCAR_hse_band KPOINTS_hse_zero_weight
    cd "$TOP/hse_band" && $MPI_RUN $VASP_EXE > vasp.out 2>&1
fi

cd "$TOP"
for d in scf band dos hse_band; do
    [ -f "$d/OUTCAR" ] || continue
    printf "%-9s %s\n" "$d" "$(grep -c 'General timing' $d/OUTCAR | sed 's/1/finished/;s/^0$/NOT finished/')"
done
echo "E-fermi used as reference (SCF / DOS):"
grep "E-fermi" scf/OUTCAR dos/OUTCAR | awk '{print $1, $4}'
