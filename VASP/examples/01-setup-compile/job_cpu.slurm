#!/bin/bash
# Slurm job script for a CPU build of VASP (pure MPI).
# Adapt partition, module names and the VASP path to your cluster.
#SBATCH --job-name=vasp_cpu
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

# ---------------- environment (site specific) ----------------
module purge
module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl
export PATH=$HOME/software/vasp.6.6.1/bin:$PATH
VASP_EXE=${VASP_EXE:-vasp_std}          # vasp_std / vasp_gam / vasp_ncl

ulimit -s unlimited
export OMP_NUM_THREADS=1

# ---------------- sanity checks ----------------
cd "$SLURM_SUBMIT_DIR" || exit 1
for f in INCAR POSCAR POTCAR; do
    [ -s "$f" ] || { echo "ERROR: $f is missing or empty"; exit 1; }
done
if [ ! -s KPOINTS ] && ! grep -qi "^ *KSPACING" INCAR; then
    echo "WARNING: no KPOINTS file and no KSPACING in INCAR (VASP will use KSPACING = 0.5)"
fi
if grep -q $'\t' INCAR || grep -U -q $'\r' INCAR; then
    echo "WARNING: INCAR contains TAB or CR characters; run dos2unix INCAR"
fi

echo "Job $SLURM_JOB_ID on $(hostname), $SLURM_NTASKS MPI ranks, start $(date)"
mpirun -np "$SLURM_NTASKS" "$VASP_EXE" > vasp.out 2>&1
echo "End $(date)"

# ---------------- quick post-run check ----------------
if grep -q "General timing and accounting" OUTCAR 2>/dev/null; then
    echo "VASP finished normally"
else
    echo "VASP did NOT finish normally, check vasp.out and OUTCAR"
fi
grep "reached required accuracy" OUTCAR | tail -1
grep "energy(sigma->0)" OUTCAR | tail -1
