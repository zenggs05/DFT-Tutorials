#!/bin/bash
#SBATCH -J omx-co-cu-cdd
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH -t 04:00:00
#SBATCH -o slurm-%j.out
#
# Charge density difference rho(AB) - rho(A) - rho(B) for CO on Cu(111).
#   1) SCF of AB, read Grid_Origin from its standard output
#   2) SCF of A and B with "scf.fixed.grid" = Grid_Origin of AB
#   3) add_gcube / diff_gcube (OpenMX tools) and planar average (cdd_planar.py)
#
# Usage:  sbatch run_cdd.sh      (or: bash run_cdd.sh on a workstation)
# Needs:  openmx in PATH; add_gcube and diff_gcube compiled from openmx4.0/source
#         (gcc add_gcube.c -lm -o add_gcube ; gcc diff_gcube.c -lm -o diff_gcube)

set -e
# module load openmx/4.0.1          # adapt to your cluster
ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}
OMX="mpirun -np ${NP} openmx"

# ---- step 1: composite system AB ----
${OMX} CO_on_Cu111_AB.dat -nt ${OMP_NUM_THREADS} > co_cu_AB.std 2>&1
ORIG=$(grep "Grid_Origin" co_cu_AB.std | tail -1 | awk '{print $2, $3, $4}')
if [ -z "${ORIG}" ]; then
  echo "Grid_Origin not found in co_cu_AB.std (level.of.stdout must be >= 1)"; exit 1
fi
echo "Grid_Origin of AB (Bohr): ${ORIG}"

# ---- step 2: fragments A (CO) and B (slab) on the same grid ----
for s in A B; do
  sed "s|^#scf.fixed.grid.*|scf.fixed.grid   ${ORIG}|" CO_on_Cu111_${s}.dat > CO_on_Cu111_${s}_run.dat
  grep "^scf.fixed.grid" CO_on_Cu111_${s}_run.dat
  ${OMX} CO_on_Cu111_${s}_run.dat -nt ${OMP_NUM_THREADS} > co_cu_${s}.std 2>&1
done

# ---- sanity check: grid origin (line 3) and grid vectors (lines 4-6) must be identical ----
for s in AB A B; do
  awk 'NR==3{print $2, $3, $4} NR>=4 && NR<=6' co_cu_${s}.tden.cube > hdr_${s}.txt
done
if cmp -s hdr_AB.txt hdr_A.txt && cmp -s hdr_AB.txt hdr_B.txt; then
  echo "grid headers identical: OK"
else
  echo "WARNING: grid headers differ -> the difference density is meaningless"; exit 1
fi
rm -f hdr_*.txt

# ---- step 3: rho(AB) - [rho(A) + rho(B)] ----
add_gcube  co_cu_A.tden.cube  co_cu_B.tden.cube   co_cu_A_B.tden.cube
diff_gcube co_cu_AB.tden.cube co_cu_A_B.tden.cube co_cu_dAB.cube

# ---- step 4: planar average along c and cumulative charge transfer ----
python3 cdd_planar.py co_cu_dAB.cube --png cdd_z.png > cdd_z.dat
echo "done: co_cu_dAB.cube, cdd_z.dat, cdd_z.png"
