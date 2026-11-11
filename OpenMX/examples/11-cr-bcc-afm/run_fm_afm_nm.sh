#!/bin/bash
#SBATCH -J omx-cr-mag
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH -t 02:00:00
#SBATCH -o slurm-%j.out
#
# Compare AFM, FM and NM states of bcc Cr with identical numerical parameters.
#   bash run_fm_afm_nm.sh          (workstation)  or  sbatch run_fm_afm_nm.sh
#   bash run_fm_afm_nm.sh collect  (only print the summary)

ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}

if [ "$1" != "collect" ]; then
  # FM: both atoms start with +3 muB
  sed -e 's/^System.Name .*/System.Name                 cr_fm/' \
      -e 's/^  2  Cr   0.5  0.5  0.5    5.5  8.5$/  2  Cr   0.5  0.5  0.5    8.5  5.5/' \
      Cr_bcc_AFM.dat > Cr_bcc_FM.dat
  # NM: spin polarization off, 7.0 + 7.0 electrons
  sed -e 's/^System.Name .*/System.Name                 cr_nm/' \
      -e 's/^scf.SpinPolarization .*/scf.SpinPolarization        off/' \
      -e 's/8\.5  5\.5$/7.0  7.0/' -e 's/5\.5  8\.5$/7.0  7.0/' \
      Cr_bcc_AFM.dat > Cr_bcc_NM.dat

  for s in AFM FM NM; do
    lc=$(echo $s | tr 'A-Z' 'a-z')
    mpirun -np ${NP} openmx Cr_bcc_${s}.dat -nt ${OMP_NUM_THREADS} > cr_${lc}.std 2>&1
  done
fi

echo "state   Utot(Hartree)        total_spin_moment(muB)"
for lc in afm fm nm; do
  e=$(grep "Utot\." cr_${lc}.out 2>/dev/null | awk '{print $2}')
  m=$(grep "Total spin moment" cr_${lc}.out 2>/dev/null | head -1 | awk '{print $NF}')
  printf "%-6s  %-20s %s\n" "$lc" "${e:-n/a}" "${m:-n/a}"
done
echo
echo "== AFM: per-atom Mulliken populations (Diff = local moment)"
grep -A 5 "Mulliken populations" cr_afm.out 2>/dev/null
