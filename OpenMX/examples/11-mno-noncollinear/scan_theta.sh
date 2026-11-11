#!/bin/bash
#SBATCH -J omx-mno-scan
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH -t 04:00:00
#SBATCH -o slurm-%j.out
#
# Constrained-DFT scan: rotate the spin of the whole MnO molecule (both atoms)
# from theta = 0 (perpendicular to the bond, along z) to theta = 90 (phi = 0,
# i.e. parallel to the bond along x), restarting every step from the previous one.
# With SOC on, Utot(theta) is the magnetic anisotropy curve of the molecule.
# The energy differences are tiny: keep scf.criterion very tight.
#
#   bash scan_theta.sh      or      sbatch scan_theta.sh

ulimit -s unlimited
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
NP=${SLURM_NTASKS:-4}
NAME=mno_nc                       # System.Name must stay the same for scf.restart on

echo "# theta(deg)  Utot(Hartree)  Ucs(Hartree)  Uele(Hartree)" > scan_theta.dat
first=1
for th in 0 15 30 45 60 75 90; do
  rs=on; [ ${first} -eq 1 ] && rs=off
  awk -v t=${th} -v rs=${rs} '
    /^<Atoms.SpeciesAndCoordinates/ {inb=1; print; next}
    /^Atoms.SpeciesAndCoordinates>/ {inb=0}
    /^scf.restart / {print "scf.restart                 " rs; next}
    { if (inb && NF>=12) {$8=t; $10=t}; print }' Mol_MnO_NC.dat > run_theta.dat
  mpirun -np ${NP} openmx run_theta.dat -nt ${OMP_NUM_THREADS} > ${NAME}.std 2>&1
  cp ${NAME}.out ${NAME}_theta${th}.out
  cp ${NAME}.std ${NAME}_theta${th}.std
  u=$(grep "Utot\." ${NAME}.out | awk '{print $2}')
  c=$(grep "Ucs\."  ${NAME}.out | awk '{print $2}')
  e=$(grep "Uele\." ${NAME}.out | awk '{print $2}')
  printf "%6s  %s  %s  %s\n" "${th}" "${u}" "${c}" "${e}" | tee -a scan_theta.dat
  first=0
done
echo "results in scan_theta.dat (1 Hartree = 27211.386 meV)"
