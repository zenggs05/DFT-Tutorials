#!/bin/bash
# C60: cluster diagonalization with ELPA1 and ELPA2, and the O(N) DC method.
#
#   NP=8 NT=1 bash run_compare.sh
#
# Runs three calculations in this directory and prints a small table:
#   c60_elpa1  (C60_cluster_elpa.dat)
#   c60_elpa2  (same input with scf.eigen.lib elpa2, written to C60_cluster_elpa2.dat)
#   c60_dc     (C60_DC.dat)
# Environment: OPENMX (default openmx), MPIRUN (default mpirun), NP (default 4), NT (default 1)
set -u
OPENMX=${OPENMX:-openmx}
MPIRUN=${MPIRUN:-mpirun}
NP=${NP:-4}
NT=${NT:-1}
export OMP_NUM_THREADS=${NT}
export OMP_STACKSIZE=${OMP_STACKSIZE:-512M}
ulimit -s unlimited 2>/dev/null || true

# ELPA2 variant of the cluster input
awk '
  tolower($1)=="system.name"   {print "System.Name                 c60_elpa2"; next}
  tolower($1)=="scf.eigen.lib" {print "scf.eigen.lib               elpa2"; next}
  {print}' C60_cluster_elpa.dat > C60_cluster_elpa2.dat

run() {   # $1 = input file, $2 = System.Name
    echo ">>> $1"
    ${MPIRUN} -np "${NP}" "${OPENMX}" "$1" -nt "${NT}" > "$2.std" 2>&1
    grep -q "The calculation was normally finished." "$2.std" || echo "    $1 did NOT finish normally, see $2.std"
}
run C60_cluster_elpa.dat  c60_elpa1
run C60_cluster_elpa2.dat c60_elpa2
run C60_DC.dat            c60_dc

ref=$(awk '/^ *Utot\./{v=$2} END{print v}' c60_elpa1.out 2>/dev/null)
printf "\n%-10s %20s %10s %10s %18s\n" name Utot_Hartree scf_steps wall_s "dE_vs_elpa1_mHa"
for n in c60_elpa1 c60_elpa2 c60_dc; do
    if [ ! -f "${n}.out" ]; then
        printf "%-10s %20s\n" "${n}" "no .out file"
        continue
    fi
    u=$(awk '/^ *Utot\./{v=$2} END{print v}' "${n}.out")
    s=$(grep -c "^ *SCF= *[0-9]" "${n}.out")
    t=$(awk '/Total Computational Time/{v=$NF} END{print v}' "${n}.out")
    d=$(awk -v u="${u}" -v r="${ref}" 'BEGIN{ if (r=="" || u=="") print "-"; else printf "%.3f", (u-r)*1000 }')
    printf "%-10s %20s %10s %10s %18s\n" "${n}" "${u}" "${s}" "${t}" "${d}"
done
echo
echo "Expected: elpa1 and elpa2 agree to ~1e-8 Hartree; DC differs by a truncation error"
echo "of the order of mHartree (cf. manual Table 5, which uses a minimal basis)."
