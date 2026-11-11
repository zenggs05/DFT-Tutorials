#!/bin/bash
#
# Summarise OpenMX molecule calculations in the current directory.
# For every <name>.out it prints: normal termination (from <name>.std),
# number of SCF steps, Utot (Hartree and eV), HOMO/LUMO (eV), gap (eV),
# |dipole| (Debye) and the Mulliken total vs. ideal electron count.
#
# Usage: bash check_results.sh [file1.out file2.out ...]
#
set -u
HA2EV=27.211386
REF_METREF=-8.033540953061   # Utot of work/Methane.dat in Sec. 5 of the OpenMX 4.0 manual

files=("$@")
if [ ${#files[@]} -eq 0 ]; then
  shopt -s nullglob
  files=(*.out)
fi
if [ ${#files[@]} -eq 0 ]; then
  echo "no .out files found"; exit 1
fi

printf "%-10s %-6s %5s %18s %14s %10s %10s %8s %9s %s\n" \
  "name" "end" "nSCF" "Utot(Ha)" "Utot(eV)" "HOMO(eV)" "LUMO(eV)" "gap(eV)" "dip(D)" "MulP/ideal"

for out in "${files[@]}"; do
  name=${out%.out}
  std=${name}.std
  ended="?"
  if [ -f "${std}" ]; then
    if grep -q "normally finished" "${std}"; then ended="ok"; else ended="NO"; fi
  fi

  nscf=$(awk '/SCF history at MD= *1[[:space:]]*$/{f=1; next} f && /^ *SCF= *[0-9]+ +NormRD/{n++} f && /Total energy \(Hartree\)/{exit} END{print n+0}' "${out}")
  utot=$(awk '$1=="Utot."{print $2; exit}' "${out}")

  read -r homo_ev lumo_ev gap_ev < <(awk -v c=${HA2EV} '
    /Eigenvalues \(Hartree\) for SCF KS-eq\./ {blk=1; next}
    blk && /HOMO =/ {homo=$3; next}
    blk && /Mulliken populations/ {blk=0}
    blk && NF==3 && $1 ~ /^[0-9]+$/ {e[$1]=$2}
    END {
      if (homo=="" || !(homo in e)) {print "NA NA NA"; exit}
      h=e[homo]*c
      if ((homo+1) in e) {l=e[homo+1]*c; printf "%.4f %.4f %.4f\n", h, l, l-h}
      else {printf "%.4f NA NA\n", h}
    }' "${out}")

  dip=$(awk '/Absolute D/{print $3; exit}' "${out}")
  mulp=$(awk '/Sum of MulP/{getline; for(i=1;i<=NF;i++){if($i ~ /^total=/){t=$(i+1)} if($i ~ /ideal/){d=$(i+1)}} print t"/"d; exit}' "${out}")

  if [ -n "${utot}" ]; then
    utot_ev=$(awk -v u="${utot}" -v c=${HA2EV} 'BEGIN{printf "%.6f", u*c}')
  else
    utot="NA"; utot_ev="NA"
  fi
  printf "%-10s %-6s %5s %18s %14s %10s %10s %8s %9s %s\n" \
    "${name}" "${ended}" "${nscf}" "${utot}" "${utot_ev}" "${homo_ev}" "${lumo_ev}" "${gap_ev}" "${dip:-NA}" "${mulp:-NA}"

  if [ "${name}" = "metref" ] && [ "${utot}" != "NA" ]; then
    awk -v u="${utot}" -v r=${REF_METREF} 'BEGIN{d=u-r; if(d<0)d=-d;
      printf "           |Utot - manual reference| = %.2e Ha  -> %s\n", d, (d<1e-6?"OK":"CHECK SETTINGS")}'
  fi
done

maxiter_note=$(for f in *.dat; do awk 'tolower($1)=="scf.maxiter"{print FILENAME": scf.maxIter="$2}' "$f"; done 2>/dev/null)
echo
echo "If nSCF equals scf.maxIter the SCF did NOT converge:"
echo "${maxiter_note}"
