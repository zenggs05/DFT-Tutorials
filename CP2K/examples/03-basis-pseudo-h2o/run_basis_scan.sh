#!/usr/bin/env bash
# run_basis_scan.sh - run the H2O basis-set ladder and print a summary table.
#
# Usage:
#     bash run_basis_scan.sh                  # uses cp2k.ssmp (or CP2K_EXE)
#     CP2K_EXE="mpirun -np 4 cp2k.psmp" bash run_basis_scan.sh
#
# For each basis the script runs h2o_szv_dzvp_tzv2p.inp with -E BASIS=... -E FAMILY=...
# and extracts: total energy (Ha), dipole moment (Debye), SCF steps, wall time (s).

set -u
CP2K_EXE=${CP2K_EXE:-cp2k.ssmp}
INP=h2o_szv_dzvp_tzv2p.inp
LADDER="SZV:MOLOPT DZVP:MOLOPT TZVP:MOLOPT TZV2P:MOLOPT TZV2PX:MOLOPT SZV:MOLOPT-SR DZVP:MOLOPT-SR"

printf "%-22s %20s %12s %8s %10s\n" "basis" "E_total[Ha]" "dipole[D]" "SCF" "time[s]"
for item in ${LADDER}; do
  b=${item%%:*}
  fam=${item##*:}
  out="h2o_${b}-${fam}.out"
  rm -f "${out}"                              # cp2k -o appends to existing files
  ${CP2K_EXE} -E BASIS="${b}" -E FAMILY="${fam}" -i "${INP}" -o "${out}" > /dev/null 2>&1
  e=$(grep "ENERGY| Total FORCE_EVAL" "${out}" | tail -1 | awk '{print $NF}')
  d=$(grep -A1 "Dipole moment \[Debye\]" "${out}" | tail -1 | awk '{print $NF}')
  n=$(grep "SCF run converged in" "${out}" | tail -1 | awk '{print $(NF-2)}')
  t=$(grep -E "^ CP2K +1 " "${out}" | tail -1 | awk '{print $NF}')
  printf "%-22s %20s %12s %8s %10s\n" "${b}-${fam}-GTH" "${e:-FAILED}" "${d:--}" "${n:--}" "${t:--}"
done
