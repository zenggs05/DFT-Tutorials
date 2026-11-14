#!/bin/bash
# Scan the polynomial window of the CWF calculation without repeating the SCF.
#
# Prerequisite: GaAs_CWF_AO.dat has been run once, so that gaas_cwf_ao_rst/
# exists. Each scan step re-uses it with "scf.restart pp" (OpenMX 4.0: one
# diagonalization, scf.maxIter=1 and MD.Type=nomd are set automatically, the
# restart files are not overwritten). System.Name must stay the same, so the
# outputs of every step are renamed before the next step starts.
#
# Usage:  NP=8 bash scan_cwf_window.sh
# Result: cwf_window_scan.txt (window, DM_func), gaas_cwf_ao_w<i>.CWF.Band and
#         gaas_cwf_ao_w<i>.out; the original outputs are restored at the end.
set -euo pipefail

NP=${NP:-4}
NT=${NT:-1}
TEMPLATE=GaAs_CWF_AO.dat
NAME=gaas_cwf_ao

# e0 e1 e2 e3 (eV relative to the chemical potential), one window per line
WINDOWS=(
  "-14.5 -13.8  1.0  8.0"
  "-14.5 -13.8  0.0  6.0"
  "-14.5 -13.8  2.0 10.0"
  "-10.0  -9.0  1.0  8.0"
)

[ -d "${NAME}_rst" ] || { echo "run ${TEMPLATE} first (missing ${NAME}_rst/)"; exit 1; }

# The scan re-uses System.Name, so it overwrites the outputs of the original
# run. Keep them as ${NAME}_orig.* and put them back at the end.
EXTS="out Band CWF.Band CWF.Dis_vs_H CWF.Hop"
for ext in ${EXTS}; do
  if [ -f "${NAME}.${ext}" ]; then cp "${NAME}.${ext}" "${NAME}_orig.${ext}"; fi
done

echo "# i  e0 e1 e2 e3   DM_func  DM_func_per_CWF_spin" > cwf_window_scan.txt
i=0
for w in "${WINDOWS[@]}"; do
  i=$((i + 1))
  inp=scan_w${i}.dat
  sed -e "s/^scf.restart .*/scf.restart                 pp/" \
      -e "s/^CWF.disentangling.Erange .*/CWF.disentangling.Erange    ${w}/" \
      -e "s/^CWF.fileout.flag .*/CWF.fileout.flag            off/" \
      "${TEMPLATE}" > "${inp}"
  mpirun -np "${NP}" openmx "${inp}" -nt "${NT}" > "scan_w${i}.std" 2>&1 || echo "step ${i}: openmx returned an error, see scan_w${i}.std"
  dm=$( (grep "DM_func\." "${NAME}.out" 2>/dev/null || true) | tail -n 1 | awk '{print $2, $3}')
  echo "${i}  ${w}   ${dm:-NA}" >> cwf_window_scan.txt
  if [ -f "${NAME}.CWF.Band" ]; then cp "${NAME}.CWF.Band" "${NAME}_w${i}.CWF.Band"; fi
  if [ -f "${NAME}.out" ]; then cp "${NAME}.out" "${NAME}_w${i}.out"; fi
done
for ext in ${EXTS}; do
  if [ -f "${NAME}_orig.${ext}" ]; then cp "${NAME}_orig.${ext}" "${NAME}.${ext}"; fi
done
cat cwf_window_scan.txt
# Note: the window "-10 -9 ..." removes the As-4s band from the window while the
# guiding functions still contain As-4s: expect a clearly larger DM_func.
