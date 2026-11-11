#!/bin/bash
# Non-interactive driver for OpenMX DosMain (total DOS + PDOS).
#
# usage: run_dosmain.sh PREFIX [tetra|gauss] [gauss_width_eV] [atom indices ...]
#   PREFIX     System.Name, i.e. PREFIX.Dos.val and PREFIX.Dos.vec must exist
#   method     tetra (default) or gauss
#   width      Gaussian width a in exp(-(E/a)^2), eV (default 0.1); also used when
#              DosMain chooses Gaussian broadening automatically (Kgrid 1 1 1, DC)
#   atoms      atoms for PDOS (default: all atoms)
# Set DOSMAIN=/path/to/DosMain if it is not in PATH. For the official Debian
# package first: export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:$LD_LIBRARY_PATH
set -e
PREFIX=$1
if [ -z "${PREFIX}" ]; then
  sed -n '2,13p' "$0"; exit 1
fi
METHOD=${2:-tetra}
WIDTH=${3:-0.1}
shift $(( $# < 3 ? $# : 3 ))
ATOMS="$*"
VAL=${PREFIX}.Dos.val
VEC=${PREFIX}.Dos.vec
DOSMAIN=${DOSMAIN:-DosMain}
[ -f "${VAL}" ] || { echo "missing ${VAL}"; exit 1; }

# header keywords of the .Dos.val file (read the same way by DosMain)
mode=$(awk 'tolower($1)=="mode"{print $2; exit}' "${VAL}")
mode=${mode:-1}
nk=$(awk 'tolower($1)=="kgrid"{print $2*$3*$4; exit}' "${VAL}")
nk=${nk:-0}
natom=$(awk 'tolower($1)=="atomnum"{print $2; exit}' "${VAL}")
[ -z "${ATOMS}" ] && ATOMS=$(seq -s ' ' 1 "${natom}")

answers() {   # $1 = 1 (DOS) or 2 (PDOS)
  if [ "${mode}" = "1" ] && [ "${nk}" != "1" ]; then
    # DosMain asks: Tetrahedron(1) or Gaussian Broadening(2)
    if [ "${METHOD}" = "gauss" ]; then echo 2; echo "${WIDTH}"; else echo 1; fi
  elif [ "${mode}" = "1" ] || [ "${mode}" = "5" ]; then
    # Kgrid 1 1 1 or divide-conquer: Gaussian is chosen automatically, width is asked
    echo "${WIDTH}"
  fi
  # modes 6 (NEGF) and 7 (DosGauss.fileout) ask nothing about the method
  echo "$1"
  if [ "$1" = "2" ] && [ "${natom}" -gt 1 ]; then echo "${ATOMS}"; fi
}

echo "== total DOS ==";  answers 1 | ${DOSMAIN} "${VAL}" "${VEC}" > ${PREFIX}.dosmain_dos.log
echo "== PDOS (atoms: ${ATOMS}) =="; answers 2 | ${DOSMAIN} "${VAL}" "${VEC}" > ${PREFIX}.dosmain_pdos.log
ls ${PREFIX}.DOS.* ${PREFIX}.PDOS.* 2>/dev/null | head -n 40
