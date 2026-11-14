#!/bin/bash
# vaspkit_batch.sh - run one VASPKIT 1.5.x task non-interactively in many calculation folders
#
# Usage:
#   bash vaspkit_batch.sh <task> <dir1> [dir2 ...]
#   bash vaspkit_batch.sh band  calc_*/band
#   bash vaspkit_batch.sh gap   strain_*/
#   bash vaspkit_batch.sh kpoints relax_*/           # writes KPOINTS (Gamma, kpr = 0.04)
#   KPR=0.03 bash vaspkit_batch.sh kpoints relax_*/
#
# Tasks (VASPKIT function numbers, see https://vaspkit.com/features.html):
#   kpoints  102  KPOINTS for SCF (Gamma-centred, resolution $KPR in 2*pi/A)
#   kpath    303  line-mode path for bulk: KPATH.in, PRIMCELL.vasp, HIGH_SYMMETRY_POINTS
#   hsekp    251  KPOINTS with zero-weight path for hybrid bands ($KPR mesh, $KPR_PATH path)
#   band     211  BAND.dat, BAND_REFORMATTED.dat, KLABELS, BAND_GAP
#   pband    213  element-projected bands PBAND_<El>.dat (needs PROCAR, LORBIT = 10/11)
#            211/213 take E-fermi from DOSCAR: after the band run, replace its DOSCAR
#            by the one of the SCF/DOS run (or set SET_FERMI_ENERGY_ZERO = .FALSE.)
#   hseband  252  bands from a zero-weight hybrid run (KPOINTS must come from 251:
#                  252-257 read the parameters 251 writes in the first KPOINTS line)
#   dos      111 + 113  TDOS.dat and PDOS_<El>.dat (needs LORBIT = 11 for 113)
#   gap      911  band gap summary
#   optics   711  linear optics from vasprun.xml (LOPTICS run), energy unit eV
#   workfunc 426  planar average of LOCPOT along z; prints the vacuum level
#   thermo-gas 502  thermal correction of a gas molecule ($TEMP K, $PRES atm, $MULT)
#   thermo-ads 501  thermal correction of an adsorbate ($TEMP K)
#
# The answers piped into VASPKIT follow the prompts of VASPKIT 1.5.1. Prompts can
# change between versions: run the task once interactively, compare, and edit the
# answer strings below if necessary (with AUTO_PLOT = .TRUE. some tasks ask an
# extra plotting question at the end; add its answer or switch AUTO_PLOT off
# for batch runs). Screen output is saved to vaspkit_<task>.log
# in each folder; the script then reads VASPKIT's data files rather than parsing
# the screen (whose layout changes between versions).
set -uo pipefail

TASK=${1:-}
shift || true
[ -n "$TASK" ] && [ $# -gt 0 ] || { sed -n '2,34p' "$0"; exit 1; }
command -v vaspkit >/dev/null 2>&1 || { echo "ERROR: vaspkit not in PATH"; exit 1; }
[ -f "$HOME/.vaspkit" ] || [ -f /etc/vaspkit ] || echo "WARNING: no ~/.vaspkit (or /etc/vaspkit) found"

KPR=${KPR:-0.04}            # k-mesh resolution, 2*pi/Angstrom
KPR_PATH=${KPR_PATH:-0.06}  # path resolution for hybrid bands
TEMP=${TEMP:-298.15}
PRES=${PRES:-1}
MULT=${MULT:-1}

answers () {    # print the stdin answers for a task
    case "$1" in
        kpoints)    printf '102\n2\n%s\n' "$KPR" ;;
        kpath)      printf '303\n' ;;
        hsekp)      printf '251\n2\n%s\n%s\n' "$KPR" "$KPR_PATH" ;;
        band)       printf '211\n' ;;
        pband)      printf '213\n' ;;
        hseband)    printf '252\n' ;;
        dos)        printf '111\n' ;;          # 113 is run as a second call below
        gap)        printf '911\n' ;;
        optics)     printf '711\n1\n' ;;
        workfunc)   printf '426\n3\n' ;;
        thermo-gas) printf '502\n%s\n%s\n%s\n' "$TEMP" "$PRES" "$MULT" ;;
        thermo-ads) printf '501\n%s\n' "$TEMP" ;;
        *) return 1 ;;
    esac
}

need_files () { # minimal input check per task
    case "$1" in
        kpoints|kpath|hsekp) echo "POSCAR" ;;
        band)                echo "EIGENVAL KPOINTS POSCAR DOSCAR" ;;
        pband)               echo "EIGENVAL KPOINTS POSCAR DOSCAR PROCAR" ;;
        hseband)             echo "EIGENVAL KPOINTS POSCAR" ;;
        dos)                 echo "DOSCAR POSCAR" ;;
        gap)                 echo "EIGENVAL POSCAR" ;;
        optics)              echo "vasprun.xml" ;;
        workfunc)            echo "LOCPOT" ;;
        thermo-gas)          echo "OUTCAR CONTCAR" ;;   # 502 reads the symmetry from CONTCAR
        thermo-ads)          echo "OUTCAR" ;;
    esac
}

answers "$TASK" >/dev/null || { echo "ERROR: unknown task '$TASK'"; exit 1; }
[ "$TASK" = "hsekp" ] && echo "NOTE: 251 reads KPATH.in in each folder; run the kpath task first."
[ "$TASK" = "hseband" ] && echo "NOTE: 252 only works with a KPOINTS file written by VASPKIT 251."
case "$TASK" in band|pband) echo "NOTE: E-fermi is read from DOSCAR; use the DOSCAR of the SCF/DOS run, not of the band run." ;; esac

printf "%-40s %-8s %s\n" "folder" "status" "result"
for d in "$@"; do
    [ -d "$d" ] || continue
    missing=""
    for f in $(need_files "$TASK"); do [ -s "$d/$f" ] || missing="$missing $f"; done
    if [ -n "$missing" ]; then
        printf "%-40s %-8s %s\n" "$d" "SKIP" "missing:$missing"
        continue
    fi
    log="vaspkit_${TASK}.log"
    ( cd "$d" && answers "$TASK" | vaspkit > "$log" 2>&1 )
    if [ "$TASK" = "dos" ]; then
        ( cd "$d" && printf '113\n' | vaspkit >> "$log" 2>&1 )
    fi

    # short result per task, read from data files where possible
    result=""
    case "$TASK" in
        kpoints)  result=$(sed -n '4p' "$d/KPOINTS" 2>/dev/null | xargs) ;;
        kpath)    result=$(ls "$d"/KPATH.in "$d"/PRIMCELL.vasp 2>/dev/null | xargs -n1 basename | xargs) ;;
        band|hseband) result=$(grep -i -m1 "band gap" "$d/BAND_GAP" 2>/dev/null | xargs) ;;
        gap)      result=$(grep -i -m1 "band gap" "$d/$log" 2>/dev/null | xargs) ;;
        dos)      result=$(ls "$d"/TDOS.dat "$d"/PDOS_*.dat 2>/dev/null | xargs -n1 basename | xargs) ;;
        workfunc) result=$(grep -i -m1 "vacuum-level" "$d/$log" 2>/dev/null | xargs) ;;
        thermo-*) result=$(grep -i -m1 "Thermal correction to G" "$d/$log" 2>/dev/null | xargs) ;;
        *)        result="see $d/$log" ;;
    esac
    if grep -qi "error" "$d/$log"; then
        printf "%-40s %-8s %s\n" "$d" "CHECK" "'error' in $log; ${result}"
    else
        printf "%-40s %-8s %s\n" "$d" "OK" "${result:-see $log}"
    fi
done
