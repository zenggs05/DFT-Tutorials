#!/usr/bin/env bash
# Gaussian 16 basis-set / integration-grid scan of a RELATIVE energy
# (HCN -> HNC isomerisation, both optimised with the same model chemistry).
#
# Usage:
#   bash run_scan.sh prepare      # write <basis>_<grid>.gjf
#   bash run_scan.sh run          # g16 < x.gjf > x.log, one after another
#   bash run_scan.sh collect      # table of E(HCN), E(HNC), dE in kJ/mol -> basis_scan.dat
#
# Optional environment variables:
#   BASES="def2SVP def2TZVP def2QZVP cc-pVDZ cc-pVTZ cc-pVQZ aug-cc-pVTZ"
#   GRIDS="FineGrid UltraFine SuperFine"      (G16 default = UltraFine, G09 default = FineGrid)
#   G16=g16
set -euo pipefail

ACTION=${1:-}
BASES=${BASES:-"def2SVP def2TZVP def2QZVP cc-pVDZ cc-pVTZ cc-pVQZ aug-cc-pVTZ"}
GRIDS=${GRIDS:-"FineGrid UltraFine SuperFine"}
G16=${G16:-g16}
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

case "$ACTION" in
  prepare|run)
    for b in $BASES; do
      for g in $GRIDS; do
        tag=${b}_${g}
        sed -e "s/@BASIS@/$b/g" -e "s/@GRID@/$g/g" -e "s/@TAG@/$tag/g" template.gjf > "$tag.gjf"
        if [ "$ACTION" = run ]; then
          echo ">>> $tag"
          $G16 < "$tag.gjf" > "$tag.log" || echo "WARNING: $tag ended with an error" >&2
        fi
      done
    done ;;
  collect)
    printf "%-14s %-10s %18s %18s %12s\n" basis grid "E_HCN(Ha)" "E_HNC(Ha)" "dE(kJ/mol)" | tee basis_scan.dat
    for b in $BASES; do
      for g in $GRIDS; do
        log=${b}_${g}.log
        [ -f "$log" ] || continue
        # one "Normal termination" per Link1 step; take the last SCF Done of each step
        awk -v b="$b" -v g="$g" '
          /SCF Done/          { e = $5 }
          /Normal termination/ { n++; en[n] = e }
          /Error termination/  { bad = 1 }
          END {
            if (bad || n < 2) { printf "%-14s %-10s   incomplete or failed\n", b, g; exit }
            printf "%-14s %-10s %18.8f %18.8f %12.3f\n", b, g, en[1], en[2], (en[2]-en[1])*2625.4996
          }' "$log" | tee -a basis_scan.dat
      done
    done ;;
  *)
    sed -n '2,13p' "$0"; exit 1 ;;
esac
