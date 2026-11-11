#!/bin/bash
# Fit two-stage RESP charges with AmberTools antechamber from the Gaussian ESP job.
# Usage: bash run_antechamber_resp.sh [gesp_or_log] [net_charge] [residue_name]
# Requires: AmberTools (antechamber, resp) in PATH.
set -euo pipefail

INPUT=${1:-paracetamol.gesp}
NETCHARGE=${2:-0}
RESNAME=${3:-LIG}
OUT=paracetamol_resp.mol2

command -v antechamber > /dev/null || { echo "antechamber not found; load AmberTools first" >&2; exit 1; }
[ -f "$INPUT" ] || { echo "input file $INPUT not found" >&2; exit 1; }

case "$INPUT" in
    *.gesp) FMT=gesp ;;
    *.log|*.out) FMT=gout ;;
    *) echo "unknown input format: $INPUT (expect .gesp or .log)" >&2; exit 1 ;;
esac

antechamber -i "$INPUT" -fi "$FMT" -o "$OUT" -fo mol2 \
            -c resp -nc "$NETCHARGE" -rn "$RESNAME" -at gaff2

echo "RESP charges written to $OUT (last column of the @<TRIPOS>ATOM block):"
awk '/@<TRIPOS>ATOM/{f=1; next} /@<TRIPOS>/{f=0} f{printf "%4s %-4s %9.4f\n", $1, $2, $NF; s+=$NF} END{printf "sum = %.4f\n", s}' "$OUT"
