#!/usr/bin/env bash
# Collect total energies, SCF convergence flags and multigrid counts
# from a CUTOFF scan (cutoff_*Ry/) or a REL_CUTOFF scan (relcutoff_*Ry/).
#
# Usage: bash collect_scan.sh cutoff|relcutoff [natoms]
# Output: cutoff_scan.csv or relcutoff_scan.csv (reference = highest scanned value)
set -euo pipefail

kind=${1:?usage: bash collect_scan.sh cutoff|relcutoff [natoms]}
natoms=${2:-1}
case $kind in
  cutoff|relcutoff) prefix=$kind ;;
  *) echo "ERROR: first argument must be 'cutoff' or 'relcutoff'" >&2; exit 1 ;;
esac

out=${prefix}_scan.csv
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT

# directories sorted numerically by the scanned value
dirs=$(ls -d ${prefix}_*Ry 2>/dev/null | sed -E "s/^${prefix}_([0-9.]+)Ry$/\1 &/" | sort -n | awk '{print $2}')
[ -n "$dirs" ] || { echo "ERROR: no ${prefix}_*Ry directories found" >&2; exit 1; }

for d in $dirs; do
  val=${d#${prefix}_}; val=${val%Ry}
  f=$d/output.out
  if [ ! -f "$f" ]; then echo "skip $d: no output.out" >&2; continue; fi
  e=$(grep "ENERGY| Total FORCE_EVAL" "$f" | tail -n 1 | awk '{print $NF}')
  if [ -z "$e" ]; then echo "skip $d: no final energy (SCF aborted or job unfinished)" >&2; continue; fi
  if grep -q "SCF run converged" "$f"; then conv=yes; else conv=NO; fi
  counts=$(grep "count for grid" "$f" | awk '{print $5}' | paste -sd ';' -)
  echo "$val $e $conv ${counts:-n/a}" >> "$tmp"
done

[ -s "$tmp" ] || { echo "ERROR: no usable output found" >&2; exit 1; }

awk -v n="$natoms" -v kind="$prefix" '
  BEGIN { ha2mev = 27211.386 }
  { v[NR] = $1; e[NR] = $2; c[NR] = $3; g[NR] = $4 }
  END {
    ref = e[NR]
    printf "%s_Ry,energy_Ha,dE_meV_per_atom,scf_converged,grid_counts\n", kind
    for (i = 1; i <= NR; i++)
      printf "%s,%.10f,%.4f,%s,%s\n", v[i], e[i], (e[i] - ref) / n * ha2mev, c[i], g[i]
  }' "$tmp" > "$out"

if command -v column >/dev/null 2>&1; then column -s, -t "$out"; else cat "$out"; fi
echo "written: $out  (natoms = $natoms, reference = last row)"
