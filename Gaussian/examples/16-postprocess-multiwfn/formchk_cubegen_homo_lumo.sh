#!/bin/bash
# formchk_cubegen_homo_lumo.sh - chk -> fchk, print HOMO/LUMO energies, write HOMO/LUMO cubes.
# Usage:  bash formchk_cubegen_homo_lumo.sh [file.chk] [npts]
#   file.chk  default: phenol.chk
#   npts      cubegen grid: 0 = 80^3 box, -2/-3/-4 = coarse/medium/fine (3/6/12 pts/Bohr). Default -3.
# Needs the Gaussian utilities formchk and cubegen in PATH (source $g16root/g16/bsd/g16.profile).
set -euo pipefail

chk=${1:-phenol.chk}
npts=${2:--3}
base=${chk%.chk}

for exe in formchk cubegen; do
    command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: $exe not found; source g16.profile first" >&2; exit 1; }
done
[ -f "$chk" ] || { echo "ERROR: $chk not found" >&2; exit 1; }

# 1) binary chk -> formatted fchk (skip if the fchk is newer than the chk)
if [ ! -f "$base.fchk" ] || [ "$chk" -nt "$base.fchk" ]; then
    formchk "$chk" "$base.fchk"
fi

# 2) HOMO/LUMO energies straight from the fchk (alpha orbitals; closed shell: HOMO = N_alpha)
awk '
/^Number of alpha electrons/ { nocc = $NF }
/^Alpha Orbital Energies/    { nmo = $NF; inside = 1; k = 0; next }
inside == 1 {
    for (i = 1; i <= NF; i++) { k++; e[k] = $i }
    if (k >= nmo) inside = 0
}
END {
    h = 27.211386
    printf "N_alpha = %d  N_MO = %d\n", nocc, nmo
    printf "HOMO (MO %d) = %10.4f eV\n", nocc, e[nocc] * h
    printf "LUMO (MO %d) = %10.4f eV\n", nocc + 1, e[nocc + 1] * h
    printf "gap         = %10.4f eV\n", (e[nocc + 1] - e[nocc]) * h
}' "$base.fchk"

# 3) orbital cubes (nprocs=0 is enough for orbitals; only ESP benefits from more cores)
cubegen 0 MO=Homo "$base.fchk" "${base}_homo.cub" "$npts" h
cubegen 0 MO=Lumo "$base.fchk" "${base}_lumo.cub" "$npts" h

# 4) quick sanity check: grid size lines (3-6) of the cube header
for c in "${base}_homo.cub" "${base}_lumo.cub"; do
    echo "== $c ($(du -h "$c" | cut -f1))"
    sed -n '3,6p' "$c"
done
echo "Open the cubes in GaussView or VMD; suggested orbital isovalue 0.02-0.05 a.u."
