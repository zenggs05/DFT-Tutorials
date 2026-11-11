#!/bin/bash
# density_esp_cube.sh - total density cube and an ESP cube on the SAME grid (for ESP-mapped surfaces).
# Usage:  bash density_esp_cube.sh [file.fchk] [nprocs] [npts]
#   file.fchk default phenol.fchk (created from phenol.chk with formchk if missing)
#   nprocs    cores used by cubegen for the ESP (default 4)
#   npts      grid of the density cube (default -3 = 6 points/Bohr)
set -euo pipefail

fchk=${1:-phenol.fchk}
np=${2:-4}
npts=${3:--3}
base=${fchk%.fchk}

if [ ! -f "$fchk" ] && [ -f "$base.chk" ]; then
    formchk "$base.chk" "$fchk"
fi
[ -f "$fchk" ] || { echo "ERROR: $fchk not found" >&2; exit 1; }

# 1) total SCF density; the header written here ('h') defines the grid
cubegen 0 Density=SCF "$fchk" "${base}_dens.cub" "$npts" h

# 2) ESP on exactly the same grid: npts must be -1 and the template cube must have a header
cubegen "$np" Potential=SCF "$fchk" "${base}_esp.cub" -1 h "${base}_dens.cub"

# 3) the grid definitions (lines 3-6) must be identical
if diff <(sed -n '3,6p' "${base}_dens.cub") <(sed -n '3,6p' "${base}_esp.cub") >/dev/null; then
    echo "OK: ${base}_dens.cub and ${base}_esp.cub share the same grid"
else
    echo "WARNING: grids differ - the ESP-mapped surface will be wrong" >&2
fi
ls -lh "${base}_dens.cub" "${base}_esp.cub"
echo "VMD: isosurface of the density at 0.001 a.u., colored by the ESP volume (e.g. -0.03..0.03 a.u.)"
