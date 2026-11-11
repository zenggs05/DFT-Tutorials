#!/bin/bash
# Run the three single points (AB, A, B) and build the charge density difference.
# Usage:  bash run_all.sh                 (local, uses $CP2K, default: mpirun -np 4 cp2k.psmp)
#         CP2K="srun cp2k.psmp" bash run_all.sh   (inside a Slurm allocation)
set -euo pipefail

CP2K=${CP2K:-"mpirun -np 4 cp2k.psmp"}

for name in full slab co; do
    if grep -q "PROGRAM ENDED AT" "${name}.out" 2>/dev/null; then
        echo "${name}: already finished, skipped"
        continue
    fi
    echo "running ${name}.inp ..."
    $CP2K -i "${name}.inp" -o "${name}.out"
    grep -q "PROGRAM ENDED AT" "${name}.out" || { echo "${name} failed" >&2; exit 1; }
done

# cube names follow <PROJECT>-ELECTRON_DENSITY-1_0.cube; check with ls if they differ
AB=$(ls full-ELECTRON_DENSITY-*.cube | head -n 1)
A=$(ls slab-ELECTRON_DENSITY-*.cube | head -n 1)
B=$(ls co-ELECTRON_DENSITY-*.cube | head -n 1)
python chgdiff.py "$AB" "$A" "$B"

# Fermi energies (hartree) of the three runs, for reference
grep -H "Fermi energy:" full.out slab.out co.out | awk '{print $1, $NF}'
