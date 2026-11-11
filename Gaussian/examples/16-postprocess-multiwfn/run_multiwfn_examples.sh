#!/bin/bash
# run_multiwfn_examples.sh - run the three Multiwfn command streams of this example non-interactively.
# Usage:  bash run_multiwfn_examples.sh
# Needs: Multiwfn in PATH (noGUI version is fine), phenol.chk / phenol_td.log / water_dimer.chk
# produced by the Gaussian jobs of this folder, and formchk for the chk -> fchk step.
# Multiwfn writes fixed file names, so every result is renamed right after each run.
set -uo pipefail

command -v Multiwfn >/dev/null 2>&1 || { echo "ERROR: Multiwfn not in PATH" >&2; exit 1; }
export OMP_STACKSIZE=${OMP_STACKSIZE:-1000M}
ulimit -s unlimited 2>/dev/null || true

for chk in phenol.chk water_dimer.chk; do
    [ -f "${chk%.chk}.fchk" ] || formchk "$chk" "${chk%.chk}.fchk"
done

# 1) ESP analysis on the rho=0.001 a.u. surface of phenol
Multiwfn phenol.fchk -silent < multiwfn_esp_surface_cmds.txt > phenol_esp_surface.out 2>&1
for f in surfanalysis.txt surfanalysis.pdb vtx.pdb; do
    [ -f "$f" ] && mv -f "$f" "phenol_$f"
done
grep -E "Global surface (minimum|maximum)|Overall surface area|Molecular polarity index" phenol_esp_surface.out

# 2) IGMH of the water dimer (fragments 1-3 and 4-6)
Multiwfn water_dimer.fchk -silent < multiwfn_igmh_dimer_cmds.txt > water_dimer_igmh.out 2>&1
mkdir -p igmh_cubes
for f in *.cub; do
    case "$f" in phenol_*|water_dimer_*) ;; *) mv -f "$f" igmh_cubes/ ;; esac
done
ls igmh_cubes

# 3) UV-Vis curve from the TD-DFT log (exports spectrum_line.txt / spectrum_curve.txt)
Multiwfn phenol_td.log -silent < multiwfn_uvvis_plot_cmds.txt > phenol_uvvis.out 2>&1
for f in spectrum_line.txt spectrum_curve.txt; do
    [ -f "$f" ] && mv -f "$f" "phenol_td_$f"
done
ls -l phenol_td_spectrum_*.txt 2>/dev/null || echo "no spectrum files: check phenol_uvvis.out"

echo "A final 'end-of-file during read' message in the .out files is harmless (end of command stream)."
