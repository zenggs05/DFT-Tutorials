#!/bin/bash
# goodvibes_run.sh - quasi-harmonic free energies with high-level single points and the
# 1 mol/L standard state, for the four stationary points of the Diels-Alder workflow.
#   G(1 M) = E(wB97X-D/def2-TZVP, SMD) + qRRHO thermal correction (B3LYP-D3(BJ)/def2-SVP, SMD)
#            + RT ln(24.46)
# GoodVibes >= 4.4 (pip install goodvibes). Older 3.x releases used -t/-c instead of
# --temp/--conc; check "goodvibes -h" for your version.
#
# Usage:  bash goodvibes_run.sh
set -uo pipefail
command -v goodvibes >/dev/null 2>&1 || { echo "ERROR: goodvibes not found (pip install goodvibes)"; exit 1; }

FREQ_LOGS=(reactant_opt.log reactant2_opt.log ts_opt.log product_opt.log)

# number of imaginary frequencies printed by the (last) Freq step; 0 if the line is absent
nimag() { grep "imaginary frequencies" "$1" | tail -n 1 | awk '{print $2}' | grep -E '^[0-9]+$' || echo 0; }

ok=1
for f in "${FREQ_LOGS[@]}"; do
    for x in "$f" "${f%.log}_TZ.log"; do
        if [ ! -f "$x" ] || ! tail -n 1 "$x" | grep -q "Normal termination"; then
            echo "ERROR: $x missing or not normally terminated"; ok=0
        fi
    done
    [ -f "$f" ] || continue
    n=$(nimag "$f")
    want=0; [ "$f" = ts_opt.log ] && want=1
    if [ "$n" -ne "$want" ]; then
        echo "ERROR: $f has $n imaginary frequencies (expected $want)"; ok=0
    else
        echo "ok: $f has $n imaginary frequencies"
    fi
done
[ "$ok" -eq 1 ] || exit 1

# per-structure table: E_SPC, ZPE, H_SPC, T.S, T.qh-S, G(T)_SPC, qh-G(T)_SPC (Hartree)
goodvibes "${FREQ_LOGS[@]}" --spc TZ --temp 298.15 --conc 1.0 -f 100 --qs grimme | tee goodvibes_table.txt

# reaction profile in kcal/mol from the YAML document (GoodVibes >= 4.4 reaction-profile format)
goodvibes "${FREQ_LOGS[@]}" --spc TZ --temp 298.15 --conc 1.0 -f 100 --qs grimme \
          --pes da_profile.yaml | tee goodvibes_profile.txt

echo "Plot: python plot_profile.py R=0 TS=<dG_TS> P=<dG_P> -o da_profile.png  (values from goodvibes_profile.txt)"
