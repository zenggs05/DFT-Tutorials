#!/bin/bash
# Collect the ENCUT scan: energy per atom and the difference to the
# largest ENCUT (the reference). Usage: bash collect_encut.sh
printf "%8s %16s %14s\n" "ENCUT" "E/atom (eV)" "dE (meV/atom)"
for d in $(ls -d encut_* 2>/dev/null | sort -t_ -k2 -n); do
    e=$(grep "free  energy   TOTEN" "$d/OUTCAR" 2>/dev/null | tail -n 1 | awk '{print $5}')
    n=$(grep "NIONS" "$d/OUTCAR" 2>/dev/null | head -n 1 | awk '{print $NF}')
    [ -n "$e" ] && [ -n "$n" ] && echo "${d#encut_} $e $n"
done | awk '{enc[NR]=$1; epa[NR]=$2/$3}
    END {ref=epa[NR];
         for (i=1;i<=NR;i++) printf "%8s %16.6f %14.2f\n", enc[i], epa[i], (epa[i]-ref)*1000}'
echo "Converged when |dE| < 1 meV/atom (for energy differences often less is needed)."
