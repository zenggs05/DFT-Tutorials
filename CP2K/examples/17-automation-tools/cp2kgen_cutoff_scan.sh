#!/bin/bash
# CUTOFF scan driven by cp2kgen (cp2k-input-tools). Chapter 17 of the CP2K tutorial.
#
#   pip install "cp2k-input-tools[yaml]"
#   bash cp2kgen_cutoff_scan.sh                 # generate, run and collect
#   DRYRUN=1 bash cp2kgen_cutoff_scan.sh        # only generate the inputs
#   CUTOFFS="300,400,500" bash cp2kgen_cutoff_scan.sh
#   CP2K_CMD="mpirun -np 4 cp2k.psmp" bash cp2kgen_cutoff_scan.sh
#
# Each generated input runs in its own directory scan_cutoff/c<value>/ so that the
# restart and wavefunction files of different runs never collide.

set -euo pipefail

TEMPLATE=${TEMPLATE:-h2o_cutoff_template.inp}
CUTOFFS=${CUTOFFS:-300,350,400,450,500,600,700,800}
CP2K_CMD=${CP2K_CMD:-cp2k.psmp}
DRYRUN=${DRYRUN:-0}
WORK=scan_cutoff

command -v cp2kgen  >/dev/null || { echo "cp2kgen not found: pip install cp2k-input-tools" >&2; exit 1; }
command -v cp2klint >/dev/null || { echo "cp2klint not found: pip install cp2k-input-tools" >&2; exit 1; }

# 1) lint the template (schema bundled with cp2k-input-tools, see chapter 17)
cp2klint "$TEMPLATE"

# 2) generate one input per cutoff; project names are zipped with the cutoffs
projects=$(echo "$CUTOFFS" | tr ',' '\n' | sed 's/^/c/' | paste -sd, -)
mkdir -p "$WORK"
cp "$TEMPLATE" "$WORK/"
(
  cd "$WORK"
  cp2kgen --zip "$TEMPLATE" \
    "force_eval/dft/mgrid/cutoff=[${CUTOFFS}]" \
    "global/project_name=[${projects}]"
)

stem=${TEMPLATE%.inp}
for c in $(echo "$CUTOFFS" | tr ',' ' '); do
  d="$WORK/c$c"
  mkdir -p "$d"
  mv "$WORK/${stem}-cutoff_${c}-project_name_c${c}.inp" "$d/run.inp"
done
echo "generated $(ls -d "$WORK"/c* | wc -l) inputs in $WORK/"
[ "$DRYRUN" = "1" ] && exit 0

# 3) syntax check with the real CP2K, then run (sequentially; use an array job on a cluster)
for c in $(echo "$CUTOFFS" | tr ',' ' '); do
  (
    cd "$WORK/c$c"
    cp2k.psmp -c run.inp > check.log 2>&1 || { echo "c$c: cp2k -c failed, see check.log" >&2; exit 1; }
    rm -f run.out
    $CP2K_CMD -i run.inp -o run.out
  )
done

# 4) collect: total energy per cutoff and difference to the largest cutoff
printf "%8s %22s %14s\n" "CUTOFF" "E_total[Ha]" "dE[meV]" > "$WORK/summary.dat"
last=$(echo "$CUTOFFS" | tr ',' '\n' | tail -1)
eref=$(grep "ENERGY| Total FORCE_EVAL" "$WORK/c$last/run.out" | tail -1 | awk '{print $NF}')
for c in $(echo "$CUTOFFS" | tr ',' ' '); do
  e=$(grep "ENERGY| Total FORCE_EVAL" "$WORK/c$c/run.out" | tail -1 | awk '{print $NF}')
  awk -v c="$c" -v e="$e" -v r="$eref" 'BEGIN{printf "%8s %22.10f %14.3f\n", c, e, (e-r)*27211.386}' >> "$WORK/summary.dat"
done
cat "$WORK/summary.dat"
