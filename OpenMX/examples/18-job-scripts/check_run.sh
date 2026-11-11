#!/bin/bash
# Quick health check of finished (or crashed) OpenMX runs.
#   bash check_run.sh                 # all *.std in the current directory
#   bash check_run.sh si64.std c60.std
# For every run it reports: normal termination, known error messages, SCF steps and the
# last NormRD / |dUele| (from System.Name.out), the requested and used cutoff, Utot,
# and the total wall time.  Exit status 1 if any run looks bad.
set -u
shopt -s nullglob
files=("$@")
[ ${#files[@]} -eq 0 ] && files=(*.std)
[ ${#files[@]} -eq 0 ] && { echo "no *.std files found"; exit 0; }

PATTERNS='Could not find|Found an undefined species|ERROR key=|errors in the inputfile|may be wrong|<ESM:Warning>|Warning: The b- or c-axis|conflicts between leads|cannot find the number of threads|not supported|input_open: can not open|You have to give a unit cell|Segmentation fault|signal 11|BAD TERMINATION|Killed|out of memory|oom-kill|error while loading shared libraries|Invalid communicator|was incorrect on entry'

bad=0
for std in "${files[@]}"; do
    name=${std%.std}
    out=${name}.out
    echo "================ ${std}"
    if grep -q "The calculation was normally finished." "${std}"; then
        echo "status      : normally finished"
    else
        echo "status      : NOT finished normally"
        bad=1
    fi
    if grep -E -n "${PATTERNS}" "${std}" | head -n 5 | sed 's|^|error/warn  : |' | grep .; then
        bad=1
    fi
    if [ -f "${out}" ]; then
        # last SCF block of the .out file: "SCF=  n  NormRD=  x  Uele=  y"
        awk '
          /SCF history at MD=/ {n=0}
          /^ *SCF= *[0-9]+ +NormRD=/ {n++; nrd=$4; u0=u1; u1=$6}
          END {
            if (n>0) {
              du = (n>1) ? u1-u0 : 0; if (du<0) du=-du
              printf "SCF         : %d steps in last block, last NormRD = %s, |dUele| = %.3e Hartree\n", n, nrd, du
            }
          }' "${out}"
        grep -m1 "scf.maxIter" "${out}" | awk '{printf "scf.maxIter : %s (in the echoed input)\n", $2}'
        grep -m1 "Required cutoff energy" "${out}" | sed 's/^ */cutoff      : /'
        grep -m1 "Num. of grids" "${out}" | sed 's/^ */grids       : /'
        awk '/^ *Utot\./{v=$2} END{if (v!="") print "Utot        : " v " Hartree"}' "${out}"
        awk '/Total Computational Time/{v=$NF} END{if (v!="") print "wall time   : " v " s (max over MPI ranks)"}' "${out}"
    else
        echo "no ${out} (run stopped before writing it)"
    fi
    grep -m1 "number of threads in each node" "${std}" | sed 's/^/threads     : /'
done
exit ${bad}
