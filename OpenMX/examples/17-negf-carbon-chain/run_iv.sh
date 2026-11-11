#!/bin/bash
# I-V curve of the carbon chain: one NEGF SCF (Step 2 + 3) per bias voltage.
# Requires lead-chain.hks from Step 1 (openmx Lead-Chain.dat) in the current directory.
#
# Environment variables (optional):
#   BIASES     list of bias voltages (V)       default "0.0 0.2 0.4 0.6 0.8 1.0"
#   OPENMX     openmx executable               default openmx
#   NP, NT     MPI processes / OpenMP threads  default 8 / 1
#   RESTART=1  start each bias from the converged density of the previous bias via
#              scf.restart on + scf.restart.filename (keep NP identical for all biases:
#              the density restart files are split per MPI rank). Verify on your build.
#   DRYRUN=1   only write the inputs ;  COLLECT_ONLY=1  only collect results
set -u
BIASES=${BIASES:-"0.0 0.2 0.4 0.6 0.8 1.0"}
OPENMX=${OPENMX:-openmx}
NP=${NP:-8}
NT=${NT:-1}
TEMPLATE=NEGF-Chain.dat
export OMP_NUM_THREADS=${NT}
ulimit -s unlimited 2>/dev/null || true

[ -f lead-chain.hks ] || [ "${DRYRUN:-0}" = "1" ] || [ "${COLLECT_ONLY:-0}" = "1" ] || {
    echo "lead-chain.hks not found: run Step 1 first (openmx Lead-Chain.dat)"; exit 1; }

prev=""
if [ "${COLLECT_ONLY:-0}" != "1" ]; then
    for V in ${BIASES}; do
        name=negf-chain-V${V}
        inp=${name}.dat
        sed -e "s#^System.Name .*#System.Name                   ${name}#" \
            -e "s#^NEGF.bias.voltage .*#NEGF.bias.voltage             ${V}#" \
            -e "s#^NEGF.tran.CurrentDensity .*#NEGF.tran.CurrentDensity      off#" \
            -e "s#^NEGF.tran.Channel .*#NEGF.tran.Channel             off#" \
            "${TEMPLATE}" > "${inp}"
        if [ "${RESTART:-0}" = "1" ] && [ -n "${prev}" ]; then
            printf 'scf.restart                   on\nscf.restart.filename          %s\n' "${prev}" >> "${inp}"
        fi
        echo "written ${inp}"
        [ "${DRYRUN:-0}" = "1" ] && continue
        mpirun -np "${NP}" "${OPENMX}" "${inp}" -nt "${NT}" > "${name}.std" 2>&1
        if grep -q "The calculation was normally finished." "${name}.std"; then
            prev=${name}
        else
            echo "WARNING: ${name} did not finish normally (see ${name}.std)"
        fi
    done
    [ "${DRYRUN:-0}" = "1" ] && exit 0
fi

# collect: '# average current for up and down spins (ampere)= I_up I_down'
echo "V,I_up_A,I_down_A,I_total_uA,SCF_steps,last_NormRD" > iv.csv
for V in ${BIASES}; do
    name=negf-chain-V${V}
    [ -f "${name}.current" ] || { echo "missing ${name}.current"; continue; }
    read -r iu id < <(awk -F'=' '/average current for up and down spins/{print $2}' "${name}.current")
    nscf=$(grep -c "SCF=" "${name}.out" 2>/dev/null || echo 0)
    nrd=$(awk '/SCF= *[0-9]+ +NormRD=/{v=$4} END{print v}' "${name}.out" 2>/dev/null)
    awk -v V="$V" -v a="$iu" -v b="$id" -v n="$nscf" -v r="$nrd" \
        'BEGIN{printf "%s,%s,%s,%.4f,%s,%s\n", V, a, b, (a+b)*1e6, n, r}' >> iv.csv
done
column -s, -t iv.csv 2>/dev/null || cat iv.csv
