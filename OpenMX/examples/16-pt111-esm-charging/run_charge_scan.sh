#!/bin/bash
# Scan the excess charge of a Pt(111) slab under ESM on3 (vacuum | slab | ideal metal).
# For every charge q (scf.system.charge, + = holes) an input is generated from
# Pt111_ESM_on3.dat, OpenMX is run, and q, sigma, mu and Utot are collected in charge_scan.csv.
#
# Environment variables (all optional):
#   CHARGES   list of charges            default "-0.02 -0.01 0.00 0.01 0.02"
#   OPENMX    openmx executable          default openmx
#   NP, NT    MPI processes / threads    default 8 / 1
#   AREA      surface-cell area in Ang^2 default 6.6540 (1x1 Pt(111), a=2.7719 Ang)
#   DRYRUN=1  only write the inputs
#   COLLECT_ONLY=1  skip running, only collect existing outputs
set -u
CHARGES=${CHARGES:-"-0.02 -0.01 0.00 0.01 0.02"}
OPENMX=${OPENMX:-openmx}
NP=${NP:-8}
NT=${NT:-1}
AREA=${AREA:-6.6540}
TEMPLATE=Pt111_ESM_on3.dat
export OMP_NUM_THREADS=${NT}
ulimit -s unlimited 2>/dev/null || true

tag() {   # -0.010 -> m0.010 ; 0.010 -> p0.010
    printf "%s%05.3f" "$(awk -v q="$1" 'BEGIN{print (q<0)?"m":"p"}')" "$(awk -v q="$1" 'BEGIN{print (q<0)?-q:q}')"
}

if [ "${COLLECT_ONLY:-0}" != "1" ]; then
    for q in ${CHARGES}; do
        name=pt111_on3_q$(tag "$q")
        inp=${name}.dat
        sed -e "s#^System.Name .*#System.Name                   ${name}#" \
            -e "s#^scf.system.charge .*#scf.system.charge             ${q}#" \
            "${TEMPLATE}" > "${inp}"
        echo "written ${inp}"
        [ "${DRYRUN:-0}" = "1" ] && continue
        mpirun -np "${NP}" "${OPENMX}" "${inp}" -nt "${NT}" > "${name}.std" 2>&1
        grep -q "The calculation was normally finished." "${name}.std" || echo "WARNING: ${name} did not finish normally"
    done
    [ "${DRYRUN:-0}" = "1" ] && exit 0
fi

# sigma (uC/cm^2) = q[e] * 1602.18 / A[Ang^2] ; positive q = positively charged slab
echo "q_e,sigma_uC_cm2,mu_Hartree,mu_eV,Utot_Hartree" > charge_scan.csv
for q in ${CHARGES}; do
    name=pt111_on3_q$(tag "$q")
    out=${name}.out
    [ -f "${out}" ] || { echo "missing ${out}"; continue; }
    mu=$(awk '/^ *Chemical potential \(Hartree\)/{v=$4} END{print v}' "${out}")
    ut=$(awk '/^ *Utot\./{v=$2} END{print v}' "${out}")
    awk -v q="$q" -v A="${AREA}" -v mu="$mu" -v ut="$ut" \
        'BEGIN{printf "%s,%.4f,%s,%.5f,%s\n", q, q*1602.18/A, mu, mu*27.211386, ut}' >> charge_scan.csv
done
column -s, -t charge_scan.csv 2>/dev/null || cat charge_scan.csv
