#!/bin/bash
# make_potcar.sh - build POTCAR from the element line of a POSCAR
#
# Usage:
#   export VASP_PP_PATH=/path/to/potentials     # contains potpaw_PBE/<name>/POTCAR
#   bash make_potcar.sh [POSCAR] [Element=Variant ...]
#
# Examples:
#   bash make_potcar.sh                         # POSCAR in the current directory
#   bash make_potcar.sh POSCAR Ti=Ti_pv Ga=Ga   # override the default choice
#   PP_SET=potpaw_PBE.54 bash make_potcar.sh    # other potential set directory
#
# Default choices = VASP "recommended PAW potentials" for standard DFT,
# PBE set potpaw.64, valid since October 2026. Calculations with many empty
# states (GW, BSE, optics) should use the _GW variants instead (override them).
# The script never ships or prints POTCAR content; it only concatenates
# files from YOUR licensed potential library and reports header fields.

POSCAR_FILE=POSCAR
if [ $# -gt 0 ] && [[ "$1" != *=* ]]; then
    POSCAR_FILE=$1
    shift
fi
PP_SET=${PP_SET:-potpaw_PBE}
OUT=${OUT:-POTCAR}

if [ -z "$VASP_PP_PATH" ]; then
    echo "ERROR: set VASP_PP_PATH to the directory that contains $PP_SET/" >&2
    exit 1
fi
if [ ! -s "$POSCAR_FILE" ]; then
    echo "ERROR: $POSCAR_FILE not found" >&2
    exit 1
fi

declare -A REC=(
 [H]=H [He]=He [Li]=Li_sv [Be]=Be_sv [B]=B [C]=C [N]=N [O]=O [F]=F [Ne]=Ne
 [Na]=Na_pv [Mg]=Mg_pv [Al]=Al [Si]=Si [P]=P [S]=S [Cl]=Cl [Ar]=Ar
 [K]=K_sv [Ca]=Ca_sv [Sc]=Sc_sv [Ti]=Ti_sv [V]=V_sv [Cr]=Cr_pv [Mn]=Mn_pv
 [Fe]=Fe_pv [Co]=Co_pv [Ni]=Ni_pv [Cu]=Cu_pv [Zn]=Zn [Ga]=Ga_d [Ge]=Ge_d
 [As]=As [Se]=Se [Br]=Br [Kr]=Kr [Rb]=Rb_sv [Sr]=Sr_sv [Y]=Y_sv [Zr]=Zr_sv
 [Nb]=Nb_sv [Mo]=Mo_sv [Tc]=Tc_pv [Ru]=Ru_pv [Rh]=Rh_pv [Pd]=Pd [Ag]=Ag
 [Cd]=Cd [In]=In_d [Sn]=Sn_d [Sb]=Sb [Te]=Te [I]=I [Xe]=Xe [Cs]=Cs_sv
 [Ba]=Ba_sv [La]=La [Ce]=Ce [Pr]=Pr_3 [Nd]=Nd_3 [Pm]=Pm_3 [Sm]=Sm_3
 [Eu]=Eu_3 [Gd]=Gd_3 [Tb]=Tb_3 [Dy]=Dy_3 [Ho]=Ho_3 [Er]=Er_3 [Tm]=Tm_3
 [Yb]=Yb_3 [Lu]=Lu_3 [Hf]=Hf_pv [Ta]=Ta_pv [W]=W_sv [Re]=Re_pv [Os]=Os_pv
 [Ir]=Ir [Pt]=Pt [Au]=Au [Hg]=Hg [Tl]=Tl_d [Pb]=Pb_d [Bi]=Bi_d [Po]=Po_d
 [At]=At [Rn]=Rn [Fr]=Fr_sv [Ra]=Ra_sv [Ac]=Ac [Th]=Th [Pa]=Pa [U]=U
 [Np]=Np [Pu]=Pu [Am]=Am [Cm]=Cm
)

# user overrides Element=Variant
for arg in "$@"; do
    el=${arg%%=*}; var=${arg#*=}
    if [ -z "$el" ] || [ -z "$var" ] || [ "$el" = "$arg" ]; then
        echo "ERROR: override must look like Element=Variant, got '$arg'" >&2
        exit 1
    fi
    REC[$el]=$var
done

# line 6 of a VASP5+ POSCAR holds the element symbols
read -r -a ELEMENTS <<< "$(sed -n '6p' "$POSCAR_FILE" | tr -d '\r')"
if [ ${#ELEMENTS[@]} -eq 0 ] || [[ "${ELEMENTS[0]}" =~ ^[0-9] ]]; then
    echo "ERROR: line 6 of $POSCAR_FILE has no element symbols (old VASP4 format?)." >&2
    echo "       Insert a line with the symbols above the atom counts." >&2
    exit 1
fi

files=()
for raw in "${ELEMENTS[@]}"; do
    el=${raw%%[_/]*}                 # Si_d or Si/... -> Si (pymatgen style labels)
    var=${REC[$el]}
    if [ -z "$var" ]; then
        echo "ERROR: no default potential for element '$el'; pass $el=Variant" >&2
        exit 1
    fi
    f="$VASP_PP_PATH/$PP_SET/$var/POTCAR"
    if [ ! -s "$f" ]; then
        echo "ERROR: $f not found (check VASP_PP_PATH, PP_SET or the variant name)" >&2
        exit 1
    fi
    files+=("$f")
done

cat "${files[@]}" > "$OUT"

echo "Wrote $OUT from $PP_SET for: ${ELEMENTS[*]}"
printf '%-6s %-12s %-34s %8s %8s\n' "elem" "variant" "TITEL" "ZVAL" "ENMAX"
i=0
maxen=0
for f in "${files[@]}"; do
    titel=$(grep -m1 "TITEL" "$f" | sed 's/.*= *//')
    zval=$(grep -m1 "ZVAL" "$f" | sed 's/.*ZVAL *= *//' | awk '{print $1}')
    enmax=$(grep -m1 "ENMAX" "$f" | sed 's/.*ENMAX *= *//' | tr -d ';' | awk '{print $1}')
    el=${ELEMENTS[$i]%%[_/]*}
    printf '%-6s %-12s %-34s %8s %8s\n' "$el" "${REC[$el]}" "$titel" "$zval" "$enmax"
    maxen=$(awk -v a="$maxen" -v b="$enmax" 'BEGIN{print (b+0 > a+0) ? b : a}')
    i=$((i + 1))
done
awk -v m="$maxen" 'BEGIN{printf "max(ENMAX) = %.1f eV  ->  1.3 x max(ENMAX) = %.0f eV\n", m, 1.3*m}'
echo "Check: element order of POTCAR must equal line 6 of $POSCAR_FILE"
grep "TITEL" "$OUT"
