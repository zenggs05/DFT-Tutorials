#!/bin/bash
# check_inputs.sh - pre-flight check of the four VASP input files
#
# Usage:  bash check_inputs.sh [directory]      (default: current directory)
#
# Checks (read-only):
#   * INCAR/POSCAR/POTCAR exist, KPOINTS exists or KSPACING is set
#   * TAB characters and Windows (CRLF) line endings
#   * POSCAR element line vs. POTCAR TITEL order, number of coordinates
#   * NELECT estimate (sum of ZVAL x atoms), ENCUT vs. max(ENMAX)
#   * number of MAGMOM values (N*value shorthand expanded)
#   * a few KPOINTS pitfalls (tetrahedron with < 4 k points, even MP meshes)
# Exit status = number of errors.

DIR=${1:-.}
cd "$DIR" || exit 99
nerr=0
err()  { printf '  [ERROR] %s\n' "$1"; nerr=$((nerr + 1)); }
warn() { printf '  [WARN ] %s\n' "$1"; }
ok()   { printf '  [ OK  ] %s\n' "$1"; }

# value of an INCAR tag (last occurrence wins), comments removed, upper case
incar_get() {
    awk -v tag="$1" '
        { line=$0; sub(/[#!].*/, "", line); n=split(line, st, ";")
          for (i=1; i<=n; i++) {
              s=st[i]; p=index(s, "="); if (p == 0) continue
              k=substr(s, 1, p-1); gsub(/[ \t]/, "", k)
              if (toupper(k) == toupper(tag)) { v=substr(s, p+1); gsub(/^[ \t]+|[ \t]+$/, "", v); val=v }
          } }
        END { print toupper(val) }' INCAR
}

echo "== files in $(pwd)"
for f in INCAR POSCAR POTCAR; do
    if [ -s "$f" ]; then ok "$f present"; else err "$f missing or empty"; fi
done
if [ -s KPOINTS ]; then
    ok "KPOINTS present"
elif [ -s INCAR ] && [ -n "$(incar_get KSPACING)" ]; then
    ok "no KPOINTS, KSPACING = $(incar_get KSPACING) in INCAR"
else
    warn "no KPOINTS and no KSPACING: VASP uses KSPACING = 0.5 (often too coarse)"
fi
[ -s INCAR ] && [ -s POSCAR ] || { echo "cannot continue"; exit "$nerr"; }

echo "== formatting"
for f in INCAR POSCAR KPOINTS; do
    [ -f "$f" ] || continue
    if grep -U -q $'\r' "$f"; then err "$f has Windows CRLF line endings (dos2unix $f)"; fi
    if grep -q $'\t' "$f"; then warn "$f contains TAB characters (replace by spaces)"; fi
done

echo "== POSCAR"
read -r -a SYMS   <<< "$(sed -n '6p' POSCAR | tr -d '\r')"
read -r -a COUNTS <<< "$(sed -n '7p' POSCAR | tr -d '\r')"
if [[ "${SYMS[0]}" =~ ^[0-9] ]]; then
    err "line 6 has no element symbols (VASP4 format); add them for safety"
    COUNTS=("${SYMS[@]}"); SYMS=(); next=7
else
    next=8
fi
nions=0
for c in "${COUNTS[@]}"; do nions=$((nions + c)); done
ok "species: ${SYMS[*]:-?}   atoms per species: ${COUNTS[*]}   NIONS = $nions"
modeline=$(sed -n "${next}p" POSCAR)
if [[ "$modeline" =~ ^[[:space:]]*[Ss] ]]; then
    ok "Selective dynamics is on"
    next=$((next + 1)); modeline=$(sed -n "${next}p" POSCAR); seldyn=1
fi
if [[ "$modeline" =~ ^[[:space:]]*[CcKk] ]]; then ok "coordinates: Cartesian"; else ok "coordinates: Direct"; fi
ncoord=$(sed -n "$((next + 1)),$((next + nions))p" POSCAR | awk 'NF >= 3' | wc -l)
if [ "$ncoord" -eq "$nions" ]; then ok "$ncoord coordinate lines"; else err "found $ncoord coordinate lines, expected $nions"; fi
if [ -n "$seldyn" ]; then
    nflag=$(sed -n "$((next + 1)),$((next + nions))p" POSCAR | awk 'NF >= 6' | wc -l)
    [ "$nflag" -eq "$nions" ] && ok "every atom has T/F flags" || err "only $nflag of $nions atoms have T/F flags"
fi

echo "== POTCAR"
if [ -s POTCAR ]; then
    mapfile -t TITELS < <(grep "TITEL" POTCAR | sed 's/.*= *//')
    mapfile -t ZVALS  < <(grep "ZVAL" POTCAR | sed 's/.*ZVAL *= *//' | awk '{print $1}')
    mapfile -t ENMAXS < <(grep "ENMAX" POTCAR | sed 's/.*ENMAX *= *//' | tr -d ';' | awk '{print $1}')
    if [ "${#TITELS[@]}" -ne "${#COUNTS[@]}" ]; then
        err "POTCAR has ${#TITELS[@]} potentials but POSCAR has ${#COUNTS[@]} species"
    fi
    nelect=0; maxen=0
    for i in "${!TITELS[@]}"; do
        pel=$(echo "${TITELS[$i]}" | awk '{print $2}'); pel=${pel%%_*}
        sel=${SYMS[$i]%%[_/]*}
        if [ -n "$sel" ] && [ "$pel" != "$sel" ]; then
            err "species $((i + 1)): POSCAR says '$sel' but POTCAR has '${TITELS[$i]}'"
        else
            ok "species $((i + 1)): ${TITELS[$i]}  ZVAL=${ZVALS[$i]}  ENMAX=${ENMAXS[$i]}"
        fi
        nelect=$(awk -v a="$nelect" -v z="${ZVALS[$i]}" -v n="${COUNTS[$i]:-0}" 'BEGIN{print a + z * n}')
        maxen=$(awk -v a="$maxen" -v b="${ENMAXS[$i]}" 'BEGIN{print (b + 0 > a + 0) ? b : a}')
    done
    ok "NELECT (neutral) = $nelect   max(ENMAX) = $maxen eV"
fi

echo "== INCAR"
encut=$(incar_get ENCUT)
if [ -z "$encut" ]; then
    warn "ENCUT not set: default = max(ENMAX); set it explicitly for comparable energies"
elif [ -n "$maxen" ]; then
    ratio=$(awk -v e="$encut" -v m="$maxen" 'BEGIN{printf "%.2f", e / m}')
    ok "ENCUT = $encut eV = $ratio x max(ENMAX)"
    awk -v r="$ratio" 'BEGIN{exit !(r < 1.0)}' && err "ENCUT is below max(ENMAX)"
fi
magmom=$(incar_get MAGMOM)
if [ -n "$magmom" ]; then
    nmag=$(echo "$magmom" | awk '{n=0; for (i=1; i<=NF; i++) { if ($i ~ /\*/) { split($i, a, "*"); n += a[1] } else n++ } print n}')
    ncl=$(incar_get LNONCOLLINEAR)$(incar_get LSORBIT)
    if [[ "$ncl" == *T* ]]; then need=$((3 * nions)); else need=$nions; fi
    if [ "$nmag" -eq "$need" ]; then ok "MAGMOM has $nmag values"; else err "MAGMOM has $nmag values, expected $need"; fi
    [ "$(incar_get ISPIN)" = "2" ] || [[ "$ncl" == *T* ]] || warn "MAGMOM set but ISPIN is not 2"
fi
ismear=$(incar_get ISMEAR)

echo "== KPOINTS"
if [ -s KPOINTS ]; then
    mode=$(sed -n '3p' KPOINTS | tr -d '\r' | sed 's/^[[:space:]]*//' | cut -c1)
    nk=$(sed -n '2p' KPOINTS | tr -d '\r' | awk '{print $1}')
    if [ "$nk" = "0" ] && [[ "$mode" =~ [GgMm] ]]; then
        read -r n1 n2 n3 _ <<< "$(sed -n '4p' KPOINTS | tr -d '\r')"
        ok "automatic mesh $n1 x $n2 x $n3 ($mode)"
        if [[ "$mode" =~ [Mm] ]] && { [ $((n1 % 2)) -eq 0 ] || [ $((n2 % 2)) -eq 0 ] || [ $((n3 % 2)) -eq 0 ]; }; then
            warn "even Monkhorst-Pack mesh: use a Gamma-centred mesh for fcc/hexagonal/2D cells"
        fi
        if [ "$ismear" = "-5" ] && [ $((n1 * n2 * n3)) -lt 4 ]; then
            err "ISMEAR = -5 needs at least 4 k points"
        fi
        if [ "$ismear" = "-5" ] && [[ "$mode" =~ [Mm] ]]; then
            warn "tetrahedron method: prefer a Gamma-centred mesh"
        fi
    elif [[ "$mode" =~ [Ll] ]]; then
        ok "line mode (band structure): needs ICHARG = 11 or use KPOINTS_OPT"
        [ "$(incar_get ICHARG)" = "11" ] || warn "line-mode KPOINTS without ICHARG = 11"
    else
        ok "explicit list or other mode (line 3 starts with '$mode')"
    fi
fi

echo
if [ "$nerr" -eq 0 ]; then echo "No errors found."; else echo "$nerr error(s) found."; fi
exit "$nerr"
