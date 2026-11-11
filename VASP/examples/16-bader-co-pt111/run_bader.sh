#!/bin/bash
# Bader analysis after a VASP run with LAECHG = .TRUE.
# Requires: chgsum.pl (VTST scripts) and the 'bader' executable in PATH, python3.
# Usage: bash run_bader.sh [group]      e.g.  bash run_bader.sh 17-18
set -e

GROUP=${1:-17-18}     # atoms 17-18 = C and O of the CO molecule in this POSCAR

for f in CHGCAR AECCAR0 AECCAR2 POSCAR POTCAR; do
    [ -s "$f" ] || { echo "missing $f - did VASP finish with LAECHG = .TRUE. ?" >&2; exit 1; }
done

grep -q "General timing" OUTCAR || echo "WARNING: OUTCAR has no timing block, the run may not have finished" >&2
grep "dimension x,y,z NGXF" OUTCAR | tail -1

# 1) reference density = core + all-electron valence  -> CHGCAR_sum
if command -v chgsum.pl >/dev/null 2>&1; then
    chgsum.pl AECCAR0 AECCAR2
else
    echo "chgsum.pl not found; you can also sum AECCAR0 + AECCAR2 with ../../scripts/chgdiff.py (see -h)" >&2
    exit 1
fi

# 2) partition with the reference, integrate the valence density of CHGCAR
bader CHGCAR -ref CHGCAR_sum > bader.log 2>&1
tail -4 ACF.dat

# 3) net charges = ZVAL (from POTCAR) - Bader population
python3 bader_net_charge.py --acf ACF.dat --poscar POSCAR --potcar POTCAR --group "$GROUP" | tee bader_net_charge.txt

# 4) sanity check against NELECT
nelect=$(grep NELECT OUTCAR | head -1 | awk '{print $3}')
echo "NELECT from OUTCAR: $nelect"
