#!/bin/bash
# crest_conformers.sh - CREST/GFN2-xTB conformer search for the diene in implicit toluene (ALPB),
# then one DFT opt+freq input per conformer (conf_NN_opt.gjf) at the level used in this workflow.
# Needs crest (>= 3.0) and xtb in PATH. CREST 3.x has no default runtype: --imtdgc selects iMTD-GC.
#
# Usage:  NP=8 bash crest_conformers.sh [input.xyz]
set -euo pipefail
XYZ=${1:-butadiene_strans.xyz}
NP=${NP:-8}
ROUTE='#p B3LYP EmpiricalDispersion=GD3BJ def2SVP SCRF=(SMD,Solvent=Toluene) Opt Freq'

command -v crest >/dev/null 2>&1 || { echo "ERROR: crest not in PATH"; exit 1; }
crest "$XYZ" --gfn2 --alpb toluene -T "$NP" --imtdgc > crest.out 2>&1
[ -f crest_conformers.xyz ] || { echo "ERROR: crest_conformers.xyz missing, see crest.out"; exit 1; }

# crest_conformers.xyz: energy-sorted multi-frame xyz (line 2 of each frame holds the xTB energy)
awk -v route="$ROUTE" '
NR == 1 { nat = $1 }
{
    i = (NR - 1) % (nat + 2)
    f = int((NR - 1) / (nat + 2)) + 1
    if (i == 0) {
        name = sprintf("conf_%02d_opt", f)
        file = name ".gjf"
        print "%chk=" name ".chk" > file
        print route > file
        print "" > file
    } else if (i == 1) {
        print "CREST conformer " f ", GFN2-xTB energy " $1 " Eh" > file
        print "" > file
        print "0 1" > file
    } else {
        printf "%-2s %12.6f %12.6f %12.6f\n", $1, $2, $3, $4 > file
        if (i == nat + 1) { print "" > file; close(file) }
    }
}' crest_conformers.xyz

n=$(ls conf_*_opt.gjf | wc -l)
echo "CREST found $n conformer(s) within its energy window; inputs: conf_*_opt.gjf"
grep -h "CREST conformer" conf_*_opt.gjf
