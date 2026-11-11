#!/bin/bash
# make_sp_inputs.sh - create <root>_TZ.gjf high-level single points from the template
# sp_highlevel.gjf (written for ts_opt). GoodVibes pairs X.log with X_TZ.log via --spc TZ.
# Usage:  bash make_sp_inputs.sh [root ...]      default roots: the four stationary points
set -euo pipefail
roots=("$@")
[ ${#roots[@]} -gt 0 ] || roots=(reactant_opt reactant2_opt ts_opt product_opt)
for r in "${roots[@]}"; do
    [ -f "$r.chk" ] || echo "note: $r.chk does not exist yet (run $r.gjf first)"
    sed "s/ts_opt/$r/g" sp_highlevel.gjf > "${r}_TZ.gjf"
    echo "wrote ${r}_TZ.gjf"
done
