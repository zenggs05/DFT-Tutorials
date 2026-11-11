#!/bin/bash
# check_vtst_patch.sh - read-only sanity check of a VTST-patched VASP source tree
#
# Usage:   bash check_vtst_patch.sh /path/to/vasp.6.6.1
#
# It checks the four edits described in vtst_patch_steps.md
#   1. VTST Fortran files copied into src/ (chain.F replaced, pyamff_fortran/)
#   2. src/.objects: VTST objects listed in the required order
#   3. src/main.F : TSIF passed to CHAIN_FORCE, chain_init called unconditionally
#   4. src/makefile: pyamff_fortran in LIB (and libs in dependencies for -j builds)
# Nothing is modified. Exit status = number of failed checks.

ROOT=${1:-.}
SRC="$ROOT/src"
nfail=0

pass() { printf '  [PASS] %s\n' "$1"; }
fail() { printf '  [FAIL] %s\n' "$1"; nfail=$((nfail + 1)); }
warn() { printf '  [WARN] %s\n' "$1"; }

if [ ! -d "$SRC" ] || [ ! -f "$SRC/main.F" ]; then
    echo "ERROR: $SRC/main.F not found. Give the VASP root directory as argument."
    exit 99
fi

echo "== 1. VTST source files in $SRC"
for f in chain.F neb.F dimer.F lanczos.F dynmat.F instanton.F bbm.F \
         sd.F cg.F qm.F lbfgs.F bfgs.F fire.F opt.F ml_pyamff.F; do
    if [ -f "$SRC/$f" ]; then pass "$f present"; else fail "$f missing"; fi
done
if [ -d "$SRC/pyamff_fortran" ]; then pass "pyamff_fortran/ present"; else fail "pyamff_fortran/ missing"; fi
if grep -qi "^ *use *dimer" "$SRC/chain.F"; then
    pass "chain.F is the VTST version (uses module dimer)"
else
    fail "chain.F does not look like the VTST version"
fi

echo "== 2. Object order in $SRC/.objects"
OBJ="$SRC/.objects"
if [ ! -f "$OBJ" ]; then
    fail ".objects not found"
else
    # all object tokens in file order, one per line
    tokens=$(grep -o '[A-Za-z0-9_/*.]*\.o' "$OBJ")
    pos() { printf '%s\n' "$tokens" | grep -n -x -F "$1" | head -1 | cut -d: -f1; }
    p_chain=$(pos chain.o)
    p_opt=$(pos opt.o)
    if [ -z "$p_chain" ]; then fail "chain.o not listed"; fi
    for o in dynmat.o neb.o dimer.o lanczos.o instanton.o; do
        p=$(pos "$o")
        if [ -z "$p" ]; then fail "$o not listed"
        elif [ -n "$p_chain" ] && [ "$p" -lt "$p_chain" ]; then pass "$o before chain.o"
        else fail "$o must appear before chain.o"; fi
    done
    if [ -z "$p_opt" ]; then
        fail "opt.o not listed"
    else
        for o in sd.o cg.o qm.o lbfgs.o bfgs.o fire.o ml_pyamff.o "pyamff_fortran/*.o"; do
            p=$(pos "$o")
            if [ -z "$p" ]; then fail "$o not listed"
            elif [ "$p" -lt "$p_opt" ]; then pass "$o before opt.o"
            else fail "$o must appear before opt.o"; fi
        done
        if [ -n "$p_chain" ] && [ "$p_opt" -lt "$p_chain" ]; then
            pass "opt.o before chain.o"
        else
            warn "opt.o is not before chain.o (VTST lists it in the block inserted before chain.o)"
        fi
    fi
fi

echo "== 3. Edits in $SRC/main.F"
call=$(grep -i -A1 "CALL *CHAIN_FORCE" "$SRC/main.F")
if [ -z "$call" ]; then
    fail "CALL CHAIN_FORCE not found in main.F"
elif printf '%s' "$call" | grep -qi "TSIF"; then
    pass "TSIF is passed to CHAIN_FORCE"
else
    fail "TSIF missing in the CHAIN_FORCE call"
fi
if grep -qi "IF *( *LCHAIN *) *CALL *chain_init" "$SRC/main.F"; then
    fail "chain_init is still guarded by IF (LCHAIN)"
elif grep -qi "CALL *chain_init" "$SRC/main.F"; then
    pass "chain_init is called unconditionally"
else
    fail "CALL chain_init not found"
fi

echo "== 4. $SRC/makefile"
if grep -E "^ *LIB *=" "$SRC/makefile" | grep -q "pyamff_fortran"; then
    pass "pyamff_fortran added to LIB"
else
    fail "pyamff_fortran not in LIB"
fi
if grep -E "^ *dependencies *:" "$SRC/makefile" | grep -q "libs"; then
    pass "dependencies rule includes libs (safe for make -j)"
else
    warn "dependencies rule has no 'libs' (needed only for parallel make)"
fi

echo
if [ "$nfail" -eq 0 ]; then
    echo "All checks passed. Now run: make veryclean && make DEPS=1 -j 16 all"
else
    echo "$nfail check(s) failed. See vtst_patch_steps.md."
fi
exit "$nfail"
