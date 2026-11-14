#!/bin/bash
# Quick text-only diagnosis of one or more CP2K output files (chapter 18, section 3.1).
#
#   bash quick_diagnose.sh run.out [more.out ...]
#
# It only greps the output (CP2K is not called) and prints:
#   normal termination, the [ABORT] message, version, MPI/OpenMP layout (and MPI-mismatch
#   symptoms), SCF convergence counts, warning summary, electron count, peak memory,
#   HFX integrals computed on the fly, and the head of the timing report.
# Hints refer to the sections of chapter 18. Output wording may differ slightly between
# CP2K versions; when a line is not found the script says so instead of guessing.

set -u
[ $# -ge 1 ] || { echo "usage: $0 cp2k_output.out [...]" >&2; exit 2; }

hint() { printf '  -> %s\n' "$*"; }

for out in "$@"; do
  echo "=================================================================="
  echo "File: $out"
  [ -r "$out" ] || { echo "  cannot read $out"; continue; }

  # ---- 1. termination -------------------------------------------------
  ver=$(grep -m1 "CP2K| version string" "$out" | sed 's/.*version string: *//')
  echo "CP2K version      : ${ver:-not found}"
  nend=$(grep -c "PROGRAM ENDED AT" "$out")
  if [ "$nend" -ge 1 ]; then
    echo "Termination       : normal (PROGRAM ENDED AT found)"
  else
    echo "Termination       : ABNORMAL (no PROGRAM ENDED AT)"
    hint "read the [ABORT] box below and the scheduler log (stderr), section 3.1"
  fi
  if grep -q "\[ABORT\]" "$out"; then
    echo "ABORT message     :"
    # The abort box is framed by two lines of asterisks. For long messages the text
    # starts ABOVE the "[ABORT]" picture line, so buffer the box from its opening line.
    box=$(awk '
      /^ \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*+ *$/ {
        if (p) { print; exit }
        nb = 0; box[++nb] = $0; next
      }
      nb > 0 && !p { box[++nb] = $0 }
      /\[ABORT\]/ && !p { p = 1; for (i = 1; i <= nb; i++) print box[i]; next }
      p { print }' "$out")
    printf '%s\n' "$box" | head -20 | sed 's/^/    /'
    # message text = columns 12-79 of the box body (left of it is the picture)
    msg=$(printf '%s\n' "$box" | awk '!/^ \*\*\*\*\*\*\*\*\*\*/ {print substr($0, 12, 67)}' \
          | tr '\n' ' ' | tr -s ' ')
    case "$msg" in
      *"not available anymore"*)      hint "keyword removed in this CP2K version: section 4.1 / 9" ;;
      *"unknown keyword"*|*"unknown subsection"*) hint "parser error: check spelling/version, section 4.1 / 9" ;;
      *"basis set"*"not found"*|*"potential"*"not found"*) hint "data file / name mismatch: section 4.2" ;;
      *"cannot be opened"*)           hint "data file missing: set CP2K_DATA_DIR, section 4.2" ;;
      *"EMAX_SPLINE"*)                hint "overlapping atoms / wrong cell: section 4.3" ;;
      *"cell volume"*|*"det(h)"*)     hint "invalid &CELL: section 4.3" ;;
      *"LSD"*)                        hint "odd electrons or wrong multiplicity: section 4.4" ;;
      *"OT not possible"*)            hint "OT + k-points: use diagonalization, section 4.4" ;;
      *"NOT converged"*)              hint "SCF not converged: section 4.4 and chapter 5" ;;
    esac
  fi

  # ---- 2. parallel layout ---------------------------------------------
  nmpi=$(grep -c "DBCSR| MPI: Number of processes" "$out")
  mpiv=$(grep -m1 "DBCSR| MPI: Number of processes" "$out" | awk '{print $NF}')
  ompv=$(grep -m1 "DBCSR| OMP: Current number of threads" "$out" | awk '{print $NF}')
  nstart=$(grep -c "PROGRAM STARTED AT" "$out")
  echo "MPI ranks         : ${mpiv:-not found}   (line found $nmpi time(s))"
  echo "OpenMP threads    : ${ompv:-not found}"
  if [ "$nmpi" -gt 1 ] || [ "$nstart" -gt 1 ]; then
    echo "  WARNING: parallel info / PROGRAM STARTED AT printed more than once."
    hint "possible MPI mismatch (several serial copies writing one file): section 6"
    hint "ignore this for multi-replica runs only if you are reading a replica log"
  fi

  # ---- 3. SCF ----------------------------------------------------------
  nconv=$(grep -c "SCF run converged" "$out")
  nnot=$(grep -c "SCF run NOT converged" "$out")
  echo "SCF converged     : $nconv time(s);  NOT converged: $nnot time(s)"
  [ "$nnot" -gt 0 ] && hint "check structure/charge/multiplicity first, then SCF settings (section 4.4)"
  grep -m1 -i "does not exist" "$out" >/dev/null && \
    hint "a restart file was requested but not found - started from ATOMIC guess (section 8.4)"
  grep -m1 -i "smearing includes the last MO" "$out" >/dev/null && \
    hint "ADDED_MOS too small for the smearing (section 4.4)"

  # ---- 4. electrons, warnings, memory ----------------------------------
  nel=$(grep -m1 "Number of electrons" "$out" | awk '{print $NF}')
  echo "Number of electrons: ${nel:-not found}   (compare with sum of valence electrons - CHARGE)"
  nw=$(grep -m1 "The number of warnings for this run is" "$out" | awk '{print $NF}')
  echo "Warnings (summary): ${nw:-not found}"
  if grep -q "\*\*\* WARNING" "$out"; then
    echo "Warning sources   :"
    grep "\*\*\* WARNING in" "$out" | sed 's/.*WARNING in \([^ ]*\).*/\1/' | sort | uniq -c | sort -rn | head -8 | sed 's/^/    /'
  fi
  mem=$(grep "MEMORY| Estimated peak process memory" "$out" | tail -1 | awk '{print $NF}')
  echo "Peak memory/rank  : ${mem:-not found} MiB"
  otf=$(grep -m1 "calculated on the fly" "$out" | awk '{print $NF}')
  if [ -n "$otf" ]; then
    echo "HFX ERIs on the fly: $otf"
    [ "$otf" != "0" ] && hint "MAX_MEMORY too small for in-core HFX (section 7.5)"
  fi

  # ---- 5. timing ---------------------------------------------------------
  if grep -q "T I M I N G" "$out"; then
    echo "Timing report (head):"
    grep -A12 "T I M I N G" "$out" | sed -n '5,13p' | sed 's/^/    /'
  fi
done
