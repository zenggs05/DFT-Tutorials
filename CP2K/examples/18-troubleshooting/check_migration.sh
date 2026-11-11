#!/bin/bash
# Heuristic scan of CP2K input files for constructs that changed in CP2K 2026.1/2026.2.
# Chapter 18 of the CP2K tutorial.
#
#   bash check_migration.sh my.inp [more.inp ...]
#
# It only reads the text (no preprocessing): @INCLUDE files are not followed and
# @IF branches are not evaluated - findings inside an @IF ... @ENDIF block are
# reported as INFO only. It does not replace "cp2k.psmp -c", which remains the
# authoritative check for the CP2K version you actually run.
#
# Exit code: 0 = nothing found, 1 = at least one ERROR (input will not parse on 2026.2).

set -u
[ $# -ge 1 ] || { echo "usage: $0 input.inp [...]" >&2; exit 2; }

status=0
for f in "$@"; do
  [ -r "$f" ] || { echo "$f: cannot read" >&2; status=1; continue; }
  awk -v FILE="$f" '
  function top(k) { return (depth - k >= 1) ? stack[depth - k] : "" }
  function report(level, msg) {
    if (in_if) level = "INFO(@IF)"
    printf "%s:%d: [%s] %s\n", FILE, NR, level, msg
    if (level == "ERROR") nerr++; else if (level == "WARN") nwarn++
  }
  {
    line = $0
    sub(/[!#].*$/, "", line)                 # strip comments
    gsub(/^[ \t]+|[ \t]+$/, "", line)
    if (toupper(line) ~ /^@IF/)    { in_if = 1; next }
    if (toupper(line) ~ /^@ENDIF/) { in_if = 0; next }
    if (line == "" || line ~ /^@/) next
    n = split(line, w, /[ \t]+/)
    word = toupper(w[1])
    if (word ~ /^&END/) {
      if (depth > 0) {
        # remember whether a closed SMEAR section had a METHOD keyword
        if (stack[depth] == "SMEAR" && !smear_method)
          report("WARN", "&SMEAR without METHOD: default changed from ENERGY_WINDOW (<= 2026.1) to GAUSSIAN (2026.2); set METHOD explicitly")
        depth--
      }
      next
    }
    if (substr(word, 1, 1) == "&") {
      name = substr(word, 2)
      if (name == "PDOS" && top(0) == "PRINT" && top(1) == "DFT")
        report("ERROR", "&DFT/&PRINT/&PDOS is the pre-2026.2 syntax; 2026.2 expects &PRINT/&DOS/&PDOS (PR #5326)")
      if (name == "SMEAR") smear_method = 0
      if (name == "OT" && has_kp)
        report("WARN", "&OT together with &KPOINTS aborts on CP2K <= 2026.2 (OT not possible with kpoint calculations)")
      if (name == "KPOINTS") has_kp = 1
      stack[++depth] = name
      next
    }
    # keywords
    if (word == "TYPE" && top(0) == "CELL_OPT" && top(1) == "MOTION")
      report("ERROR", "&MOTION/&CELL_OPT/TYPE was removed in 2026.2 (PR #5257); delete this line (DIRECT_CELL_OPT is now the only mode)")
    if ((word == "EVAL_ENERGY_FORCES" || word == "EVAL_FORCES") && top(0) == "REFTRAJ")
      report("ERROR", word " was removed in 2026.2 (PR #5401); use EVAL ENERGY or EVAL ENERGY_FORCES (keyword exists since 2024.1)")
    if (word == "DO_EWALD" && (top(0) == "XTB" || top(0) == "DFTB"))
      report("WARN", "DO_EWALD in &" top(0) " is deprecated in 2026.2 (periodicity taken from &CELL PERIODIC); still required for periodic runs on <= 2026.1")
    if (word == "METHOD" && top(0) == "SMEAR") smear_method = 1
  }
  END {
    printf "%s: %d error(s), %d warning(s)\n", FILE, nerr, nwarn
    exit (nerr > 0) ? 1 : 0
  }' "$f" || status=1
done
exit $status
