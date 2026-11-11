#!/bin/bash
#
# Check the result of "openmx -runtest" (or -runtestL).
# Reads runtest.result (default) and reports, for every input file,
# the absolute differences "diff Utot" and "diff Force" with respect to the
# reference data. According to the OpenMX manual (Sec. 6), the installation
# can be regarded as correct if all differences are within the seventh
# decimal place, i.e. smaller than 1e-7.
#
# Works with both layouts seen in practice: one record per line, or the
# number / file / Elapsed time / diff Utot / diff Force on separate lines.
#
# Usage:
#   bash check_runtest.sh [runtest.result] [threshold]
#   bash check_runtest.sh runtestL.result 1e-7
# Exit status: 0 = all tests pass, 1 = failure or nothing parsed.
#
set -u
RESULT=${1:-runtest.result}
THRESH=${2:-1e-7}

if [ ! -f "$RESULT" ]; then
  echo "file not found: $RESULT" >&2
  exit 1
fi

awk -v thr="$THRESH" '
  function absval(x) { return (x < 0) ? -x : x }
  function flush() {
    if (name != "" && u != "" && f != "") {
      n++
      bad = (absval(u + 0) >= thr + 0 || absval(f + 0) >= thr + 0)
      if (bad) nbad++
      printf "%3d  %-40s %10s %20s %20s  %s\n", n, name, t, u, f, (bad ? "FAIL" : "ok")
      name = ""; u = ""; f = ""; t = ""
    }
  }
  BEGIN {
    printf "%3s  %-40s %10s %20s %20s  %s\n", "#", "input", "time(s)", "diff Utot", "diff Force", "status"
  }
  /Total elapsed time/ { total = $NF }
  {
    for (i = 1; i <= NF; i++) {
      if ($i ~ /\.dat$/) { name = $i }
      else if ($i == "time(s)=") { t = $(i + 1) }
      else if ($i ~ /^time\(s\)=./) { t = substr($i, 10) }
      else if ($i == "Utot=") { u = $(i + 1) }
      else if ($i ~ /^Utot=./) { u = substr($i, 6) }
      else if ($i == "Force=") { f = $(i + 1); flush() }
      else if ($i ~ /^Force=./) { f = substr($i, 7); flush() }
    }
  }
  END {
    print ""
    if (n == 0) { print "no test records found in the file"; exit 1 }
    printf "tests parsed: %d   failed (>= %s): %d   total elapsed time: %s s\n", n, thr, nbad + 0, (total == "" ? "n/a" : total)
    if (nbad > 0) { print "RESULT: FAIL - see the rows marked FAIL"; exit 1 }
    print "RESULT: PASS"
  }
' "$RESULT"
