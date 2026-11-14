#!/bin/bash
# =====================================================================
# Summary of your recent Slurm jobs (run on the login node).
#   bash sacct_report.sh            # jobs of the last 2 days
#   bash sacct_report.sh 7          # jobs of the last 7 days
#   bash sacct_report.sh 7 12345    # only job 12345 (all its steps)
# Columns: id, name, state, exit code, elapsed time, allocated CPUs,
# peak memory of a step (MaxRSS), node list.
# Efficiency of one finished job: "seff <jobid>" (if your site installed
# the contrib tool seff).
# =====================================================================
DAYS=${1:-2}
JOB=${2:-}
START=$(date -d "-${DAYS} days" +%Y-%m-%d 2>/dev/null || date -v-"${DAYS}"d +%Y-%m-%d)
FMT="JobID%-16,JobName%-14,State%-12,ExitCode%-8,Elapsed%-11,AllocCPUS%-6,MaxRSS%-10,NodeList%-20"

if [ -n "${JOB}" ]; then
    sacct -j "${JOB}" --format="${FMT}"
else
    sacct -u "${USER}" -S "${START}" --format="${FMT}" -X
fi

echo
echo "--- counts by state since ${START} (main jobs only) ---"
sacct -u "${USER}" -S "${START}" -X -n -P --format=State | \
    awk '{split($1,a," "); n[a[1]]++} END {for (s in n) printf "%-14s %d\n", s, n[s]}'

echo
echo "--- failed / timed out / out of memory ---"
sacct -u "${USER}" -S "${START}" -E now -X -n -P --format=JobID,JobName,State,ExitCode \
    --state=FAILED,TIMEOUT,OUT_OF_MEMORY,NODE_FAIL | column -t -s '|'
