#!/bin/bash
# chain_opt_sp_dependency.sh - submit opt+freq, then a high-level single point that starts
# only if the first job succeeded (Slurm --dependency=afterok).
# g16_slurm.sh exits with 0 only on "Normal termination", which is what afterok tests.
#
# Usage:  bash chain_opt_sp_dependency.sh
# PBS equivalent:  jid=$(qsub -v INPUT=a.gjf g16_pbs.sh); qsub -W depend=afterok:$jid -v INPUT=b.gjf g16_pbs.sh
set -euo pipefail

jid1=$(sbatch --parsable --job-name=bz_optfreq g16_slurm.sh benzene_opt_freq.gjf)
echo "opt+freq          : job $jid1"

jid2=$(sbatch --parsable --job-name=bz_sp_tz --dependency=afterok:"$jid1" \
       --kill-on-invalid-dep=yes g16_slurm.sh benzene_sp_tz.gjf)
echo "single point (dep): job $jid2 (waits for $jid1)"

squeue -u "$USER" -o "%.10i %.12j %.8T %.20E"
