#!/bin/bash
# =====================================================================
# Submit relax -> static -> band as three dependent Slurm jobs.
# Run it on the LOGIN node (it only calls sbatch):
#   bash chain_relax_static_band.sh
# afterok:<id>          start only if job <id> ended with exit code 0
# --kill-on-invalid-dep cancel the job when its dependency can never be
#                       satisfied (otherwise it waits forever with the
#                       reason "DependencyNeverSatisfied")
#
# PBS equivalent : j1=$(qsub ...); qsub -W depend=afterok:$j1 ...
# LSF equivalent : bsub -w "done(<jobid>)" < job.sh
# =====================================================================
set -e
cd "$(dirname "$0")"

for f in POSCAR KPOINTS KPOINTS.band INCAR.relax INCAR.static INCAR.band POTCAR; do
    [ -f "$f" ] || { echo "missing $f" >&2; exit 1; }
done

# --parsable prints "jobid" or "jobid;cluster"; keep the id only
j1=$(sbatch --parsable --job-name=si-relax  --export=ALL,STEP=relax  job_step.sh)
j1=${j1%%;*}
j2=$(sbatch --parsable --job-name=si-static --export=ALL,STEP=static \
            --dependency=afterok:${j1} --kill-on-invalid-dep=yes job_step.sh)
j2=${j2%%;*}
j3=$(sbatch --parsable --job-name=si-band   --export=ALL,STEP=band \
            --dependency=afterok:${j2} --kill-on-invalid-dep=yes job_step.sh)
j3=${j3%%;*}

echo "submitted: relax=${j1}  static=${j2} (after ${j1})  band=${j3} (after ${j2})"
echo "watch with: squeue -u \$USER -o '%.10i %.12j %.2t %.10M %R'"
