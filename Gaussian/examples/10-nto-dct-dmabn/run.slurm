#!/bin/bash
#SBATCH --job-name=dmabn-nto
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out

# Runs the chapter-10 DMABN NTO/DCT input in dependency order.
# %nprocshared=8 and %mem=8GB in the inputs match --cpus-per-task=8 and ~80% of --mem.
# Adapt the module/profile lines and the scratch path to your cluster.

# module load gaussian/16.C.02
export g16root=${g16root:-/opt/gaussian}
source "$g16root/g16/bsd/g16.profile"

export GAUSS_SCRDIR=${SCRATCH_BASE:-/tmp}/$USER/g16_$SLURM_JOB_ID
mkdir -p "$GAUSS_SCRDIR"
trap 'rm -rf "$GAUSS_SCRDIR"' EXIT

cd "$SLURM_SUBMIT_DIR" || exit 1

INPUTS="nto_analysis_link1.gjf"

for inp in $INPUTS; do
    log=${inp%.gjf}.log
    echo "$(date '+%F %T')  running $inp"
    g16 < "$inp" > "$log"
    if ! tail -n 1 "$log" | grep -q "Normal termination"; then
        echo "ERROR: $inp did not end with Normal termination; later steps need its chk. Stop." >&2
        exit 1
    fi
done

for chk in *.chk; do
    formchk "$chk" "${chk%.chk}.fchk" > /dev/null
done
echo "$(date '+%F %T')  all done"
