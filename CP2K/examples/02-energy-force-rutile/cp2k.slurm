#!/bin/bash
# cp2k.slurm - DEPRECATED, kept only so that old links keep working.
#
# Despite its name, the original file was an LSF (#BSUB) script, and it ran
# "input.inp" although the input in this directory is input_file.inp.
# It has been replaced by two corrected scripts:
#
#     job_slurm.sh   Slurm version   ->  sbatch job_slurm.sh
#     job_lsf.sh     LSF version     ->  bsub < job_lsf.sh
#
# See README.md in this directory.
echo "cp2k.slurm is deprecated: use job_slurm.sh (Slurm) or job_lsf.sh (LSF)." >&2
exit 1
