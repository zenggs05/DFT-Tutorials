#!/bin/bash
# common.sh - helper functions sourced by run_encut.sh, run_kpts.sh, run_sigma.sh
#
# Functions
#   natoms POSCAR             -> total number of atoms (VASP5 POSCAR)
#   outcar_row DIR            -> "E0 TOTEN EENTRO PRESSURE MAXF NKPTS TIME" from DIR/OUTCAR
#                                (NaN for missing values; DONE=1 if the run finished)
#   run_or_submit DIR         -> run VASP in DIR (MODE=local) or sbatch a job (MODE=slurm)
#
# Environment variables used by run_or_submit
#   MODE          local | slurm | collect     (default local)
#   VASP_CMD      command for MODE=local      (default "mpirun -np ${NP:-4} vasp_std")
#   JOB_TEMPLATE  Slurm script for MODE=slurm (default ../01-setup-compile/job_cpu.slurm)

HERE_COMMON=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MODE=${MODE:-local}
VASP_CMD=${VASP_CMD:-"mpirun -np ${NP:-4} vasp_std"}
JOB_TEMPLATE=${JOB_TEMPLATE:-$HERE_COMMON/../01-setup-compile/job_cpu.slurm}

natoms() {
    sed -n '7p' "$1" | tr -d '\r' | awk '{n = 0; for (i = 1; i <= NF; i++) n += $i; print n}'
}

outcar_row() {
    local d=$1 n
    if [ ! -s "$d/OUTCAR" ]; then
        echo "NaN NaN NaN NaN NaN NaN NaN"
        return
    fi
    n=$(natoms "$d/POSCAR")
    awk -v N="$n" '
        /energy\(sigma->0\)/   { s = $0; sub(/.*energy\(sigma->0\) *= */, "", s); split(s, a, " "); e0 = a[1] }
        /free +energy +TOTEN/  { s = $0; sub(/.*TOTEN *= */, "", s); split(s, a, " "); f = a[1] }
        /EENTRO/               { s = $0; sub(/.*EENTRO *= */, "", s); split(s, a, " "); ts = a[1] }
        /external pressure/    { s = $0; sub(/.*external pressure *= */, "", s); split(s, a, " "); p = a[1] }
        /NKPTS *=/             { for (i = 1; i < NF; i++) if ($i == "NKPTS") { nk = $(i + 2); sub(/^=/, "", nk) } }
        /Elapsed time/         { t = $NF }
        /TOTAL-FORCE \(eV\/Angst\)/ {
            getline                                   # dashed line
            mx = 0
            for (i = 1; i <= N; i++) {
                getline
                ff = sqrt($4 * $4 + $5 * $5 + $6 * $6)
                if (ff > mx) mx = ff
            }
            fmax = mx
        }
        END {
            if (e0 == "") e0 = "NaN"; if (f == "") f = "NaN"; if (ts == "") ts = "NaN"
            if (p == "") p = "NaN"; if (fmax == "") fmax = "NaN"; if (nk == "") nk = "NaN"
            if (t == "") t = "NaN"
            print e0, f, ts, p, fmax, nk, t
        }' "$d/OUTCAR"
}

run_or_submit() {
    local d=$1
    case "$MODE" in
        local)
            echo "running $d"
            (cd "$d" && $VASP_CMD > vasp.out 2>&1)
            ;;
        slurm)
            if [ ! -s "$JOB_TEMPLATE" ]; then
                echo "ERROR: JOB_TEMPLATE=$JOB_TEMPLATE not found" >&2
                return 1
            fi
            (cd "$d" && cp "$JOB_TEMPLATE" job.slurm && sbatch job.slurm)
            ;;
        collect)
            : # only collect results
            ;;
        *)
            echo "ERROR: unknown MODE=$MODE (local | slurm | collect)" >&2
            return 1
            ;;
    esac
}

check_potcar() {
    if [ ! -s POTCAR ]; then
        echo "ERROR: POTCAR missing in $(pwd). Build it first, e.g." >&2
        echo "       cat \$VASP_PP_PATH/potpaw_PBE/Si/POTCAR > POTCAR" >&2
        return 1
    fi
}
