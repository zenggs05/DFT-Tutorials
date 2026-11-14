#!/usr/bin/env bash
# =============================================================================
# phonopy_workflow.sh -- VASP + phonopy phonon workflow (finite displacements or DFPT).
#                        phonopy 有限位移 / DFPT 自动化：生成位移超胞、批量建目录与提交、
#                        检查与收集 FORCE_SETS / FORCE_CONSTANTS / BORN、写 band.conf
#                        （BAND = AUTO，由 seekpath 生成路径）、作图与虚频检测。
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
#
# Steps (one sub-command each)
#   prepare  phonopy-init -d (phonopy >= 4) or phonopy -d (2.x/3.x): SPOSCAR + POSCAR-NNN;
#            fd  : one directory <prefix>-NNN per displacement (POSCAR, INCAR, KPOINTS, POTCAR, job.sh)
#            dfpt: dfpt/ with SPOSCAR as POSCAR (INCAR must contain IBRION = 7 or 8)
#            optional born/ (unit cell, LEPSILON or LCALCEPS) for the NAC correction
#   submit   submit every directory that has no complete vasprun.xml (sbatch / qsub)
#   check    vasprun.xml complete, exactly one ionic step, SCF below NELM, hessian (DFPT),
#            Born charges (born/)
#   collect  FORCE_SETS (fd) or FORCE_CONSTANTS (dfpt); BORN via phonopy-vasp-born
#   conf     write band.conf, mesh.conf (DOS/PDOS) and tprop.conf (thermal properties)
#   post     phonopy band / DOS / thermal properties with plots (-p -s), band.dat, imaginary check
#   imag     report imaginary (negative) frequencies in band.yaml / mesh.yaml
#   all      prepare + submit
# Settings given to 'prepare' are stored in phonopy_workflow.cfg and reused by the later steps.
#
# Usage examples
#   bash phonopy_workflow.sh prepare -c POSCAR --dim "2 2 2" -j job_vasp_slurm.sh
#   bash phonopy_workflow.sh submit
#   bash phonopy_workflow.sh check
#   bash phonopy_workflow.sh collect
#   bash phonopy_workflow.sh post --mesh "24 24 24" --pdos "1, 2"
#   bash phonopy_workflow.sh prepare -m dfpt -c POSCAR-unitcell -i INCAR_dfpt -k KPOINTS_dfpt \
#        --born-incar INCAR_born --born-kpoints KPOINTS_born --pa "0 1/2 1/2 1/2 0 1/2 1/2 1/2 0"
#   bash phonopy_workflow.sh post --band "0 0 0 1/2 0 1/2 1/2 1/2 1/2" --labels '$\Gamma$ X L'
#
# Dependencies: bash >= 4, awk, phonopy (2.x/3.x or >= 4.0; >= 4 uses phonopy-init),
#               seekpath (python module, for BAND = AUTO), sbatch or qsub for 'submit'.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CFG="phonopy_workflow.cfg"

usage() {
    cat <<'EOF'
Usage: phonopy_workflow.sh COMMAND [options]

Commands:
  prepare   create displaced supercells and the VASP directories
  submit    submit all directories without a complete vasprun.xml
  check     check the finished VASP runs
  collect   FORCE_SETS / FORCE_CONSTANTS (+ BORN)
  conf      write band.conf, mesh.conf, tprop.conf
  post      phonopy band structure, DOS/PDOS, thermal properties, imaginary-mode check
  imag      imaginary-mode report from band.yaml / mesh.yaml
  all       prepare + submit

Options (prepare / all):
  -m, --method fd|dfpt     finite displacements (default) or VASP DFPT (IBRION = 7/8)
  -c, --cell FILE          relaxed unit cell                 (default POSCAR)
  -d, --dim "N N N"        supercell matrix, 3 or 9 integers (default "2 2 2")
      --pa AXES            PRIMITIVE_AXES: auto (default), P, F, I, ... or 9 numbers
  -i, --incar FILE         INCAR for the supercell runs      (default INCAR)
  -k, --kpoints FILE       KPOINTS for the supercell runs    (default KPOINTS; optional with KSPACING)
  -P, --potcar FILE        POTCAR                            (default POTCAR)
  -j, --job FILE           job script copied as job.sh       (default job_vasp_slurm.sh or
                           job_vasp_pbs.sh next to this script, by the detected scheduler)
      --amplitude A        displacement distance in Angstrom (phonopy default 0.01)
      --pm                 plus and minus displacements
      --prefix NAME        directory prefix for fd           (default disp -> disp-001 ...)
      --born-incar FILE    also create born/ (unit cell, LEPSILON/LCALCEPS) with this INCAR
      --born-kpoints FILE  KPOINTS for born/ (default: KPOINTS_born if present, else -k file)
  -f, --force              overwrite existing directories / resubmit finished runs
Options (submit):
  -s, --submit-cmd CMD     e.g. "sbatch -p debug" (default: sbatch or qsub, auto-detected)
Options (conf / post / imag):
      --band "PATH"        explicit BAND string (default: BAND = AUTO, path from seekpath)
      --labels "L1 L2 .."  BAND_LABELS for --band
      --band-points N      points per segment                (default 101)
      --mesh "N N N"       q-point mesh for DOS / thermal    (default "24 24 24")
      --pdos "1 2, 3"      PDOS groups (atom indices of the primitive cell)
      --tmax T, --tstep T  temperature range of tprop.conf   (default 1000, 10 K)
      --imag-tol F         |f| threshold in THz for imaginary modes (default 0.1)
Other:
      --phonopy-old        force phonopy 2.x/3.x syntax (default: detected via phonopy-init)
  -n, --dry-run            print the commands only
  -h, --help               this help

Exit status: 0 = OK; 1 = error, unfinished runs ('check') or imaginary modes ('imag').

Notes:
  * The supercell INCAR must be a single point: IBRION = -1, NSW = 0 (fd) or IBRION = 7/8 (dfpt),
    EDIFF ~ 1E-8, PREC = Accurate, LREAL = .FALSE.; use the ENCUT/POTCAR of the relaxation.
  * KPOINTS must be scaled to the supercell (same k-point density as in the unit cell).
  * With BAND = AUTO the path refers to the seekpath standard primitive cell.
EOF
}

die()  { echo "ERROR: $*" >&2; exit 1; }
info() { echo "[phonopy_workflow] $*"; }
warn() { echo "WARNING: $*" >&2; }

# ------------------------------- defaults -----------------------------------
CMD="${1:-}"
case "$CMD" in
    -h|--help|"") usage; exit 0 ;;
    prepare|submit|check|collect|conf|post|imag|all) shift ;;
    *) usage >&2; die "unknown command '$CMD'" ;;
esac

METHOD="fd"; CELL="POSCAR"; DIM="2 2 2"; PA="auto"; INCAR="INCAR"; KPOINTS="KPOINTS"; POTCAR="POTCAR"
JOB=""; AMPLITUDE=""; PM=0; PREFIX="disp"; BORN_INCAR=""; BORN_KPOINTS=""
# read the settings stored by 'prepare' (command-line options below override them)
if [[ -f "$CFG" ]]; then
    # shellcheck disable=SC1090
    source "$CFG"
fi
FORCE=0; DRY=0; SUBMIT_CMD=""; BAND="AUTO"; LABELS=""; BAND_POINTS=101; MESH="24 24 24"; PDOS=""
TMAX=1000; TSTEP=10; IMAG_TOL=0.1; OLD=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--method)      METHOD="${2:-}"; shift 2 ;;
        -c|--cell)        CELL="${2:-}"; shift 2 ;;
        -d|--dim)         DIM="${2:-}"; shift 2 ;;
        --pa)             PA="${2:-}"; shift 2 ;;
        -i|--incar)       INCAR="${2:-}"; shift 2 ;;
        -k|--kpoints)     KPOINTS="${2:-}"; shift 2 ;;
        -P|--potcar)      POTCAR="${2:-}"; shift 2 ;;
        -j|--job)         JOB="${2:-}"; shift 2 ;;
        --amplitude)      AMPLITUDE="${2:-}"; shift 2 ;;
        --pm)             PM=1; shift ;;
        --prefix)         PREFIX="${2:-}"; shift 2 ;;
        --born-incar)     BORN_INCAR="${2:-}"; shift 2 ;;
        --born-kpoints)   BORN_KPOINTS="${2:-}"; shift 2 ;;
        -s|--submit-cmd)  SUBMIT_CMD="${2:-}"; shift 2 ;;
        --band)           BAND="${2:-}"; shift 2 ;;
        --labels)         LABELS="${2:-}"; shift 2 ;;
        --band-points)    BAND_POINTS="${2:-}"; shift 2 ;;
        --mesh)           MESH="${2:-}"; shift 2 ;;
        --pdos)           PDOS="${2:-}"; shift 2 ;;
        --tmax)           TMAX="${2:-}"; shift 2 ;;
        --tstep)          TSTEP="${2:-}"; shift 2 ;;
        --imag-tol)       IMAG_TOL="${2:-}"; shift 2 ;;
        --phonopy-old)    OLD=1; shift ;;
        -f|--force)       FORCE=1; shift ;;
        -n|--dry-run)     DRY=1; shift ;;
        -h|--help)        usage; exit 0 ;;
        *) die "unknown option '$1' (see -h)" ;;
    esac
done

[[ "$METHOD" == "fd" || "$METHOD" == "dfpt" ]] || die "--method must be fd or dfpt"
ndim=$(wc -w <<< "$DIM")
[[ "$ndim" -eq 3 || "$ndim" -eq 9 ]] || die "--dim needs 3 or 9 integers, got '$DIM'"
[[ "$DIM" =~ ^[-0-9\ ]+$ ]] || die "--dim must contain integers only"
[[ "$PREFIX" =~ ^[A-Za-z0-9_.]+$ ]] || die "--prefix may only contain letters, digits, '_' and '.'"

# phonopy >= 4.0 moved -d, -f, --fc to 'phonopy-init'
if [[ -z "$OLD" ]]; then
    if command -v phonopy-init >/dev/null 2>&1; then OLD=0; else OLD=1; fi
fi
if [[ "$OLD" == "0" ]]; then PHINIT="phonopy-init"; else PHINIT="phonopy"; fi
PYTHON="${PYTHON:-python3}"

run() {   # run a command, or only print it with --dry-run
    if [[ $DRY -eq 1 ]]; then echo "+ $*"; else "$@"; fi
}

# value of an INCAR tag (last occurrence, comments removed)
incar_get() {
    [[ -f "$2" ]] || return 0
    tr -d '\r' < "$2" | sed 's/[#!].*//' | tr ';' '\n' \
        | awk -v t="$1" '{ n = index($0, "="); if (n == 0) next
                           k = substr($0, 1, n - 1); gsub(/[ \t]/, "", k)
                           if (toupper(k) == toupper(t)) { v = substr($0, n + 1); gsub(/^[ \t]+|[ \t]+$/, "", v); last = v } }
                         END { if (last != "") print last }'
}

# run directories in phonopy order (POSCAR-001, -002, ... -> <prefix>-001, ...)
run_dirs() {
    if [[ "$METHOD" == "dfpt" ]]; then
        if [[ -d dfpt ]]; then echo "dfpt"; fi
    else
        local d
        for d in "$PREFIX"-[0-9]*; do if [[ -d "$d" ]]; then echo "$d"; fi; done | sort -t- -k2,2n
    fi
    return 0
}

vasprun_complete() { [[ -f "$1/vasprun.xml" ]] && tail -n 5 "$1/vasprun.xml" | grep -q "</modeling>"; }

save_cfg() {
    [[ $DRY -eq 1 ]] && return 0
    {
        echo "# written by phonopy_workflow.sh prepare on $(date)"
        local v
        for v in METHOD CELL DIM PA INCAR KPOINTS POTCAR JOB AMPLITUDE PM PREFIX BORN_INCAR BORN_KPOINTS; do
            printf '%s=%q\n' "$v" "${!v}"
        done
    } > "$CFG"
    info "settings stored in $CFG"
}

default_job() {
    if [[ -n "$JOB" ]]; then echo "$JOB"; return; fi
    if command -v sbatch >/dev/null 2>&1; then echo "$SCRIPT_DIR/job_vasp_slurm.sh"
    elif command -v qsub >/dev/null 2>&1; then echo "$SCRIPT_DIR/job_vasp_pbs.sh"
    else echo "$SCRIPT_DIR/job_vasp_slurm.sh"; fi
}

# ------------------------------------------------------------------ prepare
do_prepare() {
    command -v "$PHINIT" >/dev/null 2>&1 || die "$PHINIT not found (conda install -c conda-forge phonopy)"
    local f
    for f in "$CELL" "$INCAR" "$POTCAR"; do [[ -s "$f" ]] || die "file '$f' not found or empty"; done
    if [[ ! -s "$KPOINTS" && -z "$(incar_get KSPACING "$INCAR")" ]]; then
        die "no KPOINTS file '$KPOINTS' and no KSPACING in $INCAR"
    fi
    local ib nsw ncore npar
    ib=$(incar_get IBRION "$INCAR"); nsw=$(incar_get NSW "$INCAR")
    ncore=$(incar_get NCORE "$INCAR"); npar=$(incar_get NPAR "$INCAR")
    if [[ "$METHOD" == "fd" ]]; then
        [[ "${ib:--1}" == "-1" ]] || warn "IBRION = $ib in $INCAR: displaced supercells must not relax (use IBRION = -1)"
        [[ "${nsw:-0}" == "0" ]] || warn "NSW = $nsw in $INCAR: use NSW = 0 for the displaced supercells"
    else
        [[ "$ib" == "7" || "$ib" == "8" ]] || die "DFPT needs IBRION = 7 or 8 in $INCAR (found '${ib:-unset}')"
        [[ -z "$ncore" || "$ncore" == "1" ]] || warn "NCORE = $ncore: DFPT linear response requires NCORE = 1 (remove it)"
        [[ -z "$npar" ]] || warn "NPAR is set in $INCAR: remove it for DFPT"
        [[ -z "$(incar_get IVDW "$INCAR")" ]] || warn "IVDW is set: DFPT does not support vdW corrections - use finite differences"
    fi
    [[ -z "$(incar_get LREAL "$INCAR")" || "$(incar_get LREAL "$INCAR")" =~ ^\.?[Ff] ]] \
        || warn "LREAL is not .FALSE. in $INCAR: real-space projection adds force noise"
    [[ -z "$(incar_get MAGMOM "$INCAR")" ]] \
        || warn "MAGMOM in $INCAR refers to the unit cell: adapt it to the supercell (phonopy --magmom writes a MAGMOM file)"

    if [[ $FORCE -eq 0 ]]; then
        if [[ -n "$(run_dirs)" ]]; then die "run directories exist already (use -f to overwrite)"; fi
    fi
    local args=(-d --dim $DIM -c "$CELL" --pa "$PA")
    [[ -n "$AMPLITUDE" ]] && args+=(--amplitude "$AMPLITUDE")
    [[ $PM -eq 1 ]] && args+=(--pm)
    [[ $DRY -eq 1 ]] || rm -f POSCAR-[0-9]*      # stale displacements of an earlier run
    info "$PHINIT ${args[*]}"
    run "$PHINIT" "${args[@]}"

    local job
    job=$(default_job)
    [[ -f "$job" ]] || { warn "job script '$job' not found - no job.sh written (use -j)"; job=""; }
    mkd() {   # mkd DIR POSCAR-FILE INCAR KPOINTS
        run mkdir -p "$1"
        run cp "$2" "$1/POSCAR"
        run cp "$3" "$1/INCAR"
        [[ -s "$4" ]] && run cp "$4" "$1/KPOINTS"
        run cp "$POTCAR" "$1/POTCAR"
        [[ -n "$job" ]] && run cp "$job" "$1/job.sh"
        return 0
    }
    if [[ "$METHOD" == "fd" ]]; then
        local n=0 p
        if [[ $DRY -eq 1 ]]; then info "(dry run) directories $PREFIX-NNN would be created from POSCAR-NNN"; fi
        for p in $(ls POSCAR-[0-9]* 2>/dev/null | sort -t- -k2,2n); do
            mkd "$PREFIX-${p#POSCAR-}" "$p" "$INCAR" "$KPOINTS"
            n=$((n + 1))
        done
        [[ $DRY -eq 1 || $n -gt 0 ]] || die "phonopy did not write any POSCAR-NNN"
        info "created $n displacement directories ($PREFIX-*); supercell = SPOSCAR"
    else
        [[ $DRY -eq 1 || -s SPOSCAR ]] || die "SPOSCAR was not written"
        mkd dfpt SPOSCAR "$INCAR" "$KPOINTS"
        info "created dfpt/ (supercell SPOSCAR, IBRION = $ib)"
    fi
    if [[ -n "$BORN_INCAR" ]]; then
        [[ -s "$BORN_INCAR" ]] || die "born INCAR '$BORN_INCAR' not found"
        local bk="$BORN_KPOINTS"
        [[ -z "$bk" && -s KPOINTS_born ]] && bk="KPOINTS_born"
        [[ -z "$bk" ]] && bk="$KPOINTS"
        is_on() { [[ "$(incar_get "$1" "$BORN_INCAR")" =~ ^\.?[Tt] ]]; }
        is_on LEPSILON || is_on LCALCEPS || warn "$BORN_INCAR sets neither LEPSILON nor LCALCEPS"
        [[ -z "$(incar_get NCORE "$BORN_INCAR")$(incar_get NPAR "$BORN_INCAR")" ]] \
            || warn "remove NCORE/NPAR from $BORN_INCAR (LEPSILON/LCALCEPS)"
        mkd born "$CELL" "$BORN_INCAR" "$bk"
        info "created born/ (unit cell, $BORN_INCAR, k-points from $bk - use about twice the supercell density)"
    fi
    save_cfg
}

# ------------------------------------------------------------------ submit
do_submit() {
    if [[ -z "$SUBMIT_CMD" ]]; then
        if command -v sbatch >/dev/null 2>&1; then SUBMIT_CMD="sbatch"
        elif command -v qsub >/dev/null 2>&1; then SUBMIT_CMD="qsub"
        else die "neither sbatch nor qsub found (use -s CMD)"; fi
    fi
    local dirs=() d n=0
    mapfile -t dirs < <(run_dirs)
    [[ -d born ]] && dirs+=(born)
    [[ ${#dirs[@]} -gt 0 ]] || die "no run directories - run 'prepare' first"
    for d in "${dirs[@]}"; do
        if vasprun_complete "$d" && [[ $FORCE -eq 0 ]]; then
            info "$d: finished, skipped (-f to resubmit)"; continue
        fi
        [[ -f "$d/job.sh" ]] || die "$d/job.sh missing (prepare with -j JOB)"
        if [[ $DRY -eq 1 ]]; then
            echo "+ (cd $d && $SUBMIT_CMD job.sh)"
        else
            # shellcheck disable=SC2086
            out=$(cd "$d" && $SUBMIT_CMD job.sh) || die "submission failed in $d"
            info "$d: $out"
        fi
        n=$((n + 1))
    done
    info "$n job(s) submitted"
}

# ------------------------------------------------------------------ check
do_check() {
    local dirs=() d bad=0 nelm nion nscf
    mapfile -t dirs < <(run_dirs)
    [[ ${#dirs[@]} -gt 0 ]] || die "no run directories - run 'prepare' first"
    for d in "${dirs[@]}"; do
        if ! vasprun_complete "$d"; then
            echo "NOT FINISHED : $d (vasprun.xml missing or incomplete)"; bad=1; continue
        fi
        nelm=$(incar_get NELM "$d/INCAR"); nelm=${nelm:-60}
        if [[ -f "$d/OSZICAR" ]]; then
            nion=$(grep -c " F= " "$d/OSZICAR" || true)
            nscf=$(awk '/^[A-Za-z]{2,4}:/{n++} / F= /{last = n; n = 0} END{print last + 0}' "$d/OSZICAR")
            if [[ "$METHOD" == "fd" && "$nion" -ne 1 ]]; then
                echo "PROBLEM      : $d has $nion ionic steps (the supercell relaxed? IBRION = -1, NSW = 0)"; bad=1
            fi
            if [[ "$METHOD" == "fd" && "$nscf" -ge "$nelm" ]]; then
                echo "PROBLEM      : $d SCF reached NELM = $nelm (not converged)"; bad=1
            fi
        fi
        if [[ "$METHOD" == "dfpt" ]] && ! grep -q 'name="hessian"' "$d/vasprun.xml"; then
            echo "PROBLEM      : $d/vasprun.xml has no hessian (IBRION = 7/8?)"; bad=1
        fi
    done
    if [[ -d born ]]; then
        if [[ -f born/OUTCAR ]] && grep -q "BORN EFFECTIVE CHARGES" born/OUTCAR; then
            echo "born/        : Born effective charges found"
        else
            echo "NOT FINISHED : born/ (no 'BORN EFFECTIVE CHARGES' in OUTCAR)"; bad=1
        fi
    fi
    if [[ $bad -eq 0 ]]; then
        info "all ${#dirs[@]} run(s) look fine"
    else
        info "problems found - fix and resubmit (submit skips finished runs)"
    fi
    return $bad
}

# ------------------------------------------------------------------ collect
pa_numeric() { [[ $(wc -w <<< "$PA") -eq 9 ]]; }

do_collect() {
    command -v "$PHINIT" >/dev/null 2>&1 || die "$PHINIT not found"
    do_check || die "not all runs are complete - nothing collected"
    local dirs=() files=() d
    mapfile -t dirs < <(run_dirs)
    if [[ "$METHOD" == "fd" ]]; then
        for d in "${dirs[@]}"; do files+=("$d/vasprun.xml"); done
        run "$PHINIT" -f "${files[@]}"
        [[ $DRY -eq 1 || -s FORCE_SETS ]] && info "FORCE_SETS written (check 'Max drift of force constants' in post)"
    else
        run "$PHINIT" --fc dfpt/vasprun.xml
        [[ $DRY -eq 1 || -s FORCE_CONSTANTS ]] && info "FORCE_CONSTANTS written"
    fi
    if [[ -d born ]]; then
        command -v phonopy-vasp-born >/dev/null 2>&1 || die "phonopy-vasp-born not found"
        local bargs=()
        if pa_numeric; then
            bargs=(--pa "$PA")
        else
            warn "--pa is '$PA': phonopy-vasp-born accepts only a numeric matrix; without it the unit cell"
            warn "is taken as primitive. If born/POSCAR is a centred/conventional cell, rerun with --pa \"9 numbers\"."
        fi
        if [[ $DRY -eq 1 ]]; then
            echo "+ (cd born && phonopy-vasp-born ${bargs[*]}) > BORN"
        else
            (cd born && phonopy-vasp-born "${bargs[@]}") > BORN || die "phonopy-vasp-born failed"
            info "BORN written:"; cat BORN
        fi
    fi
}

# ------------------------------------------------------------------ conf
write_conf() {
    local common
    common="DIM = $DIM
PRIMITIVE_AXES = $(tr '[:lower:]' '[:upper:]' <<< "$PA")"
    [[ "$METHOD" == "dfpt" ]] && common+=$'\n'"FORCE_CONSTANTS = READ"
    if [[ "$BAND" == "AUTO" || "$BAND" == "auto" ]]; then
        if ! "$PYTHON" -c "import seekpath" >/dev/null 2>&1; then
            warn "python module seekpath not found: BAND = AUTO will fail (pip install seekpath, or give --band)"
        fi
    fi
    if [[ $DRY -eq 1 ]]; then info "(dry run) would write band.conf, mesh.conf, tprop.conf"; return 0; fi
    {
        echo "# phonopy band structure - written by phonopy_workflow.sh"
        echo "$common"
        if [[ "$BAND" == "AUTO" || "$BAND" == "auto" ]]; then
            echo "# path and labels generated by seekpath (standard primitive cell of seekpath)"
            echo "BAND = AUTO"
        else
            echo "BAND = $BAND"
            [[ -n "$LABELS" ]] && echo "BAND_LABELS = $LABELS"
        fi
        echo "BAND_POINTS = $BAND_POINTS"
        echo "BAND_CONNECTION = .TRUE."
    } > band.conf
    {
        echo "# phonopy phonon DOS (and PDOS) - written by phonopy_workflow.sh"
        echo "$common"
        echo "MESH = $MESH"
        echo "GAMMA_CENTER = .TRUE."
        if [[ -n "$PDOS" ]]; then echo "PDOS = $PDOS"; else echo "DOS = .TRUE."; fi
    } > mesh.conf
    {
        echo "# phonopy thermal properties (F, S, Cv per primitive cell) - written by phonopy_workflow.sh"
        echo "$common"
        echo "MESH = $MESH"
        echo "GAMMA_CENTER = .TRUE."
        echo "TPROP = .TRUE."
        echo "TMIN = 0"
        echo "TMAX = $TMAX"
        echo "TSTEP = $TSTEP"
    } > tprop.conf
    info "written band.conf, mesh.conf, tprop.conf (edit them before 'post' if needed)"
}

# ------------------------------------------------------------------ imag
do_imag() {
    local f found=0 rc=0
    for f in band.yaml mesh.yaml; do
        [[ -f "$f" ]] || continue
        found=1
        awk -v tol="$IMAG_TOL" -v file="$f" '
            /q-position:/ { q = $0; sub(/.*q-position:[ \t]*/, "", q); gsub(/[\[\]]/, "", q)
                            split(q, c, ","); qn = sqrt(c[1]^2 + c[2]^2 + c[3]^2) }
            /^[ \t]*frequency:/ { f = $2 + 0; n++
                            if (n == 1 || f < fmin) { fmin = f; qmin = q }
                            if (f < -tol) { nim++; if (qn < 0.05) ngam++ ; else { qs[q] = 1 } } }
            END {
                if (n == 0) { printf "%-10s: no frequencies found\n", file; exit 0 }
                printf "%-10s: %d frequencies, lowest %.4f THz at q = (%s)\n", file, n, fmin, qmin
                if (nim == 0) { printf "            no imaginary modes below -%s THz\n", tol; exit 0 }
                printf "            %d imaginary value(s) below -%s THz", nim, tol
                printf " (%d within |q| < 0.05 of Gamma)\n", ngam
                m = 0; for (k in qs) { if (m < 5) printf "            away from Gamma at q = (%s)\n", k; m++ }
                if (nim > ngam) printf "            -> imaginary modes away from Gamma: possible real instability (see chapter 12, section 10)\n"
                else            printf "            -> only near Gamma: usually numerical noise (EDIFF, supercell size, ASR)\n"
                exit 1
            }' "$f" || rc=1
    done
    [[ $found -eq 1 ]] || die "no band.yaml or mesh.yaml - run 'post' first"
    return $rc
}

# ------------------------------------------------------------------ post
ph_run() {   # ph_run CONF
    local extra=()
    if [[ "$OLD" == "1" ]]; then
        [[ -f BORN ]] && extra+=(--nac)          # phonopy 2.x/3.x: NAC must be requested
        run phonopy -c "$CELL" -p -s "${extra[@]}" "$1"
    else
        run phonopy --config "$1" -p -s          # phonopy >= 4: NAC automatic when BORN exists
    fi
}

do_post() {
    command -v phonopy >/dev/null 2>&1 || die "phonopy not found"
    [[ -f phonopy_disp.yaml ]] || die "phonopy_disp.yaml missing (run 'prepare' in this directory)"
    if [[ "$METHOD" == "fd" ]]; then
        [[ -s FORCE_SETS ]] || die "FORCE_SETS missing - run 'collect' first"
    else
        [[ -s FORCE_CONSTANTS ]] || die "FORCE_CONSTANTS missing - run 'collect' first"
    fi
    [[ -f BORN ]] && info "BORN found: non-analytical term correction (LO-TO splitting) is applied"
    [[ -f band.conf && -f mesh.conf && -f tprop.conf ]] || write_conf
    info "band structure";       ph_run band.conf
    info "DOS / PDOS";           ph_run mesh.conf
    info "thermal properties";   ph_run tprop.conf
    if command -v phonopy-bandplot >/dev/null 2>&1 && [[ $DRY -eq 1 || -f band.yaml ]]; then
        if [[ $DRY -eq 1 ]]; then echo "+ phonopy-bandplot --gnuplot band.yaml > band.dat"
        else phonopy-bandplot --gnuplot band.yaml > band.dat && info "band.dat written (gnuplot/Origin format)"; fi
    fi
    if [[ $DRY -eq 0 ]]; then
        info "files: $(ls ./*.pdf band.yaml total_dos.dat projected_dos.dat thermal_properties.yaml 2>/dev/null | tr '\n' ' ')"
        do_imag || warn "imaginary modes found (see above)"
    fi
}

# ------------------------------------------------------------------ main
case "$CMD" in
    prepare) do_prepare ;;
    submit)  do_submit ;;
    check)   do_check ;;
    collect) do_collect ;;
    conf)    write_conf ;;
    post)    do_post ;;
    imag)    do_imag ;;
    all)     do_prepare; do_submit ;;
esac
