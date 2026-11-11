#!/usr/bin/env bash
# =============================================================================
# conv_test.sh -- ENCUT / k-point / KSPACING / SIGMA (or any INCAR tag)
#                 convergence-test driver for VASP.
#                 VASP 收敛测试驱动：批量建目录、改参数、写/提交作业，配合 conv_analyze.py 分析。
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
#
# What it does
#   1. creates one directory per test value (e.g. conv_encut/encut_400/),
#   2. copies INCAR/POSCAR/KPOINTS/POTCAR (+ an optional job script) into it,
#   3. sets the tested tag in INCAR (removes old occurrences first, also in
#      'A = 1; B = 2' lines; Windows CR characters are stripped), writes a
#      Gamma/MP KPOINTS for k-point tests, or KSPACING/KGAMMA for KSPACING tests,
#   4. forces a single-point run (NSW = 0, IBRION = -1, LWAVE/LCHARG = .FALSE.)
#      unless --keep-ionic / --keep-io are given,
#   5. writes conv_info.txt (read by conv_analyze.py) and a run_all.sh helper,
#   6. optionally submits every directory (sbatch/qsub auto-detected) or runs
#      them one after another with a command you give (-r).
#
# Usage examples
#   bash conv_test.sh -t encut  -v "300 350 400 450 500 550" -j job_vasp_slurm.sh -s
#   bash conv_test.sh -t kpoints -v "4 6 8 10 12" --mp          # n -> n n n
#   bash conv_test.sh -t kpoints -v "6x6x1 9x9x1 12x12x1"       # explicit meshes (2D)
#   bash conv_test.sh -t kspacing -v "0.5 0.4 0.3 0.25 0.2 0.15"
#   bash conv_test.sh -t sigma -v "0.2 0.1 0.05 0.02 0.01" --ismear 1
#   bash conv_test.sh -t LREAL -v "Auto .FALSE."                # any INCAR tag
#   bash conv_test.sh -t encut -v "300 400 500" -r "mpirun -np 8 vasp_std"
#   python conv_analyze.py conv_encut                           # analyse afterwards
#
# Dependencies: bash >= 4, awk, coreutils; sbatch or qsub for -s.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'EOF'
Usage: conv_test.sh -t TYPE -v "VALUES" [options]

TYPE (测试类型):
  encut      ENCUT values in eV                 e.g. -v "300 350 400 450"
  kpoints    meshes: "n" (= n n n) or "AxBxC"   e.g. -v "4 6 8" or "6x6x1 9x9x1"
  kspacing   KSPACING values (1/A, KPOINTS file removed, KGAMMA = .TRUE.)
  sigma      SIGMA values in eV (use --ismear to set ISMEAR as well)
  <TAG>      any other INCAR tag, e.g. -t LREAL -v "Auto .FALSE."

Options:
  -t, --type TYPE        test type (see above)                     [required]
  -v, --values "LIST"    space-separated values                     [required]
  -i, --incar FILE       template INCAR       (default: INCAR)
  -p, --poscar FILE      structure            (default: POSCAR)
  -k, --kpoints FILE     KPOINTS              (default: KPOINTS, optional)
  -P, --potcar FILE      POTCAR               (default: POTCAR, warning if missing)
  -j, --job FILE         job script copied into every directory as 'job.sh'
  -d, --outdir DIR       parent directory     (default: conv_<type>)
  -s, --submit           submit every directory (sbatch or qsub, auto-detected)
  -c, --submit-cmd CMD   explicit submit command, e.g. "sbatch -p debug"
  -r, --run CMD          run sequentially now, e.g. "mpirun -np 8 vasp_std"
      --gamma | --mp     k-mesh type for 'kpoints' tests (default: Gamma)
      --2d               scalar k values n -> n n 1 (slabs, 2D materials)
      --ismear N         also set ISMEAR (sigma tests)
      --keep-ionic       keep NSW/IBRION of the template (default: single point)
      --keep-io          keep LWAVE/LCHARG of the template (default: .FALSE.)
  -f, --force            overwrite existing test directories
  -n, --dry-run          only print what would be done
  -h, --help             show this help

After the jobs finished:  python conv_analyze.py <outdir>
EOF
}

die()  { echo "ERROR: $*" >&2; exit 1; }
info() { echo "[conv_test] $*"; }

# ------------------------------- defaults -----------------------------------
TYPE=""; VALUES=""; INCAR="INCAR"; POSCAR="POSCAR"; KPOINTS="KPOINTS"; POTCAR="POTCAR"
JOB=""; OUTDIR=""; SUBMIT=0; SUBMIT_CMD=""; RUN_CMD=""; KMODE="Gamma"; TWO_D=0
ISMEAR=""; KEEP_IONIC=0; KEEP_IO=0; FORCE=0; DRY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--type)       TYPE="${2:-}"; shift 2 ;;
        -v|--values)     VALUES="${2:-}"; shift 2 ;;
        -i|--incar)      INCAR="${2:-}"; shift 2 ;;
        -p|--poscar)     POSCAR="${2:-}"; shift 2 ;;
        -k|--kpoints)    KPOINTS="${2:-}"; shift 2 ;;
        -P|--potcar)     POTCAR="${2:-}"; shift 2 ;;
        -j|--job)        JOB="${2:-}"; shift 2 ;;
        -d|--outdir)     OUTDIR="${2:-}"; shift 2 ;;
        -s|--submit)     SUBMIT=1; shift ;;
        -c|--submit-cmd) SUBMIT=1; SUBMIT_CMD="${2:-}"; shift 2 ;;
        -r|--run)        RUN_CMD="${2:-}"; shift 2 ;;
        --gamma)         KMODE="Gamma"; shift ;;
        --mp)            KMODE="Monkhorst-Pack"; shift ;;
        --2d)            TWO_D=1; shift ;;
        --ismear)        ISMEAR="${2:-}"; shift 2 ;;
        --keep-ionic)    KEEP_IONIC=1; shift ;;
        --keep-io)       KEEP_IO=1; shift ;;
        -f|--force)      FORCE=1; shift ;;
        -n|--dry-run)    DRY=1; shift ;;
        -h|--help)       usage; exit 0 ;;
        *) usage >&2; die "unknown option '$1'" ;;
    esac
done

[[ -n "$TYPE" ]]   || { usage >&2; die "-t/--type is required"; }
[[ -n "$VALUES" ]] || { usage >&2; die "-v/--values is required"; }
[[ -f "$INCAR" ]]  || die "template INCAR '$INCAR' not found"
[[ -f "$POSCAR" ]] || die "structure '$POSCAR' not found"
[[ -f "$POTCAR" ]] || echo "WARNING: '$POTCAR' not found - add a POTCAR to every directory before running" >&2
[[ -z "$JOB" || -f "$JOB" ]] || die "job script '$JOB' not found"
[[ $SUBMIT -eq 1 && -n "$RUN_CMD" ]] && die "use either -s/--submit or -r/--run, not both"

# resolve the submit command before anything is created
if [[ $SUBMIT -eq 1 && $DRY -eq 0 ]]; then
    [[ -n "$JOB" ]] || die "-s needs a job script (-j job_vasp_slurm.sh)"
    if [[ -z "$SUBMIT_CMD" ]]; then
        if command -v sbatch >/dev/null 2>&1; then SUBMIT_CMD="sbatch"
        elif command -v qsub >/dev/null 2>&1; then SUBMIT_CMD="qsub"
        else die "neither sbatch nor qsub found; use -c 'your submit command'"; fi
    fi
fi

TYPE_LC="$(echo "$TYPE" | tr '[:upper:]' '[:lower:]')"
case "$TYPE_LC" in
    encut)    TAG="ENCUT" ;;
    kpoints)  TAG="KPOINTS" ;;
    kspacing) TAG="KSPACING" ;;
    sigma)    TAG="SIGMA" ;;
    *)        TAG="$(echo "$TYPE" | tr '[:lower:]' '[:upper:]')"
              [[ "$TAG" =~ ^[A-Z][A-Z0-9_]*$ ]] || die "'$TYPE' is not a valid INCAR tag name" ;;
esac
[[ -n "$OUTDIR" ]] || OUTDIR="conv_${TYPE_LC}"

if [[ "$TAG" == "KPOINTS" ]]; then :; elif [[ "$TAG" != "KSPACING" && ! -f "$KPOINTS" ]] && \
     ! grep -qiE '^[[:space:]]*KSPACING[[:space:]]*=' "$INCAR"; then
    echo "WARNING: no '$KPOINTS' and no KSPACING in INCAR - VASP would use KSPACING = 0.5" >&2
fi
if grep -q $'\r' "$INCAR"; then
    echo "WARNING: '$INCAR' has Windows (CR) line endings; they are removed in the test copies" >&2
fi

# --------------------------- INCAR manipulation -----------------------------
# set_tag FILE TAG VALUE : remove TAG everywhere (case-insensitive, also inside
# 'A = 1; B = 2' lines), then append 'TAG = VALUE'. Empty VALUE = delete only.
set_tag() {
    local file="$1" tag="$2" val="${3-}"
    awk -v tag="$tag" -v val="$val" '
    BEGIN { T = toupper(tag) }
    {
        sub(/\r$/, "")
        line = $0; cmt = ""
        p = match(line, /[#!]/)
        if (p > 0) { cmt = substr(line, p); line = substr(line, 1, p - 1) }
        n = split(line, seg, ";")
        out = ""; removed = 0
        for (i = 1; i <= n; i++) {
            key = seg[i]; sub(/=.*/, "", key); gsub(/[ \t]/, "", key)
            if (index(seg[i], "=") > 0 && toupper(key) == T) { removed = 1; continue }
            out = (out == "" ? seg[i] : out ";" seg[i])
        }
        if (removed) {
            if (out ~ /^[ \t]*$/) next          # whole line was this tag
            print out; next
        }
        print $0
    }
    END { if (val != "") printf "%s = %s\n", tag, val }' "$file" > "$file.tmp_ct"
    mv "$file.tmp_ct" "$file"
}

# directory label for a value: keep [A-Za-z0-9._-], strip leading/trailing dots
# (Windows/WSL file systems do not allow names ending with '.')
label() {
    local l="${1//[^A-Za-z0-9._-]/_}"
    while [[ "$l" == .* ]]; do l="${l#.}"; done
    while [[ "$l" == *. ]]; do l="${l%.}"; done
    echo "${l:-value}"
}

# write_kpoints FILE "n1 n2 n3" MODE
write_kpoints() {
    printf '%s mesh %s (conv_test.sh)\n0\n%s\n  %s\n  0 0 0\n' "$3" "$2" "$3" "$2" > "$1"
}

# --------------------------------- main loop --------------------------------
if [[ -d "$OUTDIR" && $FORCE -eq 0 && $DRY -eq 0 ]]; then
    for v in $VALUES; do
        lab="$(label "$v")"
        [[ ! -e "$OUTDIR/${TYPE_LC}_${lab}" ]] || die "'$OUTDIR/${TYPE_LC}_${lab}' exists (use -f to overwrite)"
    done
fi
[[ $DRY -eq 1 ]] || mkdir -p "$OUTDIR"

DIRS=()
for v in $VALUES; do
    lab="$(label "$v")"
    d="$OUTDIR/${TYPE_LC}_${lab}"
    DIRS+=("$d")
    mesh=""
    if [[ "$TAG" == "KPOINTS" ]]; then
        if [[ "$v" =~ ^[0-9]+$ ]]; then
            if [[ $TWO_D -eq 1 ]]; then mesh="$v $v 1"; else mesh="$v $v $v"; fi
        elif [[ "$v" =~ ^[0-9]+[xX][0-9]+[xX][0-9]+$ ]]; then
            mesh="$(echo "$v" | tr 'xX' '  ')"
        else
            die "k-point value '$v' must be 'n' or 'AxBxC'"
        fi
    fi
    if [[ $DRY -eq 1 ]]; then
        info "would create $d with $TAG = $v${mesh:+ (mesh $mesh $KMODE)}"
        continue
    fi
    rm -rf "$d"; mkdir -p "$d"
    cp "$INCAR" "$d/INCAR"
    cp "$POSCAR" "$d/POSCAR"
    [[ -f "$POTCAR" ]] && cp "$POTCAR" "$d/POTCAR"
    [[ -n "$JOB" ]] && cp "$JOB" "$d/job.sh"
    case "$TAG" in
        KPOINTS)
            write_kpoints "$d/KPOINTS" "$mesh" "$KMODE"
            set_tag "$d/INCAR" KSPACING ""; set_tag "$d/INCAR" KGAMMA "" ;;
        KSPACING)
            set_tag "$d/INCAR" KSPACING "$v"; set_tag "$d/INCAR" KGAMMA ".TRUE." ;;
        *)
            [[ -f "$KPOINTS" ]] && cp "$KPOINTS" "$d/KPOINTS"
            set_tag "$d/INCAR" "$TAG" "$v" ;;
    esac
    [[ "$TAG" == "SIGMA" && -n "$ISMEAR" ]] && set_tag "$d/INCAR" ISMEAR "$ISMEAR"
    if [[ $KEEP_IONIC -eq 0 ]]; then
        set_tag "$d/INCAR" NSW 0; set_tag "$d/INCAR" IBRION -1
    fi
    if [[ $KEEP_IO -eq 0 ]]; then
        set_tag "$d/INCAR" LWAVE ".FALSE."; set_tag "$d/INCAR" LCHARG ".FALSE."
    fi
    printf 'type = %s\nparam = %s\nvalue = %s\n' "$TYPE_LC" "$TAG" "$v" > "$d/conv_info.txt"
    [[ -n "$mesh" ]] && printf 'mesh = %s\nkmode = %s\n' "$mesh" "$KMODE" >> "$d/conv_info.txt"
    info "created $d ($TAG = $v)"
done

[[ $DRY -eq 1 ]] && exit 0

# helper to run all tests inside ONE job (call it from your own job script)
cat > "$OUTDIR/run_all.sh" <<'EOF'
#!/usr/bin/env bash
# Run every test directory sequentially inside one job allocation.
# Usage (inside a job script):  VASP_CMD="mpirun -np 64 vasp_std" bash run_all.sh
set -uo pipefail
VASP_CMD="${VASP_CMD:-mpirun vasp_std}"
cd "$(dirname "${BASH_SOURCE[0]}")"
for d in */; do
    [[ -f "$d/conv_info.txt" ]] || continue
    if grep -q "General timing and accounting" "$d/OUTCAR" 2>/dev/null; then
        echo "skip $d (finished)"; continue
    fi
    echo "=== $d  $(date)"
    (cd "$d" && $VASP_CMD > vasp.out 2>&1) || echo "WARNING: $d returned non-zero"
done
EOF
chmod +x "$OUTDIR/run_all.sh"

# --------------------------- submit or run ----------------------------------
if [[ $SUBMIT -eq 1 ]]; then
    for d in "${DIRS[@]}"; do
        jid="$(cd "$d" && $SUBMIT_CMD job.sh)" || die "submission failed in $d"
        info "$d: $jid"
    done
elif [[ -n "$RUN_CMD" ]]; then
    for d in "${DIRS[@]}"; do
        info "running in $d: $RUN_CMD"
        (cd "$d" && eval "$RUN_CMD" > vasp.out 2>&1) || echo "WARNING: '$RUN_CMD' failed in $d (see $d/vasp.out)" >&2
    done
else
    info "inputs ready. Submit each with your job script, or call '$OUTDIR/run_all.sh' inside one job."
fi
info "analyse with: python ${SCRIPT_DIR}/conv_analyze.py $OUTDIR"
