#!/usr/bin/env bash
# =============================================================================
# env_check.sh -- self-check of a DFT working environment (login node, WSL or
#                 workstation) with a fix suggestion for every problem.
#                 计算环境自检：Python 与关键包版本、VASP_PP_PATH / PMG_VASP_PSP_DIR /
#                 MP_API_KEY、~/.vaspkit、ulimit -s、OMP_NUM_THREADS、OMP_STACKSIZE
#                 (Multiwfn)、常用程序与调度器命令，逐项给出修复建议。
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
#
# What is checked
#   python    : interpreter found, version >= 3.9, active conda environment
#   packages  : ase, pymatgen, phonopy, sumo, pyprocar, py4vasp (+ numpy, scipy,
#               matplotlib, spglib, seekpath, mp-api) - versions via
#               importlib.metadata, nothing is imported
#   VASP      : VASP_PP_PATH (directory that CONTAINS potpaw_PBE/), PMG_VASP_PSP_DIR
#               (environment or ~/.config/.pmgrc.yaml / ~/.pmgrc.yaml), MP_API_KEY /
#               PMG_MAPI_KEY (only the length is shown, never the key)
#   VASPKIT   : ~/.vaspkit exists, its PBE_PATH directory exists, vaspkit in PATH
#   limits    : ulimit -s (stack), OMP_NUM_THREADS, OMP_STACKSIZE (Multiwfn needs
#               about 200M), number of CPU cores
#   programs  : vasp_std/vasp_gam/vasp_ncl, g16/formchk (+ g16root, GAUSS_SCRDIR),
#               cp2k.psmp (+ CP2K_DATA_DIR), openmx, Multiwfn (+ Multiwfnpath),
#               mpirun/srun, sbatch/qsub/bsub, dos2unix
#
# Usage
#   bash env_check.sh                 # full report
#   bash env_check.sh -q              # only WARN/FAIL lines
#   bash env_check.sh -p ~/miniforge3/envs/dft/bin/python
#   bash env_check.sh --no-color > env_report.txt
#
# Exit status: 0 = no FAIL (warnings allowed), 1 = at least one FAIL.
# Nothing is modified: the script only reads files and environment variables.
#
# Dependencies: bash >= 4, coreutils; python (any version, to report it).
# =============================================================================
set -uo pipefail

usage() {
    cat <<'EOF'
Usage: env_check.sh [options]

Checks Python and package versions, VASP/pymatgen/Materials Project settings,
~/.vaspkit, ulimit -s, OMP_NUM_THREADS, OMP_STACKSIZE and the DFT programs in
PATH, and prints a fix for every problem found.

Options:
  -p, --python PY   Python interpreter to check (default: python3, then python)
  -q, --quiet       print only WARN and FAIL lines (and the summary)
      --no-color    plain text output (default when stdout is not a terminal)
  -h, --help        show this help

Exit status: 0 = no FAIL, 1 = at least one FAIL.
EOF
}

PY=""
QUIET=0
COLOR=1
[ -t 1 ] || COLOR=0
while [ $# -gt 0 ]; do
    case "$1" in
        -p|--python) [ $# -ge 2 ] || { echo "missing value for $1" >&2; exit 2; }; PY="$2"; shift 2 ;;
        -q|--quiet) QUIET=1; shift ;;
        --no-color) COLOR=0; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [ "${COLOR}" -eq 1 ]; then
    C_OK=$'\033[32m'; C_WARN=$'\033[33m'; C_FAIL=$'\033[31m'; C_HEAD=$'\033[1m'; C_OFF=$'\033[0m'
else
    C_OK=""; C_WARN=""; C_FAIL=""; C_HEAD=""; C_OFF=""
fi
N_OK=0; N_WARN=0; N_FAIL=0

section() { [ "${QUIET}" -eq 1 ] || printf '\n%s== %s ==%s\n' "${C_HEAD}" "$1" "${C_OFF}"; }
ok()   { N_OK=$((N_OK + 1)); [ "${QUIET}" -eq 1 ] || printf '%s[ OK ]%s %s\n' "${C_OK}" "${C_OFF}" "$1"; }
info() { [ "${QUIET}" -eq 1 ] || printf '[INFO] %s\n' "$1"; }
warn() {   # warn MESSAGE [FIX]
    N_WARN=$((N_WARN + 1))
    printf '%s[WARN]%s %s\n' "${C_WARN}" "${C_OFF}" "$1"
    [ $# -ge 2 ] && printf '       -> fix: %s\n' "$2"
    return 0
}
fail() {   # fail MESSAGE [FIX]
    N_FAIL=$((N_FAIL + 1))
    printf '%s[FAIL]%s %s\n' "${C_FAIL}" "${C_OFF}" "$1"
    [ $# -ge 2 ] && printf '       -> fix: %s\n' "$2"
    return 0
}
have() { command -v "$1" >/dev/null 2>&1; }
expand_tilde() { local p="$1"; [ "${p:0:1}" = "~" ] && p="${HOME}${p:1}"; printf '%s' "$p"; }

# ------------------------------------------------------------------ system
section "System"
info "host $(hostname 2>/dev/null || echo '?'), $(uname -sr 2>/dev/null || echo '?'), bash ${BASH_VERSION}"
NCORES=$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo 0)
info "CPU cores visible here: ${NCORES} (login nodes are NOT for production runs)"
if [ "${BASH_VERSINFO[0]}" -lt 4 ]; then
    warn "bash ${BASH_VERSION} is old; repository scripts need bash >= 4" "use a newer bash (macOS: brew install bash)"
fi

# ------------------------------------------------------------------ python
section "Python"
if [ -z "${PY}" ]; then
    if have python3; then PY=python3; elif have python; then PY=python; fi
fi
PY_OK=0
if [ -z "${PY}" ] || ! "${PY}" -c 'import sys' >/dev/null 2>&1; then
    fail "no working Python interpreter found${PY:+ (${PY})}" \
         "install Miniforge and create the env: mamba env create -f Common/examples/01-env-setup/environment.yml"
else
    PYVER=$("${PY}" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')
    PYPATH=$("${PY}" -c 'import sys; print(sys.executable)')
    if "${PY}" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)'; then
        ok "Python ${PYVER} (${PYPATH})"
        PY_OK=1
    else
        fail "Python ${PYVER} (${PYPATH}) is older than 3.9" "conda activate dft  (or create it: mamba env create -f environment.yml)"
    fi
    if [ -n "${CONDA_DEFAULT_ENV:-}" ]; then
        info "conda environment: ${CONDA_DEFAULT_ENV}"
        [ "${CONDA_DEFAULT_ENV}" = "base" ] && warn "working in the conda 'base' environment" "conda activate dft  (keep base minimal)"
    else
        info "no conda environment active"
    fi
fi

# ------------------------------------------------------------------ packages
section "Python packages"
if [ "${PY_OK}" -eq 1 ]; then
    # dist-name|required(1)/optional(0)|install hint
    PKGS="ase|1|conda install -c conda-forge ase
pymatgen|1|conda install -c conda-forge pymatgen
phonopy|1|conda install -c conda-forge phonopy
sumo|0|pip install sumo
pyprocar|0|pip install pyprocar
py4vasp|0|pip install py4vasp
numpy|1|conda install -c conda-forge numpy
scipy|0|conda install -c conda-forge scipy
matplotlib|1|conda install -c conda-forge matplotlib
spglib|1|conda install -c conda-forge spglib
seekpath|0|pip install seekpath
mp-api|0|pip install mp-api"
    while IFS='|' read -r dist req hint; do
        [ -n "${dist}" ] || continue
        ver=$("${PY}" - "${dist}" <<'PYEOF' 2>/dev/null
import sys
from importlib import metadata
try:
    print(metadata.version(sys.argv[1]))
except metadata.PackageNotFoundError:
    print("")
PYEOF
)
        if [ -n "${ver}" ]; then
            ok "$(printf '%-11s %s' "${dist}" "${ver}")"
        elif [ "${req}" = "1" ]; then
            warn "${dist} is not installed (needed by many scripts of this repository)" "${hint}"
        else
            info "$(printf '%-11s not installed (optional) -> %s' "${dist}" "${hint}")"
        fi
    done <<< "${PKGS}"
else
    info "skipped (no usable Python)"
fi

# ------------------------------------------------------------------ VASP / pymatgen / MP
section "VASP pseudopotentials, pymatgen, Materials Project"
if [ -n "${VASP_PP_PATH:-}" ]; then
    if [ -d "${VASP_PP_PATH}/potpaw_PBE" ] || [ -d "${VASP_PP_PATH}/potpaw_PBE.54" ] \
       || [ -d "${VASP_PP_PATH}/potpaw_PBE.64" ]; then
        ok "VASP_PP_PATH=${VASP_PP_PATH} (contains potpaw_PBE*)"
        [ -d "${VASP_PP_PATH}/potpaw_PBE" ] || \
            warn "ASE looks for \$VASP_PP_PATH/potpaw_PBE (no plain potpaw_PBE directory found)" \
                 "ln -s potpaw_PBE.64 \"${VASP_PP_PATH}/potpaw_PBE\"  (or name it as ASE expects)"
    elif [ -d "${VASP_PP_PATH}" ]; then
        warn "VASP_PP_PATH=${VASP_PP_PATH} exists but has no potpaw_PBE/ inside" \
             "point VASP_PP_PATH to the directory that CONTAINS potpaw_PBE/, not to potpaw_PBE itself"
    else
        warn "VASP_PP_PATH=${VASP_PP_PATH} is not a directory" "fix the path in ~/.bashrc"
    fi
else
    warn "VASP_PP_PATH is not set (ASE cannot build POTCAR files)" \
         "export VASP_PP_PATH=\$HOME/vasp_pp   # directory containing potpaw_PBE/ (see bashrc_dft_snippet.sh)"
fi

PMGRC=""
for f in "${HOME}/.config/.pmgrc.yaml" "${HOME}/.pmgrc.yaml"; do
    [ -f "$f" ] && { PMGRC="$f"; break; }
done
pmg_get() {   # pmg_get KEY : value from the pmgrc file
    [ -n "${PMGRC}" ] || return 0
    awk -v k="$1" -F':' '$1 == k { sub(/^[^:]*:[ \t]*/, ""); gsub(/["\047]/, ""); print; exit }' "${PMGRC}"
}
PSP="${PMG_VASP_PSP_DIR:-$(pmg_get PMG_VASP_PSP_DIR)}"
if [ -n "${PSP}" ]; then
    PSP=$(expand_tilde "${PSP}")
    if [ -d "${PSP}" ]; then
        ok "PMG_VASP_PSP_DIR=${PSP}${PMGRC:+ (config: ${PMGRC})}"
        ls -d "${PSP}"/POT_* >/dev/null 2>&1 || \
            warn "no POT_* sub-directories in ${PSP}" "pmg config -p <raw_potcar_dir> ${PSP}"
    else
        warn "PMG_VASP_PSP_DIR=${PSP} does not exist" "pmg config -p <raw_potcar_dir> <out_dir>; pmg config --add PMG_VASP_PSP_DIR <out_dir>"
    fi
else
    warn "PMG_VASP_PSP_DIR not set (pymatgen cannot write POTCAR)" \
         "pmg config -p ~/vasp_pp ~/pmg_potcars && pmg config --add PMG_VASP_PSP_DIR ~/pmg_potcars"
fi
[ -n "${PMGRC}" ] && [ "$(stat -c %a "${PMGRC}" 2>/dev/null || echo 600)" != "600" ] && \
    warn "${PMGRC} is readable by others (it may contain your API key)" "chmod 600 ${PMGRC}"

key_check() {   # key_check NAME VALUE
    local name="$1" val="$2"
    if [ -z "${val}" ]; then
        return 1
    elif [ "${#val}" -eq 32 ]; then
        ok "${name} is set (32 characters, new Materials Project API)"
    else
        warn "${name} has ${#val} characters; keys of the new API have 32" "create a new key at https://next-gen.materialsproject.org/api"
    fi
    return 0
}
if ! key_check MP_API_KEY "${MP_API_KEY:-}"; then
    if ! key_check PMG_MAPI_KEY "${PMG_MAPI_KEY:-$(pmg_get PMG_MAPI_KEY)}"; then
        warn "no Materials Project key (MP_API_KEY / PMG_MAPI_KEY)" \
             "export MP_API_KEY=<32-character key>  (or pmg config --add PMG_MAPI_KEY <key>)"
    fi
fi
if [ -n "${MP_API_ENDPOINT:-}${PMG_MAPI_ENDPOINT:-}" ]; then
    warn "MP_API_ENDPOINT/PMG_MAPI_ENDPOINT is set (old v2 endpoints no longer work)" "unset MP_API_ENDPOINT PMG_MAPI_ENDPOINT"
fi

# ------------------------------------------------------------------ VASPKIT
section "VASPKIT"
if [ -f "${HOME}/.vaspkit" ]; then
    ok "~/.vaspkit exists"
    pbe=$(awk '$1 == "PBE_PATH" { print $2; exit }' "${HOME}/.vaspkit")
    if [ -n "${pbe}" ]; then
        pbe=$(expand_tilde "${pbe}")
        if [ -d "${pbe}" ]; then ok "PBE_PATH=${pbe}"; else
            warn "PBE_PATH=${pbe} in ~/.vaspkit does not exist" "edit PBE_PATH in ~/.vaspkit (directory with Si/POTCAR, Fe_pv/POTCAR, ...)"
        fi
    else
        warn "no PBE_PATH line in ~/.vaspkit" "copy Common/examples/01-env-setup/vaspkit_rc.example and edit it"
    fi
    grep -q $'\r' "${HOME}/.vaspkit" 2>/dev/null && warn "~/.vaspkit has Windows (CRLF) line endings" "dos2unix ~/.vaspkit"
else
    if have vaspkit; then
        warn "vaspkit found but ~/.vaspkit is missing" "cp Common/examples/01-env-setup/vaspkit_rc.example ~/.vaspkit && vi ~/.vaspkit"
    else
        info "VASPKIT not installed (optional); template: Common/examples/01-env-setup/vaspkit_rc.example"
    fi
fi
have vaspkit && ok "vaspkit in PATH: $(command -v vaspkit)"

# ------------------------------------------------------------------ limits & OpenMP
section "Stack size and OpenMP"
STACK=$(ulimit -s 2>/dev/null || echo "?")
STACK_H=$(ulimit -H -s 2>/dev/null || echo "?")
if [ "${STACK}" = "unlimited" ]; then
    ok "ulimit -s = unlimited"
else
    warn "ulimit -s = ${STACK} kB (hard limit ${STACK_H}); Fortran codes (VASP, CP2K, Multiwfn) may crash with SIGSEGV" \
         "put 'ulimit -s unlimited' in every job script (and in ~/.bashrc for interactive tools)"
fi
if [ -n "${OMP_NUM_THREADS:-}" ]; then
    if [[ "${OMP_NUM_THREADS}" =~ ^[0-9]+$ ]]; then
        ok "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
        if [ "${NCORES}" -gt 0 ] && [ "${OMP_NUM_THREADS}" -gt "${NCORES}" ]; then
            warn "OMP_NUM_THREADS (${OMP_NUM_THREADS}) > visible cores (${NCORES})" "export OMP_NUM_THREADS=1 for pure MPI, = cpus-per-task for hybrid runs"
        fi
    else
        warn "OMP_NUM_THREADS='${OMP_NUM_THREADS}' is not a single integer" "export OMP_NUM_THREADS=1"
    fi
else
    warn "OMP_NUM_THREADS is not set: OpenMP libraries may start one thread per core in every MPI rank" \
         "export OMP_NUM_THREADS=1 in pure-MPI job scripts (= --cpus-per-task for MPI x OpenMP)"
fi
if [ -n "${OMP_STACKSIZE:-}" ]; then
    ok "OMP_STACKSIZE=${OMP_STACKSIZE}"
else
    warn "OMP_STACKSIZE is not set (Multiwfn crashes on large systems with the small default)" \
         "export OMP_STACKSIZE=200M   (Multiwfn manual recommendation; CP2K/OpenMX jobs often use 512m)"
fi

# ------------------------------------------------------------------ programs
section "Programs and schedulers in PATH"
found_any=0
for exe in vasp_std vasp_gam vasp_ncl g16 formchk cp2k.psmp openmx Multiwfn; do
    if have "${exe}"; then ok "${exe}: $(command -v "${exe}")"; found_any=1; fi
done
[ "${found_any}" -eq 1 ] || info "no DFT executable in PATH (normal on a workstation; on a cluster: module avail / module load ...)"

if have g16 || [ -n "${g16root:-}" ]; then
    [ -n "${g16root:-}" ] || warn "g16 found but g16root is not set" "export g16root=/path/containing/g16; source \$g16root/g16/bsd/g16.profile"
    if [ -n "${GAUSS_SCRDIR:-}" ]; then
        if [ -d "${GAUSS_SCRDIR}" ] && [ -w "${GAUSS_SCRDIR}" ]; then ok "GAUSS_SCRDIR=${GAUSS_SCRDIR} (writable)"
        else warn "GAUSS_SCRDIR=${GAUSS_SCRDIR} missing or not writable" "mkdir -p \"\$GAUSS_SCRDIR\" (use node-local disk in job scripts)"; fi
    else
        warn "GAUSS_SCRDIR not set (scratch goes to the working directory)" "export GAUSS_SCRDIR=/tmp/\$USER/g16_\$SLURM_JOB_ID in job scripts"
    fi
fi
if have cp2k.psmp || have cp2k; then
    if [ -n "${CP2K_DATA_DIR:-}" ] && [ -d "${CP2K_DATA_DIR}" ]; then ok "CP2K_DATA_DIR=${CP2K_DATA_DIR}"
    else info "CP2K_DATA_DIR not set: BASIS_SET_FILE_NAME/POTENTIAL_FILE_NAME need full paths unless the build has a default data dir"; fi
fi
if have Multiwfn || [ -n "${Multiwfnpath:-}" ]; then
    if [ -n "${Multiwfnpath:-}" ] && [ -f "${Multiwfnpath}/settings.ini" ]; then ok "Multiwfnpath=${Multiwfnpath} (settings.ini found)"
    else warn "Multiwfnpath not set or settings.ini missing" "export Multiwfnpath=/path/to/Multiwfn (directory with settings.ini; set nthreads there)"; fi
fi

MPI=""; for exe in mpirun srun; do have "${exe}" && MPI="${MPI} ${exe}"; done
if [ -n "${MPI}" ]; then ok "MPI launcher(s):${MPI}"; else info "no mpirun/srun in PATH (load an MPI module on the cluster)"; fi
SCHED=""; for exe in sbatch qsub bsub; do have "${exe}" && SCHED="${SCHED} ${exe}"; done
if [ -n "${SCHED}" ]; then ok "scheduler command(s):${SCHED}"; else info "no sbatch/qsub/bsub (not a cluster login node?)"; fi
if have dos2unix; then ok "dos2unix: $(command -v dos2unix)"
else warn "dos2unix not found (files edited on Windows keep CRLF endings)" "conda install -c conda-forge dos2unix  or: sed -i 's/\\r\$//' FILE"; fi

# ------------------------------------------------------------------ summary
printf '\n%sSummary:%s %d OK, %d WARN, %d FAIL\n' "${C_HEAD}" "${C_OFF}" "${N_OK}" "${N_WARN}" "${N_FAIL}"
[ "${N_FAIL}" -eq 0 ] || exit 1
exit 0
