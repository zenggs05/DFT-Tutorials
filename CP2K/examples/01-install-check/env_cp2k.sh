#!/usr/bin/env bash
# env_cp2k.sh - set up the CP2K runtime environment.
#
# Usage (must be SOURCED, not executed):
#     source env_cp2k.sh                 # uses CP2K_METHOD below
#     CP2K_METHOD=conda source env_cp2k.sh
#
# Edit the paths in the "settings" block to match your installation.
# After sourcing, the following variables are defined:
#     CP2K_EXE        name of the CP2K executable to use (cp2k.psmp / cp2k.ssmp)
#     CP2K_DATA_DIR   data directory with basis sets and pseudopotentials
#     OMP_NUM_THREADS default number of OpenMP threads (if not already set)

# ------------------------------- settings ----------------------------------
CP2K_METHOD=${CP2K_METHOD:-toolchain}   # static | conda | toolchain | spack | module

# (1) static binary from the GitHub release page (see tutorial section 4)
CP2K_STATIC_ROOT=${CP2K_STATIC_ROOT:-$HOME/opt/cp2k-2026.2}

# (2) conda-forge environment name (see tutorial section 5)
CP2K_CONDA_ENV=${CP2K_CONDA_ENV:-cp2k_env}

# (3) toolchain + build_cp2k.sh build: the CP2K_ROOT of the unpacked tarball
CP2K_SRC_ROOT=${CP2K_SRC_ROOT:-$HOME/src/cp2k-2026.2}

# (4) spack build via make_cp2k.sh: same CP2K_ROOT, uses install/bin
# (5) module name provided by your HPC centre
CP2K_MODULE=${CP2K_MODULE:-cp2k/2026.2}
# ---------------------------------------------------------------------------

# refuse to run as a script: environment changes would be lost
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "ERROR: env_cp2k.sh must be sourced:  source env_cp2k.sh" >&2
  exit 1
fi

case "${CP2K_METHOD}" in
  static)
    export PATH="${CP2K_STATIC_ROOT}/bin:${PATH}"
    export CP2K_DATA_DIR="${CP2K_STATIC_ROOT}/data"
    export CP2K_EXE=cp2k.ssmp          # symlink to cp2k-2026.2-Linux-gnu-x86_64.ssmp
    ;;
  conda)
    # make "conda activate" available in non-interactive shells
    if command -v conda > /dev/null 2>&1; then
      eval "$(conda shell.bash hook)"
      conda activate "${CP2K_CONDA_ENV}"
      export CP2K_DATA_DIR="${CONDA_PREFIX}/share/cp2k/data"
    else
      echo "ERROR: conda not found in PATH" >&2
      return 1
    fi
    export CP2K_EXE=cp2k.psmp
    ;;
  toolchain)
    if [[ -f "${CP2K_SRC_ROOT}/install/cp2k_env" ]]; then
      # shellcheck disable=SC1091
      source "${CP2K_SRC_ROOT}/install/cp2k_env"
    else
      echo "ERROR: ${CP2K_SRC_ROOT}/install/cp2k_env not found (did build_cp2k.sh finish?)" >&2
      return 1
    fi
    export CP2K_EXE=cp2k.psmp
    # The data directory is compiled in (install/share/cp2k/data); set it
    # explicitly so that the value is visible and easy to check.
    if [[ -d "${CP2K_SRC_ROOT}/install/share/cp2k/data" ]]; then
      export CP2K_DATA_DIR="${CP2K_SRC_ROOT}/install/share/cp2k/data"
    fi
    ;;
  spack)
    # make_cp2k.sh installs a launcher; adding install/bin to PATH is enough
    # for interactive use, the launcher sets LD_LIBRARY_PATH for each run:
    #     install/bin/launch cp2k.psmp -i x.inp -o x.out
    export PATH="${CP2K_SRC_ROOT}/install/bin:${PATH}"
    export CP2K_EXE=cp2k.psmp
    ;;
  module)
    module load "${CP2K_MODULE}" || return 1
    export CP2K_EXE=cp2k.psmp
    ;;
  *)
    echo "ERROR: unknown CP2K_METHOD='${CP2K_METHOD}'" >&2
    return 1
    ;;
esac

# sensible defaults for OpenMP and the Fortran stack
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
ulimit -s unlimited 2> /dev/null || true

echo "CP2K_METHOD   = ${CP2K_METHOD}"
echo "CP2K_EXE      = ${CP2K_EXE} -> $(command -v "${CP2K_EXE}" || echo 'NOT FOUND')"
echo "CP2K_DATA_DIR = ${CP2K_DATA_DIR:-<not set, compiled-in default is used>}"
echo "OMP_NUM_THREADS = ${OMP_NUM_THREADS}"
