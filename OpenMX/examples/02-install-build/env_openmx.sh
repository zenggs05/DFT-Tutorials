# OpenMX runtime environment. Use it with "source env_openmx.sh" (do not execute it).
# DFT-Tutorials, OpenMX chapter 2.
#
# Edit the defaults below or set the variables before sourcing:
#   OPENMX_ROOT       OpenMX tree (source build: $HOME/software/openmx4.0 ; deb: /opt/openmx4.0)
#   OMX_DATA          DFT_DATA19 directory used by omx_setdata (default: $OPENMX_ROOT/DFT_DATA19)
#   OMX_USE_ONEAPI=1  source the Intel oneAPI environment (Intel build)
#   ONEAPI_ROOT       oneAPI location (default: /opt/intel/oneapi)
#   OMX_DEB_RUNTIME=1 add /opt/openmx4.0/runtime/lib to LD_LIBRARY_PATH (deb package tools
#                     such as DosMain reporting "libmpifort.so.12: cannot open shared object file")
#
export OPENMX_ROOT=${OPENMX_ROOT:-$HOME/software/openmx4.0}
export OMX_DATA=${OMX_DATA:-$OPENMX_ROOT/DFT_DATA19}

# executables: prefer $OPENMX_ROOT/bin, then $OPENMX_ROOT/work
if [ -d "$OPENMX_ROOT/bin" ]; then
  export PATH=$OPENMX_ROOT/bin:$PATH
elif [ -d "$OPENMX_ROOT/work" ]; then
  export PATH=$OPENMX_ROOT/work:$PATH
fi

# Intel oneAPI runtime (MKL, Intel MPI, ifx runtime)
if [ "${OMX_USE_ONEAPI:-0}" = "1" ]; then
  if [ -f "${ONEAPI_ROOT:-/opt/intel/oneapi}/setvars.sh" ]; then
    # shellcheck disable=SC1091
    source "${ONEAPI_ROOT:-/opt/intel/oneapi}/setvars.sh" > /dev/null
  else
    echo "env_openmx.sh: setvars.sh not found under ${ONEAPI_ROOT:-/opt/intel/oneapi}" >&2
  fi
fi

# runtime libraries shipped with the official Debian package (opt-in)
if [ "${OMX_DEB_RUNTIME:-0}" = "1" ] && [ -d /opt/openmx4.0/runtime/lib ]; then
  export LD_LIBRARY_PATH=/opt/openmx4.0/runtime/lib:${LD_LIBRARY_PATH:-}
fi

# threading and stack settings (keep OMP_NUM_THREADS equal to the -nt option)
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OMP_STACKSIZE=${OMP_STACKSIZE:-512M}
export OPENBLAS_NUM_THREADS=1
ulimit -s unlimited 2>/dev/null || true

# helper: rewrite DATA.PATH in the given input files to $OMX_DATA
#   usage: omx_setdata *.dat
omx_setdata() {
  local f
  for f in "$@"; do
    if grep -qi '^[[:space:]]*DATA\.PATH' "$f"; then
      sed -i "s#^[[:space:]]*DATA\.PATH.*#DATA.PATH   ${OMX_DATA}#I" "$f"
    else
      printf 'DATA.PATH   %s\n' "${OMX_DATA}" | cat - "$f" > "$f.tmp" && mv "$f.tmp" "$f"
    fi
    echo "$f: DATA.PATH -> ${OMX_DATA}"
  done
}

if [ ! -d "$OMX_DATA/PAO" ] || [ ! -d "$OMX_DATA/VPS" ]; then
  echo "env_openmx.sh: warning: $OMX_DATA does not contain PAO/ and VPS/ (set OMX_DATA)" >&2
fi
command -v openmx >/dev/null 2>&1 || echo "env_openmx.sh: warning: openmx is not in PATH" >&2
