#!/bin/bash
#
# Build OpenMX 4.0 + patch 4.0.1 with Intel oneAPI (mpiicx/mpiifx + MKL + Intel MPI).
# DFT-Tutorials, OpenMX chapter 2 (02-install.md, section 4).
#
# Usage:
#   bash build_openmx_oneapi.sh
#   ONEAPI_ROOT=/opt/intel/oneapi PREFIX=/apps/openmx bash build_openmx_oneapi.sh
#   ARCH_FLAG="-march=core-avx2" RUN_TEST=0 bash build_openmx_oneapi.sh
#
# Environment variables (all optional):
#   PREFIX        install root, the code goes to $PREFIX/openmx4.0   (default: $HOME/software)
#   SRC_DIR       where the tarballs are stored / downloaded to       (default: $PREFIX/src)
#   ONEAPI_ROOT   oneAPI installation, setvars.sh is sourced from here (default: /opt/intel/oneapi)
#   OPT           optimization flags                                   (default: -O3)
#   ARCH_FLAG     CPU flag; default -xHOST on Intel, -march=x86-64-v3 on AMD
#   BLACS_LIB     MKL BLACS library for your MPI                       (default: -lmkl_blacs_intelmpi_lp64)
#   RUN_TEST      1 = run "openmx -runtest" after the build            (default: 1)
#   NP            MPI processes used by the runtest                    (default: 4)
#   FORCE         1 = delete an existing $PREFIX/openmx4.0 first       (default: 0)
#
set -eo pipefail     # 'set -u' is not used because setvars.sh reads unset variables

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PREFIX=${PREFIX:-$HOME/software}
SRC_DIR=${SRC_DIR:-$PREFIX/src}
ONEAPI_ROOT=${ONEAPI_ROOT:-/opt/intel/oneapi}
OPT=${OPT:--O3}
BLACS_LIB=${BLACS_LIB:--lmkl_blacs_intelmpi_lp64}
RUN_TEST=${RUN_TEST:-1}
NP=${NP:-4}
FORCE=${FORCE:-0}

OMX_URL=https://www.openmx-square.org/openmx4.0.tar.gz
PATCH_URL=https://www.openmx-square.org/bugfixed/26May08/patch4.0.1.tar.gz
OMX_DIR=$PREFIX/openmx4.0

log() { echo "[build_openmx_oneapi] $*"; }
die() { echo "[build_openmx_oneapi] ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------- 1. oneAPI environment
if ! command -v mpiicx >/dev/null 2>&1; then
  [ -f "$ONEAPI_ROOT/setvars.sh" ] || die "mpiicx not found and $ONEAPI_ROOT/setvars.sh missing"
  log "sourcing $ONEAPI_ROOT/setvars.sh"
  # shellcheck disable=SC1091
  source "$ONEAPI_ROOT/setvars.sh" > /dev/null
fi
for c in mpiicx mpiifx make tar; do
  command -v "$c" >/dev/null 2>&1 || die "'$c' not found; check your oneAPI installation"
done
command -v gcc >/dev/null 2>&1 || log "warning: system gcc not found; 'make all' needs it for bandgnu13, cif2omx, ..."
[ -n "${MKLROOT:-}" ] && [ -d "$MKLROOT" ] || die "MKLROOT is not set; source setvars.sh or load the MKL module"
if [ -d "$MKLROOT/lib/intel64" ]; then MKL_LIBDIR='$(MKLROOT)/lib/intel64'; else MKL_LIBDIR='$(MKLROOT)/lib'; fi
log "mpiicx: $(mpiicx --version 2>/dev/null | head -n 1)"
log "mpiifx: $(mpiifx --version 2>/dev/null | head -n 1)"
log "MKLROOT=$MKLROOT"

# ---------------------------------------------------------------- 2. CPU flag
if [ -z "${ARCH_FLAG:-}" ]; then
  if grep -q AuthenticAMD /proc/cpuinfo 2>/dev/null; then
    ARCH_FLAG="-march=x86-64-v3"      # -xHOST is Intel specific
  else
    ARCH_FLAG="-xHOST"
  fi
fi
log "CPU flag: $ARCH_FLAG"

# ---------------------------------------------------------------- 3. download
mkdir -p "$SRC_DIR"
cd "$SRC_DIR"
[ -f openmx4.0.tar.gz ]  || { log "downloading $OMX_URL";   wget -c "$OMX_URL"; }
[ -f patch4.0.1.tar.gz ] || { log "downloading $PATCH_URL"; wget -c "$PATCH_URL"; }

# ---------------------------------------------------------------- 4. unpack and patch
if [ -d "$OMX_DIR" ]; then
  if [ "$FORCE" = "1" ]; then
    log "removing existing $OMX_DIR (FORCE=1)"; rm -rf "$OMX_DIR"
  else
    die "$OMX_DIR already exists; remove it or rerun with FORCE=1"
  fi
fi
mkdir -p "$PREFIX"
tar -xzf "$SRC_DIR/openmx4.0.tar.gz" -C "$PREFIX"
[ -d "$OMX_DIR/source" ] || die "unexpected tarball layout: $OMX_DIR/source not found"
cp "$SRC_DIR/patch4.0.1.tar.gz" "$OMX_DIR/source/"
cd "$OMX_DIR/source"
tar -xzf patch4.0.1.tar.gz
[ -f GaAs.dat ] && mv -f GaAs.dat ../work/
log "patch 4.0.1 applied"

# ---------------------------------------------------------------- 5. rewrite the makefile
cp makefile makefile.orig
cat > makefile.frag <<EOF
# ---- MKLROOT/CC/FC/LIB inserted by build_openmx_oneapi.sh (DFT-Tutorials) ----
MKLROOT = ${MKLROOT}
CC  = mpiicx ${OPT} ${ARCH_FLAG} -fcommon -Wno-implicit-function-declaration -qopenmp -I\$(MKLROOT)/include -I\$(MKLROOT)/include/fftw
FC  = mpiifx ${OPT} ${ARCH_FLAG} -qopenmp
LIB = -L${MKL_LIBDIR} -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lifcore ${BLACS_LIB} -liomp5 -lpthread -lm -ldl
# ---- end of inserted block ----
EOF
awk '
  function ends_bs(s) { return substr(s, length(s), 1) == "\\" }
  cont { print "#" $0; cont = ends_bs($0); next }
  /^[ \t]*(MKLROOT|CC|FC|LIB)[ \t]*:?=/ { print "#" $0; cont = ends_bs($0); next }
  { print }
' makefile.orig > makefile.body
cat makefile.frag makefile.body > makefile
rm -f makefile.body
log "active compiler settings now in makefile:"
grep -nE '^[[:space:]]*(MKLROOT|CC|FC|LIB)[[:space:]]*(\+|:)?=' makefile | sed 's/^/    /'

# ---------------------------------------------------------------- 6. compile
log "make clean / make install / make all (logs: make_install.log, make_all.log)"
make clean > /dev/null 2>&1 || true
make install 2>&1 | tee make_install.log
make all     2>&1 | tee make_all.log
[ -x ../work/openmx ] || die "../work/openmx was not created; inspect make_install.log"
if ldd ../work/openmx | grep -q "not found"; then
  ldd ../work/openmx | grep "not found"
  die "unresolved shared libraries (see above); is the oneAPI environment loaded?"
fi
for t in DosMain bandgnu13 cohp polB jx cif2omx; do
  [ -x "../work/$t" ] && log "built: work/$t" || log "warning: work/$t missing (check make_all.log)"
done

cat > "$OMX_DIR/INSTALL_NOTES.txt" <<EOF
OpenMX 4.0 + patch 4.0.1 built by build_openmx_oneapi.sh on $(date '+%F %T') ($(hostname))
$(mpiicx --version 2>/dev/null | head -n 1)
$(mpiifx --version 2>/dev/null | head -n 1)
Runtime: source $ONEAPI_ROOT/setvars.sh before running openmx.
$(cat makefile.frag)
EOF

# ---------------------------------------------------------------- 7. test
if [ "$RUN_TEST" = "1" ]; then
  cd "$OMX_DIR/work"
  export OMP_NUM_THREADS=1
  log "running: mpirun -np $NP ./openmx -runtest -nt 1 (a few minutes)"
  mpirun -np "$NP" ./openmx -runtest -nt 1 > runtest.log 2>&1 || log "warning: runtest returned non-zero"
  bash "$HERE/check_runtest.sh" runtest.result || log "runtest check FAILED (see runtest.result)"
fi

log "done. Before running OpenMX: source $ONEAPI_ROOT/setvars.sh ; export PATH=$OMX_DIR/work:\$PATH"
log "and use 'DATA.PATH  $OMX_DIR/DFT_DATA19' in your input files."
