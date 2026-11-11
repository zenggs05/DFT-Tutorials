#!/bin/bash
#
# Build OpenMX 4.0 + patch 4.0.1 with GCC + Open MPI + ScaLAPACK + OpenBLAS + FFTW3.
# DFT-Tutorials, OpenMX chapter 2 (02-install.md, section 5).
#
# Usage:
#   bash build_openmx_gcc.sh                       # build into $HOME/software/openmx4.0
#   INSTALL_DEPS=1 bash build_openmx_gcc.sh        # also apt-get the dependencies (sudo)
#   PREFIX=/apps/openmx RUN_TEST=0 bash build_openmx_gcc.sh
#
# Environment variables (all optional):
#   PREFIX         install root, the code goes to $PREFIX/openmx4.0  (default: $HOME/software)
#   SRC_DIR        where the tarballs are stored / downloaded to      (default: $PREFIX/src)
#   INSTALL_DEPS   1 = apt-get install the build dependencies          (default: 0)
#   OPT            optimization flags for CC and FC                     (default: -O3)
#   SCALAPACK_LIB  ScaLAPACK link flags                                 (default: -lscalapack-openmpi)
#   BLAS_LIB       BLAS/LAPACK link flags                               (default: -lopenblas)
#   EXTRA_LIB      extra link flags, e.g. -L/opt/scalapack/lib          (default: empty)
#   RUN_TEST       1 = run "openmx -runtest" after the build            (default: 1)
#   NP             MPI processes used by the runtest                    (default: 4)
#   FORCE          1 = delete an existing $PREFIX/openmx4.0 first       (default: 0)
#
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PREFIX=${PREFIX:-$HOME/software}
SRC_DIR=${SRC_DIR:-$PREFIX/src}
INSTALL_DEPS=${INSTALL_DEPS:-0}
OPT=${OPT:--O3}
SCALAPACK_LIB=${SCALAPACK_LIB:--lscalapack-openmpi}
BLAS_LIB=${BLAS_LIB:--lopenblas}
EXTRA_LIB=${EXTRA_LIB:-}
RUN_TEST=${RUN_TEST:-1}
NP=${NP:-4}
FORCE=${FORCE:-0}

OMX_URL=https://www.openmx-square.org/openmx4.0.tar.gz
PATCH_URL=https://www.openmx-square.org/bugfixed/26May08/patch4.0.1.tar.gz
OMX_DIR=$PREFIX/openmx4.0

log() { echo "[build_openmx_gcc] $*"; }
die() { echo "[build_openmx_gcc] ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------- 1. dependencies
if [ "$INSTALL_DEPS" = "1" ]; then
  log "installing dependencies with apt-get (sudo required)"
  sudo apt-get update
  sudo apt-get install -y build-essential gfortran openmpi-bin libopenmpi-dev \
       libscalapack-openmpi-dev libopenblas-dev libfftw3-dev wget
fi
for c in mpicc mpif90 gcc make tar; do
  command -v "$c" >/dev/null 2>&1 || die "'$c' not found (install the dependencies or load modules)"
done
log "mpicc  -> $(mpicc --showme:command 2>/dev/null || echo unknown) $(mpicc -dumpversion)"
log "mpif90 -> $(mpif90 --showme:command 2>/dev/null || echo unknown) $(mpif90 -dumpversion)"

# ---------------------------------------------------------------- 2. compiler-dependent flags
GCC_MAJOR=$(mpicc -dumpversion | cut -d. -f1)
EXTRA_CFLAGS=""
if [ "$GCC_MAJOR" -ge 14 ] 2>/dev/null; then
  EXTRA_CFLAGS="$EXTRA_CFLAGS -Wno-error=implicit-int -Wno-error=incompatible-pointer-types"
fi
if [ "$GCC_MAJOR" -ge 15 ] 2>/dev/null; then
  EXTRA_CFLAGS="$EXTRA_CFLAGS -std=gnu17"
fi
case "$(uname -m)" in
  aarch64|arm64) EXTRA_CFLAGS="$EXTRA_CFLAGS -Dnosse" ;;
esac
log "GCC major version: $GCC_MAJOR ; extra CFLAGS:${EXTRA_CFLAGS:- none}"

# ---------------------------------------------------------------- 3. download
mkdir -p "$SRC_DIR"
cd "$SRC_DIR"
[ -f openmx4.0.tar.gz ]   || { log "downloading $OMX_URL";   wget -c "$OMX_URL"; }
[ -f patch4.0.1.tar.gz ]  || { log "downloading $PATCH_URL"; wget -c "$PATCH_URL"; }

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
# ---- CC/FC/LIB inserted by build_openmx_gcc.sh (DFT-Tutorials) ----
CC  = mpicc  ${OPT} -fopenmp -fcommon -Wno-implicit-function-declaration${EXTRA_CFLAGS}
FC  = mpif90 ${OPT} -fopenmp -fallow-argument-mismatch
LIB = ${SCALAPACK_LIB} ${BLAS_LIB} -lfftw3 ${EXTRA_LIB} -lmpi_mpifh -lgfortran -lpthread -lm -ldl
# ---- end of inserted block ----
EOF
# comment out the active MKLROOT/CC/FC/LIB definitions (incl. backslash continuations);
# "CC += ..." and "FC += ..." (bundled ELPA include path) are left untouched.
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
  die "unresolved shared libraries (see above)"
fi
for t in DosMain bandgnu13 cohp polB jx cif2omx; do
  [ -x "../work/$t" ] && log "built: work/$t" || log "warning: work/$t missing (check make_all.log)"
done

cat > "$OMX_DIR/INSTALL_NOTES.txt" <<EOF
OpenMX 4.0 + patch 4.0.1 built by build_openmx_gcc.sh on $(date '+%F %T') ($(hostname))
mpicc : $(mpicc --version | head -n 1)
mpif90: $(mpif90 --version | head -n 1)
$(cat makefile.frag)
EOF

# ---------------------------------------------------------------- 7. test
if [ "$RUN_TEST" = "1" ]; then
  cd "$OMX_DIR/work"
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  log "running: mpirun -np $NP ./openmx -runtest -nt 1 (a few minutes)"
  mpirun -np "$NP" ./openmx -runtest -nt 1 > runtest.log 2>&1 || log "warning: runtest returned non-zero"
  bash "$HERE/check_runtest.sh" runtest.result || log "runtest check FAILED (see runtest.result)"
fi

log "done. Add to your ~/.bashrc:  export PATH=$OMX_DIR/work:\$PATH"
log "and use 'DATA.PATH  $OMX_DIR/DFT_DATA19' in your input files."
