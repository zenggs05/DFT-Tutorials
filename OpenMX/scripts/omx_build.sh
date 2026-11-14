#!/usr/bin/env bash
#
# omx_build.sh -- download, patch, configure and build OpenMX 4.0 (+ patch 4.0.1),
#                 then run "openmx -runtest" and summarise diff Utot / diff Force.
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
#
# What it does
#   1. downloads openmx4.0.tar.gz and patch4.0.1.tar.gz (skipped if present / --offline)
#   2. unpacks into PREFIX/openmx4.0, applies the patch, moves the patched GaAs.dat to work/
#   3. comments out the active MKLROOT/CC/FC/LIB definitions of source/makefile
#      (including backslash continuation lines) and inserts new ones for the
#      chosen toolchain; "CC += ..." / "FC += ..." lines (bundled ELPA) are kept
#        gcc    : mpicc/mpif90 + ScaLAPACK + OpenBLAS + FFTW3
#                 -fcommon (GCC >= 10), -Wno-implicit-function-declaration,
#                 -fallow-argument-mismatch (gfortran >= 10),
#                 -Wno-error=implicit-int -Wno-error=incompatible-pointer-types (GCC >= 14),
#                 -std=gnu17 (GCC >= 15), -Dnosse on aarch64;
#                 Open MPI -> -lmpi_mpifh, MPICH -> -lmpifort
#        oneapi : mpiicx/mpiifx + MKL (ScaLAPACK, BLACS, FFTW3 interface) + Intel MPI;
#                 -xHOST on Intel CPUs, -march=x86-64-v3 on AMD CPUs
#   4. make clean; make install; make all  (logs: source/make_install.log, make_all.log)
#   5. builds bandgnu13 with "gcc bandgnu13.c -lm" if make all did not produce it
#   6. cd work; mpirun -np NP ./openmx -runtest -nt 1; summarises runtest.result
#      (installation is regarded as correct if all differences are below ~1e-7,
#       see manual Sec. 6)
#
# Usage
#   bash omx_build.sh -h
#   bash omx_build.sh --toolchain gcc
#   bash omx_build.sh --toolchain oneapi --prefix /apps/openmx --np 8
#   bash omx_build.sh --toolchain gcc --scalapack-lib "-L/opt/scalapack/lib -lscalapack" --no-test
#   bash omx_build.sh --toolchain gcc --offline --src-dir ~/tarballs --force
#
# Dependencies
#   bash, wget (or curl), tar, make, system gcc (for the small tools of "make all"), and
#   gcc  : mpicc/mpif90 (Open MPI or MPICH), ScaLAPACK, OpenBLAS (or other BLAS/LAPACK), FFTW3
#          (Ubuntu 24.04: build-essential gfortran openmpi-bin libopenmpi-dev
#           libscalapack-openmpi-dev libopenblas-dev libfftw3-dev)
#   oneapi: Intel oneAPI 2024+ (icx/ifx, Intel MPI, MKL); setvars.sh is sourced if needed
#
# The result is written to PREFIX/openmx4.0/work (openmx, DosMain, bandgnu13, ...).
# Afterwards: export PATH=PREFIX/openmx4.0/work:$PATH and use
#   DATA.PATH  PREFIX/openmx4.0/DFT_DATA19   in your input files.
#
set -euo pipefail

TOOLCHAIN=gcc
PREFIX=$HOME/software
SRC_DIR=""
NP=4
RUN_TEST=1
FORCE=0
OFFLINE=0
OPT="-O3"
ARCH_FLAG=""
SCALAPACK_LIB=""
BLAS_LIB="-lopenblas"
FFTW_LIB="-lfftw3"
EXTRA_LIB=""
BLACS_LIB="-lmkl_blacs_intelmpi_lp64"
ONEAPI_ROOT=${ONEAPI_ROOT:-/opt/intel/oneapi}
MPIRUN=${MPIRUN:-mpirun}

OMX_URL=https://www.openmx-square.org/openmx4.0.tar.gz
PATCH_URL=https://www.openmx-square.org/bugfixed/26May08/patch4.0.1.tar.gz

usage() {
  cat <<EOF
Usage: bash omx_build.sh [options]

Build OpenMX 4.0 + patch 4.0.1 from source and check it with -runtest.

Options:
  --toolchain gcc|oneapi   compiler/library stack (default: gcc)
  --prefix DIR             install root; code goes to DIR/openmx4.0 (default: \$HOME/software)
  --src-dir DIR            where the tarballs are stored/downloaded (default: PREFIX/src)
  --offline                do not download; tarballs must already be in --src-dir
  --force                  delete an existing PREFIX/openmx4.0 first
  --opt "FLAGS"            optimisation flags (default: -O3; avoid -ffast-math)
  --arch-flag FLAG         oneapi: CPU flag (default: -xHOST on Intel, -march=x86-64-v3 on AMD)
                           gcc: optional, e.g. -march=native (only if build and compute CPUs match)
  --scalapack-lib "LIBS"   gcc: ScaLAPACK link flags (default: -lscalapack-openmpi if that
                           library is found by ldconfig, else -lscalapack)
  --blas-lib "LIBS"        gcc: BLAS/LAPACK link flags (default: -lopenblas)
  --fftw-lib "LIBS"        gcc: FFTW3 link flags (default: -lfftw3)
  --extra-lib "LIBS"       extra link flags, e.g. -L/opt/fftw/lib
  --blacs-lib LIB          oneapi: MKL BLACS library (default: -lmkl_blacs_intelmpi_lp64;
                           use -lmkl_blacs_openmpi_lp64 with Open MPI)
  --oneapi-root DIR        oneapi: installation with setvars.sh (default: /opt/intel/oneapi)
  --np N                   MPI processes for -runtest (default: 4)
  --no-test                skip "openmx -runtest"
  -h, --help               show this help

Environment: MPIRUN (default mpirun) is used for the runtest.
EOF
}

log() { echo "[omx_build] $*"; }
die() { echo "[omx_build] ERROR: $*" >&2; exit 1; }

while [ $# -gt 0 ]; do
  case "$1" in
    --toolchain)      TOOLCHAIN=${2:?}; shift 2 ;;
    --prefix)         PREFIX=${2:?}; shift 2 ;;
    --src-dir)        SRC_DIR=${2:?}; shift 2 ;;
    --offline)        OFFLINE=1; shift ;;
    --force)          FORCE=1; shift ;;
    --opt)            OPT=${2:?}; shift 2 ;;
    --arch-flag)      ARCH_FLAG=${2:?}; shift 2 ;;
    --scalapack-lib)  SCALAPACK_LIB=${2:?}; shift 2 ;;
    --blas-lib)       BLAS_LIB=${2:?}; shift 2 ;;
    --fftw-lib)       FFTW_LIB=${2:?}; shift 2 ;;
    --extra-lib)      EXTRA_LIB=${2:?}; shift 2 ;;
    --blacs-lib)      BLACS_LIB=${2:?}; shift 2 ;;
    --oneapi-root)    ONEAPI_ROOT=${2:?}; shift 2 ;;
    --np)             NP=${2:?}; shift 2 ;;
    --no-test)        RUN_TEST=0; shift ;;
    -h|--help)        usage; exit 0 ;;
    *)                usage >&2; die "unknown option: $1" ;;
  esac
done

case "$TOOLCHAIN" in gcc|oneapi) ;; *) die "--toolchain must be gcc or oneapi" ;; esac
[[ "$NP" =~ ^[0-9]+$ ]] || die "--np must be a positive integer"
PREFIX=$(mkdir -p "$PREFIX" && cd "$PREFIX" && pwd)
SRC_DIR=${SRC_DIR:-$PREFIX/src}
OMX_DIR=$PREFIX/openmx4.0

# --------------------------------------------------------------------------- 1. toolchain
EXTRA_CFLAGS=""
EXTRA_FFLAGS=""
if [ "$TOOLCHAIN" = gcc ]; then
  for c in mpicc mpif90 gcc make tar; do
    command -v "$c" >/dev/null 2>&1 || die "'$c' not found (install the dependencies or 'module load' them)"
  done
  GCC_MAJOR=$(mpicc -dumpversion 2>/dev/null | cut -d. -f1)
  GFC_MAJOR=$(mpif90 -dumpversion 2>/dev/null | cut -d. -f1)
  [[ "$GCC_MAJOR" =~ ^[0-9]+$ ]] || die "could not get the compiler version from 'mpicc -dumpversion'"
  [[ "$GFC_MAJOR" =~ ^[0-9]+$ ]] || GFC_MAJOR=$GCC_MAJOR
  log "mpicc -> gcc $GCC_MAJOR, mpif90 -> gfortran $GFC_MAJOR"
  [ "$GCC_MAJOR" -ge 10 ] && EXTRA_CFLAGS="$EXTRA_CFLAGS -fcommon"
  EXTRA_CFLAGS="$EXTRA_CFLAGS -Wno-implicit-function-declaration"
  [ "$GCC_MAJOR" -ge 14 ] && EXTRA_CFLAGS="$EXTRA_CFLAGS -Wno-error=implicit-int -Wno-error=incompatible-pointer-types"
  [ "$GCC_MAJOR" -ge 15 ] && EXTRA_CFLAGS="$EXTRA_CFLAGS -std=gnu17"
  [ "$GFC_MAJOR" -ge 10 ] && EXTRA_FFLAGS="$EXTRA_FFLAGS -fallow-argument-mismatch"
  case "$(uname -m)" in aarch64|arm64) EXTRA_CFLAGS="$EXTRA_CFLAGS -Dnosse" ;; esac
  # Fortran MPI bindings: Open MPI -> libmpi_mpifh, MPICH/Intel MPI -> libmpifort
  if mpif90 --showme:link >/dev/null 2>&1; then
    MPI_FLIB="-lmpi_mpifh"; log "MPI flavour: Open MPI"
  else
    MPI_FLIB="-lmpifort"; log "MPI flavour: MPICH-like"
  fi
  if [ -z "$SCALAPACK_LIB" ]; then
    if ldconfig -p 2>/dev/null | grep "libscalapack-openmpi" >/dev/null; then
      SCALAPACK_LIB="-lscalapack-openmpi"
    else
      SCALAPACK_LIB="-lscalapack"
    fi
  fi
  CC_LINE="mpicc ${OPT} ${ARCH_FLAG} -fopenmp${EXTRA_CFLAGS}"
  FC_LINE="mpif90 ${OPT} ${ARCH_FLAG} -fopenmp${EXTRA_FFLAGS}"
  LIB_LINE="${SCALAPACK_LIB} ${BLAS_LIB} ${FFTW_LIB} ${EXTRA_LIB} ${MPI_FLIB} -lgfortran -lpthread -lm -ldl"
  MKL_LINE=""
else
  set +u                                   # setvars.sh reads unset variables
  if ! command -v mpiicx >/dev/null 2>&1; then
    [ -f "$ONEAPI_ROOT/setvars.sh" ] || die "mpiicx not found and $ONEAPI_ROOT/setvars.sh missing"
    log "sourcing $ONEAPI_ROOT/setvars.sh"
    # shellcheck disable=SC1091
    source "$ONEAPI_ROOT/setvars.sh" >/dev/null
  fi
  set -u
  for c in mpiicx mpiifx make tar; do
    command -v "$c" >/dev/null 2>&1 || die "'$c' not found; check the oneAPI installation"
  done
  command -v gcc >/dev/null 2>&1 || log "warning: system gcc not found; 'make all' needs it for bandgnu13 etc."
  [ -n "${MKLROOT:-}" ] && [ -d "$MKLROOT" ] || die "MKLROOT is not set; source setvars.sh or load the MKL module"
  if [ -z "$ARCH_FLAG" ]; then
    if grep -q AuthenticAMD /proc/cpuinfo 2>/dev/null; then ARCH_FLAG="-march=x86-64-v3"; else ARCH_FLAG="-xHOST"; fi
  fi
  if [ -d "$MKLROOT/lib/intel64" ]; then MKL_LIBDIR='$(MKLROOT)/lib/intel64'; else MKL_LIBDIR='$(MKLROOT)/lib'; fi
  log "mpiicx: $(mpiicx --version 2>/dev/null | head -n 1); CPU flag $ARCH_FLAG; MKLROOT=$MKLROOT"
  MKL_LINE="MKLROOT = ${MKLROOT}"
  CC_LINE="mpiicx ${OPT} ${ARCH_FLAG} -fcommon -Wno-implicit-function-declaration -qopenmp -I\$(MKLROOT)/include -I\$(MKLROOT)/include/fftw"
  FC_LINE="mpiifx ${OPT} ${ARCH_FLAG} -qopenmp"
  LIB_LINE="-L${MKL_LIBDIR} -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lifcore ${BLACS_LIB} -liomp5 ${EXTRA_LIB} -lpthread -lm -ldl"
fi

# --------------------------------------------------------------------------- 2. download
fetch() {
  local url=$1 out=$2
  if [ -f "$out" ]; then log "found $out"; return 0; fi
  [ "$OFFLINE" = 1 ] && die "--offline: $out not found"
  log "downloading $url"
  if command -v wget >/dev/null 2>&1; then
    wget -c -O "$out.part" "$url"
  elif command -v curl >/dev/null 2>&1; then
    curl -fL -o "$out.part" "$url"
  else
    die "neither wget nor curl is available"
  fi
  mv "$out.part" "$out"
}
mkdir -p "$SRC_DIR"
fetch "$OMX_URL"   "$SRC_DIR/openmx4.0.tar.gz"
fetch "$PATCH_URL" "$SRC_DIR/patch4.0.1.tar.gz"
tar -tzf "$SRC_DIR/openmx4.0.tar.gz" >/dev/null 2>&1 || die "$SRC_DIR/openmx4.0.tar.gz is not a valid tar.gz"
tar -tzf "$SRC_DIR/patch4.0.1.tar.gz" >/dev/null 2>&1 || die "$SRC_DIR/patch4.0.1.tar.gz is not a valid tar.gz"

# --------------------------------------------------------------------------- 3. unpack + patch
if [ -d "$OMX_DIR" ]; then
  if [ "$FORCE" = 1 ]; then log "removing $OMX_DIR (--force)"; rm -rf "$OMX_DIR"
  else die "$OMX_DIR exists; remove it or use --force"; fi
fi
tar -xzf "$SRC_DIR/openmx4.0.tar.gz" -C "$PREFIX"
[ -d "$OMX_DIR/source" ] || die "unexpected tarball layout: $OMX_DIR/source not found"
cp "$SRC_DIR/patch4.0.1.tar.gz" "$OMX_DIR/source/"
cd "$OMX_DIR/source"
tar -xzf patch4.0.1.tar.gz
if [ -f GaAs.dat ]; then mv -f GaAs.dat ../work/; fi
log "patch 4.0.1 applied"

# --------------------------------------------------------------------------- 4. makefile
cp makefile makefile.orig
{
  echo "# ---- inserted by omx_build.sh (DFT-Tutorials), toolchain: $TOOLCHAIN ----"
  [ -n "$MKL_LINE" ] && echo "$MKL_LINE"
  echo "CC  = $CC_LINE"
  echo "FC  = $FC_LINE"
  echo "LIB = $LIB_LINE"
  echo "# ---- end of inserted block ----"
} > makefile.frag
awk '
  function ends_bs(s) { return substr(s, length(s), 1) == "\\" }
  cont { print "#" $0; cont = ends_bs($0); next }
  /^[ \t]*(MKLROOT|CC|FC|LIB)[ \t]*:?=/ { print "#" $0; cont = ends_bs($0); next }
  { print }
' makefile.orig > makefile.body
cat makefile.frag makefile.body > makefile
rm -f makefile.body
log "active definitions in the new makefile:"
grep -nE '^[[:space:]]*(MKLROOT|CC|FC|LIB)[[:space:]]*(\+|:)?=' makefile | sed 's/^/    /'

# --------------------------------------------------------------------------- 5. compile
log "make clean / make install / make all"
make clean >/dev/null 2>&1 || true
make install 2>&1 | tee make_install.log || die "make install failed; see $OMX_DIR/source/make_install.log"
[ -x ../work/openmx ] || die "../work/openmx was not built; see $OMX_DIR/source/make_install.log"
make all 2>&1 | tee make_all.log || log "warning: 'make all' returned an error (see make_all.log)"
if [ ! -x ../work/bandgnu13 ] && [ -f bandgnu13.c ]; then
  log "building bandgnu13 separately: gcc bandgnu13.c -lm -o bandgnu13"
  if gcc bandgnu13.c -lm -o bandgnu13; then cp bandgnu13 ../work/; else log "warning: bandgnu13 failed"; fi
fi
if command -v ldd >/dev/null 2>&1 && ldd ../work/openmx 2>&1 | grep "not found" >/dev/null; then
  ldd ../work/openmx | grep "not found"
  die "unresolved shared libraries (see above)"
fi
for t in openmx DosMain bandgnu13 cohp polB jx cif2omx; do
  if [ -x "../work/$t" ]; then log "built: work/$t"; else log "warning: work/$t missing (check make_all.log)"; fi
done
{
  echo "OpenMX 4.0 + patch 4.0.1 built by omx_build.sh on $(date '+%F %T') on $(hostname)"
  echo "toolchain: $TOOLCHAIN"
  cat makefile.frag
} > "$OMX_DIR/INSTALL_NOTES.txt"

# --------------------------------------------------------------------------- 6. runtest
if [ "$RUN_TEST" = 1 ]; then
  cd "$OMX_DIR/work"
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  ulimit -s unlimited 2>/dev/null || true
  log "running: $MPIRUN -np $NP ./openmx -runtest -nt 1 (a few minutes)"
  # MPIRUN may contain options, so it is intentionally unquoted
  # shellcheck disable=SC2086
  $MPIRUN -np "$NP" ./openmx -runtest -nt 1 > runtest.log 2>&1 || log "warning: runtest returned non-zero"
  if [ -f runtest.result ]; then
    awk -v thr=1e-7 '
      function a(x) { return x < 0 ? -x : x }
      function flush() {
        if (name != "" && u != "" && f != "") {
          n++; bad = (a(u+0) >= thr+0 || a(f+0) >= thr+0); nbad += bad
          printf "%3d  %-36s %18s %18s  %s\n", n, name, u, f, bad ? "FAIL" : "ok"
          name = ""; u = ""; f = ""
        }
      }
      BEGIN { printf "%3s  %-36s %18s %18s  %s\n", "#", "input", "diff Utot", "diff Force", "status" }
      {
        for (i = 1; i <= NF; i++) {
          if ($i ~ /\.dat$/) name = $i
          else if ($i == "Utot=") u = $(i+1)
          else if ($i ~ /^Utot=./) u = substr($i, 6)
          else if ($i == "Force=") { f = $(i+1); flush() }
          else if ($i ~ /^Force=./) { f = substr($i, 7); flush() }
        }
      }
      END {
        if (n == 0) { print "no test records found in runtest.result"; exit 1 }
        printf "tests: %d, failed (>= %s): %d\n", n, thr, nbad + 0
        if (nbad > 0) { print "RUNTEST: FAIL"; exit 1 } else print "RUNTEST: PASS"
      }' runtest.result || log "runtest check FAILED; compare work/runtest.result with work/input_example/runtest.result*"
  else
    log "warning: runtest.result not written; see $OMX_DIR/work/runtest.log"
  fi
fi

log "done."
log "add to ~/.bashrc:   export PATH=$OMX_DIR/work:\$PATH"
log "in input files:     DATA.PATH  $OMX_DIR/DFT_DATA19"
[ "$TOOLCHAIN" = oneapi ] && log "before running:     source $ONEAPI_ROOT/setvars.sh"
exit 0
