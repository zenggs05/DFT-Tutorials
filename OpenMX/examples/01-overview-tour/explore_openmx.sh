#!/bin/bash
#
# Take a guided tour of an OpenMX installation (source tree or Debian package).
# Prints the directory layout, database statistics, which executables exist,
# which example inputs are available and the valence electrons of some VPS.
# Nothing is modified.
#
# DFT-Tutorials, OpenMX chapter 1 (01-overview.md).
#
# Usage:
#   bash explore_openmx.sh                      # default: $HOME/software/openmx4.0
#   bash explore_openmx.sh /opt/openmx4.0       # Debian package
#   bash explore_openmx.sh ~/openmx4.0 > tour.txt
#
set -u
ROOT=${1:-${OPENMX_ROOT:-$HOME/software/openmx4.0}}
ROOT=${ROOT%/}

hr() { printf '%s\n' "------------------------------------------------------------------"; }

if [ ! -d "$ROOT" ]; then
  echo "directory not found: $ROOT"
  echo "usage: bash explore_openmx.sh /path/to/openmx4.0"
  exit 1
fi
ROOT=$(cd "$ROOT" && pwd)

hr; echo "OpenMX tree: $ROOT"; hr
for d in source work DFT_DATA19 DFT_DATA19/PAO DFT_DATA19/VPS bin runtime/lib; do
  if [ -d "$ROOT/$d" ]; then
    printf '  %-18s %8s\n' "$d/" "$(du -sh "$ROOT/$d" 2>/dev/null | cut -f1)"
  else
    printf '  %-18s %8s\n' "$d/" "(absent)"
  fi
done

hr; echo "Version information"; hr
if [ -f "$ROOT/source/openmx.c" ]; then
  head -n 30 "$ROOT/source/openmx.c" | grep -i -m 2 "ver" | sed 's/^/  source\/openmx.c: /'
else
  echo "  source/openmx.c not available (binary installation?)"
fi
if [ -f "$ROOT/source/patch4.0.1.tar.gz" ]; then
  echo "  source/patch4.0.1.tar.gz present -> patch 4.0.1 was probably applied"
else
  echo "  no patch tarball in source/ (check whether patch 4.0.1 was applied)"
fi
[ -f "$ROOT/INSTALL_NOTES.txt" ] && { echo "  INSTALL_NOTES.txt:"; sed 's/^/    /' "$ROOT/INSTALL_NOTES.txt" | head -n 8; }

hr; echo "Database DFT_DATA19"; hr
DATA=$ROOT/DFT_DATA19
if [ -d "$DATA/PAO" ] && [ -d "$DATA/VPS" ]; then
  echo "  number of PAO files: $(ls "$DATA/PAO" | grep -c '\.pao$')"
  echo "  number of VPS files: $(ls "$DATA/VPS" | grep -c '\.vps$')"
  echo "  GGA (PBE19) VPS    : $(ls "$DATA/VPS" | grep -c '_PBE19.*\.vps$')"
  echo "  LDA (CA19)  VPS    : $(ls "$DATA/VPS" | grep -c '_CA19.*\.vps$')"
  echo "  valence electrons of a few VPS (needed for the up/down columns of atoms):"
  for v in H_PBE19 C_PBE19 O_PBE19 Si_PBE19 Fe_PBE19H Fe_PBE19S Ni_PBE19H Ni_PBE19S; do
    f="$DATA/VPS/$v.vps"
    [ -f "$f" ] && printf '    %-10s %s\n' "$v" "$(grep -m1 'valence.electron' "$f" | awk '{print $2}')"
  done
  echo "  PAO files available for carbon: $(cd "$DATA/PAO" && ls C[0-9]*.pao 2>/dev/null | tr '\n' ' ')"
  echo
  echo "  In your inputs write an ABSOLUTE path:"
  echo "    DATA.PATH   $DATA"
else
  echo "  $DATA/PAO or $DATA/VPS not found."
  echo "  Locate it with: find / -maxdepth 4 -type d -name DFT_DATA19 2>/dev/null"
fi

hr; echo "Executables"; hr
for t in openmx DosMain cohp bandgnu13 polB jx calB Z2FH kSpin cif2omx gcube2oned diff_gcube \
         add_gcube cube2xsf md2axsf bin2txt esp analysis_example; do
  loc=""
  for d in "$ROOT/bin" "$ROOT/work"; do
    [ -x "$d/$t" ] && { loc="$d/$t"; break; }
  done
  [ -z "$loc" ] && loc=$(command -v "$t" 2>/dev/null || true)
  printf '  %-17s %s\n' "$t" "${loc:-(not found: run make all in source/)}"
done

hr; echo "Example inputs in work/"; hr
if [ -d "$ROOT/work" ]; then
  n=$(ls "$ROOT/work"/*.dat 2>/dev/null | wc -l)
  echo "  $n *.dat files directly in work/, for example:"
  for f in Methane.dat H2O.dat Cdia.dat Si.dat GaAs.dat Febcc2.dat NiO-cFLL.dat NaCl.dat \
           Al111_WorkFunc_7E.dat C2H4_NEB.dat Sidia_COHP.dat DIA512_DC.dat; do
    if [ -f "$ROOT/work/$f" ]; then
      name=$(awk 'tolower($1)=="system.name"{print $2; exit}' "$ROOT/work/$f")
      solver=$(awk 'tolower($1)=="scf.eigenvaluesolver"{print $2; exit}' "$ROOT/work/$f")
      mdtype=$(awk 'tolower($1)=="md.type"{print $2; exit}' "$ROOT/work/$f")
      printf '    %-24s System.Name=%-12s solver=%-8s MD.Type=%s\n' "$f" "$name" "$solver" "$mdtype"
    fi
  done
  echo "  sub-directories of work/:"
  for d in "$ROOT"/work/*/; do [ -d "$d" ] && printf '    %s\n' "$(basename "$d")/"; done
else
  echo "  work/ not found"
fi

hr
echo "Next steps: chapter 2 (install/verify) and chapter 3 (first calculation)."
echo "Check an input file with: python list_dft_data.py check-dat your.dat"
