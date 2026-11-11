#!/bin/bash
# Create the NEB directories 00 .. NN from two relaxed endpoint runs.
#
#   bash make_images.sh IS FS 3            # IDPP path via ASE (default)
#   bash make_images.sh IS FS 3 linear     # linear path (nebmake.pl if available)
#
# IS/ and FS/ must contain the finished relaxations (CONTCAR and OUTCAR).
# Run it in the NEB root directory (where INCAR with IMAGES = N will live).
set -e
IS=${1:-IS}; FS=${2:-FS}; N=${3:-3}; METHOD=${4:-idpp}
HERE=$(cd "$(dirname "$0")" && pwd)
NN=$(printf "%02d" $((N + 1)))
PY=${PYTHON:-$(command -v python3 || command -v python)}

for d in "$IS" "$FS"; do
    [ -f "$d/CONTCAR" ] && [ -f "$d/OUTCAR" ] || { echo "ERROR: $d needs CONTCAR and OUTCAR"; exit 1; }
    grep -q "reached required accuracy" "$d/OUTCAR" || echo "WARNING: $d is not converged"
done

if [ "$METHOD" = "linear" ] && command -v nebmake.pl > /dev/null; then
    nebmake.pl "$IS/CONTCAR" "$FS/CONTCAR" "$N"
elif [ "$METHOD" = "linear" ]; then
    "$PY" "$HERE/idpp_interpolate.py" "$IS/CONTCAR" "$FS/CONTCAR" -n "$N" --method linear
else
    "$PY" "$HERE/idpp_interpolate.py" "$IS/CONTCAR" "$FS/CONTCAR" -n "$N" --method idpp
fi

# nebresults.pl / nebef.pl / nebbarrier.pl need the endpoint OUTCARs
cp "$IS/OUTCAR" 00/OUTCAR
cp "$FS/OUTCAR" "$NN"/OUTCAR

if command -v dist.pl > /dev/null; then
    echo "dist.pl IS -> FS: $(dist.pl 00/POSCAR "$NN"/POSCAR)"
    for i in $(seq 0 "$N"); do
        a=$(printf "%02d" "$i"); b=$(printf "%02d" $((i + 1)))
        echo "  $a -> $b : $(dist.pl "$a"/POSCAR "$b"/POSCAR)"
    done
fi
echo "images ready: 00 .. $NN  (set IMAGES = $N in INCAR)"
