#!/bin/bash
# Vacuum -> implicit water -> (optional) electrolyte, each in its own directory.
# Run inside an allocation (e.g. from a Slurm script) with a VASPsol-patched binary.
#
#   export VASP_PP_PATH=/path/to/paw ; bash ../10-co-pt111-adsorption/make_potcar.sh .
#   bash run_vaspsol.sh            # vacuum + solvent
#   bash run_vaspsol.sh lpb        # also the linearised Poisson-Boltzmann step
set -e
MPI=${MPI:-"mpirun -np ${SLURM_NTASKS:-32}"}
VASP_EXE=${VASP_EXE:-vasp_std}
for f in POSCAR KPOINTS POTCAR INCAR_vac INCAR_sol; do
    [ -f "$f" ] || { echo "missing $f"; exit 1; }
done

run_step () {   # $1 = directory, $2 = INCAR file, $3 = directory with WAVECAR (or "")
    mkdir -p "$1"
    cp POSCAR KPOINTS POTCAR "$1"/
    cp "$2" "$1"/INCAR
    [ -n "$3" ] && cp "$3"/WAVECAR "$1"/
    (cd "$1" && $MPI $VASP_EXE > vasp.out 2>&1)
    echo "$1: $(grep 'energy(sigma->0)' "$1"/OUTCAR | tail -1)"
}

run_step vac INCAR_vac ""
run_step sol INCAR_sol vac
if [ "$1" = "lpb" ]; then
    run_step lpb INCAR_lpb sol
    grep -i "FERMI_SHIFT" lpb/vasp.out | tail -1 || echo "no FERMI_SHIFT line found in lpb/vasp.out"
fi

ev=$(grep 'energy(sigma->0)' vac/OUTCAR | tail -1 | awk '{print $NF}')
es=$(grep 'energy(sigma->0)' sol/OUTCAR | tail -1 | awk '{print $NF}')
awk -v a="$es" -v b="$ev" 'BEGIN{printf "E(solvent) - E(vacuum) = %.3f eV\n", a-b}'
# An unpatched VASP does not apply the solvent: both energies are then (almost) identical.
echo "solvation-related lines found in sol/ (should not be empty for a VASPsol build):"
grep -i -m 5 -E "solv|EB_K|TAU|Ediel" sol/vasp.out sol/OUTCAR || echo "   none - check that the binary contains VASPsol"
