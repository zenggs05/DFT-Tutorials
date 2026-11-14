#!/usr/bin/env bash
# =============================================================================================
# batch_cube_multiwfn.sh - batch formchk + cubegen (HOMO/LUMO/density/ESP/spin) for every Gaussian
# checkpoint in a directory, plus optional non-interactive Multiwfn runs from command-stream files.
#
# For every *.chk (or *.fchk) found:
#   1. formchk  chk -> fchk (skipped when the fchk is newer, unless -f);
#   2. HOMO/LUMO energies (eV) from the fchk -> one line in batch_cube_summary.tsv;
#   3. cubegen for the requested kinds.  Open-shell fchk (Beta MOs present): orbitals are written
#      for alpha AND beta (AMO=n / BMO=n), so "homo" never puts two orbitals into one cube;
#      esp is computed on exactly the grid of the density cube (cubegen ... -1 h dens.cub);
#   4. each Multiwfn command stream (--mwfn FILE[:TAG]) is fed to "Multiwfn file.fchk -silent" inside
#      its own folder <stem>_<TAG>/, so fixed output names (surfanalysis.txt, vtx.pdb, *.cub ...)
#      never overwrite each other.
#
# Usage examples
#   bash batch_cube_multiwfn.sh                                  # *.chk in ., homo,lumo,density
#   bash batch_cube_multiwfn.sh -k homo-1,homo,lumo,lumo+1 -n -4 calc/
#   bash batch_cube_multiwfn.sh -k density,esp -p 16 phenol.chk  # ESP on the density grid
#   bash batch_cube_multiwfn.sh -k spin radical.fchk
#   bash batch_cube_multiwfn.sh --no-cube --mwfn esp_cmds.txt:esp --mwfn adch_cmds.txt:adch runs/
#   bash batch_cube_multiwfn.sh -r -o cubes/ project/            # recursive, cubes into one folder
#
# Multiwfn command streams: one answer per line, exactly what you type interactively (an empty line =
# Enter).  Run the sequence once by hand first and re-check after every Multiwfn update; examples are
# in Gaussian/examples/16-postprocess-multiwfn/ (multiwfn_esp_surface_cmds.txt, multiwfn_igmh_dimer_cmds.txt).
# A final "end-of-file during read" message in the .out file just means the stream ended.
#
# Dependencies: bash, Gaussian 16 utilities formchk and cubegen (source g16.profile first),
#               Multiwfn (optional, only for --mwfn; settings.ini is read from $Multiwfnpath).
#
# Part of DFT-Tutorials (https://github.com/zenggs05/DFT-Tutorials) - MIT License
# =============================================================================================
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bash batch_cube_multiwfn.sh [options] [DIR | file.chk | file.fchk ...]   (default: .)

Cubes
  -k, --kinds LIST     comma list (default homo,lumo,density):
                         homo, lumo, homo-N, lumo+N, mo=N   orbitals (alpha+beta if open shell)
                         density                           Density=SCF
                         spin                              Spin=SCF (open shell only)
                         esp                               Potential=SCF on the density grid
                         KIND=VALUE                        any other cubegen kind, e.g. Density=Current
  -n, --npts N         cubegen npts: 0 = 80^3, -2/-3/-4 = 3/6/12 points per Bohr (default -3)
  -p, --nproc N        cores for the ESP cube (default $SLURM_CPUS_PER_TASK or 4)
  -M, --mem SIZE       memory for formchk/cubegen, e.g. 2gb (passed as -m=SIZE, G16 C.01 or later)
      --no-cube        only formchk (+ Multiwfn)
Multiwfn
  -w, --mwfn FILE[:TAG]  run this command stream on every fchk (repeatable); TAG names the
                         output folder <stem>_TAG/ (default: file stem)
      --multiwfn PATH    Multiwfn executable (default: $MULTIWFN or Multiwfn in PATH)
Files
  -o, --outdir DIR     write cubes/results here (default: next to each fchk)
  -r, --recursive      search directories recursively
  -f, --force          redo formchk and overwrite existing cubes
      --dry-run        print the commands only
  -h, --help           this help
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }

KINDS="homo,lumo,density"
NPTS=-3
NPROC="${SLURM_CPUS_PER_TASK:-4}"
MEMOPT=""
NOCUBE=0
MWFN_CMDS=()
MWFN_BIN="${MULTIWFN:-Multiwfn}"
OUTDIR=""
RECURSIVE=0
FORCE=0
DRY=0
TARGETS=()
while [ $# -gt 0 ]; do
    case "$1" in
        -k|--kinds) KINDS="${2:?missing value for $1}"; shift 2 ;;
        -n|--npts) NPTS="${2:?missing value for $1}"; shift 2 ;;
        -p|--nproc) NPROC="${2:?missing value for $1}"; shift 2 ;;
        -M|--mem) MEMOPT="-m=${2:?missing value for $1}"; shift 2 ;;
        --no-cube) NOCUBE=1; shift ;;
        -w|--mwfn) MWFN_CMDS+=("${2:?missing value for $1}"); shift 2 ;;
        --multiwfn) MWFN_BIN="${2:?missing value for $1}"; shift 2 ;;
        -o|--outdir) OUTDIR="${2:?missing value for $1}"; shift 2 ;;
        -r|--recursive) RECURSIVE=1; shift ;;
        -f|--force) FORCE=1; shift ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) usage; exit 0 ;;
        --) shift; TARGETS+=("$@"); break ;;
        -*) usage >&2; die "unknown option $1" ;;
        *) TARGETS+=("$1"); shift ;;
    esac
done
[ ${#TARGETS[@]} -gt 0 ] || TARGETS=(.)
[[ "$NPTS" =~ ^-?[0-9]+$ ]] || die "--npts must be an integer"
[[ "$NPROC" =~ ^[0-9]+$ ]] || die "--nproc must be a non-negative integer"

run() {   # print and (unless --dry-run) execute
    echo "  + $*"
    [ "$DRY" -eq 1 ] || "$@"
}

if [ "$DRY" -eq 0 ]; then
    command -v formchk >/dev/null 2>&1 || die "formchk not found: source \$g16root/g16/bsd/g16.profile first"
    if [ "$NOCUBE" -eq 0 ]; then
        command -v cubegen >/dev/null 2>&1 || die "cubegen not found: source \$g16root/g16/bsd/g16.profile first"
    fi
fi

# check Multiwfn command files once
MW_FILES=() MW_TAGS=()
for spec in "${MWFN_CMDS[@]+"${MWFN_CMDS[@]}"}"; do
    f="${spec%%:*}"
    tag="${spec#*:}"; [ "$tag" = "$spec" ] && tag=$(basename "${f%.*}")
    [ -f "$f" ] || die "Multiwfn command file '$f' not found"
    [[ "$tag" =~ ^[A-Za-z0-9_.-]+$ ]] || die "bad TAG '$tag' (letters, digits, _ . - only)"
    MW_FILES+=("$(cd "$(dirname "$f")" && pwd)/$(basename "$f")")
    MW_TAGS+=("$tag")
done
if [ ${#MW_FILES[@]} -gt 0 ] && [ "$DRY" -eq 0 ]; then
    command -v "$MWFN_BIN" >/dev/null 2>&1 || die "Multiwfn ('$MWFN_BIN') not found; use --multiwfn PATH"
    export OMP_STACKSIZE="${OMP_STACKSIZE:-200M}"
    ulimit -s unlimited 2>/dev/null || true
fi

# collect checkpoints: prefer the chk when both chk and fchk exist
FILES=()
for t in "${TARGETS[@]}"; do
    if [ -d "$t" ]; then
        depth=(-maxdepth 1); [ "$RECURSIVE" -eq 1 ] && depth=()
        while IFS= read -r f; do FILES+=("$f"); done < <(
            find "$t" "${depth[@]+"${depth[@]}"}" -type f \( -name '*.chk' -o -name '*.fchk' -o -name '*.fch' \) | sort |
            awk '{ s = $0; sub(/\.(chk|fchk|fch)$/, "", s); if (!(s in seen)) { seen[s] = 1; print } }')
    elif [ -f "$t" ]; then
        FILES+=("$t")
    else
        warn "'$t' not found - skipped"
    fi
done
[ ${#FILES[@]} -gt 0 ] || die "no .chk/.fchk files found"

if [ -n "$OUTDIR" ]; then
    run mkdir -p "$OUTDIR"
    [ "$DRY" -eq 1 ] || OUTDIR=$(cd "$OUTDIR" && pwd)
fi
SUMMARY="${OUTDIR:-.}/batch_cube_summary.tsv"
if [ "$DRY" -eq 0 ]; then
    printf 'file\tn_alpha\tn_beta\topen_shell\thomo_eV\tlumo_eV\tgap_eV\n' > "$SUMMARY"
fi

# fchk helpers ----------------------------------------------------------------------------
fchk_int() { awk -v k="$2" 'index($0, k) == 1 { print $NF; exit }' "$1"; }
open_shell() { grep -q '^Beta Orbital Energies' "$1"; }
orbital_energies() {   # prints "homo_eV lumo_eV" for the alpha set
    awk '
    /^Number of alpha electrons/ { nocc = $NF }
    /^Alpha Orbital Energies/    { nmo = $NF; inside = 1; k = 0; next }
    inside == 1 { for (i = 1; i <= NF; i++) { k++; e[k] = $i }; if (k >= nmo) inside = 0 }
    END {
        h = 27.211386
        if (nocc > 0 && nocc + 1 <= nmo) printf "%.4f %.4f\n", e[nocc] * h, e[nocc + 1] * h
        else if (nocc > 0) printf "%.4f NA\n", e[nocc] * h
        else print "NA NA"
    }' "$1"
}

# cubegen wrapper: cube kind file npts [template]
cube() {
    local np="$1" kind="$2" fchk="$3" out="$4" npts="$5" tmpl="${6:-}"
    if [ -f "$out" ] && [ "$FORCE" -eq 0 ]; then echo "  = $out exists (use -f to redo)"; return 0; fi
    local cmd=(cubegen)
    [ -n "$MEMOPT" ] && cmd+=("$MEMOPT")
    cmd+=("$np" "$kind" "$fchk" "$out" "$npts" h)
    [ -n "$tmpl" ] && cmd+=("$tmpl")
    run "${cmd[@]}" || { warn "cubegen $kind failed for $fchk"; return 1; }
}

orbital_cube() {   # label index fchk base open
    local label="$1" idx="$2" fchk="$3" base="$4" open="$5" na="$6" nb="$7" rc=0
    if [ "$open" -eq 1 ]; then
        local ia ib
        ia=$(resolve "$idx" "$na"); ib=$(resolve "$idx" "$nb")
        if [ "$ia" -ge 1 ]; then cube 0 "AMO=$ia" "$fchk" "${base}_${label}_alpha.cub" "$NPTS" || rc=1; fi
        if [ "$ib" -ge 1 ]; then cube 0 "BMO=$ib" "$fchk" "${base}_${label}_beta.cub" "$NPTS" || rc=1
        else echo "  - beta $label skipped: no beta electrons"; fi
    else
        local i
        i=$(resolve "$idx" "$na")
        if [ "$i" -ge 1 ]; then cube 0 "MO=$i" "$fchk" "${base}_${label}.cub" "$NPTS" || rc=1
        else warn "$label: orbital index $i < 1"; rc=1; fi
    fi
    return "$rc"
}

resolve() {   # "homo-1" / "lumo+2" / "25" -> absolute MO index for nocc=$2
    local s="$1" nocc="$2"
    case "$s" in
        homo) echo "$nocc" ;;
        lumo) echo $(( nocc + 1 )) ;;
        homo-*) echo $(( nocc - ${s#homo-} )) ;;
        lumo+*) echo $(( nocc + 1 + ${s#lumo+} )) ;;
        *) echo "$s" ;;
    esac
}

# main loop -------------------------------------------------------------------------------
NFAIL=0
IFS=',' read -r -a KLIST <<< "$KINDS"
for f in "${FILES[@]}"; do
    dir=$(cd "$(dirname "$f")" && pwd)
    name=$(basename "$f")
    stem="${name%.*}"
    echo "== $dir/$name"
    case "$name" in
        *.chk)
            fchk="$dir/$stem.fchk"
            if [ "$FORCE" -eq 1 ] || [ ! -f "$fchk" ] || [ "$dir/$name" -nt "$fchk" ]; then
                fc=(formchk); [ -n "$MEMOPT" ] && fc+=("$MEMOPT")
                run "${fc[@]}" "$dir/$name" "$fchk" || { warn "formchk failed for $name"; NFAIL=$((NFAIL + 1)); continue; }
            fi ;;
        *) fchk="$dir/$name" ;;
    esac
    out="${OUTDIR:-$dir}"
    base="$out/$stem"
    if [ "$DRY" -eq 1 ]; then
        na=1; nb=1; open=0
    else
        [ -s "$fchk" ] || { warn "$fchk missing or empty"; NFAIL=$((NFAIL + 1)); continue; }
        na=$(fchk_int "$fchk" "Number of alpha electrons"); nb=$(fchk_int "$fchk" "Number of beta electrons")
        [[ "$na" =~ ^[0-9]+$ ]] || { warn "$fchk: cannot read the number of electrons (not an fchk?)"; NFAIL=$((NFAIL + 1)); continue; }
        open=0; open_shell "$fchk" && open=1
        read -r eh el <<< "$(orbital_energies "$fchk")"
        gap=NA
        [ "$eh" != NA ] && [ "$el" != NA ] && gap=$(awk -v a="$eh" -v b="$el" 'BEGIN{printf "%.4f", b - a}')
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$fchk" "$na" "$nb" "$open" "$eh" "$el" "$gap" >> "$SUMMARY"
        echo "  N_alpha=$na N_beta=$nb open_shell=$open  HOMO(alpha)=$eh eV  LUMO(alpha)=$el eV  gap=$gap eV"
    fi

    if [ "$NOCUBE" -eq 0 ]; then
        ok=1
        for k in "${KLIST[@]}"; do
            k=$(echo "$k" | tr -d ' ')
            lk=$(echo "$k" | tr 'A-Z' 'a-z')
            case "$lk" in
                "") ;;
                homo|lumo|homo-[0-9]*|lumo+[0-9]*)
                    orbital_cube "$lk" "$lk" "$fchk" "$base" "$open" "$na" "$nb" || ok=0 ;;
                mo=[0-9]*)
                    orbital_cube "mo${lk#mo=}" "${lk#mo=}" "$fchk" "$base" "$open" "$na" "$nb" || ok=0 ;;
                density)
                    cube 0 Density=SCF "$fchk" "${base}_dens.cub" "$NPTS" || ok=0 ;;
                spin)
                    if [ "$open" -eq 1 ]; then cube 0 Spin=SCF "$fchk" "${base}_spin.cub" "$NPTS" || ok=0
                    else echo "  - spin skipped: closed-shell wavefunction"; fi ;;
                esp)
                    cube 0 Density=SCF "$fchk" "${base}_dens.cub" "$NPTS" || ok=0
                    cube "$NPROC" Potential=SCF "$fchk" "${base}_esp.cub" -1 "${base}_dens.cub" || ok=0 ;;
                *=*)
                    tagk=$(echo "$k" | tr -c 'A-Za-z0-9_.-' '_' | sed 's/_*$//')
                    cube "$NPROC" "$k" "$fchk" "${base}_${tagk}.cub" "$NPTS" || ok=0 ;;
                *)
                    warn "unknown kind '$k' (see -h)"; ok=0 ;;
            esac
        done
        [ "$ok" -eq 1 ] || NFAIL=$((NFAIL + 1))
    fi

    for i in "${!MW_FILES[@]}"; do
        wdir="${base}_${MW_TAGS[$i]}"
        echo "  Multiwfn ${MW_TAGS[$i]} -> $wdir/"
        if [ "$DRY" -eq 1 ]; then
            echo "  + (cd $wdir && $MWFN_BIN $fchk -silent < ${MW_FILES[$i]} > ${MW_TAGS[$i]}.out)"
            continue
        fi
        mkdir -p "$wdir"
        ( cd "$wdir" && "$MWFN_BIN" "$fchk" -silent < "${MW_FILES[$i]}" > "${MW_TAGS[$i]}.out" 2>&1 ) || true
        if [ -s "$wdir/${MW_TAGS[$i]}.out" ]; then
            echo "    output: $wdir/${MW_TAGS[$i]}.out ($(find "$wdir" -type f | wc -l) file(s))"
        else
            warn "Multiwfn produced no output for $name (${MW_TAGS[$i]})"; NFAIL=$((NFAIL + 1))
        fi
    done
done

[ "$DRY" -eq 1 ] || echo "summary: $SUMMARY"
if [ "$NFAIL" -gt 0 ]; then
    warn "$NFAIL problem(s) - see the messages above"
    exit 1
fi
echo "done: ${#FILES[@]} file(s). Suggested isovalues: orbitals 0.02-0.05 a.u., density surface 0.001 a.u. for ESP mapping."
