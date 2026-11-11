# =====================================================================
# DFT environment snippet - append to ~/.bashrc (or source it from there)
#   cat bashrc_dft_snippet.sh >> ~/.bashrc     # after editing the paths
#   source ~/.bashrc
# Every path below is an EXAMPLE: change it to your own installation.
# Keep secrets (API keys) out of git: chmod 600 ~/.bashrc
# =====================================================================

# Variables needed by batch jobs come first; aliases and conda activation
# (interactive shells only) are at the end of the file.

# ---------- VASP pseudopotentials (directory that CONTAINS potpaw_PBE) ---
# $HOME/vasp_pp/potpaw_PBE/Si/POTCAR, $HOME/vasp_pp/potpaw_PBE/Fe_pv/POTCAR ...
export VASP_PP_PATH="$HOME/vasp_pp"
# ASE builds the POTCAR from VASP_PP_PATH; how ASE starts VASP:
export ASE_VASP_COMMAND="mpirun -np 8 vasp_std"
# (ASE >= 3.23 prefers ~/.config/ase/config.ini, see ase_config.ini.example;
#  ASE 3.29 still reads these variables but warns ASEEnvDeprecationWarning)

# ---------- pymatgen ---------------------------------------------------
# Preferred: run "pmg config -p <raw> <out>" and "pmg config --add ..." once
# (writes ~/.config/.pmgrc.yaml). Environment variables override the file:
# export PMG_VASP_PSP_DIR="$HOME/pmg_potcars"
# export PMG_DEFAULT_FUNCTIONAL="PBE_54"

# ---------- Materials Project (new API, keys are 32 characters) ---------
# export MP_API_KEY="put-your-32-character-key-here"     # used by mp-api
# export PMG_MAPI_KEY="$MP_API_KEY"                      # used by pymatgen
unset MP_API_ENDPOINT PMG_MAPI_ENDPOINT   # old v2 endpoints no longer work

# ---------- Gaussian 16 ------------------------------------------------
# export g16root="/opt/gaussian"                 # directory containing g16/
# [ -f "$g16root/g16/bsd/g16.profile" ] && source "$g16root/g16/bsd/g16.profile"
# export GAUSS_SCRDIR="/tmp/$USER/gauss_scr"     # node-local, must exist
# mkdir -p "$GAUSS_SCRDIR" 2>/dev/null

# ---------- CP2K -------------------------------------------------------
# export CP2K_DATA_DIR="$HOME/cp2k/data"         # BASIS_MOLOPT_UZH, POTENTIAL_UZH ...
# export ASE_CP2K_COMMAND="mpirun -np 4 cp2k.psmp -s"

# ---------- OpenMX -----------------------------------------------------
# DATA.PATH in the .dat input must point to the DFT_DATA19 directory;
# ASE's OpenMX calculator reads it from:
# export OPENMX_DFT_DATA_PATH="$HOME/openmx4.0/DFT_DATA19"

# ---------- VASPKIT / Multiwfn -----------------------------------------
# export PATH="$HOME/software/vaspkit.1.5.1/bin:$PATH"
# export Multiwfnpath="$HOME/software/Multiwfn"
# export PATH="$Multiwfnpath:$PATH"
export OMP_STACKSIZE=200M          # Multiwfn/OpenMP codes crash with small stacks

# ---------- stack size for Fortran codes (also put it in job scripts) ---
ulimit -s unlimited 2>/dev/null || ulimit -s "$(ulimit -H -s)" 2>/dev/null

# ---------------------- interactive-only below -------------------------
case $- in *i*) ;; *) return 0 2>/dev/null ;; esac

# conda/mamba (Miniforge). "conda init bash" writes a similar block itself.
# [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ] && source "$HOME/miniforge3/etc/profile.d/conda.sh"
# conda activate dft

alias sq='squeue -u $USER -o "%.10i %.9P %.20j %.2t %.10M %.6D %R"'
alias lt='ls -lhtr'
alias d2u='dos2unix'
# last TOTEN of every OUTCAR below the current directory
toten() {
    find . -name OUTCAR | sort | while read -r f; do
        printf "%-50s %s\n" "$f" \
            "$(grep 'free  energy   TOTEN' "$f" | tail -n 1 | awk '{print $5}')"
    done
}
