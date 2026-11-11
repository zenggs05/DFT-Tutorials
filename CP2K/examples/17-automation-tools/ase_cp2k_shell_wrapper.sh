#!/bin/bash
# Wrapper that starts CP2K in interactive shell mode for the ASE CP2K calculator.
# Chapter 17 of the CP2K tutorial.
#
# Why: ASE <= 3.28 checks that its CP2K command contains the string "cp2k_shell",
# but CMake builds of CP2K (2026.1+) only install cp2k.psmp / cp2k.ssmp. The file
# name of this wrapper contains "cp2k_shell", so pointing ASE_CP2K_COMMAND at it
# satisfies the check; inside, cp2k.psmp is started with --shell (-s).
# ASE >= 3.29 accepts "cp2k.psmp -s" directly and does not need this wrapper.
#
# Usage:
#   chmod +x ase_cp2k_shell_wrapper.sh
#   export ASE_CP2K_COMMAND="$PWD/ase_cp2k_shell_wrapper.sh"
#   export CP2K_SHELL_NP=4          # MPI ranks (default 1)
#   export OMP_NUM_THREADS=1
#   python ase_h2o_opt.py
#
# Alternative without a script: a symbolic link whose name starts with cp2k_shell
# also switches CP2K into shell mode, e.g.
#   ln -s "$(command -v cp2k.psmp)" ~/bin/cp2k_shell.psmp
#
# Notes
#   - ASE talks to CP2K through stdin/stdout, so nothing else may be printed to
#     stdout here (no echo before exec).
#   - With MPI, stdin must reach rank 0. mpirun/mpiexec forward stdin to rank 0 by
#     default; if you use srun, check the --input option of your Slurm version.
#   - Load your CP2K environment (module load / source cp2k_env / conda activate)
#     before starting Python, or add it below (redirect its output to stderr).

set -e
NP=${CP2K_SHELL_NP:-1}
CP2K_EXE=${CP2K_EXE:-cp2k.psmp}

if [ "$NP" -gt 1 ]; then
  exec mpirun -np "$NP" "$CP2K_EXE" --shell
else
  exec "$CP2K_EXE" --shell
fi
