#!/bin/bash
# Quasi-harmonic thermochemistry with GoodVibes (option names as documented for GoodVibes 4.4).
# Install:  pip install goodvibes     (or: conda install -c conda-forge goodvibes)
# Run after the Gaussian jobs of this directory have finished:  bash goodvibes_run.sh
# Full option list:  goodvibes -h
set -e

echo "== 1) H2O: Grimme qRRHO entropy + Head-Gordon enthalpy, 100 cm-1 cutoff, 1 atm"
goodvibes h2o_opt_freq_raman.log -q -f 100

echo "== 2) H2O: same, but 1 mol/L standard state (G rises by ~1.89 kcal/mol at 298.15 K)"
goodvibes h2o_opt_freq_raman.log -q -f 100 --conc 1.0

echo "== 3) H2O: no frequency scaling, to compare with Gaussian's own thermochemistry"
goodvibes h2o_opt_freq_raman.log -q -v 1.0

echo "== 4) Ethanol: high-level SP (Link1 step) + low-level thermal corrections, 1 mol/L"
goodvibes highlevel_sp_plus_thermo.log -q --spc link --conc 1.0

echo "== 5) Methanol: Truhlar-type quasi-harmonic treatment at 373.15 K"
goodvibes freq_temperature_pressure.log -q --qs truhlar -f 100 --temp 373.15
