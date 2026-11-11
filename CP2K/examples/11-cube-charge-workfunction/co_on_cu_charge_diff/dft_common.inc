# Shared &DFT settings for the charge-density-difference workflow.
# Included by full.inp, slab.inp and co.inp. The three runs MUST use the same
# cell, CUTOFF, REL_CUTOFF and STRIDE so that the cube files share one grid.
# ${NADD} (number of added MOs) is set with @SET in each input file.
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    SURFACE_DIPOLE_CORRECTION T
    SURF_DIP_DIR Z
    &MGRID
      CUTOFF     500
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &KPOINTS
      SCHEME MONKHORST-PACK 4 4 1
    &END KPOINTS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF   1.0E-7
      MAX_SCF   300
      ADDED_MOS ${NADD}
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING
        METHOD  BROYDEN_MIXING
        ALPHA   0.1
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
    &PRINT
      &E_DENSITY_CUBE
        STRIDE 1 1 1
      &END E_DENSITY_CUBE
      &V_HARTREE_CUBE
        STRIDE 1 1 1
      &END V_HARTREE_CUBE
      &HIRSHFELD ON
      &END HIRSHFELD
    &END PRINT
