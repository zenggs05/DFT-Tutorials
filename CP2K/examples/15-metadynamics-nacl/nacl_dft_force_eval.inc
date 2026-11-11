! nacl_dft_force_eval.inc - PBE-D3(BJ) FORCE_EVAL for Na+ Cl- + 62 H2O (chapter 15)
! Production level for AIMD-based metadynamics / blue-moon sampling.
! The including file must define COORD.  Atom 1 = Na, atom 2 = Cl, then O H H ...
! 512 valence electrons (Na q9, Cl q7, O q6, H q1); OT + ASPC as in chapter 12.

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF 400
      REL_CUTOFF 50
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-6
      MAX_SCF 30
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
      &PRINT
        &RESTART
          BACKUP_COPIES 0
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        DISPERSION_FUNCTIONAL PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 12.42 12.42 12.42
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME ${COORD}
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND Na
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q9
      POTENTIAL GTH-PBE-q9
    &END KIND
    &KIND Cl
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q7
      POTENTIAL GTH-PBE-q7
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
    ! COLVAR 1: Na-Cl distance (internal unit: bohr)
    &COLVAR
      &DISTANCE
        ATOMS 1 2
      &END DISTANCE
    &END COLVAR
    ! COLVAR 2: coordination number of Na by water oxygens (dimensionless)
    &COLVAR
      &COORDINATION
        ATOMS_FROM 1
        KINDS_TO O
        R0 [angstrom] 3.0
        NN 12
        ND 24
      &END COORDINATION
    &END COLVAR
  &END SUBSYS
&END FORCE_EVAL
