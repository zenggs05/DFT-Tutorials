! h2o_force_eval.inc - shared FORCE_EVAL for the isolated H2O molecule (chapter 14)
! Included by h2o_opt_tight.inp and h2o_vib_ir.inp, so that the geometry
! optimization and the frequency calculation use EXACTLY the same settings.
! The including file must define COORD.
!
! PBE / TZV2P-MOLOPT-PBE-GTH (UZH), 600 Ry, tight SCF (finite-difference
! Hessians amplify force noise), non-periodic Poisson solver (MT).

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    &MGRID
      CUTOFF 600
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-14
    &END QS
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER MT
    &END POISSON
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.0E-8
      MAX_SCF 30
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_ALL
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-8
        MAX_SCF 20
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
    &END XC
    &PRINT
      ! dipole moment of every displaced geometry -> IR intensities
      ! (needed by INTENSITIES in h2o_vib_ir.inp, harmless in GEO_OPT)
      &MOMENTS ON
        PERIODIC F
        REFERENCE COM
      &END MOMENTS
    &END PRINT
  &END DFT
  &SUBSYS
    &CELL
      ABC 12.0 12.0 12.0
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME ${COORD}
      COORD_FILE_FORMAT XYZ
      &CENTER_COORDINATES
      &END CENTER_COORDINATES
    &END TOPOLOGY
    &KIND O
      BASIS_SET TZV2P-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
    &KIND H
      BASIS_SET TZV2P-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
      ! isotope test (D2O): uncomment the next line
      ! MASS 2.014102
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
