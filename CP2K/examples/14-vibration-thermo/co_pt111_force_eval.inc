! co_pt111_force_eval.inc - shared FORCE_EVAL for CO atop on a Pt(111) 3x3x3 slab
! Included by co_pt111_opt.inp and co_on_pt_partial_hessian.inp (chapter 14).
! The including file must define COORD.
!
! System : 27 Pt (q18) + C (q4) + O (q6) = 496 valence electrons, closed shell
! Method : PBE + D3(BJ), DZVP-MOLOPT-PBE-GTH (UZH), Gamma point only,
!          diagonalization + Broyden mixing + Fermi-Dirac smearing (metal)
! Cell   : hexagonal 3x3 surface cell, a(Pt) = 3.97 A (about the PBE value),
!          ~16 A vacuum.  Teaching setup: add &KPOINTS (3 3 1) for production.

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    SURFACE_DIPOLE_CORRECTION T
    SURF_DIP_DIR Z
    &MGRID
      CUTOFF 500
      REL_CUTOFF 60
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      ! tight: the CO stretch is computed by finite differences of forces
      EPS_SCF 1.0E-8
      MAX_SCF 200
      ! 248 occupied MOs
      ADDED_MOS 100
      &DIAGONALIZATION ON
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &MIXING ON
        METHOD BROYDEN_MIXING
        ALPHA 0.2
        NBUFFER 8
      &END MIXING
      &SMEAR ON
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 500
      &END SMEAR
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
      ABC 8.42164 8.42164 24.0
      ALPHA_BETA_GAMMA 90.0 90.0 60.0
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME ${COORD}
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND Pt
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q18
      POTENTIAL GTH-PBE-q18
    &END KIND
    &KIND C
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q4
      POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND O
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q6
      POTENTIAL GTH-PBE-q6
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
