! cu111_h_force_eval.inc - shared FORCE_EVAL for H on a Cu(111) 3x3x3 slab
! Included by geo_opt_endpoint.inp, h_diffusion_cineb.inp, ts_guess_vib.inp
! and dimer_ts.inp (chapter 13) and by ../14-vibration-thermo/ts_freq_check.inp.
! The including file must define COORD; EPS_SCF is optional (default 1.0E-7).
!
! System : 27 Cu (q11) + 1 H = 298 valence electrons, closed shell
! Method : PBE, DZVP-MOLOPT-PBE-GTH (UZH), Gamma point only, diagonalization +
!          Broyden mixing + Fermi-Dirac smearing (Cu is a metal: no OT here)
! Cell   : hexagonal 3x3 surface cell, a(Cu) = 3.63 A (about the PBE value),
!          ~16 A vacuum along z.  Use your own converged lattice constant.
! NOTE   : a 3x3 metal slab at Gamma only is a teaching setup; for production
!          add &KPOINTS (e.g. MONKHORST-PACK 3 3 1) and check convergence.

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT_UZH
    POTENTIAL_FILE_NAME POTENTIAL_UZH
    ! asymmetric slab (H on one side only)
    SURFACE_DIPOLE_CORRECTION T
    SURF_DIP_DIR Z
    &MGRID
      CUTOFF 400
      REL_CUTOFF 50
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF ${EPS_SCF-1.0E-7}
      MAX_SCF 200
      ! 149 occupied MOs; extra MOs for the Fermi-Dirac tail
      ADDED_MOS 60
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
    &END XC
  &END DFT
  &SUBSYS
    &CELL
      ABC 7.70039 7.70039 22.0
      ALPHA_BETA_GAMMA 90.0 90.0 60.0
      PERIODIC XYZ
    &END CELL
    &TOPOLOGY
      COORD_FILE_NAME ${COORD}
      COORD_FILE_FORMAT XYZ
    &END TOPOLOGY
    &KIND Cu
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q11
      POTENTIAL GTH-PBE-q11
    &END KIND
    &KIND H
      BASIS_SET DZVP-MOLOPT-PBE-GTH-q1
      POTENTIAL GTH-PBE-q1
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
