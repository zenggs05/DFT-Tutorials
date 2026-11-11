! kind_templates.inc - reusable KIND blocks for H, C, N, O with UZH basis sets
! Tutorial: CP2K chapter 3, section 10.  Include inside &SUBSYS:
!     @INCLUDE 'kind_templates.inc'
! The main input must use the files BASIS_MOLOPT_UZH and POTENTIAL_UZH
! (and BASIS_ADMM_UZH as an additional BASIS_SET_FILE_NAME if ADMM is 1).
!
! Variables (all optional, set them in the main input or with -E):
!   Q      orbital basis quality: DZVP | TZVP | TZV2P        (default DZVP)
!   XCTAG  functional family of the UZH pair: PBE | SCAN | PBE0 (default PBE)
!          equivalent aliases: GGA | MGGA | HYB
!          The tag must match the functional used in &XC.
!   ADMM   1 adds the AUX_FIT basis needed by ADMM hybrid calculations (default 0)
!
! Valence charges: H q1, C q4, N q5, O q6. Basis and potential suffixes agree.
! KIND blocks for elements that are not present in the structure are normally ignored.

&KIND H
  BASIS_SET ${Q-DZVP}-MOLOPT-${XCTAG-PBE}-GTH-q1
  POTENTIAL GTH-${XCTAG-PBE}-q1
  @IF ${ADMM-0} == 1
  BASIS_SET AUX_FIT admm-dzp-q1
  @ENDIF
&END KIND

&KIND C
  BASIS_SET ${Q-DZVP}-MOLOPT-${XCTAG-PBE}-GTH-q4
  POTENTIAL GTH-${XCTAG-PBE}-q4
  @IF ${ADMM-0} == 1
  BASIS_SET AUX_FIT admm-dzp-q4
  @ENDIF
&END KIND

&KIND N
  BASIS_SET ${Q-DZVP}-MOLOPT-${XCTAG-PBE}-GTH-q5
  POTENTIAL GTH-${XCTAG-PBE}-q5
  @IF ${ADMM-0} == 1
  BASIS_SET AUX_FIT admm-dzp-q5
  @ENDIF
&END KIND

&KIND O
  BASIS_SET ${Q-DZVP}-MOLOPT-${XCTAG-PBE}-GTH-q6
  POTENTIAL GTH-${XCTAG-PBE}-q6
  @IF ${ADMM-0} == 1
  BASIS_SET AUX_FIT admm-dzp-q6
  @ENDIF
&END KIND
