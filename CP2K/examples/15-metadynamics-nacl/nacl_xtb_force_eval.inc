! nacl_xtb_force_eval.inc - GFN1-xTB FORCE_EVAL for Na+ Cl- + 62 H2O (chapter 15)
! Cheap level for equilibration and for a first exploration of the free-energy
! profile ("pre-sampling"); refine with nacl_dft_force_eval.inc afterwards.
! The including file must define COORD.  Atom 1 = Na, atom 2 = Cl, then O H H ...
! NOTE: GFN1-xTB describes liquid water only qualitatively.

&FORCE_EVAL
  METHOD QUICKSTEP
  &DFT
    &QS
      METHOD XTB
      EXTRAPOLATION ASPC
      EXTRAPOLATION_ORDER 3
      &XTB
        ! native CP2K GFN1-xTB (default GFN_TYPE)
        GFN_TYPE 1
        ! Ewald for the periodic cell: required up to CP2K 2026.1; in 2026.2 it is
        ! deprecated (periodicity is taken from CELL/PERIODIC) but still accepted
        DO_EWALD T
      &END XTB
    &END QS
    &SCF
      SCF_GUESS MOPAC
      EPS_SCF 1.0E-6
      MAX_SCF 50
      &OT ON
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF ON
        EPS_SCF 1.0E-6
        MAX_SCF 20
      &END OUTER_SCF
      &PRINT
        &RESTART
          BACKUP_COPIES 0
        &END RESTART
      &END PRINT
    &END SCF
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
