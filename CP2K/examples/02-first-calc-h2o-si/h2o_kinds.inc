! h2o_kinds.inc - KIND blocks for O and H, included by preprocessor_demo.inp
! The variables BASIS (optional, default DZVP) and XC (set in the main input)
! are expanded here, i.e. included files see the variables of the parent file.
&KIND O
  BASIS_SET ${BASIS-DZVP}-MOLOPT-GTH
  POTENTIAL GTH-${XC}-q6
&END KIND
&KIND H
  BASIS_SET ${BASIS-DZVP}-MOLOPT-GTH
  POTENTIAL GTH-${XC}-q1
&END KIND
