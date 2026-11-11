SYSTEM   = NiO AFM-II 2x2x2, plain DFT ground state for linear-response U

PREC     = Accurate
ENCUT    = 520
EDIFF    = 1E-6
NELM     = 150
LASPH    = .TRUE.
LREAL    = .FALSE.
ISMEAR   = 0
SIGMA    = 0.05

ISPIN    = 2
# species order in POSCAR: Ni(1, perturbed)  Ni(15)  O(16)
MAGMOM   = 2 -2 2 -2 2 -2 2 -2 2 -2 2 -2 2 -2 2 -2 16*0
LORBIT   = 11          # prints s/p/d occupations per site ("total charge" table)
LMAXMIX  = 4           # CHGCAR must carry d one-centre terms for ICHARG = 11

# symmetry: atom 1 is a different species, so VASP already lowers the symmetry.
NSW      = 0
LWAVE    = .TRUE.
LCHARG   = .TRUE.
