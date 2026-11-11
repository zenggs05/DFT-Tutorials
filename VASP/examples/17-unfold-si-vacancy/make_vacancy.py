#!/usr/bin/env python3
"""Build a 3x3x3 supercell of the Si primitive cell and remove one atom (vacancy).

Writes POSCAR_sc (53 atoms). The supercell matrix with respect to the primitive cell
is diag(3, 3, 3), i.e. use  easyunfold generate ... --matrix "3 3 3".
Requires ASE.
"""
from ase.io import read, write

prim = read("POSCAR_prim", format="vasp")
sc = prim.repeat((3, 3, 3))
del sc[0]                                   # vacancy at the origin
write("POSCAR_sc", sc, format="vasp", direct=True, sort=True)
print("POSCAR_sc: %d atoms, cell lengths %s A" % (len(sc), sc.cell.lengths().round(4)))
