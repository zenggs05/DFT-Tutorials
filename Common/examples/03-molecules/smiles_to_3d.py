#!/usr/bin/env python3
"""SMILES -> 3D structure with RDKit (ETKDG conformers + MMFF94 cleanup).

  python smiles_to_3d.py "c1ccccc1" -o benzene_rdkit.xyz
  python smiles_to_3d.py "CC(=O)Nc1ccc(O)cc1" -n 30 -o paracetamol.xyz
  python smiles_to_3d.py "[NH4+]" -o nh4.xyz          # charge read from SMILES

Requirements: pip install rdkit   (or: conda install -c conda-forge rdkit)

What it does
  1. parse SMILES, add explicit hydrogens
  2. embed N conformers with ETKDGv3 (fixed random seed -> reproducible)
  3. optimize every conformer with MMFF94 (UFF if MMFF has no parameters)
  4. keep the lowest-energy conformer, write xyz; the comment line holds
     SMILES, total charge and spin multiplicity for the QM input.
The force-field geometry is only a starting point for DFT optimization.
"""
import argparse
import sys

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
except ImportError:
    sys.exit("RDKit not found: pip install rdkit")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("smiles")
    p.add_argument("-n", "--nconf", type=int, default=10)
    p.add_argument("-o", "--output", default="mol.xyz")
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    mol = Chem.MolFromSmiles(args.smiles)
    if mol is None:
        sys.exit("Invalid SMILES: %s" % args.smiles)
    mol = Chem.AddHs(mol)
    charge = Chem.GetFormalCharge(mol)
    mult = Descriptors.NumRadicalElectrons(mol) + 1

    params = AllChem.ETKDGv3()
    params.randomSeed = args.seed
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=args.nconf,
                                           params=params))
    if not cids:
        sys.exit("Embedding failed; try more conformers or another seed")

    if AllChem.MMFFHasAllMoleculeParams(mol):
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
        ff = "MMFF94"
    else:
        res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=2000)
        ff = "UFF"
    energies = [e for (_, e) in res]
    best = min(range(len(energies)), key=energies.__getitem__)
    conf = mol.GetConformer(cids[best])

    with open(args.output, "w", newline="\n") as f:
        f.write("%d\n" % mol.GetNumAtoms())
        f.write("SMILES=%s charge=%d mult=%d %s_E=%.3f_kcal/mol\n"
                % (args.smiles, charge, mult, ff, energies[best]))
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            f.write("%-2s %12.6f %12.6f %12.6f\n"
                    % (atom.GetSymbol(), pos.x, pos.y, pos.z))

    spread = max(energies) - min(energies)
    print("%s: %d atoms, charge %d, multiplicity %d"
          % (args.smiles, mol.GetNumAtoms(), charge, mult))
    print("%d conformers optimized with %s; lowest E = %.3f kcal/mol, "
          "spread = %.3f kcal/mol" % (len(energies), ff, energies[best],
                                      spread))
    print("written:", args.output)


if __name__ == "__main__":
    main()
