import ast
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]

class TestFragmentSyntax(unittest.TestCase):
    def test_embedded_method_compiles(self):
        # README states it is a *method fragment*, not a standalone executable.
        src=(ROOT/'Scripts/dataset/dpdata-system-perturb.py').read_text()
        compile('class SystemExample:\n'+src,'<fragment>','exec')
        ast.parse('class SystemExample:\n'+src)
    def test_index_selection_branches_present(self):
        src=(ROOT/'Scripts/dataset/dpdata-system-perturb.py').read_text()
        self.assertIn('pert_indices = list(range(',src)
        self.assertIn('pert_indices = [idx for idx',src)
        self.assertIn('pert_natoms = int(atom_pert_prob * len(pert_indices))',src)
if __name__=='__main__':unittest.main()
