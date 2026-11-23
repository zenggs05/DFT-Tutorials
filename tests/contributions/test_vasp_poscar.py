import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
def load(rel):
    path = ROOT / rel
    spec = importlib.util.spec_from_file_location('module_' + path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

ELASTIC = 'VASP/examples/13-elastic-mos2-2d/fit_2d_elastic.py'
EFFMASS = 'VASP/examples/17-effmass-gaas/fit_effmass.py'

class TestPoscarScaling(unittest.TestCase):
    def setUp(self):
        self.a = load(ELASTIC)
        self.b = load(EFFMASS)
    def check(self, scale, area, lattice):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / 'POSCAR'
            p.write_text('synthetic\n'+scale+'\n1 0 0\n0 1 0\n0 0 1\nX\n1\nDirect\n0 0 0\n')
            self.assertAlmostEqual(self.a.area_from_poscar(p), area)
            np.testing.assert_allclose(self.b.read_lattice(p), np.diag(lattice))
    def test_normal_scale(self): self.check('2', 4, [2,2,2])
    def test_negative_volume(self): self.check('-8', 4, [2,2,2])
    def test_three_scales(self): self.check('2 3 4', 6, [2,3,4])
    def test_invalid_scales(self):
        for scale in ('0', 'nan', '2 -3 4', '2 3'):
            with self.subTest(scale=scale), tempfile.TemporaryDirectory() as td:
                p = Path(td) / 'POSCAR'
                p.write_text('test\n'+scale+'\n1 0 0\n0 1 0\n0 0 1\n')
                with self.assertRaises(ValueError): self.a.area_from_poscar(p)
                with self.assertRaises(ValueError): self.b.read_lattice(p)
if __name__ == '__main__': unittest.main()
