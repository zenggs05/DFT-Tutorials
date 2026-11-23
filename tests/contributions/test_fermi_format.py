import importlib.util
from pathlib import Path
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]

def load(name,rel):
    spec=importlib.util.spec_from_file_location(name,ROOT/rel)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

class TestFermiNumbers(unittest.TestCase):
    def test_formats(self):
        a=load('wf','CP2K/examples/11-cube-charge-workfunction/workfunction.py')
        b=load('ct','CP2K/scripts/cube_tools.py')
        for number in ('-0.125','-1.25E-1','-1.25D-1','-1.25d-1','-.125'):
            with self.subTest(number=number), tempfile.TemporaryDirectory() as td:
                p=Path(td)/'cp2k.out';p.write_text(' Fermi energy: '+number+'\n')
                self.assertAlmostEqual(a.fermi_energy(p),-.125*a.HA2EV)
                self.assertAlmostEqual(b.read_fermi_ev(p),-.125*b.HA2EV)
if __name__=='__main__':unittest.main()
