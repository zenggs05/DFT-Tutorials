import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
ROOT = Path(__file__).resolve().parents[2]
VIB = ROOT/'CP2K/scripts/vib_thermo.py'
HARM = ROOT/'CP2K/examples/14-vibration-thermo/harmonic_thermo.py'

def load(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

class TestThermochemistry(unittest.TestCase):
    def test_low_temperature_limit(self):
        m=load(VIB)
        z,u,s=m.vib_terms([4000.0], 1.0)
        self.assertTrue(all(math.isfinite(x) for x in (z,u,s)))
        self.assertAlmostEqual(z,u,places=5)
        self.assertAlmostEqual(s,0.0,places=10)
    def test_finite_typical_temperature(self):
        m=load(VIB)
        z,u,s=m.vib_terms([250,1000],298.15)
        self.assertGreater(u,z)
        self.assertGreater(s,0.0)
    def test_reject_invalid_temperature_vib_terms(self):
        m=load(VIB)
        for t in [0,-1,float('nan'),float('inf'),float('-inf')]:
            with self.subTest(t=t), self.assertRaises(ValueError):m.vib_terms([100],t)
    def test_harmonic_cli_low_temp_and_invalid(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'freq.mol'
            p.write_text('[FREQ]\n4000\n')
            ok=subprocess.run([sys.executable,str(HARM),str(p),'-T','1'],capture_output=True,text=True)
            self.assertEqual(ok.returncode,0,ok.stderr)
            self.assertIn('U_vib(T)',ok.stdout)
            self.assertNotIn('nan',ok.stdout.lower())
            for t in ['0','nan','inf','-1']:
                res=subprocess.run([sys.executable,str(HARM),str(p),'-T',t],capture_output=True,text=True)
                self.assertNotEqual(res.returncode,0,t)
                self.assertIn('temperature',res.stderr.lower())
    def test_vib_cli_nonfinite(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'dummy.out';p.write_text('dummy')
            for t in ['nan','inf','-1']:
                res=subprocess.run([sys.executable,str(VIB),str(p),'-T',t],capture_output=True,text=True)
                self.assertNotEqual(res.returncode,0,t)
                self.assertIn('temperature',res.stderr.lower())
if __name__ == '__main__': unittest.main()
