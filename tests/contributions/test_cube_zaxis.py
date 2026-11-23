import importlib.util
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2]
FILE=ROOT/'Scripts/analysis/extract_vacuum_level_cube.py'

class TestCubeAxis(unittest.TestCase):
    def test_origin_and_sample_spacing(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td);(t/'sample.cube').write_text('unused by mocked reader')
            cube={'data':np.arange(4).reshape(1,1,4).astype(float),
                 'origin':np.array([0.,0.,1.5]), 'spacing':np.diag([2.,2.,2.])}
            ase=types.ModuleType('ase');io=types.ModuleType('ase.io');sub=types.ModuleType('ase.io.cube')
            sub.read_cube=lambda handle:cube
            before=Path.cwd()
            try:
                os.chdir(t)
                # Script selects TkAgg for interactive use; ignore this in headless tests.
                with mock.patch.dict(sys.modules,{'ase':ase,'ase.io':io,'ase.io.cube':sub}):
                    with mock.patch.object(matplotlib,'use'):
                        spec=importlib.util.spec_from_file_location('extract_tested',FILE)
                        mod=importlib.util.module_from_spec(spec)
                        spec.loader.exec_module(mod)
                    with mock.patch.object(plt,'show'):
                        mod.get_vacuum_level_and_export('sample.cube')
                    arr=np.loadtxt('Electrostatic_Potential_1D.txt',skiprows=1)
                    np.testing.assert_allclose(arr[:,0],[1.5,3.5,5.5,7.5])
                    self.assertTrue((t/'Electrostatic_Potential.png').exists())
                    cube['spacing'][0,2]=0.25
                    with self.assertRaises(ValueError):mod.get_vacuum_level_and_export('sample.cube')
            finally:os.chdir(before);plt.close('all')
if __name__=='__main__':unittest.main()
