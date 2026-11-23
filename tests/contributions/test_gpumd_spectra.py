import ast
import os
import pathlib
import runpy
import sys
import tempfile
import types
import unittest
from unittest import mock
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=pathlib.Path(__file__).resolve().parents[2]

class TestGPUMDSpectra(unittest.TestCase):
    def test_float_omega_parser_in_three_scripts(self):
        # Evaluate actual expression from the AST, rather than searching for a literal text snippet.
        paths = ['plot_kappa_shc.py','plot_shc_quantum_correction.py','plot_compute_shc.py']
        for name in paths:
            path=ROOT/'Scripts/gpumd'/name
            tree=ast.parse(path.read_text())
            found=[]
            for node in ast.walk(tree):
                if isinstance(node,ast.Assign):
                    target_names=[]
                    for tar in node.targets:
                        if isinstance(tar,ast.Name): target_names.append(tar.id)
                        if isinstance(tar,ast.Tuple): target_names += [el.id for el in tar.elts if isinstance(el,ast.Name)]
                    if 'num_omega' in target_names: found.append(ast.get_source_segment(path.read_text(),node))
            with self.subTest(name=name):
                self.assertTrue(found)
                self.assertTrue(any('float(line.split()[5])' in x for x in found),found)
                self.assertEqual(float('400.0'),400.0)
    def test_float_omega_in_complete_kappa_script(self):
        # Execute the full script with 400.0 and a lightweight ASE interface stub.
        class Cell:
            def lengths(self): return np.array([10., 10., 10.])
        class Atoms:
            cell=Cell()
            pbc=None
            def __mul__(self, factor): return self
        with tempfile.TemporaryDirectory() as td:
            t=pathlib.Path(td)
            (t/'run.in').write_text('time_step 0.5\nnvt 1000 300\ncompute_shc 2 2 0 5 400.0\ncompute_hnemd 100 0.0001 0.0 0.0\n')
            np.savetxt(t/'kappa.out',np.ones((10,5))*2)
            s=np.zeros((8,3));s[:3,0]=[-.002,0,.002];s[:3,1]=1
            s[-5:,0]=[10,20,30,40,50];s[-5:,1:]=1
            np.savetxt(t/'shc.out',s)
            np.savetxt(t/'thermo.out',np.ones((12,11))*10)
            np.save(t/'Gc.npy',np.ones(5))
            ase=types.ModuleType('ase');io=types.ModuleType('ase.io')
            io.read=lambda *_:Atoms();io.write=lambda *args:None
            before=pathlib.Path.cwd()
            try:
                os.chdir(t)
                with mock.patch.dict(sys.modules,{'ase':ase,'ase.io':io}):
                    runpy.run_path(str(ROOT/'Scripts/gpumd/plot_kappa_shc.py'),run_name='__main__')
                self.assertTrue((t/'hnemd.png').exists())
                self.assertTrue((t/'shc.png').exists())
            finally:
                os.chdir(before);plt.close('all')
    def test_masked_integration_full_script(self):
        with tempfile.TemporaryDirectory() as td:
            t=pathlib.Path(td)
            (t/'run.in').write_text('nvt 10000 300\ncompute_shc 2 2 0 5 400.0\ncompute_hnemd 100 0.0001 0 0\n')
            block=np.zeros((8,3))
            block[:3,0]=[-0.002,0,0.002]
            block[:3,1]=1
            block[-5:,0]=np.array([10,20,30,40,60])*2*np.pi
            block[-5:,1:]=np.array([[1,0.5],[2,0.5],[3,0.5],[4,0.5],[5,0.5]])
            np.savetxt(t/'shc.out',np.vstack([block,block]))
            thermo=np.ones((12,11))*10
            np.savetxt(t/'thermo.out',thermo)
            np.save(t/'Gc.npy',np.ones(5))
            ase=types.ModuleType('ase');io=types.ModuleType('ase.io');io.read=lambda *_:None
            before=pathlib.Path.cwd();argv=sys.argv
            try:
                os.chdir(t);sys.argv=['plot_shc_multiple.py','2']
                with mock.patch.dict(sys.modules,{'ase':ase,'ase.io':io}):
                    runpy.run_path(str(ROOT/'Scripts/gpumd/plot_shc_multiple.py'),run_name='__main__')
            finally:
                os.chdir(before);sys.argv=argv;plt.close('all')
            self.assertTrue((t/'shc-multiple.png').exists())
            self.assertGreater((t/'shc-multiple.png').stat().st_size,0)
if __name__=='__main__':unittest.main()
