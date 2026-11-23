import os
from pathlib import Path
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
ROOT=Path(__file__).resolve().parents[2]

class FakeCell:
    def lengths(self):return np.array([10.,10.,10.])
class FakeAtoms:
    cell=FakeCell()
    pbc=None
    def __mul__(self,coords):return self

def run_case(direction):
    with tempfile.TemporaryDirectory() as td:
        t=Path(td)
        (t/'run.in').write_text(f'time_step 0.5\nnvt 100 300\ncompute_shc 2 2 {direction} 5 400\ncompute_hnemd 100 0.0001 0.0001 0.0001\n')
        kappa=np.tile([1.,2.,9.,10.,5.],(10,1))
        np.savetxt(t/'kappa.out',kappa)
        shc=np.zeros((8,3));shc[:3,0]=[-.002,0,.002];shc[:3,1]=1
        shc[-5:,0]=[10,20,30,40,50];shc[-5:,1]=np.array([1,2,3,4,5]);shc[-5:,2]=1
        np.savetxt(t/'shc.out',shc)
        np.savetxt(t/'thermo.out',np.ones((12,11))*10)
        np.save(t/'Gc.npy',np.ones(5))
        ase=types.ModuleType('ase');io=types.ModuleType('ase.io')
        io.read=lambda *_:FakeAtoms();io.write=lambda *args:None
        before=Path.cwd()
        try:
            os.chdir(t)
            with mock.patch.dict(sys.modules,{'ase':ase,'ase.io':io}):
                runpy.run_path(str(ROOT/'Scripts/gpumd/plot_kappa_shc.py'),run_name='__main__')
            fig=plt.figure(1)
            line=fig.axes[0].lines[0]
            x=np.array(line.get_xdata())
            y=np.array(line.get_ydata())
            assert abs(x[-1] - 0.0005)<1e-12, x
            expected=1.0 if direction==0 else 9.0
            assert np.all(y==expected),y
            assert (t/'hnemd.png').exists()
            assert (t/'shc.png').exists()
        finally:os.chdir(before);plt.close('all')

class TestPlot(unittest.TestCase):
    def test_x_plot_time_and_columns(self):run_case(0)
    def test_y_plot_time_and_columns(self):run_case(1)
if __name__=='__main__':unittest.main()
