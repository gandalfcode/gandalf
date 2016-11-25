from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import unittest

class SoundWaveTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/hydro_tests/soundwave.dat")
        self.sim.SetParam("Nhydro",64)
        self.run_id="SOUNDWAVE_SPH"
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error = 1e-4
    

    def test_error(self):
        p=run_async()
        p.wait()
        snap(-1)
        errnorm=L1errornorm("soundwave","x","rho",0.01,0.99)
        self.assertLess(errnorm,self.expected_l1error)

class SoundWaveMeshlessTest(SoundWaveTest):
    def setUp(self):
        self.sim=newsim(paramfile="tests/hydro_tests/soundwave.dat",sim='meshlessfv',ndim=1)
        self.sim.SetParam("Nhydro",64)
        self.sim.SetParam("kernel","m4")
        self.sim.SetParam("riemann_solver","hllc")
        self.run_id="SOUNDWAVE_MESHLESS"
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error = 2e-3