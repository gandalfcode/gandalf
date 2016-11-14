from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import os
import unittest

class SoundWaveTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/hydro_tests/soundwave.dat")
        self.sim.SetParam("Nhydro",64)
        self.expected_l1error = 1e-4
    

    def test_error(self):
        run()
        self.errnorm=L1errornorm("soundwave","x","rho",0.01,0.99)
        self.assertLess(self.errnorm,self.expected_l1error)

class SoundWaveMeshlessTest(SoundWaveTest):
    def setUp(self):
        self.sim=newsim(paramfile="tests/hydro_tests/soundwave.dat",sim='meshlessfv',ndim=1)
        self.sim.SetParam("Nhydro",64)
        self.sim.SetParam("kernel","m4")
        self.sim.SetParam("run_id","test")
        #self.sim.SetParam("sim","meshlessfv")
        setupsim()
        self.expected_l1error = 2e-3