from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import unittest

class SoundWaveTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/nbody_tests/nbodyjeans.dat")
        self.expected_l1error = 8e-5
    

    def test_error(self):
        p=run_async()
        p.wait()
        snap(-1)
        errnorm=L1errornorm("jeans","x","vx",0.01,0.99,type='star')
        self.assertLess(errnorm,self.expected_l1error)