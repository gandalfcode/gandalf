from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import unittest
import numpy as np


class BinaryTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/nbody_tests/binary.dat")
    

    def test_error(self):
        p=run_async()
        p.wait()
        x0=get_data('x')
        y0=get_data('y')
        z0=get_data('z')
        snap(-1)
        x1=get_data('x')
        y1=get_data('y')
        z1=get_data('z')
        np.testing.assert_allclose(z1, np.zeros(2))
        diag=self.sim.GetDiagnostics()
        self.assertLess(diag.Eerror,1e-8)