from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import unittest

class SedovTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/hydro_tests/sedov.dat")
        self.run_id="SEDOV_SPH"
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error = 4e-2
    

    def test_error(self):
        p=run_async()
        p.wait()
        snap(-1)
        errnorm=L1errornorm("sedov","R","press",xmin=0.4,xmax=1)
        self.assertLess(errnorm,self.expected_l1error)

class SedovMeshlessTest(SedovTest):
    def setUp(self):
        self.sim=newsim(paramfile="tests/hydro_tests/sedov.dat",sim='meshlessfv')
        self.run_id="SEDOV_MESHLESS"
        self.sim.SetParam("run_id",self.run_id)
        self.sim.SetParam("riemann_solver",'hllc')
        self.expected_l1error = 3e-2