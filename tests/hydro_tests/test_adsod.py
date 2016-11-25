from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import os
import unittest

class AdSodTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/hydro_tests/adsod.dat")
        self.run_id="ADSOD_SPH"
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error = 9e-3
    

    def test_error(self):
        p=run_async()
        p.wait()
        loadsim(self.run_id)
        snap(-1)
        errnorm=L1errornorm("shocktube","x","vx",-1.0,1.0)
        self.assertLess(errnorm,self.expected_l1error)

class AdSodMeshlessTest(AdSodTest):
    def setUp(self):
        self.sim=newsim(paramfile="tests/hydro_tests/adsod.dat",sim='meshlessfv',ndim=1)
        self.run_id="ADSOD_MESHLESS"
        self.sim.SetParam("run_id",self.run_id)
        self.sim.SetParam("riemann_solver",'hllc')
        self.expected_l1error = 7e-3