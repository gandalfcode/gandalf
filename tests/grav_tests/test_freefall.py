from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm, lagrangian_radii
from freefall import timeratiofreefall,freefall_analytical_radius
import unittest
import numpy as np



class FreeFallTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/grav_tests/freefall.dat")
        self.run_id="FREEFALL_SPH"
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error = 9e-3
    

    def test_error(self):
        p=run_async()
        fetcher_5=CreateTimeData('lr',lagrangian_radii,mfrac=0.5)
        fetcher_t=CreateTimeData('tr',timeratiofreefall)
        p.wait()
        time=fetcher_t.fetch()[1]
        lr_5=fetcher_5.fetch()[1]
        analytical_r=np.empty_like(lr_5)
        for i,t in enumerate(time):
            analytical_r[i]=freefall_analytical_radius(t)
        errnorm= np.linalg.norm(analytical_r*lr_5[0] - lr_5, ord=1)/time.size
        self.assertLess(errnorm,self.expected_l1error)

class FreeFallMeshlessTest(FreeFallTest):
    def setUp(self):
        self.sim=newsim(paramfile="tests/grav_tests/freefall.dat",sim='meshlessfv')
        self.run_id="FREEFALL_MESHLESS"
        self.sim.SetParam("run_id",self.run_id)
        self.sim.SetParam("riemann_solver",'hllc')
        self.expected_l1error = 1e-2