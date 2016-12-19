from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import unittest

class JeansTest(unittest.TestCase):
    sim = 'sph'
    run_id = 'JEANS_SPH'
    expected_l1error = 6e-3
    params = {}
    def setUp(self):
        print self.sim, self.run_id
        self.sim=newsim("tests/gravhydro_tests/jeans.dat", sim=self.sim)
        self.sim.SetParam("run_id",self.run_id)
        for param in self.params:
            self.sim.SetParam(param, self.params[param])

    def test_error(self):
        p=run_async()
        p.wait()
        snap(-1)
        errnorm=L1errornorm("jeans","x","vx",0.01,0.99)
        print self.run_id, errnorm, self.expected_l1error
        self.assertLess(errnorm,self.expected_l1error)



class JeansTest_SPHRelative(JeansTest):
    run_id = "JEANS_SPH_RELATIVE"
    expected_l1error = 5.3e-3
    params = { 'relative_opening_criterion' : 1 }
    

class JeansTest_SPHEigen(JeansTest):
    run_id = "JEANS_SPH_EIGEN"
    expected_l1error = 5.3e-3
    params = { 'gravity_mac' : 'eigenmac' }

class JeansTest_Meshless(JeansTest):
    sim = 'mfvmuscl'
    run_id = "JEANS_MFM"
    expected_l1error = 1.3e-3
    params = { 'riemann_solver' : 'hllc',
               'relative_opening_criterion' : 1,
               'zero_mass_flux' : 1,
               'h_fac' : 1.0 }
