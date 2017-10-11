from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import unittest

class RadwsTest(unittest.TestCase):
    rho_c = 3.19347096432e-10
    U_c   = 8102481.56317
    tmax  = 116623.901728
    tol = 0.1

    # Just output a single snapshot
    params = { 'radws_table' : 'eos.bell.cc.dat',
               'dt_snap' : '18560' }
    def setUp(self):
        self.sim=newsim("tests/gravhydro_tests/radws_test.dat")
 
        for param in self.params:
            self.sim.SetParam(param, self.params[param])

    def test_error(self):
        p=run_async()
        p.wait()

        # Compute central density / internal energy
        tmax = snap(-1).t
        rho = get_data('rho')
        U   = get_data('u')
        
        args = rho.argsort()[-10:]
        rho_c = rho[args].mean()
        U_c   = U[args].mean()

        # Check we're at least close:
        self.assertLess( abs(1 - tmax/self.tmax)   , self.tol)
        self.assertLess( abs(1 - rho_c/self.rho_c) , self.tol)
        self.assertLess( abs(1 - U_c/self.U_c)     , self.tol)
