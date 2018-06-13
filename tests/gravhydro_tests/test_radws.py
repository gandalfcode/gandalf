from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import unittest

class _RadwsTest(unittest.TestCase):
    rho_c = 8.88283823757e-14
    U_c   = 148076.49794
    tmax  = 78931.3355517
    tol = 0.1

    # Just output a single snapshot
    sim = 'gradhsph'
    params = { 'radws_table' : 'eos.bell.cc.dat',
               'tend' : '12555', 'dt_snap' : '12555' }
    def setUp(self):
        self.sim=newsim("tests/gravhydro_tests/radws_test.dat", sim=self.sim)

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

        #print tmax, rho_c, U_c

        # Check we're at least close:
        self.assertLess( abs(1 - tmax/self.tmax)   , self.tol)
        self.assertLess( abs(1 - rho_c/self.rho_c) , self.tol)
        self.assertLess( abs(1 - U_c/self.U_c)     , self.tol)

class _Radws_mfm(_RadwsTest):
    rho_c = 2.02672009908e-13 
    U_c   = 155461.06611
    tmax  = 76101.3956448


    # Just output a single snapshot
    sim = 'meshlessfv'
    params = { 'radws_table' : 'eos.bell.cc.dat', 'courant_mult' : '0.25',
               'tend' : '12100', 'dt_snap' : '12100' }
    
