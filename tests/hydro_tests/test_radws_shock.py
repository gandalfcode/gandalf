from gandalf.analysis.facade import *
import numpy as np
import unittest

## For this test there is no analytical solution, so we will just compare
## the meshless to SPH...
class _RadWS_Shock(unittest.TestCase):
    def setUp(self):
        self.expected_l1error_1 = 3e-3
        self.expected_l1error_2 = 3e-3


    def run_test(self, sim, params={}):
        sim = newsim(paramfile="tests/hydro_tests/sod_radws.dat", sim=sim)
        sim.SetParam('radws_table', 'eos.bell.cc.dat')
        for k, p in params.items():
            sim.SetParam(k, p)

        run_async().wait()
        snap(-1)

        x   = get_data('x')
        rho = get_data('rho')

        # Re-scale the density to get sensible error-norm estimates
        rho /= sim.simparams.floatparams['rhofluid1']
        
        return x, rho

    def interpolate(self, x1, y1, x2):
        args = np.argsort(x1)
        return np.interp(x2, x1[args], y1[args], period=4.)

    def test_error(self):
        x_sph, rho_sph = self.run_test('sph')
        x_mfm, rho_mfm = self.run_test('mfvmuscl')

        L1_1 = abs(self.interpolate(x_sph, rho_sph, x_mfm) - rho_mfm).mean()
        L1_2 = abs(self.interpolate(x_mfm, rho_mfm, x_sph) - rho_sph).mean()

        self.assertLess(L1_1, self.expected_l1error_1)
        self.assertLess(L1_2, self.expected_l1error_2)
