from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import TimeData
from gandalf.analysis.compute import L1errornorm, lagrangian_radii
import unittest
from plot_dustybox import DriftVelocitySolutionFactory
import numpy as np



class DustyBoxTest(unittest.TestCase):
    params = {}
    run_id="DUSTYBOX_SPH"
    expected_l1error_gas = 8e-4
    expected_l1error_dust = 8e-4
    energy_error = 4.1e-4
    param_file = "tests/dust_tests/dustybox.dat"
    def setUp(self):
        self.sim=newsim(self.param_file)
        self.sim.SetParam("run_id",self.run_id)
        for p in self.params:
            self.sim.SetParam(p, self.params[p])

    def check_energy_conservation(self):
        '''Check the energy error is never too large'''

        # Skip for test-particle, which doesn't conserve energy
        if self.energy_error is None:
            return
        
        t  = []
        Ek = []
        U  = []
        s = snap(0)
        while True:
            t.append(s.t)
            
            m = get_data('m', type='sph')
            U.append((m*get_data('u', type='sph')).sum())

            _E = 0
            for tp in ['sph', 'dust']:
                m = get_data('m', type=tp)
                for j in range(self.sim.simparams.intparams['ndim']):
                    _E += (0.5*m*get_data('v'+chr(ord('x')+j),type=tp)**2).sum()
            Ek.append(_E)

            try:
                s = next()
            except:
                break
        
        Etot = np.array(Ek) + np.array(U)
        
        self.assertLess(max(abs(1 - Etot/Etot[0])), self.energy_error)


    def test_error(self):
        p=run_async()
        fetcher_t=TimeData('t')
        vx_fetcher=TimeData('vx',id=0)
        sol = DriftVelocitySolutionFactory()
        p.wait()
        time=fetcher_t.fetch()[1]
        vx_gas=vx_fetcher.fetch(type='sph')[1]
        vx_dust=vx_fetcher.fetch(type='dust')[1]
        
        analytical_gas=np.empty_like(vx_gas)
        analytical_dust=np.empty_like(vx_dust)
        for i,t in enumerate(time):
            analytical_gas[i]=sol.vg(t)
            analytical_dust[i]=sol.vd(t)
        errnorm_gas= np.linalg.norm(analytical_gas-vx_gas, ord=1)/time.size
        errnorm_dust= np.linalg.norm(analytical_dust-vx_dust, ord=1)/time.size
        print(errnorm_gas,self.expected_l1error_gas)
        print(errnorm_dust,self.expected_l1error_dust)
        self.assertLess(errnorm_gas,self.expected_l1error_gas)
        self.assertLess(errnorm_dust,self.expected_l1error_dust)
        self.check_energy_conservation()

        
class DustyBoxTestParticle(DustyBoxTest):
    run_id = "DUSTYBOX_SPH2"
    params = {"dust_forces" : "test_particle"}
    expected_l1error_gas  = 1e-12
    expected_l1error_dust = 9e-4
    energy_error = None


class DustyBoxTestStrong(DustyBoxTest):
    run_id = "DUSTYBOX_SPH_STRONG"
    expected_l1error_gas  = 6e-4
    expected_l1error_dust = 6e-4
    energy_error = 0.05
    params = {"drag_coeff" : 100 }

class DustyBoxTestStrongTS(DustyBoxTest):
    run_id = "DUSTYBOX_SPH_STRONG2"
    expected_l1error_gas  = 1e-12
    expected_l1error_dust = 1.4e-3
    energy_error = None
    params = {"drag_coeff" : 100,
              "dust_forces" : "test_particle" }
    



### Test the meshless
class DustyBoxMfv(DustyBoxTest):
    param_file = "tests/dust_tests/dustybox_meshless.dat"
    run_id = "DUSTYBOX_MFV"
    expected_l1error_gas  = 5e-4
    expected_l1error_dust = 5e-4
    energy_error = 1e-15


class DustyBoxMfvTestParticle(DustyBoxTest):
    param_file = "tests/dust_tests/dustybox_meshless.dat"
    run_id = "DUSTYBOX_MFV2"
    expected_l1error_gas  = 1e-10
    expected_l1error_dust = 1e-13
    energy_error = None
    params = {"dust_forces" : "test_particle"}


class DustyBoxMfvStrong(DustyBoxTest):
    param_file = "tests/dust_tests/dustybox_meshless.dat"
    run_id = "DUSTYBOX_MFV"
    expected_l1error_gas  = 2e-4
    expected_l1error_dust = 4e-5
    energy_error = 1e-15
    params = {"drag_coeff" : 100}

class DustyBoxMfvTestParticleStrong(DustyBoxTest):
    param_file = "tests/dust_tests/dustybox_meshless.dat"
    run_id = "DUSTYBOX_MFV2"
    expected_l1error_gas  = 3e-10
    expected_l1error_dust = 3e-11
    energy_error = None
    params = {"dust_forces" : "test_particle",
              "drag_coeff" : 100 }
