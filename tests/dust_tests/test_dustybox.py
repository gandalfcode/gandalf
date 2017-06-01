from gandalf.analysis.facade import *
from gandalf.analysis.data_fetcher import TimeData
from gandalf.analysis.compute import L1errornorm, lagrangian_radii
import unittest
from plot_dustybox import DriftVelocitySolutionFactory
import numpy as np



class DustyBoxTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/dust_tests/dustybox.dat")
        self.run_id="DUSTYBOX_SPH"
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error_gas = 8e-4
        self.expected_l1error_dust = 8e-4
        self.energy_error = 4.1e-4

    def check_energy_conservation(self):
        '''Check the energy error is never too large'''
        
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
        self.assertLess(errnorm_gas,self.expected_l1error_gas)
        self.assertLess(errnorm_dust,self.expected_l1error_dust)
        
        self.check_energy_conservation()
