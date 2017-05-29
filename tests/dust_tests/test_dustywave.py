from gandalf.analysis.facade import *
import unittest
from dustywave_sol import DustyWave
import numpy as np
import matplotlib.pyplot as plt


class DustyWaveTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/dust_tests/dustywave.dat")
        self.run_id="DUSTYWAVE_SPH"
        self.sim.SetParam("drag_law", "epstein")
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error_gas  = 4.1e-6
        self.expected_l1error_dust = 2.3e-5
    

    def test_error(self):
        p=run_async()
        p.wait()
        s = snap(-1)
        
        x_gas   = get_data('x',  type='sph')
        vx_gas  = get_data('vx', type='sph')
        x_dust  = get_data('x',  type='dust')
        vx_dust = get_data('vx', type='dust')

        sol = DustyWave(self.sim, s.t)

        Ng, Nd = len(x_gas), len(x_dust)
        errnorm_gas  = np.linalg.norm(sol.v_gas(x_gas)  -vx_gas,  ord=1) / Ng
        errnorm_dust = np.linalg.norm(sol.v_dust(x_dust)-vx_dust, ord=1) / Nd

        self.assertLess(errnorm_gas,self.expected_l1error_gas)
        self.assertLess(errnorm_dust,self.expected_l1error_dust)


class DustyWaveTestParticleTest(DustyWaveTest):
    def setUp(self):
        self.sim=newsim("tests/dust_tests/dustywave.dat")
        self.sim.SetParam("dust_forces", "test_particle")
        self.sim.SetParam("drag_law", "epstein")
        self.run_id="DUSTYWAVE2_SPH"
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error_gas  = 2.4e-6
        self.expected_l1error_dust = 5.4e-7


class DustyWaveTestMeshless(DustyWaveTest):
    def setUp(self):
        self.sim=newsim("tests/dust_tests/dustywave_meshless.dat")
        self.sim.SetParam("drag_law", "epstein")
        self.run_id="DUSTYWAVE_MFV"
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error_gas  = 1.6e-5
        self.expected_l1error_dust = 5.7e-7

