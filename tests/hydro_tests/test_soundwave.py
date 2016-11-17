from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import os
import unittest

class SoundWaveTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/hydro_tests/soundwave.dat")
        self.sim.SetParam("Nhydro",64)
        self.run_id="SOUNDWAVE_SPH"
        self.sim.SetParam("run_id",self.run_id)
        self.expected_l1error = 1e-4
        setupsim()
    

    def test_error(self):
        if self.sim.MPI:
            from mpi4py import MPI
            comm=MPI.COMM_SELF.Spawn('bin/gandalf', args=self.run_id+'.param', maxprocs=4)
            comm.Barrier()
            loadsim(self.run_id)
            snap(-1)
        else:
            run()
        errnorm=L1errornorm("soundwave","x","rho",0.01,0.99)
        self.assertLess(errnorm,self.expected_l1error)

class SoundWaveMeshlessTest(SoundWaveTest):
    def setUp(self):
        self.sim=newsim(paramfile="tests/hydro_tests/soundwave.dat",sim='meshlessfv',ndim=1)
        self.sim.SetParam("Nhydro",64)
        self.sim.SetParam("kernel","m4")
        self.run_id="SOUNDWAVE_MESHLESS"
        self.sim.SetParam("run_id",self.run_id)
        #self.sim.SetParam("sim","meshlessfv")
        setupsim()
        self.expected_l1error = 2e-3