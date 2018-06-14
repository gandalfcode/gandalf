from gandalf.analysis.facade import *
from gandalf.analysis.compute import L1errornorm
import unittest

####
# We use these tests to check the tree stocking, so set ntreebuildstep = 1024.
# 
# Check this for both KD and Brute force trees.

class AdSodTest(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/hydro_tests/adsod.dat")
        self.run_id="ADSOD_SPH"
        self.sim.SetParam("run_id",self.run_id)
        self.sim.SetParam("ntreebuildstep", 1024)
        self.sim.SetParam("pruning_level_min", 1)
        self.sim.SetParam("neib_search", "kdtree")

        self.expected_l1error = 9e-3
    

    def test_error(self):
        p=run_async()
        p.wait()
        snap(-1)
        errnorm=L1errornorm("shocktube","x","vx",-1.0,1.0)
        self.assertLess(errnorm,self.expected_l1error)

class AdSodMeshlessTest(AdSodTest):
    def setUp(self):
        self.sim=newsim(paramfile="tests/hydro_tests/adsod.dat",sim="meshlessfv",ndim=1)
        self.run_id="ADSOD_MESHLESS"
        self.sim.SetParam("run_id",self.run_id)
        self.sim.SetParam("ntreebuildstep", 1024)
        self.sim.SetParam("pruning_level_min", 1)
        self.sim.SetParam("neib_search", "kdtree")
        self.sim.SetParam("riemann_solver","hllc")
        self.expected_l1error = 7e-3


class AdSodTest_BFTree(unittest.TestCase):
    def setUp(self):
        self.sim=newsim("tests/hydro_tests/adsod.dat")
        self.run_id="ADSOD_SPH"
        self.sim.SetParam("run_id",self.run_id)
        self.sim.SetParam("ntreebuildstep", 1024)
        self.sim.SetParam("pruning_level_min", 1)
        self.sim.SetParam("neib_search", "bruteforce")
        self.expected_l1error = 9e-3
    

    def test_error(self):
        p=run_async()
        p.wait()
        snap(-1)
        errnorm=L1errornorm("shocktube","x","vx",-1.0,1.0)
        self.assertLess(errnorm,self.expected_l1error)

class AdSodMeshlessTest_BFTree(AdSodTest):
    def setUp(self):
        self.sim=newsim(paramfile="tests/hydro_tests/adsod.dat",sim="meshlessfv",ndim=1)
        self.run_id="ADSOD_MESHLESS"
        self.sim.SetParam("run_id",self.run_id)
        self.sim.SetParam("ntreebuildstep", 1024)
        self.sim.SetParam("pruning_level_min", 1)
        self.sim.SetParam("neib_search", "bruteforce")
        self.sim.SetParam("riemann_solver","hllc")
        self.expected_l1error = 7e-3
