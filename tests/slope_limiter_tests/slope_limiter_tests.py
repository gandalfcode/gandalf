#==============================================================================
#!/usr/bin/python
# slope_limiter_tests.py
# Run a few of the standard tests with the meshless methods and different 
# slope limiters in order to determine their performance
#==============================================================================
import os
from gandalf.analysis.facade import *

tests = [ 'ADSOD', 'SEDOV', 'GRESHO' ]
limiters = ['zeroslope', 'tvdscalar', 'balsara2004', 
            'springel2009', 'gizmo' ]


def run_slope_limiter_test(test, slope_limiter, mfm=False):
    sim = newsim(test.lower() + '.dat')

    # Setup the run-specific paramters:
    run_id = os.path.join('results',
                          test + '-' + slope_limiter + '-')
    if mfm:
        run_id += 'mfm'
        zero_mass_flux = 1
    else:
        run_id += 'mfv'
        zero_mass_flux = 0
        
    sim.SetParam('run_id',         run_id)
    sim.SetParam('slope_limiter',  slope_limiter)
    sim.SetParam('zero_mass_flux', zero_mass_flux)
    sim.SetParam('riemann_solver', 'hllc')
    
    print 'About to run sim:', run_id

    setupsim()
    run()


for test in tests:
    for slope_limiter in limiters:
        for mfm in [True, False]:
            run_slope_limiter_test(test, slope_limiter, mfm)
