from gandalf.analysis.facade import *
import numpy as np
import matplotlib.pyplot as plt

def run_test(sim, **params):
    sim = newsim(paramfile="sod_radws.dat", sim=sim)
    for k, p in params.items():
        sim.SetParam(k, p)

    run_async().wait()
    snap(-1)
    
    x   = get_data('x')
    rho = get_data('rho')

    # Re-scale the density to get sensible error-norm estimates
    rho /= sim.simparams.floatparams['rhofluid1']
        
    return x, rho

def interpolate(x1, y1, x2):
    args = np.argsort(x1)
    return np.interp(x2, x1[args], y1[args], period=4.)



x_SPH, rho_SPH = run_test('sph')
x_MFM, rho_MFM = run_test('mfvmuscl')
x_AD , rho_AD  = run_test('sph', gas_eos='energy_eqn')

xi = np.linspace(-2, 2, 10**4)

l, = plt.plot(x_SPH, rho_SPH, '.', label='SPH')
plt.plot(xi, interpolate(x_SPH, rho_SPH, xi), '--', c=l.get_color())

l, = plt.plot(x_MFM, rho_MFM, '.', label='MFM')
plt.plot(xi, interpolate(x_MFM, rho_MFM, xi), '--', c=l.get_color())

l, = plt.plot(x_AD, rho_AD, 'k.', label='ideal gas, SPH')
plt.plot(xi, interpolate(x_AD, rho_AD, xi), '--', c=l.get_color())



plt.legend(loc='best')
plt.show()
