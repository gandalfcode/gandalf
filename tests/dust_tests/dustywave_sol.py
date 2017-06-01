import numpy as np
from scipy.integrate import ode

class _DustyWaveSolution(object):
    def __init__(self, t, k,
                 rho_g, rho_d, drho_g, drho_d, v_g, v_d):
        self._t = t
        self._k = k
        self._rho_g = rho_g
        self._rho_d = rho_d
        self._drho_g = drho_g
        self._drho_d = drho_d
        self._v_g = v_g
        self._v_d = v_d


    def v_gas(self, x):
        return (self._v_g * np.exp(1j*self._k*x)).real
    def v_dust(self, x):
        return (self._v_d * np.exp(1j*self._k*x)).real

    def rho_gas(self, x):
        return (self._rho_g + self._drho_g * np.exp(1j*self._k*x)).real
    def rho_dust(self, x):
        return (self._rho_d + self._drho_d * np.exp(1j*self._k*x)).real

    @property
    def time(self):
        return self._t

class DustyWaveSolver(object):

    def __init__(self, rho_g=1., rho_d=1., cs=1., K=1., delta=1e-3,
                 wavelength=1., feedback=True):

        self.rho_g = rho_g
        self.rho_d = rho_d
        self.cs = cs
        self.K = K
        self.delta = delta
        self.wavelength=wavelength
        self.feedback=feedback

    def _solve_system(self, times):
        '''Solve the dusty-wave problem up to specified times'''
        k = 2*np.pi / self.wavelength
        cs2 = self.cs**2

        ts_inv = self.K * self.cs * self.rho_g
        if self.feedback:
            Kg = (self.rho_d/self.rho_g) * ts_inv
            Kd = ts_inv
        else:
            Kg = 0
            Kd = ts_inv

        def f(t, y,_=None):
            drho_g, v_g = y[:2]
            drho_d, v_d = y[2:]
            
            dydt = np.empty([4], dtype='c8')
            dydt[0] = - 1j*k*self.rho_g * v_g
            dydt[2] = - 1j*k*self.rho_d * v_d
            dydt[1] = - Kg*(v_g - v_d) - 1j*k*drho_g * cs2
            dydt[3] =   Kd*(v_g - v_d)
            
            return dydt

        _jac = np.zeros([4,4], dtype='c8')
        _jac[0,1] = -1j*k*self.rho_g 
        _jac[2,3] = -1j*k*self.rho_d
        _jac[1,:] = [-1j*k * cs2, -Kg, 0,  Kg]
        _jac[3,:] = [0          ,  Kd, 0, -Kd]
        def jac(t, y,_=None):
            return _jac

        # Do the right going part of the wave (we can get the left going wave
        # for free by taking the real part).
        e = -1j * self.delta
        IC = np.array([self.rho_g*e, self.cs*e, 
                       self.rho_d*e, self.cs*e], dtype='c8')
        integ = ode(f, jac).set_integrator('zvode', method='bdf', 
                                           rtol=1e-12,atol=1e-12)
        integ.set_initial_value(IC, 0).set_f_params(None).set_jac_params(None)

        sol = []
        for ti in times:
            if ti > integ.t: integ.integrate(ti)
            sol.append(integ.y)
            if not integ.successful(): raise RuntimeError('Integ failed')

        return np.array(sol)
    

    def __call__(self, time):
        '''Solve the dustwave problem at a given time or list of times'''
        try:
            iter(time)
        except TypeError:
            time = [time,]
        time = sorted(time)
        sol = self._solve_system(list(time))

        def _create_solution(t, drho_g, v_g, drho_d, v_d):
            return _DustyWaveSolution(t, 2*np.pi/self.wavelength,
                                      self.rho_g, self.rho_d,
                                      drho_g, drho_d, v_g, v_d)

        return [_create_solution(t, *s) for t, s in zip(time,sol)]

def DustyWave(sim, time):
    '''(Semi)-Analytical solution for the DustyWave problem'''
    # Extract the parameters
    simfloatparams = sim.simparams.floatparams
    simstringparams = sim.simparams.stringparams
    rho_g  = simfloatparams["rhofluid1"]
    delta  = simfloatparams["amp"]
    xL     = simfloatparams["boxmin[0]"]
    xR     = simfloatparams["boxmax[0]"]
    if sim.simparams.stringparams["gas_eos"] == "isothermal":
        temp0  = simfloatparams["temp0"]
        mu_bar = simfloatparams["mu_bar"]
        csound = np.sqrt(temp0/mu_bar)
    else:
        gamma = simfloatparams["gamma_eos"]
        press  = simfloatparams["press1"]
        csound = np.sqrt(gamma*press/rho_g)
    wlambda = xR - xL

    eps = simfloatparams["dust_mass_factor"]
    K   = simfloatparams["drag_coeff"]
    feedback = simstringparams["dust_forces"] == "full_twofluid"
    assert(simstringparams["drag_law"] == "epstein" or
           simstringparams["drag_law"] == "LP2012")
    
    # Approximate the epstein drag coefficient
    if simstringparams["drag_law"] == "LP2012":
        K *= csound / (rho_g*rho_g*eps)

    solver = DustyWaveSolver(rho_g, rho_g*eps, csound, K, delta,
                             wlambda, feedback)
    return solver(time)[0]
        

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    
    dw = DustyWaveSolver()
    
    N = 101
    times = np.linspace(0, 5, N)
    x     = np.linspace(0, 1, 1001)
    
    sol = dw(times)

    f, ax = plt.subplots(2,1)
    l1, = ax[0].plot(x, sol[0].rho_gas(x),  'k-')     
    l2, = ax[0].plot(x, sol[0].rho_dust(x), 'k:')
    l3, = ax[1].plot(x, sol[0].v_gas(x),  'k-') 
    l4, = ax[1].plot(x, sol[0].v_dust(x), 'k:')


    ax[0].set_ylabel(r'$\rho$')
    ax[0].set_ylim(0.99,1.01)

    ax[1].set_xlabel(r'$x$')
    ax[1].set_ylabel(r'$v$')
    ax[1].set_ylim(-0.001,0.001)

    def animate(i):
        l1.set_ydata(sol[i].rho_gas(x))
        l2.set_ydata(sol[i].rho_dust(x))
        l3.set_ydata(sol[i].v_gas(x))
        l4.set_ydata(sol[i].v_dust(x))
        return l1,l2,l3,l4

    def init():
        return l1,l2,l3,l4

    ani = animation.FuncAnimation(f, animate, np.arange(0, N), 
                                  interval=50, blit=True)


    plt.show()
    
