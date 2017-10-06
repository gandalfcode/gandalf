#==============================================================================
#  plot_dustywave.py
#  Plot the gas and dust velocity for the final snapshot in the DUSTWAVE
#  problem.
#
#  This file is part of GANDALF :
#  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
#  https://github.com/gandalfcode/gandalf
#  Contact : gandalfcode@gmail.com
#
#  Copyright (C) 2013  D. A. Hubber, G. Rosotti
#
#  GANDALF is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  GANDALF is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License (http://www.gnu.org/licenses) for more details.
#==============================================================================
from gandalf.analysis.facade import *
import sys ; sys.path.append("../dust_tests/")
from dustywave_sol import DustyWave
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter



def power_formatter(x, p):
    tstr = '%e'%x
    e = tstr.find('e')
    p = tstr.find('.')

    pw = int(tstr[e+1:])
    if pw:               
        t_str = r'${}\times 10^{{{}}}$'.format(tstr[:p],pw)
    else:
        t_str = r'${}$'.format(tstr[:p])
        
    return t_str
                          
# Set up sims with both SPH and MFV

Ks = [1.0, 10.0, 100.0, 1000.0]
schemes = ['sph', 'mfvmuscl']
h_fac   = [1.2, 1.0]

f, subs = plt.subplots(len(Ks), len(schemes), sharex=True, sharey=True)

for j, [s, h] in enumerate(zip(schemes, h_fac)):
    for i, K in enumerate(Ks):
        print  K
        sim = newsim("dustywave.dat", sim=s)
        sim.SetParam("drag_coeff", K)
        sim.SetParam("h_fac", h)
        setupsim()
        run()

        t = snap(-1).t

        # Set up the analytical solution
        sol = DustyWave(sim, t)
    
        # Load the data
        x_gas   = get_data('x',  type='sph')
        vx_gas  = get_data('vx', type='sph')
        x_dust  = get_data('x',  type='dust')
        vx_dust = get_data('vx', type='dust')
 
        Ng, Nd = len(x_gas), len(x_dust)


        # Plot the numerical and analytical solutions
        sub = subs[i][j]
        sub.plot(x_dust, vx_dust,            'rx', markerfacecolor='None')
        sub.plot(x_gas, vx_gas,           'k.', markerfacecolor='None')

        sub.plot(x_gas, sol.v_gas(x_gas), 'k' )
        sub.plot(x_dust, sol.v_dust(x_dust), 'r--' )

        if i+1 == len(Ks): sub.set_xlabel('$x$')
        if j == 0:         sub.set_ylabel('$v_x$')

        sub.set_ylim(-1.2e-4, 1.2e-4)
        sub.yaxis.set_major_formatter(FuncFormatter(power_formatter))

        sub.text(0.0, -1e-4,'$K={}$'.format(K),
                 verticalalignment='center', horizontalalignment='left')

subs[0][0].text(1.0, 1e-4, 'SPH',
                verticalalignment='center', horizontalalignment='right')

subs[0][1].text(1.0, 1e-4, 'MFM',
                verticalalignment='center', horizontalalignment='right')

plt.subplots_adjust(top=0.99, right=0.98,
                    left=0.165, bottom=0.1,
                    hspace=0.02, wspace=0.05)
plt.savefig("dustywave.pdf")
plt.show()
