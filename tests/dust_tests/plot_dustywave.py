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
from dustywave_sol import DustyWave
import numpy as np
import matplotlib.pyplot as plt
import sys


    
if __name__=="__main__":
    # Load the sim
    if len(sys.argv) > 1:
        sim = loadsim(sys.argv[1])
    else:
        sim = loadsim('DUSTYWAVE')

    s = snap(-1)

    # Set up the analytical solution
    sol = DustyWave(sim, s.t)
    
    # Load the data
    x_gas   = get_data('x',  type='sph')
    vx_gas  = get_data('vx', type='sph')
    x_dust  = get_data('x',  type='dust')
    vx_dust = get_data('vx', type='dust')
 
    Ng, Nd = len(x_gas), len(x_dust)
    errnorm_gas  = np.linalg.norm(sol.v_gas(x_gas)  -vx_gas,  ord=1) / Ng
    errnorm_dust = np.linalg.norm(sol.v_dust(x_dust)-vx_dust, ord=1) / Nd
    print 'Error norms (gas, dust):', errnorm_gas, errnorm_dust


    # Plot the numerical and analytical solutions
    plt.plot(x_gas, vx_gas,           'k.')
    plt.plot(x_gas, sol.v_gas(x_gas), 'k' )

    plt.plot(x_dust, vx_dust,            'r.')
    plt.plot(x_dust, sol.v_dust(x_dust), 'r' )
    
    plt.xlabel('x')
    plt.ylabel('v_x')

    plt.show()
