#==============================================================================
# plummertest.py
# Run the Plummer sphere test using initial conditions specified in the file 
# specified from the command line (default is plummer.dat) and then plot the 
# evolution of lagrangian radii
#==============================================================================
from gandalf.analysis.facade import *
from gandalf.analysis.compute import lagrangian_radii
import time
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("param_file",default='plummer.dat',nargs='?')
args=parser.parse_args()

# Create new shocktube simulation from 'adshock.dat' file
newsim(args.param_file)
setupsim()


run()

CreateTimeData('lr1',lagrangian_radii,mfrac=0.1)
CreateTimeData('lr2',lagrangian_radii,mfrac=0.5)
CreateTimeData('lr3',lagrangian_radii,mfrac=0.9)

time_plot("t","lr3",linestyle='-')
limit('lr3',0,4)
time_plot("t","lr2",overplot=True,linestyle='-')
time_plot("t","lr1",overplot=True,linestyle='-')

block()
