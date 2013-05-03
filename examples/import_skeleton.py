from seren.analysis.facade import *

paramfile='testimport.dat' #change this to your needs
sim=newsim(paramfile)

#you can still change some parameters by calling SetParam if you want

#before uploading the initial condition, you need to call PreSetupForPython
#if you forget, you will get an error
#After calling this function, you can no longer change the parameters
#(if you try, you will get an error)
sim.PreSetupForPython()


#---------------------------------
# DO YOUR INITIALIZATION HERE

# To import an array, you have to do like that:
# sim.ImportArray(array, string),
# where array is a numpy array and string defines the quantity
# that you are importing (e.g., 'x' or 'vy').
# At minimum, you need to import the coordinates, the mass arrays
# and the internal energies
# (quantities not imported are set to zero)
#--------------------------------

#Once you are finished, call setup
#If you forget, run will do it for you, but you can't
#do plots before calling run (if you try, you will get
#an error)
setupsim()

#You can now do plots, that will be updated as the simulation runs
plot('x','y')


#Now you can call run
run()
#Block does not exit when the script ends
block()

