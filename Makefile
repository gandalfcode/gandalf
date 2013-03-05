# -------------------------------------------------
# SEREN Makefile
# Here you can set the compile-time options
# -------------------------------------------------

#-------------------------------------------------
# Compiler options
#-------------------------------------------------

CC = g++ # C++ compiler
F2PY = f2py2.7
#OPT = -pg -O3 -fPIC
OPT = -O3 -ffast-math -fPIC #-g -Wall
#OPT = -g -pg -fprofile-arcs -ftest-coverage -fPIC
#OPT = -g -Wall -fbounds-check

PYTHON = python2.7 # Name of the python interpreter


#-------------------------------------------------
# Compile time options
#-------------------------------------------------
OUTPUT_LEVEL              = 1
PRECISION                 = SINGLE
# If set to 0, the number of dimensions can be set at runtime
NDIM                      = 0
DEBUG                     = 1
# Turn on expensive verifications (only needed to debug)
VERIFY_ALL                = 0


# Select location of python and numpy libraries.  If blank, make will try to 
# find the location of the libraries automatically using installed python 
# utilities.  If you have multiple versions of python installed on your 
# computer, then select the prefered version with the PYTHON variable above.
# -----------------------------------------------------------------------------
PYLIB =
NUMPY = 



# Don't delete this command! Makes sure that the variables defined here are passed
# to the others make
export

all:
	@+$(MAKE) -C src

clean:
	@+$(MAKE) clean -C src
