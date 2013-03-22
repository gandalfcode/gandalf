# =============================================================================
# Makefile
# =============================================================================

CC                 = g++-4
PYTHON             = python2.7
F2PY               = f2py2.7
COMPILER_MODE      = FAST
OPENMP             = 1

NDIM               = 0
PRECISION          = DOUBLE

OUTPUT_LEVEL       = 2
DEBUG              = 2
VERIFY_ALL         = 0


# Select location of python and numpy libraries.  If blank, make will try to 
# find the location of the libraries automatically using installed python 
# utilities.  If you have multiple versions of python installed on your 
# computer, then select the prefered version with the PYTHON variable above.
# -----------------------------------------------------------------------------
PYLIB = 
NUMPY = 


# Don't delete this command! Makes sure that defined variables here are passed
# to the others make
export

all:
	@+$(MAKE) -C src

executable:
	@+$(MAKE) executable -C src

clean:
	@+$(MAKE) clean -C src
