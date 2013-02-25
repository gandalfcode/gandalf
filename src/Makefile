# =============================================================================
# MAKEFILE
# =============================================================================


ifneq ($(PYTHON),)
ifeq ($(NUMPY),)
NUMPY = $(shell $(PYTHON) -c "import numpy; print numpy.get_include()")
endif
ifeq ($(PYLIB),)
PYLIB = $(shell $(PYTHON) -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc()")
endif
endif

#PYLIB = /sw/include/python2.7
#NUMPY = /sw/lib/python2.7/site-packages/numpy/core/include

#location of the python files produced by swig
PYFOLDER = ../analysis/swig_generated

# Dimensionality of the code
# ----------------------------------------------------------------------------
ifeq ($(NDIM),1)
CFLAGS += -DNDIM=1 -DVDIM=1 -DFIXED_DIMENSIONS
else ifeq ($(NDIM),2)
CFLAGS += -DNDIM=2 -DVDIM=2 -DFIXED_DIMENSIONS
else ifeq ($(NDIM),3)
CFLAGS += -DNDIM=3 -DVDIM=3 -DFIXED_DIMENSIONS
else ifeq ($(NDIM),0)
CFLAGS += -DNDIM=3 -DVDIM=3
else
ERROR += "Invalid value for NDIM : "$(NDIM)"\n"
endif


# Precision options
# ----------------------------------------------------------------------------
ifeq ($(PRECISION),SINGLE)
CFLAGS += -DSINGLE_PRECISION
else ifeq ($(PRECISION),DOUBLE)
CFLAGS += -DDOUBLE_PRECISION
endif


# Debug flags
# ----------------------------------------------------------------------------
ifeq ($(DEBUG),1)
CFLAGS += -DDEBUG1
else ifeq ($(DEBUG),2)
CFLAGS += -DDEBUG1 -DDEBUG2
endif


# Include expensive verification code
# ----------------------------------------------------------------------------
ifeq ($(VERIFY_ALL),1)
CFLAGS += -DVERIFY_ALL
endif



# Object files to be compiled
# ----------------------------------------------------------------------------
SWIG_HEADERS = Parameters.i SimUnits.i Sph.i SphSnapshot.i SphSimulation.i
WRAP_OBJ = Parameters_wrap.o SimUnits_wrap.o Sph_wrap.o SphSnapshot_wrap.o SphSimulation_wrap.o
OBJ = Parameters.o SimUnits.o SphSnapshot.o SphSimulation.o
OBJ += SphSimulationIC.o SphSimulationIO.o SphSimulationTimesteps.o
OBJ += SphAnalysis.o
OBJ += M4Kernel.o QuinticKernel.o
OBJ += Sph.o GradhSph.o
OBJ += EnergyPEC.o
OBJ += SphIntegration.o SphLeapfrogKDK.o
#OBJ += SphNeighbourSearch.o 
OBJ += BruteForceSearch.o GridSearch.o
OBJ += AdiabaticEOS.o IsothermalEOS.o
OBJ += SimGhostParticles.o
OBJ += toymain.o
OBJ += Exception.o


.SUFFIXES: .cpp .i .o


%_wrap.cxx: %.i
	swig -c++ -python -outdir $(PYFOLDER) $(CFLAGS) $<

%.o: %.cxx
	$(CC) $(OPT) $(CFLAGS) -c $< -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY)

%.o: %.cpp
	$(CC) $(OPT) $(CFLAGS) -c $<

# Depending if we are on Linux or on Mac OS, sets the correct
# flags for building a dynamic library
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SHARED_OPTIONS = -bundle -flat_namespace -undefined suppress
else ifeq ($(UNAME_S),Linux)
SHARED_OPTIONS = -shared
endif

# =============================================================================
toy : $(WRAP_OBJ) $(OBJ)
	$(CC) $(CFLAGS) $(OPT) $(SHARED_OPTIONS) $(OBJ) $(WRAP_OBJ) -o _SphSim.so
	$(F2PY) --quiet -m shocktub -c shocktub.f 
	$(CC) $(CFLAGS) $(OPT) -o ../bin/toymain $(OBJ)
	mv *.so ../analysis


# =============================================================================
clean ::
	\rm -f *_wrap.cxx
	\rm -f *.o
	\rm -f *.so
	\rm -f ../analysis/*.so
	\rm -f *.pyc
