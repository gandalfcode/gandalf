# =============================================================================
# MAKEFILE
# =============================================================================


CC = g++
PYTHON = python2.7

OPT = -fPIC -g

OUTPUT_LEVEL              = 1
PRECISION                 = SINGLE
NDIM                      = 3
DEBUG                     = 2

#PYLIB = /sw/include/python2.7
#NUMPY = /sw/lib/python2.7/site-packages/numpy/core/include

ifneq ($(PYTHON),)
ifeq ($(NUMPY),)
NUMPY = $(shell $(PYTHON) -c "import numpy; print numpy.get_include()")
endif
ifeq ($(PYLIB),)
PYLIB = $(shell $(PYTHON) -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc()")
endif
endif


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


ifeq ($(DEBUG),1)
CFLAGS += -DDEBUG1
else ifeq ($(DEBUG),2)
CFLAGS += -DDEBUG1 -DDEBUG2
endif


SWIG_HEADERS = Parameters.i SimUnits.i Sph.i SphSnapshot.i SphSimulation.i

WRAP_OBJ = Parameters_wrap.o SimUnits_wrap.o Sph_wrap.o SphSnapshot_wrap.o SphSimulation_wrap.o

OBJ = Parameters.o SimUnits.o SphSnapshot.o SphSimulation.o SphKernel.o 
OBJ += Sph.o GradhSph.o
OBJ += SphIntegration.o SphLeapfrogKDK.o
OBJ += SphNeighbourSearch.o BruteForceSearch.o
OBJ += IsothermalEOS.o
OBJ += toymain.o


.SUFFIXES: .cpp .i .o


%_wrap.cxx: %.i
	swig -c++ -python $(CFLAGS) $<

%.o: %.cxx
	$(CC) $(OPT) $(CFLAGS) -c -fPIC $< -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY)

%.o: %.cpp
	$(CC) $(OPT) $(CFLAGS) -c -fPIC $<


# =============================================================================
toy2 : $(WRAP_OBJ) $(OBJ)
	@echo -e $(PYLIB)
	ld -bundle -flat_namespace -undefined suppress $(OBJ) $(WRAP_OBJ) -o _SphSim.so
	ld -bundle -flat_namespace -undefined suppress SphSnapshot.o SphSnapshot_wrap.o -o _SphSnap.so
	$(CC) $(CFLAGS) -o toymain $(OBJ)




# =============================================================================
toy ::
	swig -c++ -python $(CFLAGS) Parameters.i
	swig -c++ -python $(CFLAGS) SimUnits.i
	swig -c++ -python $(CFLAGS) Sph.i
	swig -c++ -python $(CFLAGS) SphSnapshot.i
	swig -c++ -python $(CFLAGS) SphSimulation.i
	$(CC) $(OPT) $(CFLAGS) SphKernel.cpp
	$(CC) $(OPT) $(CFLAGS) SphIntegration.cpp
	$(CC) $(OPT) $(CFLAGS) SphNeighbourSearch.cpp
	$(CC) $(OPT) $(CFLAGS) EOS.cpp
	$(CC) $(OPT) $(CFLAGS) toymain.cpp
	$(CC) $(OPT) $(CFLAGS) Parameters.cpp Parameters_wrap.cxx -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY)
	$(CC) $(OPT) $(CFLAGS) SimUnits.cpp SimUnits_wrap.cxx -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY)
	$(CC) $(OPT) $(CFLAGS) Sph.cpp Sph_wrap.cxx -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY)
	$(CC) $(OPT) $(CFLAGS) SphSnapshot.cpp SphSnapshot_wrap.cxx -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY)
	$(CC) $(OPT) $(CFLAGS) SphSimulation.cpp SphSimulation_wrap.cxx -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY)
	ld -bundle -flat_namespace -undefined suppress SphKernel.o EOS.o SphIntegration.o SphNeighbourSearch.o Sph.o Sph_wrap.o Parameters.o Parameters_wrap.o SimUnits.o SimUnits_wrap.o SphSnapshot.o SphSnapshot_wrap.o SphSimulation.o SphSimulation_wrap.o -o _SphSim.so
	ld -bundle -flat_namespace -undefined suppress SphSnapshot.o SphSnapshot_wrap.o -o _SphSnap.so
	$(CC) $(CFLAGS) -o toymain SphKernel.o EOS.o SphIntegration.o SphNeighbourSearch.o Sph.o Parameters.o SimUnits.o SphSnapshot.o SphSimulation.o toymain.o
#\rm -f *_wrap.cxx
#\rm -f *.o


# =============================================================================
clean ::
	\rm -f *_wrap.cxx
	\rm -f *.o
	\rm -f *.so
	\rm -f *.pyc
