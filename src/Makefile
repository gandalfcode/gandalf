#==============================================================================
#  GANDALF v0.1.0 Makefile
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


ifneq ($(PYTHON),)
ifeq ($(NUMPY),)
NUMPY = $(shell $(PYTHON) -c "import numpy; print numpy.get_include()")
endif
ifeq ($(PYLIB),)
PYLIB = $(shell $(PYTHON) -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc()")
endif
endif


#location of the python files produced by swig
PYFOLDER = ../analysis/swig_generated
GTEST = /Users#/david/astro/code/gtest-1.6.0/include


# Compiler mode flags
# ----------------------------------------------------------------------------
ifeq ($(CPP),g++)
ifeq ($(COMPILER_MODE),FAST)
OPT += -O3 -ffast-math -fPIC -fno-exceptions -fno-rtti
else ifeq ($(COMPILER_MODE),STANDARD)
OPT += -O3 -fPIC -fno-exceptions -fno-rtti
else ifeq ($(COMPILER_MODE),PROFILE)
OPT += -O3 -fPIC -fno-exceptions -fno-rtti -pg
else ifeq ($(COMPILER_MODE),DEBUG)
OPT += -g -Wall -fbounds-check -fPIC -O3  #-fprofile-arcs -ftest-coverage
endif
ifeq ($(OPENMP),1)
OPT += -fopenmp 
endif
endif

ifeq ($(CPP),icpc)
ifeq ($(COMPILER_MODE),FAST)
OPT += -O3 -fast -fPIC -fno-exceptions
else ifeq ($(COMPILER_MODE),DEBUG)
OPT += -g -fbounds-check -fPIC -O3  #-fprofile-arcs -ftest-coverage
endif
ifeq ($(OPENMP),1)
OPT += -openmp
endif
endif

ifeq ($(CPP),mpic++)
MPI = 1
CFLAGS += -DMPI_PARALLEL
ifeq ($(COMPILER_MODE),FAST)
OPT += -O3 -ffast-math -fPIC -fno-exceptions -fno-rtti
else ifeq ($(COMPILER_MODE),STANDARD)
OPT += -O3 -fPIC -fno-exceptions -fno-rtti
else ifeq ($(COMPILER_MODE),PROFILE)
OPT += -O3 -fPIC -fno-exceptions -fno-rtti -pg
else ifeq ($(COMPILER_MODE),DEBUG)
OPT += -g -Wall -fbounds-check -fPIC -O3  #-fprofile-arcs -ftest-coverage
endif
ifeq ($(OPENMP),1)
OPT += -fopenmp 
endif
endif

ifeq ($(CPP),g++-4)
ifeq ($(COMPILER_MODE),FAST)
OPT += -O3 -ffast-math -fPIC -fno-exceptions -fno-rtti
else ifeq ($(COMPILER_MODE),STANDARD)
OPT += -O3 -fPIC -fno-exceptions -fno-rtti
else ifeq ($(COMPILER_MODE),PROFILE)
OPT += -O3 -fPIC -fno-exceptions -fno-rtti -pg
else ifeq ($(COMPILER_MODE),DEBUG)
OPT += -g -Wall -fbounds-check -fPIC -O3 #-fprofile-arcs -ftest-coverage
endif
ifeq ($(OPENMP),1)
OPT += -fopenmp -DOPENMP
endif
endif


ifeq ($(CPP),clang)
ifeq ($(COMPILER_MODE),FAST)
OPT += -O3 -ffast-math -fPIC
else ifeq ($(COMPILER_MODE),DEBUG)
OPT += -g -pg -Wall -fbounds-check
endif
ifeq ($(OPENMP),1)
OPT += -fopenmp -DOPENMP
endif
endif


# Precision options
# ----------------------------------------------------------------------------
ifeq ($(PRECISION),SINGLE)
CFLAGS += -DGANDALF_SINGLE_PRECISION
else ifeq ($(PRECISION),DOUBLE)
CFLAGS += -DGANDALF_DOUBLE_PRECISION
endif


# FFTW linkage
#-----------------------------------------------------------------------------
ifeq ($(FFTW),1)
CFLAGS += -DFFTW_TURBULENCE
#OPT += -I/Users/davidhubber/astro/code/gandalf/src -L/Users/davidhubber/astro/code/gandalf/src -lm -lfftw3f
#OPT += -I/usr/local/include -L/usr/local/lib -lfftw3
OPT += -I/sw/include -L/sw/lib -lfftw3
endif


# Debug output flags
# ----------------------------------------------------------------------------
ifeq ($(OUTPUT_LEVEL),1)
CFLAGS += -DDEBUG1
else ifeq ($(OUTPUT_LEVEL),2)
CFLAGS += -DDEBUG1 -DDEBUG2
else ifeq ($(OUTPUT_LEVEL),3)
CFLAGS += -DDEBUG1 -DDEBUG2 -DDEBUG3
endif


# Include expensive verification code
# ----------------------------------------------------------------------------
ifeq ($(DEBUG_LEVEL),1)
CFLAGS +=
else ifeq ($(DEBUG_LEVEL),2)
CFLAGS += -DVERIFY_ALL
else
CFLAGS += -DNDEBUG
endif


# Reorder particles flag
ifeq ($(REORDER_PARTICLES),1)
CFLAGS += -DREORDER_PARTICLES
endif


#CFLAGS += -DVERIFY_BINARY_TREE


# Object files to be compiled
# ----------------------------------------------------------------------------
SWIG_HEADERS = Simulation.i
WRAP_OBJ = Simulation_wrap.o
OBJ = Parameters.o SimUnits.o 
OBJ += Simulation.o 
OBJ += SphSimulation.o Sph.o
OBJ += GradhSphSimulation.o GradhSph.o
OBJ += GodunovSphSimulation.o GodunovSph.o
OBJ += SM2012SphSimulation.o SM2012Sph.o
OBJ += NbodySimulation.o
OBJ += M4Kernel.o QuinticKernel.o GaussianKernel.o TabulatedKernel.o
OBJ += SphIntegration.o SphLeapfrogKDK.o SphLeapfrogDKD.o 
OBJ += SphGodunovIntegration.o EnergyGodunovIntegration.o RiemannSolver.o
#OBJ += SphNeighbourSearch.o 
OBJ += SphNeighbourSearch.o BruteForceSearch.o BinaryTree.o
OBJ += AdiabaticEOS.o IsothermalEOS.o BarotropicEOS.o Barotropic2EOS.o
OBJ += EnergyPEC.o EnergyLeapfrogDKD.o
#OBJ += SimGhostParticles.o
#OBJ += Render.o
OBJ += Nbody.o NbodyLeapfrogKDK.o NbodyLeapfrogDKD.o
OBJ += NbodyHermite4.o NbodyHermite4TS.o
OBJ += NbodySystemTree.o
OBJ += Sinks.o
OBJ += Ghosts.o
OBJ += SphSnapshot.o
OBJ += CodeTiming.o

ifeq ($(MPI),1)
OBJ += MpiNode.o MpiControl.o
endif


TEST_OBJ = #TestScaling.o Parameters.o SimUnits.o Exception.o

.SUFFIXES: .cpp .i .o


%_wrap.cxx: %.i
	swig -c++ -python -outdir $(PYFOLDER) $(CFLAGS) $<

%_wrap.o: %_wrap.cxx
	$(CPP) $(OPT) $(CFLAGS) -fexceptions -frtti -c $< -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY)

%.o: %.cxx
	$(CPP) $(OPT) $(CFLAGS) -c $<

%.o: %.cpp
	$(CPP) $(OPT) $(CFLAGS) -c $< -I$(GTEST)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SHARED_OPTIONS = -bundle -flat_namespace -undefined suppress
else ifeq ($(UNAME_S),Linux)
SHARED_OPTIONS = -shared
endif
.SECONDARY:

# =============================================================================
all : $(TEST_OBJ) $(WRAP_OBJ) $(OBJ) Exception.o shocktub.so _SphSim.so executable
#uses the interface file rather than directly parsing the fortran file
#if it's needed, interface file can be generated by running
#	f2py -m shocktub -h shocktube.pyf shocktub.f 
#$(CPP) $(CFLAGS) $(OPT) -o ../bin/testing $(TEST_OBJ) /Users/david/astro/code/gtest-1.6.0/make/gtest_main.a
	cp shocktub.so ../analysis
	cp _SphSim.so ../analysis/swig_generated

# =============================================================================

executable : $(OBJ) gandalf.o Exception.o
	$(CPP) $(CFLAGS) $(OPT) -o gandalf $(OBJ) Exception.o gandalf.o
	cp gandalf ../bin/gandalf

_SphSim.so : $(WRAP_OBJ) $(OBJ) Exception.o Render.o
	$(CPP) $(CFLAGS) $(OPT) $(SHARED_OPTIONS) $(WRAP_OBJ) $(OBJ) Exception.o Render.o -o _SphSim.so

shocktub.so : shocktub.f shocktub.pyf
	$(F2PY) --quiet -c shocktub.f shocktub.pyf

Exception.o : Exception.cpp
	$(CPP) $(OPT) $(CFLAGS) -fexceptions -c $<

SphSimulation.o : SimulationIO.hpp SimulationIC.hpp SimAnalysis.hpp


force: ;

clean ::
	\rm -f *_wrap.cxx
	\rm -f *.o
	\rm -f *.so
	\rm -f ../analysis/*.so
	\rm -f ../analysis/swig_generated/*.so
	\rm -f *.pyc
	\rm -f ../analysis/swig_generated/SphSim.py
