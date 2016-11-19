#==================================================================================================
#  GANDALF v0.4.0 Makefile
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
#==================================================================================================


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
#GTEST = $(GTEST_DIR)
#GTEST = /Users#/david/astro/code/gtest-1.6.0/include


VPATH = $(PWD)/src:$(PWD)/src/Common:$(PWD)/src/Feedback
VPATH +=$(PWD)/src/GradhSph:$(PWD)/src/Headers:$(PWD)/src/Hydrodynamics
VPATH +=$(PWD)/src/Ic:$(PWD)/src/MeshlessFV:$(PWD)/src/Mpi:$(PWD)/src/Nbody
VPATH +=$(PWD)/src/Radiation:$(PWD)/src/SM2013:$(PWD)/src/Thermal:$(PWD)/src/Tree
VPATH +=$(PWD)/src/UnitTesting
INCLUDE += -I$(PWD)/src/Headers

ifeq ($(strip $(CPP)),)
  CPP=$(CXX)
endif

ifneq ($(strip $(MPICPP)),)
  CPP=$(MPICPP)
endif

# Compiler mode flags
# -------------------------------------------------------------------------------------------------
ifeq ($(COMPILER_MODE),FAST)
GCCOPT   += -O3 -ffast-math -fPIC -fno-exceptions -fno-rtti
ICCOPT   += -O3 -fast -fPIC -fno-exceptions
CLANGOPT += -O3 -ffast-math -fPIC
else ifeq ($(COMPILER_MODE),DEBUG)
GCCOPT   += -O3 -g -Wall -Wno-unknown-pragmas -Wno-reorder -fbounds-check -fPIC  -fno-rtti #-fprofile-arcs -ftest-coverage
ICCOPT   += -O3 -g -fPIC #-fprofile-arcs -ftest-coverage
CLANGOPT += -O3 -g -pg -Wall
else ifeq ($(COMPILER_MODE),STANDARD)
GCCOPT   += -O3 -fPIC -fno-exceptions -fno-rtti
else ifeq ($(COMPILER_MODE),PROFILE)
GCCOPT   += -O3 -fPIC -fno-exceptions -fno-rtti -pg
endif

ifeq ($(OPENMP),1)
GCCOPT   += -fopenmp
ICCOPT   += -openmp
CLANGOPT += -fopenmp -DOPENMP
endif

# Select the compiler
# -------------------------------------------------------------------------------------------------
ifeq ($(findstring g++,$(CPP)),g++)
OPT += $(GCCOPT)
else ifeq ($(CPP),icpc)
OPT += $(ICCOPT)
else ifeq ($(CPP),clang)
OPT += $(CLANGOPT)
else ifeq ($(CPP),mpic++)
MPI = 1
OPT += $(GCCOPT)
endif


ifeq ($(MPI), 1)
CFLAGS += -DMPI_PARALLEL
ifneq ($(CPP), mpic++)
LIB += -lmpi -lmpi++
endif
endif


ifeq ($(CXX),clang++)
 OPT += -D__extern_always_inline=inline
endif

# Precision options
# -------------------------------------------------------------------------------------------------
ifeq ($(PRECISION),SINGLE)
CFLAGS += -DGANDALF_SINGLE_PRECISION
else ifeq ($(PRECISION),DOUBLE)
CFLAGS += -DGANDALF_DOUBLE_PRECISION
endif


# FFTW linkage
#--------------------------------------------------------------------------------------------------
ifeq ($(FFTW),1)
CFLAGS += -DFFTW_TURBULENCE
ifeq ($(FFTW_INCLUDE),)
INCLUDE += -I/sw/include
else
INCLUDE += $(FFTW_INCLUDE)
endif
ifeq ($(FFTW_LIBRARY),)
LIB += -lfftw3
else
LIB += $(FFTW_LIBRARY)
endif
endif


# GNU scientific library linkage
#--------------------------------------------------------------------------------------------------
ifeq ($(GSL),1)
CFLAGS += -DGANDALF_GSL
ifeq ($(GSL_INCLUDE),)
INCLUDE += -I/usr/local/include
else
INCLUDE += $(GSL_INCLUDE)
endif
ifeq ($(GSL_LIBRARY),)
LIB += -L/usr/local/lib -lgsl -lgslcblas
else
LIB += $(GSL_LIBRARY)
endif
endif


# Debug output flags
# -------------------------------------------------------------------------------------------------
ifeq ($(OUTPUT_LEVEL),1)
CFLAGS += -DDEBUG1
else ifeq ($(OUTPUT_LEVEL),2)
CFLAGS += -DDEBUG1 -DDEBUG2
else ifeq ($(OUTPUT_LEVEL),3)
CFLAGS += -DDEBUG1 -DDEBUG2 -DDEBUG3
endif


# Include expensive verification code
# -------------------------------------------------------------------------------------------------
ifeq ($(DEBUG_LEVEL),1)
CFLAGS +=
else ifeq ($(DEBUG_LEVEL),2)
CFLAGS += -DVERIFY_ALL
else
CFLAGS += -DNDEBUG
endif


# Reorder particles flag
#ifeq ($(REORDER_PARTICLES),1)
#CFLAGS += -DREORDER_PARTICLES
#endif



# Object files to be compiled
# -------------------------------------------------------------------------------------------------
SWIG_HEADERS = Simulation.i
WRAP_OBJ = Simulation_wrap.o
OBJ = Parameters.o SimUnits.o
OBJ += Simulation.o Hydrodynamics.o
OBJ += SphSimulation.o Sph.o
OBJ += GradhSphSimulation.o GradhSph.o
OBJ += SM2012SphSimulation.o SM2012Sph.o
OBJ += MeshlessFVSimulation.o FV.o MeshlessFV.o MfvCommon.o
OBJ += MfvMusclSimulation.o MfvMuscl.o
OBJ += MfvRungeKuttaSimulation.o MfvRungeKutta.o
OBJ += NbodySimulation.o
OBJ += M4Kernel.o QuinticKernel.o GaussianKernel.o TabulatedKernel.o
OBJ += SphIntegration.o SphLeapfrogKDK.o SphLeapfrogDKD.o
OBJ += RiemannSolver.o
OBJ += SphNeighbourSearch.o
OBJ += HydroTree.o Tree.o KDTree.o OctTree.o BruteForceTree.o
OBJ += MeshlessFVTree.o
OBJ += GradhSphTree.o
OBJ += SM2012SphTree.o
OBJ += Ewald.o
OBJ += AdiabaticEOS.o BarotropicEOS.o Barotropic2EOS.o
OBJ += PolytropicEOS.o IsothermalEOS.o RadwsEOS.o
OBJ += IonisingRadiationEOS.o MCRadiationEOS.o
OBJ += MultipleSourceIonisation.o
OBJ += KDRadiationTree.o TreeMonteCarlo.o
OBJ += MonochromaticIonisationMonteCarlo.o
OBJ += chealpix.o TreeRay.o TreeRayOnTheSpot.o
OBJ += EnergyEquation.o EnergyRadws.o
OBJ += Nbody.o NbodyLeapfrogKDK.o NbodyLeapfrogDKD.o
OBJ += NbodyHermite4.o NbodyHermite4TS.o NbodyHermite6TS.o
OBJ += NbodySystemTree.o
OBJ += Sinks.o
OBJ += Ghosts.o
OBJ += SphSnapshot.o
OBJ += CodeTiming.o
OBJ += Dust.o
OBJ += Particle.o RandomNumber.o
OBJ += Supernova.o SupernovaDriver.o
OBJ += Ic.o BinaryAccretionIc.o BondiAccretionIc.o
OBJ += BossBodenheimerIc.o ContactDiscontinuityIc.o DustyBoxIc.o
OBJ += EvrardCollapseIc.o EwaldIc.o FilamentIc.o GreshoVortexIc.o
OBJ += HierarchicalSystemIc.o IsothermalSphereIc.o KhiIc.o NohIc.o PlummerSphereIc.o
OBJ += PolytropeIc.o RtiIc.o SedovBlastwaveIc.o ShearflowIc.o ShocktubeIc.o
OBJ += SilccIc.o SoundwaveIc.o SpitzerExpansionIc.o TurbulentCoreIc.o

ifeq ($(MPI),1)
OBJ += MpiNode.o MpiTree.o MpiControl.o MpiKDTreeDecomposition.o
endif


TEST_OBJ = TestScaling.o TestTree.o TestOctTree.o TestSinks.o

.SUFFIXES: .cpp .i .o


%_wrap.cxx: %.i
	swig -c++ -python -outdir $(PYFOLDER) $(CFLAGS) $(INCLUDE) $<

%_wrap.o: %_wrap.cxx
	$(CPP) $(OPT) $(CFLAGS) -fexceptions -frtti -c $< -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY) $(INCLUDE)

%.o: %.cxx
	$(CPP) $(CFLAGS) $(INCLUDE) -c $< $(OPT)

%.o: %.cpp
	$(CPP) $(OPT) $(CFLAGS) -I$(GTEST)/include $(INCLUDE) -c $<

%.o: %.c
	$(CPP) $(OPT) $(CFLAGS) -I$(GTEST)/include $(INCLUDE) -c $<

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SHARED_OPTIONS = -bundle -flat_namespace -undefined suppress
else ifeq ($(UNAME_S),Linux)
SHARED_OPTIONS = -shared
endif
.SECONDARY:

#==================================================================================================
all : $(WRAP_OBJ) $(OBJ) Exception.o _SphSim.so executable
	cp _SphSim.so ../analysis/swig_generated

#==================================================================================================

executable : $(OBJ) gandalf.o Exception.o
	$(CPP) $(CFLAGS) $(INCLUDE) -o gandalf $(OBJ) Exception.o gandalf.o $(OPT) $(LIB)
	cp gandalf ../bin/gandalf

unittests : $(OBJ) $(TEST_OBJ) Exception.o
	$(CPP) $(CFLAGS) $(INCLUDE) $(OPT) $(LIB) -o testgandalf $(OBJ) $(TEST_OBJ) Exception.o $(GTEST)/make/gtest_main.a
	cp testgandalf ../bin/testgandalf

_SphSim.so : $(WRAP_OBJ) $(OBJ) Exception.o Render.o
	$(CPP) $(CFLAGS) $(INCLUDE) $(OPT) $(LIB) $(SHARED_OPTIONS) $(WRAP_OBJ) $(OBJ) Exception.o Render.o -o _SphSim.so

Exception.o : Exception.cpp
	$(CPP) $(OPT) $(INCLUDE) $(CFLAGS) -fexceptions -c $<

SphSimulation.o : SimulationIO.hpp SimulationIC.hpp SimAnalysis.hpp


force: ;

clean ::
	\rm -f *_wrap.cxx
	\rm -f */*_wrap.cxx
	\rm -f *.o
	\rm -f *.so
	\rm -f ../analysis/*.so
	\rm -f ../analysis/swig_generated/*.so
	\rm -f *.pyc
	\rm -f ../analysis/swig_generated/SphSim.py
	\rm -f *.d



ifeq ($(BUILD_DEPENDENCIES), 1)
# Handle the dependencies using black magic
get_dir_str = $(patsubst %/,%,$(dir $1))
WHICH_GPP   = which g++

ifeq ($(findstring g++,$(WHICH_GPP)),g++)
BUILD_DEPS = @set -e; g++ -MM -MG $(INCLUDE) $< |\
	 sed 's;\($*\)\.o[ :]*;$(call get_dir_str, $@)/\1.o $@ : ;g' > $@
else
BUILD_DEPS = @echo "Warning: Cannot build dependencies for $(notdir $<)"
endif

else
# Ignore dependencies
BUILD_DEPS =
endif

# Compute the list of .cpp files from the OBJ files
#    Note: chealpix is a .c file, so handle the special case
CPP_TGT := $(patsubst %.o,  %.cpp, $(OBJ))
CPP_TGT := $(patsubst chealpix.cpp, chealpix.c, $(CPP_TGT))
REL_DEPS := $(patsubst %.cpp,  %.d, $(CPP_TGT))
REL_DEPS := $(patsubst %.c,  %.d, $(REL_DEPS))



%.o : %.d

.PHONY : depends tidy

depends: $(REL_DEPS)
	@echo "Building Dependencies..."

%.d: %.cpp
	$(BUILD_DEPS)

%.d: %.c
	$(BUILD_DEPS)


-include $(REL_DEPS)
