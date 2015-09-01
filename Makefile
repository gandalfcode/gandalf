#==================================================================================================
#  GANDALF v0.3.0 Makefile frontend
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


CPP                = g++-5
PYTHON             = python2.7
F2PY               = f2py2.7
FFTW               = 1
GSL                = 0
COMPILER_MODE      = DEBUG
PRECISION          = DOUBLE
OPENMP             = 0
OUTPUT_LEVEL       = 1
DEBUG_LEVEL        = 0
REORDER_PARTICLES  = 0


# Select location of python and numpy libraries.  If blank, make will try to
# find the location of the libraries automatically using installed python
# utilities.  If you have multiple versions of python installed on your
# computer, then select the prefered version with the PYTHON variable above.
#--------------------------------------------------------------------------------------------------
PYLIB =
NUMPY =
GTEST = $(GTEST_DIR)


# Don't delete this command! Makes sure that defined variables here are passed
# to the others make
export

gandalf:
	@+$(MAKE) -C src

all:
	@+$(MAKE) -C src

executable:
	@+$(MAKE) executable -C src

unittests:
	@+$(MAKE) unittests -C src

clean:
	@+$(MAKE) clean -C src
