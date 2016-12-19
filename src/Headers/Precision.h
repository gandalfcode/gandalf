//=============================================================================
//  Precision.h
//  Contains macro definitions for floating point precision in all routines.
//  Also contains locally-defined MPI data types.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=============================================================================


#ifndef _PRECISION_H_
#define _PRECISION_H_


// MPI communication floating point data types
//-----------------------------------------------------------------------------
#if defined(MPI_PARALLEL)
#include "mpi.h"
#if defined(GANDALF_SINGLE_PRECISION)
#define GANDALF_MPI_FLOAT MPI_FLOAT
#define GANDALF_MPI_DOUBLE MPI_DOUBLE
#elif defined(GANDALF_DOUBLE_PRECISION)
#define GANDALF_MPI_FLOAT MPI_DOUBLE
#define GANDALF_MPI_DOUBLE MPI_DOUBLE
#endif
#endif

// Floating point data types
//-----------------------------------------------------------------------------
#if defined(GANDALF_SINGLE_PRECISION)
typedef float FLOAT;
typedef double DOUBLE;
#elif defined(GANDALF_DOUBLE_PRECISION)
typedef double FLOAT;
typedef double DOUBLE;
#endif

#if defined(GANDALF_SNAPSHOT_SINGLE_PRECISION)
typedef float SNAPFLOAT;
#else
typedef double SNAPFLOAT;
#endif

#endif
