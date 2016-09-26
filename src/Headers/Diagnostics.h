//=================================================================================================
//  Diagnostics.h
//  ..
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
//=================================================================================================


#ifndef _DIAGNOSTICS__H
#define _DIAGNOSTICS__H

#include "Precision.h"
#ifdef MPI_PARALLEL
#include <stddef.h>
#include "mpi.h"
#include "Exception.h"
#endif


//=================================================================================================
//  Structure Diagnostics
/// \brief  Structure containing snapshot of current diagnostic quantities.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=================================================================================================
template <int ndim>
struct Diagnostics
{
  int Nhydro;                          ///< Total no. of SPH particles
  int Nstar;                           ///< Total no. of star particles
  int Ndead;                           ///< Total no. of dead SPH particles
  DOUBLE Eerror;                       ///< Total energy error
  DOUBLE Etot;                         ///< Total energy
  DOUBLE utot;                         ///< Total thermal energy
  DOUBLE ketot;                        ///< Total kinetic energy
  DOUBLE gpetot;                       ///< Total grav. potential energy
  DOUBLE mtot;                         ///< Total mass in simulation
  DOUBLE mom[ndim];                    ///< Total momentum vector
  DOUBLE angmom[3];                    ///< Total angular momentum vector
  DOUBLE force[ndim];                  ///< Net force
  DOUBLE rcom[ndim];                   ///< Position of centre of mass
  DOUBLE vcom[ndim];                   ///< Velocity of centre of mass

#ifdef MPI_PARALLEL
  static MPI_Datatype CreateMpiDataType() {
    MPI_Datatype diagnostic_type;
    MPI_Datatype types[1] = {MPI_BYTE};
    MPI_Aint offsets[1] = {0};
    int blocklen[1] = {sizeof(Diagnostics<ndim>)};
    MPI_Type_create_struct(1,blocklen,offsets,types,&diagnostic_type);
    return diagnostic_type;
  }
#endif

};

#endif
