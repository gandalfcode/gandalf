//=================================================================================================
//  MpiNode.h
//  Contains MPI node data class definition.
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


#ifndef _MPI_NODE_H_
#define _MPI_NODE_H_


#include <string>
#include "Precision.h"
#include "Hydrodynamics.h"
#include "Particle.h"
#include "SmoothingKernel.h"
using namespace std;



//=================================================================================================
//  Class MpiNode
/// \brief   MPI node data class
/// \details MPI node data class
/// \author  D. A. Hubber, G. Rosotti
/// \date    09/10/2013
//=================================================================================================
template <int ndim>
class MpiNode
{
 public:

  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  MpiNode();
  ~MpiNode();


  // Other functions
  //-----------------------------------------------------------------------------------------------
  void PackNodeData(void);
  void UnpackNodeData(void);
  void UpdateBoundingBoxData(int, Hydrodynamics<ndim> *, SmoothingKernel<ndim> *);


  // MPI node variables
  //-----------------------------------------------------------------------------------------------
  int ifirst;                          ///< i.d. of first ghost from node
  int ilast;                           ///< i.d. of last ghost from node
  int Nhydro;                          ///< No. of hydro particles on node
  int Ntot;                            ///< Total no. of particles on node
  int Ntotmax;                         ///< Max. no. of particles on node
  int Nghost;                          ///< No. of ghost particles originally
                                       ///< from node exported to current node
  FLOAT hmax;                          ///< Maximum smoothing length on node
  FLOAT worktot;                       ///< Total 'work' on each node
  Box<ndim> domain;                    ///< Box representing node volume
  Box<ndim> rbox;                      ///< Bounding box of particle positions in node
  Box<ndim> hbox;                      ///< Bounding box of h-extent in node

  int *ids;                            ///< List of particles in tree

};
#endif
