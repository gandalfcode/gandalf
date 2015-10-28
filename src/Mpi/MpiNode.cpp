//=============================================================================
//  MpiNode.cpp
//  Contains functions relating to the MPI node class.
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Constants.h"
#include "Precision.h"
#include "SmoothingKernel.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
#include "MpiNode.h"
using namespace std;



//=============================================================================
//  MpiNode::MpiNode()
/// MPI node class constructor.
//=============================================================================
template <int ndim>
MpiNode<ndim>::MpiNode()
{
  ifirst = -1;
  ilast = -1;
  Nhydro = 0;
  Ntot = 0;
  Ntotmax = 0;
  Nghost = 0;
  hmax = 0.0;
  worktot = 0.0;
}



//=============================================================================
//  MpiNode::~MpiNode()
/// MPI node class destructor.
//=============================================================================
template <int ndim>
MpiNode<ndim>::~MpiNode()
{
}



//=============================================================================
//  MpiNode::PackNodeData
/// Pack all node data into contiguous arrays for broadcasting by MPI.
//=============================================================================
template <int ndim>
void MpiNode<ndim>::PackNodeData(void)
{
  return;
}



//=============================================================================
//  MpiNode::UnpackNodeData
/// Unpack all node data broadcast from other domains.
//=============================================================================
template <int ndim>
void MpiNode<ndim>::UnpackNodeData(void)
{
  return;
}



//=============================================================================
//  MpiNode::UpdateBoundingBoxData
/// Update all local domain bounding box data.
//=============================================================================
template <int ndim>
void MpiNode<ndim>::UpdateBoundingBoxData
(int Npart,                            ///< No. of SPH particles
 Hydrodynamics<ndim> *hydro,           ///< Pointer to SPH object
 SmoothingKernel<ndim> *kernptr)       ///< Pointer to kernel object
{
  int i;                               // Particle counter
  int k;                               // Dimension counter
  FLOAT hrange;                        // ..
  const FLOAT kernange = kernptr->kernrange;

  // Initialise bounding box values
  for (k=0; k<ndim; k++) rbox.boxmin[k] = big_number;
  for (k=0; k<ndim; k++) rbox.boxmax[k] = -big_number;
  for (k=0; k<ndim; k++) hbox.boxmin[k] = big_number;
  for (k=0; k<ndim; k++) hbox.boxmax[k] = -big_number;


  // Loop over all particles and compute new bounding boxes
  //---------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);

    hrange = 2.0*kernange*part.h;
    for (k=0; k<ndim; k++) {
      rbox.boxmin[k] = min(rbox.boxmin[k], part.r[k]);
      rbox.boxmax[k] = max(rbox.boxmax[k], part.r[k]);
      hbox.boxmin[k] = min(hbox.boxmin[k], part.r[k] - hrange);
      hbox.boxmax[k] = max(hbox.boxmax[k], part.r[k] + hrange);
    }
  }

#ifdef OUTPUT_ALL
  cout << "RBOX : " << rbox.boxmin[0] << "    " << rbox.boxmax[0] << endl;
  cout << "HBOX : " << hbox.boxmin[0] << "    " << hbox.boxmax[0] << endl;
  cout << "DOMAIN : " << domain.boxmin[0] << "    " << domain.boxmax[0] << endl;
#endif

  return;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class MpiNode<1>;
template class MpiNode<2>;
template class MpiNode<3>;
