//=============================================================================
//  DomainBox.h
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
//=============================================================================


#ifndef _DOMAIN_BOX__H
#define _DOMAIN_BOX__H


#include <string>
#include <vector>
#include <algorithm>
#include "Precision.h"
using namespace std;


//=============================================================================
//  Structure Box
/// \brief  Simplified bounding box data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   15/11/2013
//=============================================================================
template <int ndim>
struct Box {
  FLOAT boxmin[ndim];                   ///< Minimum bounding box extent
  FLOAT boxmax[ndim];                   ///< Maximum bounding box extent
};

#ifdef MPI_PARALLEL
template <int ndim>
MPI_Datatype CreateBoxType (Box<ndim> dummy) {
  MPI_Datatype box_type;
  MPI_Datatype types[1] = {MPI_BYTE};
  MPI_Aint offsets[1] = {0};
  int blocklen[1] = {sizeof(Box<ndim>)};

  MPI_Type_create_struct(1,blocklen,offsets,types,&box_type);

  return box_type;
}
#endif


//=============================================================================
//  Structure DomainBox
/// \brief  Bounding box data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
template <int ndim>
struct DomainBox {
  string x_boundary_lhs;                ///< x-dimension LHS boundary condition
  string x_boundary_rhs;                ///< x-dimension RHS boundary condition
  string y_boundary_lhs;                ///< y-dimension LHS boundary condition
  string y_boundary_rhs;                ///< y-dimension RHS boundary condition
  string z_boundary_lhs;                ///< z-dimension LHS boundary condition
  string z_boundary_rhs;                ///< z-dimension RHS boundary condition
  FLOAT boxmin[ndim];                   ///< Minimum bounding box extent
  FLOAT boxmax[ndim];                   ///< Maximum bounding box extent
  FLOAT boxsize[ndim];                  ///< Side-lengths of bounding box
  FLOAT boxhalf[ndim];                  ///< Half side-lengths of bounding box
};


//=============================================================================
/// \brief  Helper function to find if any of the boundaries is "special" (that is, mirror or periodic)
/// \author D. A. Hubber, G. Rosotti
/// \date   28/10/2013
/// \return A boolean saying whether any "special" boundary was found
//=============================================================================
template <int ndim>
bool IsAnyBoundarySpecial(const DomainBox<ndim>& box) {
  vector<string> special;
  special.push_back("mirror");
  special.push_back("periodic");

  if (ndim >= 1) {
    if (std::find(special.begin(), special.end(), box.x_boundary_lhs) != special.end() ) return true;
    if (std::find(special.begin(), special.end(), box.x_boundary_rhs) != special.end() ) return true;
  }
  if (ndim >=2) {
    if (std::find(special.begin(), special.end(), box.y_boundary_lhs) != special.end() ) return true;
    if (std::find(special.begin(), special.end(), box.y_boundary_rhs) != special.end() ) return true;
  }
  if (ndim ==3) {
    if (std::find(special.begin(), special.end(), box.z_boundary_lhs) != special.end() ) return true;
    if (std::find(special.begin(), special.end(), box.z_boundary_rhs) != special.end() ) return true;
  }

  return false;
}


//=============================================================================
/// \brief  Helper function to say if a value is contained inside an interval
//=============================================================================
inline bool valueInRange(FLOAT value, FLOAT min, FLOAT max)
{ return (value >= min) && (value <= max); }


//=============================================================================
/// \brief  Helper function to find if two boxes overlap
/// \author D. A. Hubber, G. Rosotti
/// \date   12/11/2013
/// \return A boolean saying whether the boxes overlap
//=============================================================================
template <int ndim>
inline bool BoxesOverlap (Box<ndim>& A, Box<ndim>& B) {
  bool coord_overlap[ndim];

  for (int i=0; i<ndim; i++) {
    coord_overlap[i] = valueInRange(A.boxmin[i], B.boxmin[i], B.boxmax[i]) ||
                       valueInRange(B.boxmin[i], A.boxmin[i], A.boxmax[i]);
  }

  bool result=true;
  for (int i=0; i<ndim; i++) {
    result = result && coord_overlap[i];
  }

  return result;

}


#endif
