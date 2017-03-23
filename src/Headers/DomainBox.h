//=================================================================================================
//  DomainBox.h
//  Contans basic box and Domain box data structures.  Also contains various
//  helper routines related to the domain box, boundaries and for
//  MPI communication.
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


#ifndef _DOMAIN_BOX_H_
#define _DOMAIN_BOX_H_


#include <string>
#include <vector>
#include <algorithm>
#include "Exception.h"
#include "Precision.h"
#include "Constants.h"
using namespace std;


enum boundaryEnum{openBoundary, periodicBoundary, mirrorBoundary, wallBoundary, Nboundarytypes};


//=================================================================================================
//  Structure Box
/// \brief  Simplified bounding box data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   15/11/2013
//=================================================================================================
template <int ndim>
struct Box {
  FLOAT min[ndim];                     ///< Minimum bounding box extent
  FLOAT max[ndim];                     ///< Maximum bounding box extent

  Box () {
    for (int k=0; k<ndim; k++) {
      min[k]=-big_number;
      max[k]=big_number;
    }
  }
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



//=================================================================================================
//  Structure DomainBox
/// \brief  Bounding box data structure. Everything is 3D as the Ewald sum needs all 3 dimensions.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=================================================================================================
template <int ndim>
struct DomainBox : public Box<3> {
  boundaryEnum boundary_lhs[3];        ///< LHS boundary types
  boundaryEnum boundary_rhs[3];        ///< RHS boundary types
  FLOAT size[3];                       ///< Side-lengths of bounding box
  FLOAT half[3];                       ///< Half side-lengths of bounding box
  bool PeriodicGravity;                ///< Whether the domain is using periodic gravity.
};



//=================================================================================================
/// \brief  Helper function to find if any of the boundaries is "special" (mirror or periodic)
/// \author D. A. Hubber, G. Rosotti
/// \date   28/10/2013
/// \return A boolean saying whether any "special" boundary was found
//=================================================================================================
template <int ndim>
bool IsAnyBoundarySpecial(const DomainBox<ndim>& box)
{
  vector<boundaryEnum> special;
  special.push_back(mirrorBoundary);
  special.push_back(periodicBoundary);

  if (ndim >= 1) {
    if (std::find(special.begin(), special.end(), box.boundary_lhs[0]) != special.end() ) return true;
    if (std::find(special.begin(), special.end(), box.boundary_rhs[0]) != special.end() ) return true;
  }
  if (ndim >=2) {
    if (std::find(special.begin(), special.end(), box.boundary_lhs[1]) != special.end() ) return true;
    if (std::find(special.begin(), special.end(), box.boundary_rhs[1]) != special.end() ) return true;
  }
  if (ndim ==3) {
    if (std::find(special.begin(), special.end(), box.boundary_lhs[2]) != special.end() ) return true;
    if (std::find(special.begin(), special.end(), box.boundary_rhs[2]) != special.end() ) return true;
  }

  return false;
}

//=================================================================================================
/// \brief  Helper function to find if any of the boundaries is periodic
/// \author R. A. Booth
/// \date   28/10/2015
/// \return A boolean saying whether any periodic boundary was found
//=================================================================================================
template <int ndim>
bool IsAnyBoundaryPeriodic(const DomainBox<ndim>& box)
{
  for (int k=0; k < ndim; k++){
	if(box.boundary_lhs[k] == periodicBoundary || box.boundary_rhs[k] == periodicBoundary)
		return true ;
  }

  return false;
}


//=================================================================================================
/// \brief  Helper function to find if any of the boundaries is periodic
/// \author R. A. Booth
/// \date   28/10/2015
/// \return A boolean saying whether any periodic boundary was found
//=================================================================================================
template <int ndim>
bool IsAnyBoundaryReflecting(const DomainBox<ndim>& box)
{
  for (int k=0; k < ndim; k++){
	if(box.boundary_lhs[k] == mirrorBoundary || box.boundary_rhs[k] == mirrorBoundary)
		return true ;
  }

  return false;
}


//=================================================================================================
/// \brief  Helper function to set the boundary type enum from the string parameter.
//=================================================================================================
inline boundaryEnum setBoundaryType(string boundaryString)
{
  if (boundaryString == "open") return openBoundary;
  else if (boundaryString == "periodic") return periodicBoundary;
  else if (boundaryString == "mirror") return mirrorBoundary;
  else if (boundaryString == "wall") return wallBoundary;
  else {
    ExceptionHandler::getIstance().raise("Invalid boundary type chosen");
    return Nboundarytypes;
  }
}



//=================================================================================================
/// \brief  Helper function to say if a value is contained inside an interval
//=================================================================================================
inline bool valueInRange(const FLOAT value, const FLOAT min, const FLOAT max)
{
  return (value >= min) && (value <= max);
}



//=================================================================================================
/// \brief  Helper function to find if two boxes overlap
/// \author D. A. Hubber, G. Rosotti
/// \date   12/11/2013
/// \return A boolean saying whether the boxes overlap
//=================================================================================================
template <int ndim>
inline bool BoxesOverlap (const Box<ndim>& A, const Box<ndim>& B)
{
  bool coord_overlap[ndim];

  for (int i=0; i<ndim; i++) {
    coord_overlap[i] = valueInRange(A.min[i], B.min[i], B.max[i]) ||
                       valueInRange(B.min[i], A.min[i], A.max[i]);
  }

  bool result=true;
  for (int i=0; i<ndim; i++) {
    result = result && coord_overlap[i];
  }

  return result;

}



//=================================================================================================
/// \brief  Helper function to find if two boxes overlap
/// \author D. A. Hubber, G. Rosotti
/// \date   12/11/2013
/// \return A boolean saying whether the boxes overlap
//=================================================================================================
template <int ndim>
inline void NearestPeriodicVector
 (const DomainBox<ndim> &box,
  FLOAT dr[ndim],
  FLOAT dr_corr[ndim])
{
  for (int k=0; k<ndim; k++) dr_corr[k] = 0.0;
  for (int k=0; k<ndim; k++) {
    if (box.boundary_lhs[k] == periodicBoundary && box.boundary_rhs[k] == periodicBoundary) {
      if (dr[k] > box.half[k]) dr_corr[k] = -box.size[k];
      else if (dr[k] < -box.half[k]) dr_corr[k] = box.size[k];
    }
  }
  for (int k=0; k<ndim; k++) dr[k] += dr_corr[k];

  return;
}



#endif
