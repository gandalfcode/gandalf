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
#include "Precision.h"
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
  FLOAT boxmin[3];                     ///< Minimum bounding box extent
  FLOAT boxmax[3];                     ///< Maximum bounding box extent
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
/// \brief  Bounding box data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=================================================================================================
template <int ndim>
struct DomainBox {
  boundaryEnum boundary_lhs[3];        ///< LHS boundary types
  boundaryEnum boundary_rhs[3];        ///< RHS boundary types
  FLOAT boxmin[3];                     ///< Minimum bounding box extent
  FLOAT boxmax[3];                     ///< Maximum bounding box extent
  FLOAT boxsize[3];                    ///< Side-lengths of bounding box
  FLOAT boxhalf[3];                    ///< Half side-lengths of bounding box
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
/// \brief  Helper function to set the boundary type enum from the string parameter.
//=================================================================================================
inline boundaryEnum setBoundaryType(string boundaryString)
{
  if (boundaryString == "open") return openBoundary;
  else if (boundaryString == "periodic") return periodicBoundary;
  else if (boundaryString == "mirror") return mirrorBoundary;
  else if (boundaryString == "wall") return wallBoundary;
  else {
    exit(0);
  }
}



//=================================================================================================
/// \brief  Helper function to say if a value is contained inside an interval
//=================================================================================================
static inline bool valueInRange(const FLOAT value, const FLOAT min, const FLOAT max)
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
static inline bool BoxesOverlap (const Box<ndim>& A, const Box<ndim>& B)
{
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



//=================================================================================================
/// \brief  Helper function to find if two boxes overlap
/// \author D. A. Hubber, G. Rosotti
/// \date   12/11/2013
/// \return A boolean saying whether the boxes overlap
//=================================================================================================
template <int ndim>
static inline void NearestPeriodicVector
 (const DomainBox<ndim> &box,
  FLOAT dr[ndim],
  FLOAT dr_corr[ndim])
{
  for (int k=0; k<ndim; k++) dr_corr[k] = 0.0;
  for (int k=0; k<ndim; k++) {
    if (box.boundary_lhs[k] == periodicBoundary && box.boundary_rhs[k] == periodicBoundary) {
      if (dr[k] > box.boxhalf[k]) dr_corr[k] = -box.boxsize[k];
      else if (dr[k] < -box.boxhalf[k]) dr_corr[k] = box.boxsize[k];
    }
  }
  for (int k=0; k<ndim; k++) dr[k] += dr_corr[k];

  return;
}


//=================================================================================================
/// \brief  Helper function to find if two boxes overlap
/// \author D. A. Hubber, G. Rosotti
/// \date   12/11/2013
/// \return A boolean saying whether the boxes overlap
//=================================================================================================
template <int ndim>
inline bool FractionalBoxOverlap
 (Box<ndim> &box1,                     ///< ..
  Box<ndim> &box2,                     ///< ..
  DomainBox<ndim> &simbox,             ///< ..
  FLOAT &overlapfrac)                  ///< ..
{
  int k;
  FLOAT dr[ndim];

  // First, calculate relative position vector between boxes
  for (k=0; k<ndim; k++) dr[k] = box2.r[k] - box1.r[k];

  // Calculate closest position vector to nearest


  return false;
}


//=================================================================================================
/// \brief  Make the ghost particles based upon the boundary conditions
/// \author R. A. Booth
/// \date   26/10/2015
/// \return The number of neighbours found
//=================================================================================================
template<int ndim>
class MirrorNeighbourFinder
{
public:
	MirrorNeighbourFinder(const DomainBox<ndim>& simbox)
	: _domain(simbox),
	  _any_periodic(false), _any_mirror(false), _need_mirrors(false), _any_special(false)
	{
		for (int k=0; k <ndim; k++)
		{
			_periodic_bound[k] = _domain.boundary_lhs[k] == periodicBoundary ;
			if (_periodic_bound[k])
			 assert(_domain.boundary_rhs[k] == periodicBoundary) ;

			_mirror_bound[k][0] = _domain.boundary_lhs[k] == mirrorBoundary ;
			_mirror_bound[k][1] = _domain.boundary_rhs[k] == mirrorBoundary ;

			_any_periodic |= _periodic_bound[k] ;
			_any_mirror   |= (_mirror_bound[k][0] | _mirror_bound[k][1]) ;
		}
	}

	//=================================================================================================
	/// \brief  Set a target cell to compute mirrors for.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return A boolean saying whether the boxes overlap
	//=================================================================================================
	template<template <int> class TreeCell>
	void SetTargetCell(const TreeCell<ndim>& cell)
	{
	  _need_mirrors = false ;
	  for (int k=0; k < ndim; k++){
		_centre[k] = cell.rcell[k] ;

		if (_any_mirror){
		  _cell.boxmin[k] = cell.hboxmin[k] ;
	      _cell.boxmax[k] = cell.hboxmax[k] ;

		  // Compute whether we need a mirror in the gather case
		  _need_mirror[k][0] = (_mirror_bound[k][0] & _cell.boxmin[k] < _domain.boxmin[k]);
		  _need_mirror[k][1] = (_mirror_bound[k][1] & _cell.boxmax[k] > _domain.boxmax[k]);
		  _need_mirrors     |= (_need_mirror[k][0] | _need_mirror[k][1]) ;
	    }
	  }
	  _any_special = _need_mirrors | _any_periodic ;
	}

	//=================================================================================================
	/// \brief  Find the smallest separation vector dr[] in a periodic box
	/// \author D. A. Hubber, G. Rosotti
	/// \date   12/11/2013
	/// \return A boolean saying whether the boxes overlap
	//=================================================================================================
	void NearestPeriodicVector(FLOAT dr[ndim]) const
	{
   	  if (_any_periodic)
		{
   		  for (int k=0; k<ndim; k++) {
   		    if (_periodic_bound[k]) {
   		      if (dr[k] > _domain.boxhalf[k])
   		        dr[k] -=_domain.boxsize[k];
   		      else if (dr[k] < -_domain.boxhalf[k])
   		        dr[k] += _domain.boxsize[k];
   		    }
   		  }
		}
	}

	//=================================================================================================
	/// \brief Construct the centres and reflection signs of a cell. This list will include the
	///        original cell or the periodic neighbour of the original cell.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class TreeCell>
	int ConstructGhostVectorsGather(const TreeCell<ndim>& cell, FLOAT * r, int * sign) const
	{
	  // First find the nearest periodic mirror
	  if (_any_periodic){
	    FLOAT dr[ndim] ;

	    for (int k=0; k <ndim; k++){
		  dr[k] = cell.rcell[k] - _centre[k] ;
	    }

	    NearestPeriodicVector(dr) ;
	    for (int k=0; k <ndim; k++){
		  r[k] = _centre[k] + dr[k];
		  sign[k] = 1;
	    }
	  }
	  else {
		for (int k=0; k <ndim; k++){
		  r[k] = cell.rcell[k] ;
		  sign[k] = 1;
		}
	  }

	  // Number of Ghost cells
	  int Nghost = 1 ;

	  // Now recursively reflect the cells
	  if (_need_mirrors){
	    int nc = Nghost ;
	    // Loop over the possible directions for reflections
	    for (int k = 0; k < ndim; k++){
	      // Save the current number of images
		  Nghost = nc ;

		  // Do reflections on the left edge
		  if (_need_mirror[k][0]){
			for (int n=0; n < Nghost; n++){
			   for (int l=0; l < ndim; ++l){
			     r[nc*ndim + l] = r[n*ndim +l] ;
			     sign[nc*ndim + l] = sign[n*ndim +l] ;
			   }
			   r[nc*ndim + k] = 2*_domain.boxmin[k] - r[n*ndim +k] ;
			   sign[nc*ndim +k] *= -1 ;
			   nc++ ;
			}
		  }
		 // Do reflections on the right edge
	     if (_need_mirror[k][1]){
		   for (int n=0; n < Nghost; n++){
			 for (int l=0; l < ndim; ++l){
			   r[nc*ndim + l] = r[n*ndim +l] ;
			   sign[nc*ndim + l] = sign[n*ndim +l] ;
			 }
			 r[nc*ndim + k] = 2*_domain.boxmax[k] - r[n*ndim +k] ;
			 sign[nc*ndim +k] *= -1 ;
	         nc++ ;
		    }
		  }
	    }
	    // Update to the new total number of images
	    Nghost = nc ;
	  }
	  return Nghost ;
	}

	//=================================================================================================
	/// \brief Compute all interacting neighbours in the target cell using scatter-gather routines
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class TreeCell>
	int ConstructGhostVectorsScatterGather(const TreeCell<ndim>& cell, FLOAT * r, int * sign) const
	{
	  // First find the nearest periodic mirror
	  if (_any_periodic){
	    FLOAT dr[ndim] ;

	    for (int k=0; k <ndim; k++){
		  dr[k] = cell.rcell[k] - _centre[k] ;
	    }

	    NearestPeriodicVector(dr) ;

	    for (int k=0; k <ndim; k++){
		  r[k] = _centre[k] + dr[k];
		  sign[k] = 1;
	    }
	  }
	  else {
		for (int k=0; k <ndim; k++){
		  r[k] = cell.rcell[k] ;
		  sign[k] = 1;
		}
	  }

	  // Number of Ghost cells
	  int Nghost = 1 ;

	  // Now recursively reflect the cells
	  if (_any_mirror){
	    int nc = Nghost ;
	    // Loop over the possible directions for reflections
	    for (int k = 0; k < ndim; k++){
	      // Save the current number of images
		  Nghost = nc ;

		  // Do reflections on the left edge
		  if (_mirror_bound[k][0]){
		    if ((2*_domain.boxmin[k] - cell.hboxmin[k]) > _cell.boxmin[k]){
			  for (int n=0; n < Nghost; n++){
			     for (int l=0; l < ndim; ++l){
			       r[nc*ndim + l] = r[n*ndim +l] ;
			       sign[nc*ndim + l] = sign[n*ndim +l] ;
			     }
			     r[nc*ndim + k] = 2*_domain.boxmin[k] - r[n*ndim +k] ;
			     sign[nc*ndim +k] *= -1 ;
			     nc++ ;
			  }
		    }
		  }
		 // Do reflections on the right edge
	     if (_mirror_bound[k][1]){
		   if ((2*_domain.boxmax[k] - cell.hboxmax[k]) < _cell.boxmax[k]){
		     for (int n=0; n < Nghost; n++){
			   for (int l=0; l < ndim; ++l){
			     r[nc*ndim + l] = r[n*ndim +l] ;
			     sign[nc*ndim + l] = sign[n*ndim +l] ;
			   }
			   r[nc*ndim + k] = 2*_domain.boxmax[k] - r[n*ndim +k] ;
			   sign[nc*ndim +k] *= -1 ;
	           nc++ ;
		      }
		   }
		  }
	    }
	    // Update to the new total number of images
	    Nghost = nc ;
	  }
	  return Nghost ;
	}

	template< template<int> class TreeCell, template<int> class ParticleType>
	void CorrectGhostParticlePosition(const TreeCell<ndim>& cell, const FLOAT * r, const int * sign,
			                          ParticleType<ndim>& p) const
	{
		if (_any_special){
		  for (int k=0; k < ndim; k++){
		    FLOAT dr_k = p.r[k] - cell.rcell[k] ;
		    p.r[k] = r[k] + dr_k * sign[k] ;
		    p.v[k] *= sign[k] ;
		  }
		}
	}

private:
	DomainBox<ndim> _domain ;
	Box<ndim> _cell ;
	FLOAT _centre[ndim] ;
	bool _mirror_bound[ndim][2] ;
	bool _need_mirror[ndim][2] ;
	bool _periodic_bound[ndim] ;
	bool _any_periodic, _any_mirror, _need_mirrors, _any_special;
};



#endif
