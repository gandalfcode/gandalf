//=================================================================================================
//  MirrorNeighbours.hpp
//  Contains the definitions of a class used for constructing periodic and mirror ghosts
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

#ifndef _MIRROR_NEIGHBOURS_H_
#define _MIRROR_NEIGHBOURS_H_

#include "DomainBox.h"
#include "Precision.h"

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
	/// \brief  Find 'periodic' correction vector.
	/// \author D. A. Hubber, G. Rosotti
	/// \date   12/11/2013
	/// \return A boolean saying whether the boxes overlap
	//=================================================================================================
	void PeriodicDistanceCorrection(const FLOAT dr[ndim], FLOAT dr_corr[ndim]) const
	{
	  if (_any_periodic)
		{
	   	  for (int k=0; k<ndim; k++) {
	   	    if (_periodic_bound[k]) {
	   	      if (dr[k] > _domain.boxhalf[k])
	   	        dr_corr[k] =- _domain.boxsize[k];
	   	      else if (dr[k] < -_domain.boxhalf[k])
	   	        dr_corr[k] = _domain.boxsize[k];
	   	      else
	   	    	dr_corr[k] = 0 ;
	   	    }
	   	  }
		}
	}

	//=================================================================================================
	/// \brief  Find the maximum number of mirrors required
	/// \author R. A. Booth
	/// \date   28/10/2015
	/// \return An integer specifying the max number of neighbours (>=1)
	//=================================================================================================
	int MaxNumMirrors() const
	{
	  int NumMirrors = 1;
		if (_any_mirror){
		  for (int k=0; k < ndim; k++){
			NumMirrors *= 1 + _mirror_bound[k][0] + _mirror_bound[k][1] ;
		}
	  }
	  return NumMirrors ;
	}
	/*
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
	*/
	//=================================================================================================
	/// \brief Construct the centres and reflection signs of a cell. This list will include the
	///        original cell or the periodic neighbour of the original cell.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class Particle>
	int ConstructGhostsScatterGather(const Particle<ndim>& p, Particle<ndim>* ngbs) const
	{
	  // First find the nearest periodic mirror
	  ngbs[0] = p ;
	  if (_any_periodic){
	    FLOAT dr[ndim] ;

	    for (int k=0; k <ndim; k++){
		  dr[k] = p.r[k] - _centre[k] ;
	    }

	    NearestPeriodicVector(dr) ;
	    for (int k=0; k <ndim; k++){
		  ngbs[0].r[k] = _centre[k] + dr[k];
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
			FLOAT dx = 2*_domain.boxmin[k] - p.r[k] - _cell.boxmin[k] ;
			if (dx*dx < p.hrangesqd){
			  for (int n=0; n < Nghost; n++){
			    ngbs[nc] = ngbs[n] ;
			    ngbs[nc].r[k] = 2*_domain.boxmin[k] - ngbs[n].r[k] ;
			    ngbs[nc].v[k] *= -1 ;
			    nc++ ;
			 }
		   }
		 }
		 // Do reflections on the right edge
	     if (_mirror_bound[k][1]){
		   FLOAT dx = 2*_domain.boxmax[k] - p.r[k] - _cell.boxmax[k] ;
		   if (dx*dx < p.hrangesqd){
		     for (int n=0; n < Nghost; n++){
		       ngbs[nc] = ngbs[n] ;
			   ngbs[nc].r[k] = 2*_domain.boxmax[k] - ngbs[n].r[k] ;
			   ngbs[nc].v[k] *= -1 ;
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

    //=================================================================================================
	/// \brief Construct the centres and reflection signs of a cell. This list will include the
	///        original cell or the periodic neighbour of the original cell.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class Particle>
	int ConstructGhostsGather(const Particle<ndim>& p, Particle<ndim>* ngbs) const
	{
	  // First find the nearest periodic mirror
	  ngbs[0] = p ;
	  if (_any_periodic){
	    FLOAT dr[ndim] ;

	    for (int k=0; k <ndim; k++){
		  dr[k] = p.r[k] - _centre[k] ;
	    }

	    NearestPeriodicVector(dr) ;
	    for (int k=0; k <ndim; k++){
		  ngbs[0].r[k] = _centre[k] + dr[k];
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
			   ngbs[nc] = ngbs[n] ;
			   ngbs[nc].r[k] = 2*_domain.boxmin[k] - ngbs[n].r[k] ;
			   ngbs[nc].v[k] *= -1 ;
			   nc++ ;
			}
		  }
		 // Do reflections on the right edge
	     if (_need_mirror[k][1]){
		   for (int n=0; n < Nghost; n++){
		     ngbs[nc] = ngbs[n] ;
			 ngbs[nc].r[k] = 2*_domain.boxmax[k] - ngbs[n].r[k] ;
			 ngbs[nc].v[k] *= -1 ;
			 nc++ ;
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



#endif//_MIRROR_NEIGHBOURS_H_

