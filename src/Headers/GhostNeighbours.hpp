//=================================================================================================
//  GhostNeighbours.hpp
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

#ifndef _GHOST_NEIGHBOURS_H_
#define _GHOST_NEIGHBOURS_H_

#include "DomainBox.h"
#include "Precision.h"
#include "Particle.h"
#include "TreeCell.h"

//=================================================================================================
/// \brief  Make the ghost particles based upon the boundary conditions
/// \author R. A. Booth
/// \date   26/10/2015
/// \return The number of neighbours found
//=================================================================================================
template<int ndim>
class GhostNeighbourFinder
{
private:

  const DomainBox<ndim>& _domain ;
  Box<ndim> _cell ;
  Box<ndim> _hbox ;
  Box<ndim> _vbox ;
  FLOAT _centre[ndim] ;
  bool _mirror_bound[ndim][2] ;
  bool _need_mirror[ndim][2] ;
  bool _periodic_bound[ndim] ;
  bool _any_periodic, _any_mirror, _need_mirrors, _any_special;

  void 	SetBoundaryFlags() {
    for (int k=0; k <ndim; k++) {
      _periodic_bound[k] = _domain.boundary_lhs[k] == periodicBoundary ;
      if (_periodic_bound[k])
    	  assert(_domain.boundary_rhs[k] == periodicBoundary) ;

      _mirror_bound[k][0] = _domain.boundary_lhs[k] == mirrorBoundary ;
      _mirror_bound[k][1] = _domain.boundary_rhs[k] == mirrorBoundary ;

      _any_periodic |= _periodic_bound[k] ;
      _any_mirror   |= (_mirror_bound[k][0] | _mirror_bound[k][1]) ;

      _need_mirror[k][0] = _need_mirror[k][1] = false ;
    }
  }

public:

	GhostNeighbourFinder(const DomainBox<ndim>& simbox)
	: _domain(simbox),
	  _any_periodic(false), _any_mirror(false), _need_mirrors(false), _any_special(false)
	{
	  SetBoundaryFlags() ;
	}

	GhostNeighbourFinder(const DomainBox<ndim>& simbox, const TreeCellBase<ndim>& cell)
	: _domain(simbox),
	  _any_periodic(false), _any_mirror(false), _need_mirrors(false), _any_special(false)
	{
	  SetBoundaryFlags() ;
	  SetTargetCell(cell) ;
	}


	//=================================================================================================
	/// \brief  Set a target cell to compute mirrors for.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return A boolean saying whether the boxes overlap
	//=================================================================================================
	void SetTargetCell(const TreeCellBase<ndim>& cell)
	{
	  _need_mirrors = false ;
	  for (int k=0; k < ndim; k++){
//		_centre[k] = cell.rcell[k] ;
    cell.ComputeCellCentre(_centre);
		_cell.min[k] = cell.bb.min[k] ;
		_cell.max[k] = cell.bb.max[k] ;
        _hbox.min[k] = cell.hbox.min[k] ;
        _hbox.max[k] = cell.hbox.max[k] ;
        _vbox.min[k] = cell.vbox.min[k] ;
        _vbox.max[k] = cell.vbox.max[k] ;
		if (_any_mirror){
		  // Compute whether we need a mirror in the gather case
		  _need_mirror[k][0] = (_mirror_bound[k][0] & (_hbox.min[k] < _domain.min[k]));
		  _need_mirror[k][1] = (_mirror_bound[k][1] & (_hbox.max[k] > _domain.max[k]));
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
	type_flag NearestPeriodicVector(FLOAT dr[ndim]) const
	{
	  type_flag bound ;
   	  if (_any_periodic)
		{
   		  for (int k=0; k<ndim; k++) {
   		    if (_periodic_bound[k]) {
   		      if (dr[k] > _domain.half[k]) {
   		        dr[k] -=_domain.size[k];
   		        bound.set(periodic_bound_flags[k][1]) ;
   		      }
   		      else if (dr[k] < -_domain.half[k]) {
   		        dr[k] += _domain.size[k];
   		        bound.set(periodic_bound_flags[k][0]) ;
   		      }
   		    }
   		  }
		}
   	  return bound ;
	}

	//=================================================================================================
	/// \brief  Apply 'periodic' correction.
	/// \author D. A. Hubber, G. Rosotti
	/// \date   12/11/2013
	/// \return A boolean saying whether the boxes overlap
	//=================================================================================================
	type_flag ApplyPeriodicDistanceCorrection(FLOAT r[ndim], FLOAT dr[ndim]) const
	{
	  type_flag bound ;
	  if (_any_periodic)
		{
	   	  for (int k=0; k<ndim; k++)
	   	    if (_periodic_bound[k]) {
	   	      if (dr[k] > _domain.half[k]) {
	   	        dr[k] += - _domain.size[k];
	   	        r[k]  += - _domain.size[k];
	   	        bound.set(periodic_bound_flags[k][1]) ;
	   	      }
	   	      else if (dr[k] < -_domain.half[k]) {
	   	        dr[k] += _domain.size[k];
	   	        r[k]  += _domain.size[k] ;
	   	        bound.set(periodic_bound_flags[k][0]) ;
	   	      }
	   	    }
		}
	  return bound ;
	}

    //=================================================================================================
    /// \brief  Find 'periodic' correction vector.
    /// \author D. A. Hubber, G. Rosotti
    /// \date   12/11/2013
    /// \return A boolean saying whether the boxes overlap
    //=================================================================================================
    type_flag PeriodicDistanceCorrection(const FLOAT dr[ndim], FLOAT dr_corr[ndim]) const
    {
      type_flag bound ;
      if (_any_periodic)
        {
          for (int k=0; k<ndim; k++) {
            if (_periodic_bound[k]) {
              if (dr[k] > _domain.half[k]) {
                dr_corr[k] =- _domain.size[k];
                bound.set(periodic_bound_flags[k][1]) ;
              }
              else if (dr[k] < -_domain.half[k]) {
                dr_corr[k] = _domain.size[k];
                bound.set(periodic_bound_flags[k][0]) ;
              }
              else
                dr_corr[k] = 0 ;
            }
          }
        }
      return bound ;
    }
	//=================================================================================================
	/// \brief  Find out whether two boxes overlap in a domain that might be periodic
	/// \author R. A. Booth
	/// \date   09/10/2016
	/// \return A boolean saying whether the boxes overlap
	//=================================================================================================
	bool PeriodicBoxOverlap(const Box<ndim>& box1, const Box<ndim>& box2) const
	{
	  if (!_any_periodic)
	    return BoxOverlap(ndim, box1.min, box1.max, box2.min, box2.max) ;
	  else {
	    // Find the smallest possible distance between box centres
        FLOAT dr[ndim];
        FLOAT dr_corr[ndim];
	    for (int k=0; k<ndim; k++)
	      dr[k] = 0.5*((box2.max[k] + box2.min[k]) - (box1.max[k] + box1.min[k])) ;

	    // Find whether the boxes can overlap.
	    if (PeriodicDistanceCorrection(dr, dr_corr).is_periodic()) {
	      Box<ndim> periodic_box ;
	      for (int k=0; k<ndim; k++) {
	        periodic_box.min[k] = box2.min[k] + dr_corr[k];
	        periodic_box.max[k] = box2.max[k] + dr_corr[k];
	      }
	      return BoxOverlap(ndim, box1.min, box1.max, periodic_box.min, periodic_box.max) ;
	    }
	    else {
	      return BoxOverlap(ndim, box1.min, box1.max, box2.min, box2.max) ;
	    }
	  }
	}

	//=================================================================================================
	/// \brief  Find the maximum number of mirrors required
	/// \author R. A. Booth
	/// \date   28/10/2015
	/// \return An integer specifying the max number of neighbours (>=1)
	//=================================================================================================
	int MaxNumGhosts() const
	{
	  int NumGhosts = 1;
		if (_any_mirror){
		  for (int k=0; k < ndim; k++){
			NumGhosts *= 1 + _mirror_bound[k][0] + _mirror_bound[k][1] ;
		}
	  }
	  return NumGhosts ;
	}

	//=================================================================================================
	/// \brief Construct the centres and reflection signs of a cell. This list will include the
	///        original cell or the periodic neighbour of the original cell.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class InParticleType, class OutParticleType>
	int ConstructGhostsScatterGather(const InParticleType<ndim>& p, vector<OutParticleType>& ngbs) const
	{
	  // First find the nearest periodic mirror
	  ngbs.push_back(p);
	  if (_any_periodic)
		_MakePeriodicGhost(ngbs.back()) ;


	  // Number of Ghost cells
	  int Nghost = 1 ;

	  if (_any_mirror)
		Nghost = _MakeReflectedScatterGatherGhosts(ngbs) ;

	  return Nghost ;
	}

    //=================================================================================================
	/// \brief Construct the centres and reflection signs of a cell. This list will include the
	///        original cell or the periodic neighbour of the original cell.
	/// \author R. A. Booth
	/// \date   13/10/2016
	/// \return The number of neighbours found
	//===============================================================================================
	template<template <int> class ParticleType>
	int ConstructAllGhosts(const ParticleType<ndim>& p, ParticleType<ndim>* ngbs) const
	{
	  // First find the nearest periodic mirror
	  ngbs[0] = p ;
	  if (_any_periodic)
		_MakePeriodicGhost(ngbs[0]) ;


	  // Number of Ghost cells
	  int Nghost = 1 ;
	  // Now recursively reflect the cells
	  if (_need_mirrors)
		Nghost = _MakeReflectedGhostsAll(ngbs) ;

	  return Nghost ;
	}

    //=================================================================================================
    /// \brief Construct the centres and reflection signs of a cell. This list will include the
    ///        original cell or the periodic neighbour of the original cell.
    /// \author R. A. Booth
    /// \date   27/10/2015
    /// \return The number of neighbours found
    //===============================================================================================
    template<template <int> class ParticleType, class OutParticleType>
    int ConstructGhostsGather(const ParticleType<ndim>& p, vector<OutParticleType>& ngbs) const
    {
      // First find the nearest periodic mirror
      ngbs.push_back(p);
      if (_any_periodic)
        _MakePeriodicGhost(ngbs.back()) ;


      // Number of Ghost cells
      int Nghost = 1 ;
      // Now recursively reflect the cells
      if (_need_mirrors)
        Nghost = _MakeReflectedGhostsGather(ngbs) ;

      return Nghost ;
    }


	template< template<int> class TreeCell, template<int> class ParticleType>
	void CorrectGhostParticlePosition(const TreeCell<ndim>& cell, const FLOAT * r, const int * sign,
			                          ParticleType<ndim>& p) const
	{
		if (_any_special){
		  for (int k=0; k < ndim; k++){
		    FLOAT dr_k = p.r[k] - cell.rcell(k) ;
		    p.r[k] = r[k] + dr_k * sign[k] ;
		    p.v[k] *= sign[k] ;
		  }
		}
	}

private:
    //=================================================================================================
	//  _MakePeriodicGhost
	/// \brief Do the actual construction of the nearest periodic ghost.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<class ParticleType>
	void _MakePeriodicGhost(ParticleType& p) const {
	  FLOAT dr[ndim] ;

	  for (int k=0; k <ndim; k++)
	    dr[k] = p.r[k] - _centre[k] ;

	  type_flag bound_flag = NearestPeriodicVector(dr) ;

	  if (bound_flag.is_periodic())
	    for (int k=0; k <ndim; k++)
	      p.r[k] = _centre[k] + dr[k];

	  p.flags.set(bound_flag.get()) ;

	}

    //=================================================================================================
	//  _MakeReflectedGhostsGather
	/// \brief Do the actual construction of the mirror ghosts within the gather range. Assumes the
	/// first particle is already saved in ngbs.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<class ParticleType>
	int _MakeReflectedGhostsGather(vector<ParticleType>& ngbs) const {
	  int nc = 1 ;
	  const int old_size = ngbs.size()-1;

	  // Loop over the possible directions for reflections
	  for (int k = 0; k < ndim; k++){
		// Save the current number of images
		int Nghost = nc ;
		// Do reflections on the left edge
		if (_need_mirror[k][0]){
		  for (int n=0; n < Nghost; n++){
			double rk = 2*_domain.min[k] - ngbs[n].r[k] ;
			if (rk > _hbox.min[k]) {
              ngbs.push_back(ngbs[n+old_size]);
			  reflect<ParticleType::NDIM>(ngbs.back(), k, _domain.min[k]) ;
			  ngbs[nc].flags.set(mirror_bound_flags[k][0]) ;
			  nc++ ;
			}
		  }
		}
		// Do reflections on the right edge
		if (_need_mirror[k][1]){
		  for (int n=0; n < Nghost; n++){
			double rk = 2*_domain.max[k] - ngbs[n].r[k] ;
			if (rk < _hbox.max[k]) {
              ngbs.push_back(ngbs[n+old_size]);
			  reflect<ParticleType::NDIM>(ngbs.back(), k, _domain.max[k]) ;
			  ngbs[nc].flags.set(mirror_bound_flags[k][1]) ;
			  nc++ ;
			}
		  }
		}
	  }
	  return nc ;
	}

	//=================================================================================================
	//  _MakeReflectedGhostsScatterGather
	/// \brief Do the actual construction of the mirror ghosts within the scatter-gather range. Assumes
	/// the  first particle is already saved in ngbs.
	/// \author R. A. Booth
	/// \date   27/10/2015
	/// \return The number of neighbours found
	//===============================================================================================
	template<class ParticleType>
	int _MakeReflectedScatterGatherGhosts(vector<ParticleType>& ngbs) const {
	  int nc = 1 ;
	  const int old_size = ngbs.size()-1;
	  const ParticleType& real_particle = ngbs.back();
	  FLOAT h2 = real_particle.hrangesqd ;

	  // Loop over the possible directions for reflections
	  for (int k = 0; k < ndim; k++){
		// Save the current number of images
		const int Nghost = nc ;

		// Do reflections on the left edge
		if (_mirror_bound[k][0]){
		  FLOAT x  = 2*_domain.min[k] - real_particle.r[k];
		  FLOAT dx = x - _cell.min[k] ;
		  if (dx*dx < h2 || x > _hbox.min[k]){
			for (int n=0; n < Nghost; n++){
			  ngbs.push_back(ngbs[n+old_size]);
			  reflect<ParticleType::NDIM>(ngbs.back(), k, _domain.min[k]) ;
			  ngbs.back().flags.set(mirror_bound_flags[k][0]) ;
			  nc++;
			}
		  }
		}
		// Do reflections on the right edge
		if (_mirror_bound[k][1]){
		  FLOAT x  = 2*_domain.max[k] - real_particle.r[k];
		  FLOAT dx = x - _cell.max[k];
		  if (dx*dx < h2 || x < _hbox.max[k]){
			for (int n=0; n < Nghost; n++){
			 ngbs.push_back(ngbs[n+old_size]);
			 reflect<ParticleType::NDIM>(ngbs.back(), k, _domain.max[k]) ;
			 ngbs.back().flags.set(mirror_bound_flags[k][1]) ;
			 nc++ ;
			}
		  }
		}
	  }
	  return nc ;
	}


    //=================================================================================================
    //  _MakeReflectedGhostsAll
    /// \brief Do the actual construction of all of the mirror ghosts. Assumes that the first
    /// particle is already saved in ngbs.
    /// \author R. A. Booth
    /// \date   27/10/2015
    /// \return The number of neighbours found
    //===============================================================================================
    template<template <int> class ParticleType>
    int _MakeReflectedGhostsAll(ParticleType<ndim>* ngbs) const {
      int nc = 1 ;
      // Loop over the possible directions for reflections
      for (int k = 0; k < ndim; k++){
        // Save the current number of images
        int Nghost = nc ;

        // Do reflections on the left edge
        if (_mirror_bound[k][0]){
          for (int n=0; n < Nghost; n++){
            ngbs[nc] = ngbs[n] ;
            reflect(ngbs[nc], k, _domain.min[k]) ;
            ngbs[nc].flags.set(mirror_bound_flags[k][0]) ;
            nc++;
          }
        }

        // Do reflections on the right edge
        if (_mirror_bound[k][1]){
          for (int n=0; n < Nghost; n++){
            ngbs[nc] = ngbs[n] ;
            reflect(ngbs[nc], k, _domain.max[k]) ;
            ngbs[nc].flags.set(mirror_bound_flags[k][1]) ;
            nc++ ;
          }
        }
      }
      return nc ;
    }
};

#endif//_GHOST_NEIGHBOURS_H_
