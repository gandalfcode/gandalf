//=============================================================================
//  Ghosts.h
//  Contains definitions for ghost particle class.
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


#ifndef _GHOSTS_H_
#define _GHOSTS_H_


#include <map>
#include <string>
#include <list>
#include "Diagnostics.h"
#include "DomainBox.h"
#include "Precision.h"
#include "Parameters.h"
#include "SimUnits.h"
#include "SphKernel.h"
#include "Sph.h"
#include "Nbody.h"
#if defined MPI_PARALLEL
#include "MpiControl.h"
#endif
using namespace std;



//=============================================================================
//  Class Ghosts
/// \brief   Main ghost particle class.
/// \details Class for creating and updating ghost particles for periodic
///          boundary conditions.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class Ghosts
{
 public:

  // Main ghost particle functions
  //---------------------------------------------------------------------------
  virtual void SearchGhostParticles(FLOAT, DomainBox<ndim>, Sph<ndim> *)=0;
  virtual void CopySphDataToGhosts(DomainBox<ndim>, Sph<ndim> *)=0;
  virtual void CheckBoundaries(DomainBox<ndim>, Sph<ndim> *)=0;

//  DomainBox<ndim> simbox;               ///< Simulation boundary data
//  Sph<ndim> *sph;                       ///< SPH algorithm pointer

  static const FLOAT ghost_range = 1.8;

};

template <int ndim>
class PeriodicGhosts : public Ghosts<ndim>
{
  void CreateGhostParticle(int, int, FLOAT, FLOAT, Sph<ndim> *,int);
public:
  using Ghosts<ndim>::ghost_range;

  virtual void SearchGhostParticles(FLOAT, DomainBox<ndim>, Sph<ndim> *);
  virtual void CopySphDataToGhosts(DomainBox<ndim>, Sph<ndim> *);
  virtual void CheckBoundaries(DomainBox<ndim>, Sph<ndim> *);
};

template <int ndim>
class NullGhosts : public Ghosts<ndim>
{
public:
  using Ghosts<ndim>::ghost_range;

  virtual void SearchGhostParticles(FLOAT, DomainBox<ndim>, Sph<ndim> *);
  virtual void CopySphDataToGhosts(DomainBox<ndim>, Sph<ndim> *);
  virtual void CheckBoundaries(DomainBox<ndim>, Sph<ndim> *);
};

#if defined MPI_PARALLEL
template <int ndim>
class MPIGhosts : public Ghosts<ndim>
{
private:
  MpiControl<ndim>* mpicontrol;
public:
  using Ghosts<ndim>::ghost_range;

  MPIGhosts(MpiControl<ndim>* mpicontrol_aux): mpicontrol(mpicontrol_aux) {};

  virtual void SearchGhostParticles(FLOAT, DomainBox<ndim>, Sph<ndim> *);
  virtual void CopySphDataToGhosts(DomainBox<ndim>, Sph<ndim> *);
  virtual void CheckBoundaries(DomainBox<ndim>, Sph<ndim> *);
};
#endif

#endif
