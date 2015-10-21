//=================================================================================================
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
//=================================================================================================


#ifndef _GHOSTS_H_
#define _GHOSTS_H_


#include <map>
#include <string>
#include <list>
#include "Diagnostics.h"
#include "DomainBox.h"
#include "Hydrodynamics.h"
#include "Precision.h"
#include "Parameters.h"
#include "SimUnits.h"
#include "SmoothingKernel.h"
#include "Nbody.h"
#if defined MPI_PARALLEL
#include "MpiControl.h"
#endif
using namespace std;



//=================================================================================================
//  Class Ghosts
/// \brief   Main ghost particle class.
/// \details Class for creating and updating ghost particles for periodic boundary conditions.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class Ghosts
{
 public:

  // Main ghost particle functions
  //-----------------------------------------------------------------------------------------------
  virtual void SearchGhostParticles(FLOAT, DomainBox<ndim>, Hydrodynamics<ndim> *)=0;
  virtual void CopyHydroDataToGhosts(DomainBox<ndim>, Hydrodynamics<ndim> *)=0;
  virtual void CheckBoundaries(DomainBox<ndim>, Hydrodynamics<ndim> *)=0;

//  DomainBox<ndim> simbox;               ///< Simulation boundary data
//  Hydrodynamics<ndim> *sph;                       ///< SPH algorithm pointer

  static const FLOAT ghost_range; //= 1.6;

};



//=================================================================================================
//  Class PeriodicGhosts
/// \brief   ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class PeriodicGhosts : public Ghosts<ndim>
{
public:
  using Ghosts<ndim>::ghost_range;

  virtual void CheckBoundaries(DomainBox<ndim>, Hydrodynamics<ndim> *);
};



//=================================================================================================
//  Class PeriodicGhostsSpecific
/// \brief   ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template <int> class ParticleType>
class PeriodicGhostsSpecific : public PeriodicGhosts<ndim>
{
public:
  using Ghosts<ndim>::ghost_range;
  
  virtual void CopyHydroDataToGhosts(DomainBox<ndim>, Hydrodynamics<ndim> *);
  virtual void SearchGhostParticles(FLOAT, DomainBox<ndim>, Hydrodynamics<ndim> *) {};
};



//=================================================================================================
//  Class NullGhosts
/// \brief   ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class NullGhosts : public Ghosts<ndim>
{
public:
  using Ghosts<ndim>::ghost_range;

  virtual void SearchGhostParticles(FLOAT, DomainBox<ndim>, Hydrodynamics<ndim> *) {};
  virtual void CopyHydroDataToGhosts(DomainBox<ndim>, Hydrodynamics<ndim> *);
  virtual void CheckBoundaries(DomainBox<ndim>, Hydrodynamics<ndim> *);
};



#if defined MPI_PARALLEL
//=================================================================================================
//  Class MpiGhosts
/// \brief   ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class MpiGhosts : public Ghosts<ndim>
{
public:
  using Ghosts<ndim>::ghost_range;

  //MpiGhosts(MpiControl<ndim>* mpicontrol_aux): mpicontrol(mpicontrol_aux) {};
  MpiGhosts(){};

  virtual void CheckBoundaries(DomainBox<ndim>, Hydrodynamics<ndim> *);
};



//=================================================================================================
//  Class MpiGhostsSpecific
/// \brief   ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim, template <int> class ParticleType>
class MpiGhostsSpecific : public MpiGhosts<ndim> {

  MpiControlType<ndim,ParticleType>* mpicontrol;
  //using MpiGhosts<ndim>::mpicontrol;

public:
  virtual void SearchGhostParticles(FLOAT, DomainBox<ndim>, Hydrodynamics<ndim> *);
  virtual void CopyHydroDataToGhosts(DomainBox<ndim>, Hydrodynamics<ndim> *);

  MpiGhostsSpecific(MpiControl<ndim>* mpicontrol_aux): mpicontrol(static_cast<MpiControlType<ndim,ParticleType> *> (mpicontrol_aux)) {};

};
#endif

#endif
