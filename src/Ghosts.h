//=============================================================================
//  Ghosts.h
//  Contains definitions for ghost particle class.
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

  Ghosts();
  ~Ghosts();

  // Main ghost particle functions
  // --------------------------------------------------------------------------
  void SearchGhostParticles(DomainBox<ndim>, Sph<ndim> *);
  void CreateGhostParticle(int, int, FLOAT, FLOAT, Sph<ndim> *);
  void CopySphDataToGhosts(Sph<ndim> *);
  void CheckBoundaries(DomainBox<ndim>, Sph<ndim> *);

  DomainBox<ndim> simbox;               ///< Simulation boundary data
  Sph<ndim> *sph;                       ///< SPH algorithm pointer

  static const FLOAT ghost_range = 1.1;

};
#endif
