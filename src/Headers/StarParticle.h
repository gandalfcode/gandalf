//=============================================================================
//  StarParticle.h
//  Main star particle data structure
//=============================================================================


#ifndef _STAR_PARTICLE_H_
#define _STAR_PARTICLE_H_


#include "Precision.h"
#include "Constants.h"
#include "NbodyParticle.h"


//=============================================================================
//  Structure StarParticle
/// \brief  Star particle data structure
/// \author D. A. Hubber
/// \date   15/04/2013
//=============================================================================
template <int ndim>
class StarParticle: public NbodyParticle<ndim>
{
  //DOUBLE Tsurface;                  ///< Surface temperature
  //DOUBLE L;                         ///< Luminosity
};
#endif
