//=============================================================================
//  SystemParticle.h
//  Main system particle data structure
//=============================================================================


#ifndef _SYSTEM_PARTICLE_H_
#define _SYSTEM_PARTICLE_H_


#include "Precision.h"
#include "Constants.h"
#include "NbodyParticle.h"


static const int Ncompmax = 4;



//=============================================================================
//  Structure SystemParticle
/// \brief  System particle data structure
/// \author D. A. Hubber
/// \date   10/05/2013
//=============================================================================
template <int ndim>
class SystemParticle: public NbodyParticle<ndim>
{
public:
  int Nchildren; ///Number of nbody children
  NbodyParticle<ndim>* children[Ncompmax]; ///Array of pointers to children

};
#endif
