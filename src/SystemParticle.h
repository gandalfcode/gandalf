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
static const int Npertmax = 4;



//=============================================================================
//  Structure SystemParticle
/// \brief  System particle data structure
/// \author D. A. Hubber, G. Rosotti
/// \date   10/05/2013
//=============================================================================
template <int ndim>
class SystemParticle: public NbodyParticle<ndim>
{
public:
  int inode;                                ///< NN-tree node id
  int Nchildren;                            ///< Number of nbody children
  int Npert;                                ///< Number of perturbers
  NbodyParticle<ndim>* children[Ncompmax];  ///< Array of ptrs to children
  NbodyParticle<ndim>* perturber[Npertmax]; ///< Array of ptrs to perturbers

};
#endif
