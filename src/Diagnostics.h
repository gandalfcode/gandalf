#ifndef _PRECISION__H
#define _PRECISION__H

#include "Precision.h"

//=============================================================================
//  Structure Diagnostics
/// \brief  Structure containing snapshot of current diagnostic quantities.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
template <int ndim>
struct Diagnostics {
  DOUBLE Eerror;                    ///< Total energy error
  DOUBLE Etot;                      ///< Total energy
  DOUBLE utot;                      ///< Total thermal energy
  DOUBLE ketot;                     ///< Total kinetic energy
  DOUBLE gpetot;                    ///< Total grav. potential energy
  DOUBLE mom[ndim];                 ///< Total momentum vector
  DOUBLE angmom[3];                 ///< Total angular momentum vector
  DOUBLE force[ndim];               ///< Net force
  DOUBLE force_grav[ndim];          ///< Net gravitational force
  DOUBLE rcom[ndim];                ///< Position of centre of mass
  DOUBLE vcom[ndim];                ///< Velocity of centre of mass
};

#endif
