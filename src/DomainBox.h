#ifndef _DOMAIN_BOX__H
#define _DOMAIN_BOX__H


#include <string>
using namespace std;

#include "Precision.h"

//=============================================================================
//  Structure DomainBox
/// \brief  Bounding box data structure.
/// \author D. A. Hubber, G. Rosotti
/// \date   03/04/2013
//=============================================================================
template <int ndim>
struct DomainBox {
  string x_boundary_lhs;                ///< x-dimension LHS boundary condition
  string x_boundary_rhs;                ///< x-dimension RHS boundary condition
  string y_boundary_lhs;                ///< y-dimension LHS boundary condition
  string y_boundary_rhs;                ///< y-dimension RHS boundary condition
  string z_boundary_lhs;                ///< z-dimension LHS boundary condition
  string z_boundary_rhs;                ///< z-dimension RHS boundary condition
  FLOAT boxmin[ndim];                   ///< Minimum bounding box extent
  FLOAT boxmax[ndim];                   ///< Maximum bounding box extent
  FLOAT boxsize[ndim];                  ///< Side-lengths of bounding box
  FLOAT boxhalf[ndim];                  ///< Half side-lengths of bounding box
};

#endif
