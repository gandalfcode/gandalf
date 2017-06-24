//=================================================================================================
//  ICRegularization.cpp
//  Implementation of the initial conditions regularization
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

#include <vector>

#include "Debug.h"
#include "Nbody.h"
#include "Parameters.h"
#include "Particle.h"
#include "Ic.h"

namespace Regularization {



template<int ndim>
ParticleRegularizer<ndim>::ParticleRegularizer
(Parameters* simparams,
 const DomainBox<ndim>& icbox,
 DomainBox<ndim>& simbox_aux)
: Nreg(simparams->intparams["Nreg"]),
  localBox(icbox),
  simbox(simbox_aux)
{
}



template <int ndim>
void ParticleRegularizer<ndim>::operator()
 (Hydrodynamics<ndim>* hydro,
  NeighbourSearch<ndim> *neib,
  Nbody<ndim>* nbody,
  const RegularizerFunction<ndim>& regularizer) {
  using std::vector ;
  using std::min ;
  using std::max ;

  vector<FLOAT> rreg(ndim*hydro->Nhydromax);    // Array of particle positions

  debug1("[ICRegularization::operator()]");

  //===============================================================================================
  for (int ireg=0; ireg<Nreg; ireg++) {

    // Buid/re-build tree, create ghosts and update particle properties
    for (int i=0; i<hydro->Nhydro; i++) hydro->GetParticlePointer(i).flags.set(active);
    neib->BuildTree(true, 0, 1, 1, 0.0, hydro);
    neib->SearchBoundaryGhostParticles(0, localBox, hydro);
    neib->BuildGhostTree(true, 0, 1, 1, 0.0, hydro);
    neib->UpdateAllProperties(hydro, nbody, simbox);

    //=============================================================================================
#pragma omp parallel default(none) shared(rreg, hydro, neib, regularizer, cout)
    {
      FLOAT dr[ndim];
      FLOAT drsqd;
      vector<int> neiblist(hydro->Nhydromax);


      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (int i=0; i<hydro->Nhydro; i++) {
        Particle<ndim> &part = hydro->GetParticlePointer(i);
        const FLOAT invhsqd = (FLOAT) 1.0/(part.h*part.h);
        for (int k=0; k<ndim; k++) rreg[ndim*i + k] = (FLOAT) 0.0;

        // Find list of gather neighbours
        int Nneib = neib->GetGatherNeighbourList(part.r, hydro->kernrange*part.h,
                                                 hydro->GetParticleArrayUnsafe(), hydro->Ntot,
                                                 hydro->Nhydromax, &(neiblist[0]));


        // Loop over all neighbours and calculate position correction for regularisation
        //-----------------------------------------------------------------------------------------
        for (int jj=0; jj<Nneib; jj++) {
          const int j = neiblist[jj];
          const Particle<ndim> &neibpart = hydro->GetParticlePointer(j);

          for (int k=0; k<ndim; k++) dr[k] = neibpart.r[k] - part.r[k];
          drsqd = DotProduct(dr, dr, ndim);
          if (drsqd >= part.hrangesqd) continue;

          for (int k=0; k<ndim; k++) {
            rreg[ndim*i + k] -=
              dr[k]*hydro->kernp->w0_s2(drsqd*invhsqd)*regularizer(part, neibpart);
          }

        }
        //-----------------------------------------------------------------------------------------

        for (int k=0; k<ndim; k++) {
          rreg[ndim*i + k] = min(rreg[ndim*i + k], localBox.size[k]);
          rreg[ndim*i + k] = max(rreg[ndim*i + k], -localBox.size[k]);
        }

      }
      //-------------------------------------------------------------------------------------------


      // Apply all regularisation corrections to particle positions
      //-------------------------------------------------------------------------------------------
#pragma omp for
      for (int i=0; i<hydro->Nhydro; i++) {
        Particle<ndim> &part = hydro->GetParticlePointer(i);

        // Safety condition to enforce that a particle cannot move more than a single smoothing
        // length (i.e. roughly the interparticle separation) in one regularisation step.
        FLOAT rdiff = (FLOAT) 0.0;
        for (int k=0; k<ndim; k++) rdiff += rreg[ndim*i + k]*rreg[ndim*i + k];
        rdiff = sqrt(rdiff);
        if (rdiff > 0.5*part.h) {
          FLOAT runit[ndim];
          for (int k=0; k<ndim; k++) runit[k] = rreg[ndim*i + k]/rdiff;
          for (int k=0; k<ndim; k++) rreg[ndim*i + k] = 0.5*part.h*runit[k];
        }

        // Apply regularisation step and safely wrap particles around periodic boundaries
        for (int k=0; k<ndim; k++) {
          part.r[k] += rreg[ndim*i + k];

          // Wrap the particle positions
          while (part.r[k] > localBox.max[k]) part.r[k] -= localBox.size[k];
          while (part.r[k] < localBox.min[k]) part.r[k] += localBox.size[k];
        }
      }
    }
    //=============================================================================================

  }
  //================================================================================================

  return;
}


template class ParticleRegularizer<1>;
template class ParticleRegularizer<2>;
template class ParticleRegularizer<3>;

} // namespace Regularization
