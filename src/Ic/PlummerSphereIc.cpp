//=================================================================================================
//  PlummerSphereIc.cpp
//  Class for generating initial conditions for ...
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


#include <fstream>
#include <sstream>
#include "Precision.h"
#include "Debug.h"
#include "Ic.h"
using namespace std;



//=================================================================================================
//  PlummerSphereIc::PlummerSphereIc
/// Set-up SILCC-type simulation initial conditions.
//=================================================================================================
template <int ndim>
PlummerSphereIc<ndim>::PlummerSphereIc(Simulation<ndim>* _sim, FLOAT _invndim) :
  Ic<ndim>(_sim, _invndim)
{
  if (simparams->intparams["ndim"] != 3) {
    ExceptionHandler::getIstance().raise("Plummer sphere can only be generated in 3d");
  }
}



//=================================================================================================
//  PlummerSphereIc::Generate
/// ...
//=================================================================================================
template <int ndim>
void PlummerSphereIc<ndim>::Generate(void)
{
  Nbody<ndim>* nbody = sim->nbody;

  // Only compile for 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    bool flag;                           // Aux. flag
    int i,j,k;                           // Particle and dimension counter
    FLOAT raux;                          // Aux. float variable
    //FLOAT rcentre[ndim];                 // Position of centre of Plummer sphere
    FLOAT vplummer;                      // Velocity of Plummer components
    FLOAT x1,x2,x3,x4,x5,x6,x7;          // Aux. random number variables
    FLOAT rad,vm,ve,t1,t2,w,z;           // Other variables

    // Local copies of important parameters
    int Nhydro      = simparams->intparams["Nhydro"];
    int Nstar       = simparams->intparams["Nstar"];
    FLOAT gamma_eos = simparams->floatparams["gamma_eos"];
    FLOAT gasfrac   = simparams->floatparams["gasfrac"];
    FLOAT starfrac  = simparams->floatparams["starfrac"];
    FLOAT mplummer  = simparams->floatparams["mplummer"];
    FLOAT rplummer  = simparams->floatparams["rplummer"];
    FLOAT radius    = simparams->floatparams["radius"];
    FLOAT rstar     = simparams->floatparams["rstar"];

    debug1("[Ic::PlummerSphere]");

    hydro->Nhydro = Nhydro;
    hydro->Ntot = Nhydro;
    nbody->Nstar = Nstar;
    sim->AllocateParticleMemory();

    //for (k=0; k<ndim; k++) rcentre[k] = 0.0;
    raux = gasfrac + starfrac;
    gasfrac /= raux;
    starfrac /= raux;


    // Loop over all particles (gas and stars)
    //=============================================================================================
    for (j=0; j<Nhydro+Nstar; j++) {

      do {
        flag = false;
        x1 = sim->randnumb->floatrand();
        x2 = sim->randnumb->floatrand();
        x3 = sim->randnumb->floatrand();

        if (x1 == (FLOAT) 0.0 && x2 == (FLOAT) 0.0 && x3 == (FLOAT) 0.0) flag = true;
        rad = (FLOAT) 1.0 / sqrt(pow(x1,-(FLOAT) 2.0/ (FLOAT) 3.0) - (FLOAT) 1.0);
        if (rad > radius/rplummer) flag = true;
      } while (flag);

      z = ((FLOAT) 1.0 - (FLOAT) 2.0*x2)*rad;


      // Set position depending on particle type
      //-------------------------------------------------------------------------------------------
      if (j >= Nstar && j < Nstar + Nhydro) {
        i = j - Nstar;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        part.r[2] = z;
        part.r[0] = sqrt(rad*rad - z*z)*cos(twopi*x3);
        part.r[1] = sqrt(rad*rad - z*z)*sin(twopi*x3);
        part.m = gasfrac / (FLOAT) Nhydro;
      }
      else {
        i = j;
        nbody->stardata[i].r[0] = sqrt(rad*rad - z*z)*cos(twopi*x3);
        nbody->stardata[i].r[1] = sqrt(rad*rad - z*z)*sin(twopi*x3);
        nbody->stardata[i].r[2] = z;
        nbody->stardata[i].m    = starfrac / (FLOAT) Nstar;
      }

      // Maximum velocity for this distance
      ve = sqrt((FLOAT) 2.0 / sqrt((FLOAT) 1.0 + rad*rad));


      // Velocity of particle
      //-------------------------------------------------------------------------------------------
      do {
        x4 = sim->randnumb->floatrand();
        x5 = sim->randnumb->floatrand();
        t1 = (FLOAT) 0.1*x5;
        t2 = x4*x4*pow((FLOAT) 1.0 - x4*x4,(FLOAT) 3.5);
      } while (t1 > t2);

      vm = ve*x4;
      x6 = sim->randnumb->floatrand();
      x7 = sim->randnumb->floatrand();
      w = ((FLOAT) 1.0 - (FLOAT) 2.0*x6)*vm;


      // Set velocity depending on particle type
      //-------------------------------------------------------------------------------------------
      if (j >= Nstar && j < Nstar + Nhydro) {
        i = j - Nstar;
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (k=0; k<ndim; k++) part.v[k] = (FLOAT) 0.0;
        FLOAT sound = sqrt((FLOAT) 0.16666666666666666 / sqrt((FLOAT) 1.0 + rad*rad));
        part.rho   = (FLOAT) 1.0;
        part.u     = sound*sound/(gamma_eos - (FLOAT) 1.0);
      }
      else {
        i = j;
        nbody->stardata[i].v[0] = sqrt(vm*vm - w*w)*cos(twopi*x7);
        nbody->stardata[i].v[1] = sqrt(vm*vm - w*w)*sin(twopi*x7);
        nbody->stardata[i].v[2] = w;
      }

    }
    //=============================================================================================


    // Instanly move to COM
    //ConvertToComFrame();
    vplummer = sqrt(mplummer/rplummer);

    // Now scale variables to required physical size
    for (i=0; i<Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (k=0; k<ndim; k++) {
        part.r[k] = part.r[k]*rplummer;
        part.v[k] = part.v[k]*vplummer;
      }
      part.m = part.m*mplummer;
      if (i < Nhydro) part.u = part.u*(mplummer/rplummer);
    }

    for (i=0; i<Nstar; i++) {
      for (k=0; k<ndim; k++) {
        nbody->stardata[i].r[k] = nbody->stardata[i].r[k]*rplummer;
        nbody->stardata[i].v[k] = nbody->stardata[i].v[k]*vplummer;
      }
      nbody->stardata[i].m      = nbody->stardata[i].m*mplummer;
      nbody->stardata[i].radius = rstar;
      nbody->stardata[i].h      = nbody->kernp->invkernrange*rstar;
      nbody->stardata[i].invh   = 1.0 / (nbody->stardata[i].h + small_number_dp);
    }

  }

  return;
}



template class PlummerSphereIc<1>;
template class PlummerSphereIc<2>;
template class PlummerSphereIc<3>;
