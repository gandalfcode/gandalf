//=============================================================================
//  SM2012SphSimulation.cpp
//  Contains all main functions controlling Saitoh & Makino (2012) SPH
//  simulation work-flow.
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


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <time.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
#include "CodeTiming.h"
#include "Exception.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "Simulation.h"
#include "Parameters.h"
#include "Nbody.h"
#include "Sph.h"
#include "RiemannSolver.h"
#include "Ghosts.h"
#include "Sinks.h"
using namespace std;



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class SM2012SphSimulation<1>;
template class SM2012SphSimulation<2>;
template class SM2012SphSimulation<3>;




//=============================================================================
//  SM2012SphSimulation::ProcessSphParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
//=============================================================================
template <int ndim>
void SM2012SphSimulation<ndim>::ProcessSphParameters(void)
{
  aviscenum avisc = noav;              // Artificial viscosity enum
  acondenum acond = noac;              // Artificial conductivity enum
  eosenum gas_eos = noeos;             // Gas EOS enum
  tdaviscenum tdavisc = notdav;        // Time-dependent viscosity enum

  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;
  string KernelName = stringparams["kernel"];


  // Set the enum for artificial viscosity
  if (stringparams["avisc"] == "none")
    avisc = noav;
  else if (stringparams["avisc"] == "mon97" &&
           stringparams["time_dependent_avisc"] == "mm97")
    avisc = mon97mm97;
  else if (stringparams["avisc"] == "mon97" &&
           stringparams["time_dependent_avisc"] == "cd2010")
    avisc = mon97cd2010;
  else if (stringparams["avisc"] == "mon97")
    avisc = mon97;
  else {
    string message = "Unrecognised parameter : avisc = " +
      simparams->stringparams["avisc"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the enum for artificial viscosity
  if (stringparams["time_dependent_avisc"] == "none")
    tdavisc = notdav;
  else if (stringparams["time_dependent_avisc"] == "mm97")
    tdavisc = mm97;
  else if (stringparams["time_dependent_avisc"] == "cd2010")
    tdavisc = cd2010;
  else {
    string message = "Unrecognised parameter : time_dependent_avisc = " +
      simparams->stringparams["time_dependent_avisc"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the enum for artificial conductivity
  if (stringparams["acond"] == "none")
    acond = noac;
  else if (stringparams["acond"] == "wadsley2008")
    acond = wadsley2008;
  else if (stringparams["acond"] == "price2008")
    acond = price2008;
  else {
    string message = "Unrecognised parameter : acond = " +
        simparams->stringparams["acond"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set gas EOS values
  if (stringparams["gas_eos"] == "isothermal")
    gas_eos = isothermal;
  else if (stringparams["gas_eos"] == "barotropic")
    gas_eos = barotropic;
  else if (stringparams["gas_eos"] == "energy_eqn")
    gas_eos = energy_eqn;
  else if (stringparams["gas_eos"] == "constant_temp")
    gas_eos = constant_temp;
  else {
    string message = "Unrecognised parameter : gas_eos = " +
        simparams->stringparams["gas_eos"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Create Saitoh-Makino (2012) SPH object
  //===========================================================================
  if (stringparams["sph"] == "sm2012") {
    if (intparams["tabulated_kernel"] == 1) {
      sph = new SM2012Sph<ndim, TabulatedKernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["alpha_visc"],
         floatparams["beta_visc"], floatparams["h_fac"], floatparams["h_converge"],
         avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
    }
    else if (intparams["tabulated_kernel"] == 0){
      // Depending on the kernel, instantiate a different SM2012 object
      if (KernelName == "m4") {
	sph = new SM2012Sph<ndim, M4Kernel>
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
      }
      else if (KernelName == "quintic") {
	sph = new SM2012Sph<ndim, QuinticKernel>
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
      }
      else if (KernelName == "gaussian") {
	sph = new SM2012Sph<ndim, GaussianKernel>
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
      }
      else {
	string message = "Unrecognised parameter : kernel = " +
	  simparams->stringparams["kernel"];
	ExceptionHandler::getIstance().raise(message);
      }
    }
    else {
      string message = "Invalid option for the tabulated_kernel parameter: " +
	stringparams["tabulated_kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }

  //===========================================================================
  else {
    string message = "Invalid or unrecognised parameter : sph = "
      + simparams->stringparams["sph"];
    ExceptionHandler::getIstance().raise(message);
  }
  //===========================================================================


  // Create SPH particle integration object
  //---------------------------------------------------------------------------
  if (stringparams["sph_integration"] == "lfkdk") {
    hydroint = new SphLeapfrogKDK<ndim, SM2012SphParticle>(floatparams["accel_mult"],
			              floatparams["courant_mult"],
			              floatparams["energy_mult"],
				      gas_eos, tdavisc);
  }
  else if (stringparams["sph_integration"] == "lfdkd") {
    hydroint = new SphLeapfrogDKD<ndim, SM2012SphParticle>(floatparams["accel_mult"],
			              floatparams["courant_mult"],
			              floatparams["energy_mult"],
				      gas_eos, tdavisc);
    integration_step = max(integration_step,2);
  }
  else {
    string message = "Unrecognised parameter : sph_integration = "
      + simparams->stringparams["sph_integration"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Energy integration object
  //---------------------------------------------------------------------------
  if (stringparams["energy_integration"] == "null" ||
      stringparams["energy_integration"] == "none") {
    uint = new NullEnergy<ndim>(floatparams["energy_mult"]);
  }
  else {
    string message = "Unrecognised parameter : energy_integration = "
      + simparams->stringparams["energy_integration"];
    ExceptionHandler::getIstance().raise(message);
  }

#if defined MPI_PARALLEL

  //mpicontrol = new MpiControlType<ndim, SM2012SphParticle>;

#endif

  // Create neighbour searching object based on chosen method in params file
  //-------------------------------------------------------------------------
  string tree_type = stringparams["neib_search"];

  multipole_method multipole = monopole ;
  if (stringparams["multipole"] == "monopole")
    multipole = monopole ;
  else if (stringparams["multipole"] == "quadrupole")
    multipole = quadrupole ;
  else if (stringparams["multipole"] == "fast_monopole")
    multipole = fast_monopole ;
  else if (stringparams["multipole"] == "fast_quadrupole")
    multipole = fast_quadrupole ;
  else {
    string message = "Multipole type not recognised.";
    ExceptionHandler::getIstance().raise(message);
  }

  sphneib = new SM2012SphTree<ndim,SM2012SphParticle>
     (tree_type, intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
      floatparams["thetamaxsqd"], sph->kernp->kernrange, floatparams["macerror"],
      stringparams["gravity_mac"], multipole, &simbox, sph->kernp, timing, sph->types);


#if defined MPI_PARALLEL
  mpicontrol->SetNeibSearch(sphneib);
#endif

  // Depending on the dimensionality, calculate expected neighbour number
  //---------------------------------------------------------------------------
  if (ndim == 1)
    sph->Ngather = (int) (2.0*sph->kernp->kernrange*sph->h_fac);
  else if (ndim == 2)
    sph->Ngather = (int) (pi*pow(sph->kernp->kernrange*sph->h_fac,2));
  else if (ndim == 3)
    sph->Ngather = (int) (4.0*pi*pow(sph->kernp->kernrange*sph->h_fac,3)/3.0);

  // Create ghosts object
  if (IsAnyBoundarySpecial(simbox))
    LocalGhosts = new PeriodicGhostsSpecific<ndim,SM2012SphParticle >();
  else
    LocalGhosts = new NullGhosts<ndim>();
#ifdef MPI_PARALLEL
  MpiGhosts = new MpiGhostsSpecific<ndim, SM2012SphParticle>(mpicontrol);
#endif


  return;
}
