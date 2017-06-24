//=================================================================================================
//  GradhSphSimulation.cpp
//  Contains all main functions controlling grad-h SPH simulation work-flow.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body statics And Lagrangian Fluids
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



//=================================================================================================
//  Simulation::ProcessSphParameters
/// Process all the options chosen in the parameters file, setting various
/// simulation variables and creating important simulation objects.
//=================================================================================================
template <int ndim>
void GradhSphSimulation<ndim>::ProcessSphParameters(void)
{
  aviscenum avisc = noav;              // Artificial viscosity enum
  acondenum acond = noac;              // Artificial conductivity enum
  eosenum eos_type = noeos;            // Gas EOS enum
  tdaviscenum tdavisc = notdav;        // Time-dependent viscosity enum

  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, double> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;
  string KernelName = stringparams["kernel"];
  string gas_eos = stringparams["gas_eos"];
  string gas_radiation = stringparams["radiation"];

  debug2("[GradhSphSimulation::ProcessSphParameters]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("GRADH_PROCESS_PARAMETERS");


  // Set the enum for artificial viscosity
  if (stringparams["avisc"] == "none") {
    avisc = noav;
    tdavisc = notdav;
  }
  else if (stringparams["avisc"] == "mon97" && stringparams["time_dependent_avisc"] == "mm97") {
    avisc = mon97mm97;
    tdavisc = mm97;
  }
  else if (stringparams["avisc"] == "mon97" && stringparams["time_dependent_avisc"] == "cd2010") {
    avisc = mon97cd2010;
    tdavisc = cd2010;
  }
  else if (stringparams["avisc"] == "mon97") {
    avisc = mon97;
    tdavisc = notdav;
  }
  else {
    string message = "Unrecognised parameter : avisc = " + simparams->stringparams["avisc"] +
      "   or time_dependent_avisc : " + simparams->stringparams["time_dependent_avisc"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the enum for artificial conductivity
  if (stringparams["acond"] == "none") {
    acond = noac;
  }
  else if (stringparams["acond"] == "wadsley2008") {
    acond = wadsley2008;
  }
  else if (stringparams["acond"] == "price2008") {
    acond = price2008;
  }
  else {
    string message = "Unrecognised parameter : acond = " + simparams->stringparams["acond"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set gas EOS values
  if (stringparams["gas_eos"] == "isothermal") {
    eos_type = isothermal;
  }
  else if (stringparams["gas_eos"] == "locally_isothermal") {
    eos_type = locally_isothermal;
  }
  else if (stringparams["gas_eos"] == "polytropic") {
    eos_type = polytropic;
  }
  else if (stringparams["gas_eos"] == "barotropic") {
    eos_type = barotropic;
  }
  else if (stringparams["gas_eos"] == "barotropic2") {
    eos_type = barotropic2;
  }
  else if (stringparams["gas_eos"] == "energy_eqn") {
    eos_type = energy_eqn;
  }
  else if (stringparams["gas_eos"] == "constant_temp") {
    eos_type = constant_temp;
  }
  else if (stringparams["gas_eos"] == "rad_ws" || stringparams["gas_eos"] == "radws") {
    eos_type = radws;
  }
  else if (stringparams["gas_eos"] == "disc_locally_isothermal") {
    eos_type = disc_locally_isothermal;
  }
  else {
    string message = "Unrecognised eos parameter : gas_eos = " + simparams->stringparams["gas_eos"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Create 'grad-h' SPH object depending on choice of kernel
  //===============================================================================================
  if (intparams["tabulated_kernel"] == 1) {
    sph = new GradhSph<ndim, TabulatedKernel>
      (intparams["hydro_forces"], intparams["self_gravity"], floatparams["alpha_visc"],
       floatparams["beta_visc"], floatparams["h_fac"], floatparams["h_converge"],
       avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
  }
  else if (intparams["tabulated_kernel"] == 0) {
    // Depending on the kernel, instantiate a different GradSph object
    if (KernelName == "m4") {
      sph = new GradhSph<ndim, M4Kernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["alpha_visc"],
         floatparams["beta_visc"], floatparams["h_fac"], floatparams["h_converge"],
         avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
    }
    else if (KernelName == "quintic") {
      sph = new GradhSph<ndim, QuinticKernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["alpha_visc"],
         floatparams["beta_visc"], floatparams["h_fac"], floatparams["h_converge"],
         avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
    }
    else if (KernelName == "gaussian") {
      sph = new GradhSph<ndim, GaussianKernel>
        (intparams["hydro_forces"], intparams["self_gravity"], floatparams["alpha_visc"],
         floatparams["beta_visc"], floatparams["h_fac"], floatparams["h_converge"],
         avisc, acond, tdavisc, stringparams["gas_eos"], KernelName, simunits, simparams);
    }
    else {
      string message = "Unrecognised parameter : kernel = " + simparams->stringparams["kernel"];
      ExceptionHandler::getIstance().raise(message);
    }
  }
  else {
    string message = "Invalid option for the tabulated_kernel parameter: " +
      stringparams["tabulated_kernel"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Create SPH particle integration object
  //-----------------------------------------------------------------------------------------------
  if (stringparams["sph_integration"] == "lfkdk") {
    hydroint = new SphLeapfrogKDK<ndim, GradhSphParticle>
      (floatparams["accel_mult"], floatparams["courant_mult"],
       floatparams["energy_mult"], eos_type, tdavisc);
  }
  else if (stringparams["sph_integration"] == "lfdkd") {
    hydroint = new SphLeapfrogDKD<ndim, GradhSphParticle>
      (floatparams["accel_mult"], floatparams["courant_mult"],
       floatparams["energy_mult"], eos_type, tdavisc);
    integration_step = max(integration_step,2);
  }
  else {
    string message = "Unrecognised parameter : sph_integration = "
      + simparams->stringparams["sph_integration"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Energy integration object
  //-----------------------------------------------------------------------------------------------
  if (stringparams["energy_integration"] == "Radws" ||
      stringparams["energy_integration"] == "radws"||
      stringparams["energy_integration"] == "rad_ws") {
    uint = new EnergyRadws<ndim, GradhSphParticle>
      (floatparams["energy_mult"], stringparams["radws_table"],
       floatparams["temp_ambient"], &simunits, sph->eos);
  }
  else if (stringparams["energy_integration"] == "null" ||
           stringparams["energy_integration"] == "none") {
    uint = new NullEnergy<ndim>(floatparams["energy_mult"]);
  }
  else {
    string message = "Unrecognised parameter : energy_integration = "
      + simparams->stringparams["energy_integration"];
    ExceptionHandler::getIstance().raise(message);
  }


  //-----------------------------------------------------------------------------------------------
#if defined MPI_PARALLEL
  if (stringparams["mpi_decomposition"] == "kdtree") {
    mpicontrol = new MpiKDTreeDecomposition<ndim,GradhSphParticle>(simbox);
  }
  else {
    string message = "Unrecognised parameter : mpi_decomposition = "
      + simparams->stringparams["mpi_decomposition"];
    ExceptionHandler::getIstance().raise(message);
  }
  mpicontrol->timing = timing;
  rank = mpicontrol->rank;
  Nmpi = mpicontrol->Nmpi;
#endif



  // Create neighbour searching object based on chosen method in params file
  //-----------------------------------------------------------------------------------------------
  string tree_type = stringparams["neib_search"] ;
  if (gas_radiation == "tree_ray") {
    if (tree_type == "octtree" && ndim == 3) {
      tree_type = "treeray" ;
    }
    else {
      string message = "Error: Tree Ray needs an Oct Tree and ndim=3";
      ExceptionHandler::getIstance().raise(message);
    }
  }

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


  sphneib = new GradhSphTree<ndim,GradhSphParticle>
    (tree_type, intparams["Nleafmax"], Nmpi, intparams["pruning_level_min"], intparams["pruning_level_max"],
     floatparams["thetamaxsqd"], sph->kernp->kernrange, floatparams["macerror"],
     stringparams["gravity_mac"], multipole, &simbox, sph->kernp, timing, sph->types);

  // Here I do a horrible hack to get at the underlying tree, needed for the dust.
  TreeBase<ndim> * t = NULL, * gt = NULL, *mpit = NULL ;
  typedef GradhSphTree<ndim,GradhSphParticle> TreeType ;
  TreeType *pTree = reinterpret_cast<TreeType*>(sphneib) ;
  t = pTree->tree ; gt = pTree->ghosttree ;

#ifdef MPI_PARALLEL
  mpit = pTree->mpighosttree ;
#endif


#if defined MPI_PARALLEL
  mpicontrol->SetNeibSearch(sphneib);
#endif


  // Radiation transport object
  //-----------------------------------------------------------------------------------------------
  if (gas_radiation == "treeray" && ndim == 3) {
    radiation = new TreeRay<ndim,1,GradhSphParticle,TreeRayCell>
      (Nmpi, intparams["on_the_spot"], intparams["nside"], intparams["ilNR"], intparams["ilNTheta"],
       intparams["ilNPhi"], intparams["ilNNS"], intparams["ilFinePix"], floatparams["maxDist"],
       floatparams["rayRadRes"], floatparams["relErr"], stringparams["errControl"],
       simbox, &simunits, simparams, sphneib);
  }
  else if (gas_radiation == "ionisation") {
    radiation = new MultipleSourceIonisation<ndim,GradhSphParticle>
      (sphneib, floatparams["mu_bar"], floatparams["mu_ion"], floatparams["temp0"],
       floatparams["temp_ion"], floatparams["Ndotmin"], floatparams["gamma_eos"],
       pow(simunits.r.outscale*simunits.r.outcgs, 3.)/
       pow(simunits.m.outscale*simunits.m.outcgs, 2.),
       simunits.temp.outscale, pow(simunits.r.outscale*simunits.r.outcgs,-4)*
       pow(simunits.t.outscale*simunits.t.outcgs,+2)/simunits.m.outscale*simunits.m.outcgs);
  }
  else if (gas_radiation == "none") {
    radiation = new NullRadiation<ndim>();
  }
  else {
    string message = "Unrecognised parameter : radiation = " + gas_radiation;
    ExceptionHandler::getIstance().raise(message);
  }


  // Create ghost particle object
  //-----------------------------------------------------------------------------------------------
  if (IsAnyBoundarySpecial(simbox)) {
    LocalGhosts = new PeriodicGhostsSpecific<ndim,GradhSphParticle >();
  }
  else {
    LocalGhosts = new NullGhosts<ndim>();
  }
#ifdef MPI_PARALLEL
  MpiGhosts = new MpiGhostsSpecific<ndim, GradhSphParticle>(mpicontrol);
#endif


  // Depending on the dimensionality, calculate expected neighbour number
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    sph->Ngather = (int) (2.0*sph->kernp->kernrange*sph->h_fac);
  }
  else if (ndim == 2) {
    sph->Ngather = (int) (pi*pow(sph->kernp->kernrange*sph->h_fac, 2));
  }
  else if (ndim == 3) {
    sph->Ngather = (int) (4.0*pi*pow(sph->kernp->kernrange*sph->h_fac, 3)/3.0);
  }

  // Setup the dust force object
  //-----------------------------------------------------------------------------------------------
  sphdust = DustFactory<ndim, GradhSphParticle>::ProcessParameters(simparams, timing, sph->types,
                                                                   simbox, t, gt, mpit) ;

#if defined MPI_PARALLEL
  if (sphdust != NULL) sphdust->SetMpiControl(mpicontrol);
#endif

  return;
}



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class GradhSphSimulation<1>;
template class GradhSphSimulation<2>;
template class GradhSphSimulation<3>;
