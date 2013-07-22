//=============================================================================
//  NbodySimulation.cpp
//  Contains all main functions controlling SPH simulation work-flow.
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
#include "Exception.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "Simulation.h"
#include "Parameters.h"
#include "Nbody.h"
#include "Sph.h"
#include "Ghosts.h"
using namespace std;



// Create template class instances of the main NbodySimulation object for
// each dimension used (1, 2 and 3)
template class NbodySimulation<1>;
template class NbodySimulation<2>;
template class NbodySimulation<3>;



//=============================================================================
//  NbodySimulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various 
/// simulation variables and creating important simulation objects.
//=============================================================================
template <int ndim>
void NbodySimulation<ndim>::ProcessParameters(void)
{
  aviscenum avisc;                  // Artificial viscosity enum
  acondenum acond;                  // Artificial conductivity enum

  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, float> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;

  debug2("[NbodySimulation::ProcessParameters]");

  // Sanity check for valid dimensionality
  if (ndim < 1 || ndim > 3) {
    string message = "Invalid dimensionality chosen : ndim = " + ndim;
    ExceptionHandler::getIstance().raise(message);
  }

  // Set-up all output units for scaling parameters
  if (intparams["dimensionless"] == 0) {
    simunits.r.outunit = stringparams["routunit"];
    simunits.m.outunit = stringparams["moutunit"];
    simunits.t.outunit = stringparams["toutunit"];
    simunits.v.outunit = stringparams["voutunit"];
    simunits.a.outunit = stringparams["aoutunit"];
    simunits.rho.outunit = stringparams["rhooutunit"];
    simunits.press.outunit = stringparams["pressoutunit"];
    simunits.f.outunit = stringparams["foutunit"];
    simunits.E.outunit = stringparams["Eoutunit"];
    simunits.mom.outunit = stringparams["momoutunit"];
    simunits.angmom.outunit = stringparams["angmomoutunit"];
    simunits.angvel.outunit = stringparams["angveloutunit"];
    simunits.dmdt.outunit = stringparams["dmdtoutunit"];
    simunits.u.outunit = stringparams["uoutunit"];
    simunits.dudt.outunit = stringparams["dudtoutunit"];
    simunits.temp.outunit = stringparams["tempoutunit"];
    simunits.SetupUnits(simparams);
  }

  // Create SPH object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["sph"] == "gradh") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      sph = new GradhSph<ndim, TabulatedKernel> 
	(intparams["hydro_forces"], intparams["self_gravity"],
	 floatparams["alpha_visc"], floatparams["beta_visc"],
	 floatparams["h_fac"], floatparams["h_converge"], 
	 avisc, acond, stringparams["gas_eos"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	sph = new GradhSph<ndim, M4Kernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "quintic") {
	sph = new GradhSph<ndim, QuinticKernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "gaussian") {
	sph = new GradhSph<ndim, GaussianKernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
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
  // --------------------------------------------------------------------------
  else {
    string message = "Invalid or unrecognised parameter : sph = " 
      + simparams->stringparams["sph"];
    ExceptionHandler::getIstance().raise(message);
  }
  // --------------------------------------------------------------------------


  // Create N-body object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["nbody"] == "lfkdk") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyLeapfrogKDK<ndim, TabulatedKernel> 
	(intparams["nbody_softening"], intparams["sub_systems"],
	 floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	nbody = new NbodyLeapfrogKDK<ndim, M4Kernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
	nbody = new NbodyLeapfrogKDK<ndim, QuinticKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
	nbody = new NbodyLeapfrogKDK<ndim, GaussianKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
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
  // --------------------------------------------------------------------------
  else if (stringparams["nbody"] == "hermite4") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite4<ndim, TabulatedKernel> 
	(intparams["nbody_softening"], intparams["sub_systems"],
	 floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	nbody = new NbodyHermite4<ndim, M4Kernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
	nbody = new NbodyHermite4<ndim, QuinticKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
	nbody = new NbodyHermite4<ndim, GaussianKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
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
  // --------------------------------------------------------------------------
  else if (stringparams["nbody"] == "hermite4ts") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyHermite4TS<ndim, TabulatedKernel>
	(intparams["nbody_softening"], intparams["sub_systems"],
	 floatparams["nbody_mult"], KernelName, intparams["Npec"]);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	nbody = new NbodyHermite4TS<ndim, M4Kernel>
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "quintic") {
	nbody = new NbodyHermite4TS<ndim, QuinticKernel>
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (KernelName == "gaussian") {
	nbody = new NbodyHermite4TS<ndim, GaussianKernel>
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName, intparams["Npec"]);
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
  // --------------------------------------------------------------------------
  else {
    string message = "Unrecognised parameter : nbody = " 
      + simparams->stringparams["nbody"];
    ExceptionHandler::getIstance().raise(message);
  }
  // --------------------------------------------------------------------------



  // Create sub-system object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (intparams["sub_systems"] == 1) {

    // ------------------------------------------------------------------------
    if (stringparams["sub_system_integration"] == "lfkdk") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
	subsystem = new NbodyLeapfrogKDK<ndim, TabulatedKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["subsys_mult"], KernelName);
      }
      else if (intparams["tabulated_kernel"] == 0) {
	// Depending on the kernel, instantiate a different GradSph object
	if (KernelName == "m4") {
	  subsystem = new NbodyLeapfrogKDK<ndim, M4Kernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["subsys_mult"], KernelName);
	}
	else if (KernelName == "quintic") {
	  subsystem = new NbodyLeapfrogKDK<ndim, QuinticKernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["subsys_mult"], KernelName);
	}
	else if (KernelName == "gaussian") {
	  subsystem = new NbodyLeapfrogKDK<ndim, GaussianKernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["subsys_mult"], KernelName);
	}
	else {
	  string message = "Unrecognised parameter : kernel = " +
	    simparams->stringparams["kernel"];
	  ExceptionHandler::getIstance().raise(message);
	}
      }
      else {
	string message = "Invalid option for the tabulated_kernel parameter: "
	  + stringparams["tabulated_kernel"];
	ExceptionHandler::getIstance().raise(message);
      }
    }
    // ------------------------------------------------------------------------
    else if (stringparams["sub_system_integration"] == "hermite4") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
	subsystem = new NbodyHermite4<ndim, TabulatedKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["subsys_mult"], KernelName);
      }
      else if (intparams["tabulated_kernel"] == 0) {
	// Depending on the kernel, instantiate a different GradSph object
	if (KernelName == "m4") {
	  subsystem = new NbodyHermite4<ndim, M4Kernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["subsys_mult"], KernelName);
	}
	else if (KernelName == "quintic") {
	  subsystem = new NbodyHermite4<ndim, QuinticKernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["subsys_mult"], KernelName);
	}
	else if (KernelName == "gaussian") {
	  subsystem = new NbodyHermite4<ndim, GaussianKernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["subsys_mult"], KernelName);
	}
	else {
	string message = "Unrecognised parameter : kernel = " +
	  simparams->stringparams["kernel"];
	ExceptionHandler::getIstance().raise(message);
	}
      }
      else {
	string message = "Invalid option for the tabulated_kernel parameter: "
	  + stringparams["tabulated_kernel"];
	ExceptionHandler::getIstance().raise(message);
      }
    }
    // ------------------------------------------------------------------------
    else if (stringparams["sub_system_integration"] == "hermite4ts") {
      string KernelName = stringparams["kernel"];
      if (intparams["tabulated_kernel"] == 1) {
	subsystem = new NbodyHermite4TS<ndim, TabulatedKernel>
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["subsys_mult"], KernelName, intparams["Npec"]);
      }
      else if (intparams["tabulated_kernel"] == 0) {
	// Depending on the kernel, instantiate a different GradSph object
	if (KernelName == "m4") {
	  subsystem = new NbodyHermite4TS<ndim, M4Kernel>
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["subsys_mult"], KernelName, intparams["Npec"]);
	}
	else if (KernelName == "quintic") {
	  subsystem = new NbodyHermite4TS<ndim, QuinticKernel>
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["subsys_mult"], KernelName, intparams["Npec"]);
	}
	else if (KernelName == "gaussian") {
	  subsystem = new NbodyHermite4TS<ndim, GaussianKernel>
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["subsys_mult"], KernelName, intparams["Npec"]);
	}
	else {
	  string message = "Unrecognised parameter : kernel = " +
	    simparams->stringparams["kernel"];
	  ExceptionHandler::getIstance().raise(message);
	}
      }
      else {
	string message = "Invalid option for the tabulated_kernel parameter: "
	  + stringparams["tabulated_kernel"];
	ExceptionHandler::getIstance().raise(message);
      }
    }
    // ------------------------------------------------------------------------
    else {
      string message = "Unrecognised parameter : sub_system_integration = " 
	+ simparams->stringparams["sub_system_integration"];
      ExceptionHandler::getIstance().raise(message);
    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------


  // Boundary condition variables
  // --------------------------------------------------------------------------
  simbox.x_boundary_lhs = stringparams["x_boundary_lhs"];
  simbox.x_boundary_rhs = stringparams["x_boundary_rhs"];
  simbox.y_boundary_lhs = stringparams["y_boundary_lhs"];
  simbox.y_boundary_rhs = stringparams["y_boundary_rhs"];
  simbox.z_boundary_lhs = stringparams["z_boundary_lhs"];
  simbox.z_boundary_rhs = stringparams["z_boundary_rhs"];
  simbox.boxmin[0] = floatparams["boxmin[0]"]/simunits.r.outscale;
  simbox.boxmin[1] = floatparams["boxmin[1]"]/simunits.r.outscale;
  simbox.boxmin[2] = floatparams["boxmin[2]"]/simunits.r.outscale;
  simbox.boxmax[0] = floatparams["boxmax[0]"]/simunits.r.outscale;
  simbox.boxmax[1] = floatparams["boxmax[1]"]/simunits.r.outscale;
  simbox.boxmax[2] = floatparams["boxmax[2]"]/simunits.r.outscale;
  for (int k=0; k<3; k++) {
    simbox.boxsize[k] = simbox.boxmax[k] - simbox.boxmin[k];
    simbox.boxhalf[k] = 0.5*simbox.boxsize[k];
  }


  // Set all other parameter variables
  // --------------------------------------------------------------------------
  nbody_single_timestep = intparams["nbody_single_timestep"];
  nbody->Nstar          = intparams["Nstar"];
  nbodytree.gpefrac     = floatparams["gpefrac"];

  dt_snap               = floatparams["dt_snap"]/simunits.t.outscale;
  dt_python             = floatparams["dt_python"];
  Nlevels               = intparams["Nlevels"];
  noutputstep           = intparams["noutputstep"];
  Nstepsmax             = intparams["Nstepsmax"];
  out_file_form         = stringparams["out_file_form"];
  run_id                = stringparams["run_id"];
  tend                  = floatparams["tend"]/simunits.t.outscale;
  tsnapnext             = floatparams["tsnapfirst"]/simunits.t.outscale;

  // Flag that we've processed all parameters already
  ParametersProcessed = true;

  return;
}



//=============================================================================
//  NbodySimulation::PostInitialConditionsSetup
/// ..
//=============================================================================
template <int ndim>
void NbodySimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[NbodySimulation::PostInitialConditionsSetup]");

  // Set time variables here (for now)
  Noutsnap = 0;
  nresync = 0;
  //tsnapnext = dt_snap;

  // Compute all initial N-body terms
  // --------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Zero all acceleration terms
    for (i=0; i<nbody->Nstar; i++) {
      for (k=0; k<ndim; k++) nbody->stardata[i].a[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].adot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = 0.0;
      for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = 0.0;
      nbody->stardata[i].gpot = 0.0;
      nbody->stardata[i].active = true;
      nbody->stardata[i].level = level_step;
      nbody->stardata[i].nstep = 0;
      nbody->stardata[i].nlast = 0;
      nbody->nbodydata[i] = &(nbody->stardata[i]);
    }

    nbody->Nnbody = nbody->Nstar;
    nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

  }

  // Set particle values for initial step (e.g. r0, v0, a0)
  nbody->EndTimestep(n,nbody->Nstar,nbody->nbodydata);

  this->CalculateDiagnostics();
  this->diag0 = this->diag;

  this->setup = true;

  return;
}



//=============================================================================
//  NbodySimulation::MainLoop
/// Main SPH simulation integration loop.
//=============================================================================
template <int ndim>
void NbodySimulation<ndim>::MainLoop(void)
{
  int i;                            // Particle loop counter
  int it;                           // Time-symmetric iteration counter
  int k;                            // Dimension counter

  debug2("[NbodySimulation::MainLoop]");

  // If we are using sub-systems, create N-body tree here
  // --------------------------------------------------------------------------
  if (nbody->sub_systems == 1) {
    nbody->reset_tree = 1;

    // If we are obliged to re-build the tree, then recompute the grav.
    // potential for all star particles (could be optimised in the future).
    if (nbody->reset_tree == 1) {
      
      // Zero all acceleration terms
      for (i=0; i<nbody->Nstar; i++) {
        for (k=0; k<ndim; k++) nbody->stardata[i].a[k] = 0.0;
        for (k=0; k<ndim; k++) nbody->stardata[i].adot[k] = 0.0;
        for (k=0; k<ndim; k++) nbody->stardata[i].a2dot[k] = 0.0;
        for (k=0; k<ndim; k++) nbody->stardata[i].a3dot[k] = 0.0;
        nbody->stardata[i].gpot = 0.0;
        nbody->stardata[i].active = true;
        nbody->nbodydata[i] = &(nbody->stardata[i]);
      }
      nbody->Nnbody = nbody->Nstar;
     
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
      nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

      nbodytree.CreateNbodySystemTree(nbody);
      nbodytree.BuildSubSystems(nbody);
    }

    nbodytree.FindPerturberLists(nbody);

  }
  // --------------------------------------------------------------------------


  // Compute timesteps for all particles
  if (Nlevels == 1)
    this->ComputeGlobalTimestep();
  else 
    this->ComputeBlockTimesteps();

  // Advance time variables
  n = n + 1;
  Nsteps = Nsteps + 1;
  t = t + timestep;

  // Advance SPH particles positions and velocities
  nbody->AdvanceParticles(n,nbody->Nnbody,nbody->nbodydata,timestep);

  // Compute N-body forces
  // --------------------------------------------------------------------------
  if (nbody->Nstar > 0) {

    // Iterate end-of-step
    // ------------------------------------------------------------------------
    for (it=0; it<nbody->Npec; it++) {

      // Zero all acceleration terms
      for (i=0; i<nbody->Nnbody; i++) {
        if (nbody->nbodydata[i]->active) {
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->adot[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a2dot[k] = 0.0;
          for (k=0; k<ndim; k++) nbody->nbodydata[i]->a3dot[k] = 0.0;
          nbody->nbodydata[i]->gpot = 0.0;
        }
      }

      // Calculate forces, force derivatives etc.., for active stars/systems
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
      
      // Calculate correction step for all stars at end of step
      nbody->CorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);

    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------


  // Now loop over children and, if they are systems, integrate
  // their internal motion
  // --------------------------------------------------------------------------
  for (i=0; i<nbody->Nnbody; i++) {
    if (nbody->nbodydata[i]->Ncomp > 1) {
      // The cast is needed because the function is defined only in
      // SystemParticle, not in NbodyParticle.  
      // The safety of the cast relies on the correctness of the Ncomp value
      subsystem->IntegrateInternalMotion(static_cast<SystemParticle<ndim>* > 
					 (nbody->nbodydata[i]), timestep);
    }
  }


  // Set all end-of-step variables
  nbody->EndTimestep(n,nbody->Nnbody,nbody->nbodydata);

  return;
}



//=============================================================================
//  NbodySimulation::ComputeGlobalTimestep
/// Computes global timestep for SPH simulation.  Calculates the minimum 
/// timestep for all SPH and N-body particles in the simulation.
//=============================================================================
template <int ndim>
void NbodySimulation<ndim>::ComputeGlobalTimestep(void)
{
  int i;                            // Particle counter
  DOUBLE dt = big_number_dp;        // Particle timestep
  DOUBLE dt_min = big_number_dp;    // Local copy of minimum timestep

  debug2("[NbodySimulation::ComputeGlobalTimestep]");

  // --------------------------------------------------------------------------
  if (n == nresync) {

    n = 0;
    level_max = 0;
    level_step = level_max + integration_step - 1;
    nresync = integration_step;

    // Compute minimum timestep due to stars/systems
    for (i=0; i<nbody->Nnbody; i++)
      dt_min = min(dt_min,nbody->Timestep(nbody->nbodydata[i]));

    // Set all particles to same timestep
    timestep = dt_min;
    for (i=0; i<nbody->Nnbody; i++) {
      nbody->nbodydata[i]->level = 0;
      nbody->nbodydata[i]->nstep = 
        pow(2,level_step - nbody->nbodydata[i]->level);
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->dt = timestep;
    }

  }
  // --------------------------------------------------------------------------

  return;
}



//=============================================================================
//  NbodySimulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=============================================================================
template <int ndim>
void NbodySimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                // Particle counter
  int istep;                            // ??
  int level;                            // Particle timestep level
  int last_level;                       // Previous timestep level
  int level_max_old;                    // Old level_max
  int level_max_nbody = 0;              // level_max for SPH particles only
  int nstep;                            // ??
  int nfactor;                          // ??
  DOUBLE dt;                            // Aux. timestep variable
  DOUBLE dt_min_nbody = big_number_dp;  // Maximum N-body particle timestep

  int *ninlevel;

  debug2("[NbodySimulation::ComputeBlockTimesteps]");

  // Synchronise all timesteps and reconstruct block timestep structure.
  // ==========================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;
    for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt = big_number_dp;

    // Compute minimum timestep due to stars/systems
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->Timestep(nbody->nbodydata[i]);
      if (dt < timestep) timestep = dt;
      if (dt < dt_min_nbody) dt_min_nbody = dt;
      nbody->nbodydata[i]->dt = dt;
    }

    // Calculate new block timestep levels
    level_max = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max = timestep*powf(2.0,level_max);
    
    // Calculate the maximum level occupied by all SPH particles
    level_max_nbody = min((int) (invlogetwo*log(dt_max/dt_min_nbody)) + 1, 
			  level_max);

    // Populate timestep levels with N-body particles
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->Timestep(nbody->nbodydata[i]);
      level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
      level = max(level,0);
      nbody->nbodydata[i]->level = level;
      nbody->nbodydata[i]->nlast = n;
      nbody->nbodydata[i]->nstep = 
	pow(2,level_step - nbody->nbodydata[i]->level);
    }

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

  }

  // If not resynchronising, check if any SPH/N-body particles need to move  
  // up or down timestep levels.
  // ==========================================================================
  else {

    level_max_old = level_max;
    level_max = 0;


    // Find all N-body particles at the beginning of a new timestep
    // ------------------------------------------------------------------------
    for (i=0; i<nbody->Nnbody; i++) {

      // Skip particles that are not at end of step
      if (nbody->nbodydata[i]->nlast == n) {
	nstep = nbody->nbodydata[i]->nstep;
	last_level = nbody->nbodydata[i]->level;

	// Compute new timestep value and level number
	dt = nbody->Timestep(nbody->nbodydata[i]);
	nbody->nbodydata[i]->dt = dt;
	level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);

	// Move up one level (if levels are correctly synchronised) or 
	// down several levels if required
	if (level < last_level && last_level > 1 && n%(2*nstep) == 0) 
	  nbody->nbodydata[i]->level = last_level - 1;
	else if (level > last_level)
	  nbody->nbodydata[i]->level = level;
	else
	  nbody->nbodydata[i]->level = last_level;

	nbody->nbodydata[i]->nlast = n;
	nbody->nbodydata[i]->nstep = 
	  pow(2,level_step - nbody->nbodydata[i]->level);

	/*if (last_level != nbody->nbodydata[i]->level) {
	  cout << "Changing level : " << i << "    " << last_level << "    "
	       << level << "     " << nbody->nbodydata[i]->level << "     "
	       << invlogetwo*log(dt_max/dt) << "    " << dt << "    " 
	       << timestep*(FLOAT) nbody->nbodydata[i]->nstep << "    " 
	       << dt_max << endl;
	       }*/

      }

      // Find maximum level of all N-body particles
      level_max_nbody = max(level_max_nbody,nbody->nbodydata[i]->level);
      level_max = max(level_max,nbody->nbodydata[i]->level);
    }
    // ------------------------------------------------------------------------
      

    // For now, don't allow levels to be removed
    level_max = max(level_max,level_max_old);
    level_step = level_max + integration_step - 1;

    for (i=0; i<nbody->Nnbody; i++) {
      if (nbody->nbodydata[i]->nlast == n) nbody->nbodydata[i]->nstep = 
	pow(2,level_step - nbody->nbodydata[i]->level);
    }

    // Update all timestep variables if we have removed or added any levels
    // ------------------------------------------------------------------------
    if (level_max != level_max_old) {

      // Increase maximum timestep level if correctly synchronised
      istep = pow(2,level_step - level_max_old + 1);
      if (level_max <= level_max_old - 1 && level_max_old > 1 && n%istep == 0)
	level_max = level_max_old - 1;
      else if (level_max == level_max_old)
	level_max = level_max_old;

      // Adjust integer time if levels added or removed
      if (level_max > level_max_old) {
	nfactor = pow(2,level_max - level_max_old);
	n *= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
      }
      else if (level_max < level_max_old) {
	nfactor = pow(2,level_max_old - level_max);
	n /= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast /= nfactor;
	for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep /= nfactor;
      }

    }
    // ------------------------------------------------------------------------

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

  }
  // ==========================================================================


#if defined(VERIFY_ALL)
  //VerifyBlockTimesteps();
#endif

  return;

  // Some validations
  // --------------------------------------------------------------------------
  ninlevel = new int[level_max+1];
  for (int l=0; l<level_max+1; l++) ninlevel[l] = 0;
  for (i=0; i<nbody->Nnbody; i++) ninlevel[nbody->nbodydata[i]->level]++;
  cout << "Level occupancy; timestep : " << timestep << endl;
  for (int l=0; l<level_max+1; l++) 
    cout << "level : " << l << "     N : " << ninlevel[l] << endl;
  delete[] ninlevel;
  cin >> i;
  return;
}
