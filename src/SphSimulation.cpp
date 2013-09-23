//=============================================================================
//  SphSimulation.cpp
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
#include "RiemannSolver.h"
#include "Ghosts.h"
#include "Sinks.h"
using namespace std;



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class SphSimulation<1>;
template class SphSimulation<2>;
template class SphSimulation<3>;



//=============================================================================
//  SphSimulation::ProcessParameters
/// Process all the options chosen in the parameters file, setting various 
/// simulation variables and creating important simulation objects.
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::ProcessParameters(void)
{
  aviscenum avisc;                  // Artificial viscosity enum
  acondenum acond;                  // Artificial conductivity enum

  // Local references to parameter variables for brevity
  map<string, int> &intparams = simparams->intparams;
  map<string, float> &floatparams = simparams->floatparams;
  map<string, string> &stringparams = simparams->stringparams;

  debug2("[SphSimulation::ProcessParameters]");

  // Sanity check for valid dimensionality
  if (ndim < 1 || ndim > 3) {
    string message = "Invalid dimensionality chosen : ndim = " + ndim;
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the enum for artificial viscosity
  if (stringparams["avisc"] == "none")
    avisc = noneav;
  else if (stringparams["avisc"] == "mon97")
    avisc = mon97;
  else {
    string message = "Unrecognised parameter : avisc = " +
      simparams->stringparams["avisc"];
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the enum for artificial conductivity
  if (stringparams["acond"] == "none")
    acond = noneac;
  else if (stringparams["acond"] == "wadsley2008")
    acond = wadsley2008;
  else if (stringparams["acond"] == "price2008")
    acond = price2008;
  else {
    string message = "Unrecognised parameter : acond = " +
        simparams->stringparams["acond"];
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
  else if (stringparams["sph"] == "sm2012") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      sph = new SM2012Sph<ndim, TabulatedKernel> 
        (intparams["hydro_forces"], intparams["self_gravity"],
	 floatparams["alpha_visc"], floatparams["beta_visc"],
	 floatparams["h_fac"], floatparams["h_converge"],
	 avisc, acond, stringparams["gas_eos"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0){
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	sph = new SM2012Sph<ndim, M4Kernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "quintic") {
	sph = new SM2012Sph<ndim, QuinticKernel> 
	  (intparams["hydro_forces"], intparams["self_gravity"],
	   floatparams["alpha_visc"], floatparams["beta_visc"],
	   floatparams["h_fac"], floatparams["h_converge"],
	   avisc, acond, stringparams["gas_eos"], KernelName);
      }
      else if (KernelName == "gaussian") {
	sph = new SM2012Sph<ndim, GaussianKernel> 
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


  // Thermal physics object.  If energy equation is chosen, also initiate
  // the energy integration object.
  // --------------------------------------------------------------------------
  string gas_eos = stringparams["gas_eos"];
  if (gas_eos == "energy_eqn" || gas_eos == "constant_temp") {
    sph->eos = new Adiabatic<ndim>(floatparams["temp0"],
				   floatparams["mu_bar"],
				   floatparams["gamma_eos"]);
  }
  else if (gas_eos == "isothermal")
    sph->eos = new Isothermal<ndim>(floatparams["temp0"],
				    floatparams["mu_bar"],
				    floatparams["gamma_eos"],
                                    &simunits);
  else if (gas_eos == "barotropic")
    sph->eos = new Barotropic<ndim>(floatparams["temp0"],
				    floatparams["mu_bar"],
				    floatparams["gamma_eos"],
				    floatparams["rho_bary"],
                                    &simunits);
  else {
    string message = "Unrecognised parameter : gas_eos = " + gas_eos;
    ExceptionHandler::getIstance().raise(message);
  }



  // Create neighbour searching object based on chosen method in params file
  // --------------------------------------------------------------------------
  if (stringparams["neib_search"] == "bruteforce")
    sphneib = new BruteForceSearch<ndim>;
  else if (stringparams["neib_search"] == "grid")
    sphneib = new GridSearch<ndim>;
  else if (stringparams["neib_search"] == "tree") {
    sphneib = new BinaryTree<ndim>(intparams["Nleafmax"],
                                   floatparams["thetamaxsqd"],
                                   sph->kernp->kernrange,
                                   stringparams["gravity_mac"],
                                   stringparams["multipole"]);
  }
  else {
    string message = "Unrecognised parameter : neib_search = " 
      + simparams->stringparams["neib_search"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Create SPH particle integration object
  // --------------------------------------------------------------------------
  if (stringparams["sph_integration"] == "lfkdk") {
    sphint = new SphLeapfrogKDK<ndim>(floatparams["accel_mult"],
			              floatparams["courant_mult"]);
  }
  else if (stringparams["sph_integration"] == "lfdkd") {
    sphint = new SphLeapfrogDKD<ndim>(floatparams["accel_mult"],
			              floatparams["courant_mult"]);
    integration_step = max(integration_step,2);
  }
  else {
    string message = "Unrecognised parameter : sph_integration = " 
      + simparams->stringparams["sph_integration"];
    ExceptionHandler::getIstance().raise(message);
  }


  // Energy integration object
  // --------------------------------------------------------------------------
  if (stringparams["energy_integration"] == "PEC") {
    uint = new EnergyPEC<ndim>(floatparams["energy_mult"]);
  }
  else if (stringparams["energy_integration"] == "lfdkd") {
    uint = new EnergyLeapfrogDKD<ndim>(floatparams["energy_mult"]);
  }
  else {
    string message = "Unrecognised parameter : energy_integration = "
      + simparams->stringparams["energy_integration"];
    ExceptionHandler::getIstance().raise(message);
  }

  
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
  // Create N-body object based on chosen method in params file
  // --------------------------------------------------------------------------
  else if (stringparams["nbody"] == "lfdkd") {
    string KernelName = stringparams["kernel"];
    if (intparams["tabulated_kernel"] == 1) {
      nbody = new NbodyLeapfrogDKD<ndim, TabulatedKernel> 
	(intparams["nbody_softening"], intparams["sub_systems"],
	 floatparams["nbody_mult"], KernelName);
    }
    else if (intparams["tabulated_kernel"] == 0) {
      // Depending on the kernel, instantiate a different GradSph object
      if (KernelName == "m4") {
	nbody = new NbodyLeapfrogDKD<ndim, M4Kernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "quintic") {
	nbody = new NbodyLeapfrogDKD<ndim, QuinticKernel> 
	  (intparams["nbody_softening"], intparams["sub_systems"],
	   floatparams["nbody_mult"], KernelName);
      }
      else if (KernelName == "gaussian") {
	nbody = new NbodyLeapfrogDKD<ndim, GaussianKernel> 
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
    integration_step = max(integration_step,2);
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
	   floatparams["nbody_mult"], KernelName);
      }
      else if (intparams["tabulated_kernel"] == 0) {
	// Depending on the kernel, instantiate a different GradSph object
	if (KernelName == "m4") {
	  subsystem = new NbodyLeapfrogKDK<ndim, M4Kernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["nbody_mult"], KernelName);
	}
	else if (KernelName == "quintic") {
	  subsystem = new NbodyLeapfrogKDK<ndim, QuinticKernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["nbody_mult"], KernelName);
	}
	else if (KernelName == "gaussian") {
	  subsystem = new NbodyLeapfrogKDK<ndim, GaussianKernel> 
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
	   floatparams["nbody_mult"], KernelName);
      }
      else if (intparams["tabulated_kernel"] == 0) {
	// Depending on the kernel, instantiate a different GradSph object
	if (KernelName == "m4") {
	  subsystem = new NbodyHermite4<ndim, M4Kernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["nbody_mult"], KernelName);
	}
	else if (KernelName == "quintic") {
	  subsystem = new NbodyHermite4<ndim, QuinticKernel> 
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["nbody_mult"], KernelName);
	}
	else if (KernelName == "gaussian") {
	  subsystem = new NbodyHermite4<ndim, GaussianKernel> 
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
	   floatparams["nbody_mult"], KernelName, intparams["Npec"]);
      }
      else if (intparams["tabulated_kernel"] == 0) {
	// Depending on the kernel, instantiate a different GradSph object
	if (KernelName == "m4") {
	  subsystem = new NbodyHermite4TS<ndim, M4Kernel>
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["nbody_mult"], KernelName, intparams["Npec"]);
	}
	else if (KernelName == "quintic") {
	  subsystem = new NbodyHermite4TS<ndim, QuinticKernel>
	    (intparams["nbody_softening"], intparams["sub_systems"],
	     floatparams["nbody_mult"], KernelName, intparams["Npec"]);
	}
	else if (KernelName == "gaussian") {
	  subsystem = new NbodyHermite4TS<ndim, GaussianKernel>
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


  // Sink particles
  // --------------------------------------------------------------------------
  sink_particles = intparams["sink_particles"];
  sinks.sink_particles = intparams["sink_particles"];
  sinks.create_sinks = intparams["create_sinks"];
  sinks.smooth_accretion = intparams["smooth_accretion"];
  sinks.rho_sink = floatparams["rho_sink"]
    /simunits.rho.outscale/simunits.rho.outcgs;
  sinks.alpha_ss = floatparams["alpha_ss"];
  sinks.smooth_accrete_frac = floatparams["smooth_accrete_frac"];
  sinks.smooth_accrete_dt = floatparams["smooth_accrete_dt"];
  sinks.sink_radius_mode = stringparams["sink_radius_mode"];

  if (sinks.sink_radius_mode == "fixed")
    sinks.sink_radius = floatparams["sink_radius"]/simunits.r.outscale;
  else
    sinks.sink_radius = floatparams["sink_radius"];


  // Set all other parameter variables
  // --------------------------------------------------------------------------
  sph->Nsph             = intparams["Nsph"];
  sph->riemann_solver   = stringparams["riemann_solver"];
  sph->slope_limiter    = stringparams["slope_limiter"];
  sph->riemann_order    = intparams["riemann_order"];
  sph->create_sinks     = intparams["create_sinks"];

  nbody->Nstar          = intparams["Nstar"];
  nbody_single_timestep = intparams["nbody_single_timestep"];

  dt_python             = floatparams["dt_python"];
  dt_snap               = floatparams["dt_snap"]/simunits.t.outscale;
  level_diff_max        = intparams["level_diff_max"];
  Nlevels               = intparams["Nlevels"];
  noutputstep           = intparams["noutputstep"];
  Nstepsmax             = intparams["Nstepsmax"];
  out_file_form         = stringparams["out_file_form"];
  run_id                = stringparams["run_id"];
  sph_single_timestep   = intparams["sph_single_timestep"];
  tend                  = floatparams["tend"]/simunits.t.outscale;
  tsnapnext             = floatparams["tsnapfirst"]/simunits.t.outscale;

  // Flag that we've processed all parameters already
  ParametersProcessed = true;

  return;
}



//TODO: make this mess more modular (note: initial h computation
//should be done inside the neighbour search)
//=============================================================================
//  SphSimulation::PostInitialConditionsSetup
/// ..
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::PostInitialConditionsSetup(void)
{
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[SphSimulation::PostInitialConditionsSetup]");

  // Set time variables here (for now)
  Noutsnap = 0;
  nresync = 0;
  //tsnapnext = dt_snap;

  // Set initial smoothing lengths and create initial ghost particles
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Set all relevant particle counters
    sph->Nghost = 0;
    sph->Nghostmax = sph->Nsphmax - sph->Nsph;
    sph->Ntot = sph->Nsph;
    for (int i=0; i<sph->Nsph; i++) sph->sphdata[i].active = true;

    // Compute mean mass
    sph->mmean = 0.0;
    for (i=0; i<sph->Nsph; i++) sph->mmean += sph->sphdata[i].m;
    sph->mmean /= (FLOAT) sph->Nsph;

    sph->InitialSmoothingLengthGuess();
    sphneib->UpdateTree(sph,*simparams);

    sphneib->neibcheck = false;
    sphneib->UpdateAllSphProperties(sph,nbody);

    // Search ghost particles
    ghosts.SearchGhostParticles(simbox,sph);

    // Update neighbour tree
    sphneib->UpdateTree(sph,*simparams);

    level_step = 1;

    // Zero accelerations (perhaps here)
    for (i=0; i<sph->Ntot; i++) sph->sphdata[i].active = true;

    // Calculate all SPH properties
    sphneib->UpdateAllSphProperties(sph,nbody);

    // Search ghost particles
    ghosts.SearchGhostParticles(simbox,sph);

    // Update neighbour tre
    sphneib->UpdateTree(sph,*simparams);
    sphneib->neibcheck = true;
    sphneib->UpdateAllSphProperties(sph,nbody);

  }


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
      nbody->stardata[i].level = 0; //level_step;
      nbody->stardata[i].nstep = 0;
      nbody->stardata[i].nlast = 0;
      nbody->nbodydata[i] = &(nbody->stardata[i]);
    }

    nbody->Nnbody = nbody->Nstar;
    nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);
    if (sph->self_gravity == 1)
      nbody->CalculateDirectSPHForces(nbody->Nnbody,sph->Nsph,
				      sph->sphdata,nbody->nbodydata);
    nbody->CalculateAllStartupQuantities(nbody->Nnbody,nbody->nbodydata);

  }


  // Compute all initial SPH force terms
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // Zero accelerations (here for now)
    for (i=0; i<sph->Ntot; i++) {
      for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
      for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
      sph->sphdata[i].gpot = (FLOAT) 0.0;
      sph->sphdata[i].dudt = (FLOAT) 0.0;
      sph->sphdata[i].active = true;
      sph->sphdata[i].level = 0;
      sph->sphdata[i].nstep = 0;
      sph->sphdata[i].nlast = 0;
    }

    ghosts.CopySphDataToGhosts(sph);
    sphneib->UpdateTree(sph, *simparams);

    // Calculate SPH gravity and hydro forces, depending on which are activated
    if (sph->hydro_forces == 1 && sph->self_gravity == 1)
      sphneib->UpdateAllSphForces(sph);
    else if (sph->hydro_forces == 1)
      sphneib->UpdateAllSphHydroForces(sph);
    else if (sph->self_gravity == 1)
      sphneib->UpdateAllSphGravForces(sph);

    // Compute contribution to grav. accel from stars
    for (i=0; i<sph->Nsph; i++)
      if (sph->sphdata[i].active)
        sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,
                                   sph->sphdata[i]);

    // Add accelerations
    for (i=0; i<sph->Nsph; i++) {
      sph->sphdata[i].active = false;
      for (k=0; k<ndim; k++)
        sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
    }

    ghosts.CopySphDataToGhosts(sph);

  }

  // Set particle values for initial step (e.g. r0, v0, a0)
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EndTimestep(n,sph->Nsph,sph->sphdata);
  sphint->EndTimestep(n,sph->Nsph,sph->sphdata);
  nbody->EndTimestep(n,nbody->Nstar,nbody->nbodydata);

  this->CalculateDiagnostics();
  this->OutputDiagnostics();
  this->diag0 = this->diag;
  this->setup = true;

  return;
}



//=============================================================================
//  SphSimulation::MainLoop
/// Main SPH simulation integration loop.
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::MainLoop(void)
{
  int activecount;                  // ..
  int i;                            // Particle loop counter
  int it;                           // Time-symmetric iteration counter
  int k;                            // Dimension counter

  debug2("[SphSimulation::MainLoop]");


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
  sphint->AdvanceParticles(n,sph->Nsph,sph->sphdata,(FLOAT) timestep);
  if (simparams->stringparams["gas_eos"] == "energy_eqn")
    uint->EnergyIntegration(n,sph->Nsph,sph->sphdata,(FLOAT) timestep);
  nbody->AdvanceParticles(n,nbody->Nnbody,nbody->nbodydata,timestep);

  // Check all boundary conditions
  ghosts.CheckBoundaries(simbox,sph);


  // Compute all SPH quantities
  // --------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    
    // Reorder particles

    // Search ghost particles
    ghosts.SearchGhostParticles(simbox,sph);

    // Update neighbour tree
    sphneib->UpdateTree(sph,*simparams);

    // Iterate if we need to immediately change SPH particle timesteps
    // (e.g. due to feedback, or sudden change in neighbour timesteps)
    // ------------------------------------------------------------------------
    do {

      // Update cells containing active particles
      sphneib->UpdateActiveParticleCounters(sph);

      // Calculate all SPH properties
      sphneib->UpdateAllSphProperties(sph,nbody);
      
      // Copy properties from original particles to ghost particles
      ghosts.CopySphDataToGhosts(sph);
      
      // Zero accelerations
#pragma parallel for default(none) private(k) shared(sph)
      for (i=0; i<sph->Ntot; i++) {
        if (sph->sphdata[i].active) {
          for (k=0; k<ndim; k++) sph->sphdata[i].a[k] = (FLOAT) 0.0;
          for (k=0; k<ndim; k++) sph->sphdata[i].agrav[k] = (FLOAT) 0.0;
          sph->sphdata[i].gpot = (FLOAT) 0.0;
          sph->sphdata[i].gpe = (FLOAT) 0.0;
          sph->sphdata[i].dudt = (FLOAT) 0.0;
          sph->sphdata[i].levelneib = 0;
        }
      }
      
      // Calculate SPH gravity and hydro forces, depending on which are activated
      if (sph->hydro_forces == 1 && sph->self_gravity == 1)
        sphneib->UpdateAllSphForces(sph);
      else if (sph->hydro_forces == 1)
        sphneib->UpdateAllSphHydroForces(sph);
      else if (sph->self_gravity == 1)
        sphneib->UpdateAllSphGravForces(sph);
      
      // Compute contribution to grav. accel from stars
      for (i=0; i<sph->Nsph; i++)
        if (sph->sphdata[i].active)
          sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,
                                     sph->sphdata[i]);
      
      // Add accelerations
      for (i=0; i<sph->Nsph; i++) {
        if (sph->sphdata[i].active) {
          for (k=0; k<ndim; k++)
            sph->sphdata[i].a[k] += sph->sphdata[i].agrav[k];
        }
      }

      // Check if all neighbouring timesteps are acceptable
      if (Nlevels > 1)
        activecount = sphint->CheckTimesteps(level_diff_max,n,
                                             sph->Nsph,sph->sphdata);
      else activecount = 0;
      activecount = 0;

    } while (activecount > 0);
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------


  // Compute N-body forces
  // --------------------------------------------------------------------------
  if (nbody->Nnbody > 0) {

    // Iterate for P(EC)^n schemes
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
          nbody->nbodydata[i]->gpe = 0.0;
        }
      }

      // Calculate forces, force derivatives etc.., for active stars/systems
      nbody->CalculateDirectGravForces(nbody->Nnbody,nbody->nbodydata);

      if (sph->self_gravity == 1)
        nbody->CalculateDirectSPHForces(nbody->Nnbody,sph->Nsph,
                                        sph->sphdata,nbody->nbodydata);

      // Calculate correction step for all stars at end of step
      nbody->CorrectionTerms(n,nbody->Nnbody,nbody->nbodydata,timestep);

    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------


  // Compute correction steps for all SPH particles
  if (sph->Nsph > 0) {
    sphint->CorrectionTerms(n,sph->Nsph,sph->sphdata,(FLOAT) timestep);
    if (simparams->stringparams["gas_eos"] == "energy_eqn")
      uint->EnergyCorrectionTerms(n,sph->Nsph,sph->sphdata,(FLOAT) timestep);
  }


  // Search for new sink particles (if activated)
  if (sink_particles == 1) {
    if (sinks.create_sinks == 1) sinks.SearchForNewSinkParticles(n,sph,nbody);
    if (sinks.Nsink > 0) {
      for (int s=0; s<sinks.Nsink; s++) 
        sinks.AccreteMassToSinks(sph,nbody,n,timestep,s);
    }
  }


  // End-step terms for all SPH particles
  if (sph->Nsph > 0) {
    if (simparams->stringparams["gas_eos"] == "energy_eqn")
      uint->EndTimestep(n,sph->Nsph,sph->sphdata);
    sphint->EndTimestep(n,sph->Nsph,sph->sphdata);
  }

  // End-step terms for all star particles
  if (nbody->Nstar > 0)
    nbody->EndTimestep(n,nbody->Nnbody,nbody->nbodydata);

  return;
}



//=============================================================================
//  SphSimulation::ComputeGlobalTimestep
/// Computes global timestep for SPH simulation.  Calculates the minimum 
/// timestep for all SPH and N-body particles in the simulation.
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::ComputeGlobalTimestep(void)
{
  int i;                            // Particle counter
  DOUBLE dt;                        // Particle timestep
  DOUBLE dt_min = big_number_dp;    // Local copy of minimum timestep

  debug2("[SphSimulation::ComputeGlobalTimestep]");

  // --------------------------------------------------------------------------
  if (n == nresync) {

    n = 0;
    level_max = 0;
    level_step = level_max + integration_step - 1;
    nresync = integration_step;

    // Find minimum timestep from all SPH particles
    // ------------------------------------------------------------------------
#pragma omp parallel default(none) private(i,dt) \
  shared(dt_min)
    {
      dt = big_number_dp;
#pragma omp for
      for (i=0; i<sph->Nsph; i++) {
        sph->sphdata[i].dt = sphint->Timestep(sph->sphdata[i],
                                              sph->hydro_forces);
        dt = min(dt,sph->sphdata[i].dt);
      }
      
      // If integrating energy equation, include energy timestep
      if (simparams->stringparams["gas_eos"] == "energy_eqn") {
#pragma omp for
        for (i=0; i<sph->Nsph; i++) {
          sph->sphdata[i].dt = min(sph->sphdata[i].dt,
                                   uint->Timestep(sph->sphdata[i]));
          dt = min(dt,sph->sphdata[i].dt);
        }

      }

#pragma omp critical
      {
        if (dt < dt_min) dt_min = dt;
      }
    }
    // ------------------------------------------------------------------------


    // Now compute minimum timestep due to stars/systems
    for (i=0; i<nbody->Nnbody; i++)
      dt_min = min(dt_min,nbody->Timestep(nbody->nbodydata[i]));

    
    // Set all particles to same timestep
    timestep = dt_min;
    for (i=0; i<sph->Nsph; i++) {
      sph->sphdata[i].level = 0;
      sph->sphdata[i].levelneib = 0;
      sph->sphdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
      sph->sphdata[i].nlast = n;
      sph->sphdata[i].dt = timestep;
    }
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
//  SphSimulation::ComputeBlockTimesteps
/// Compute timesteps for all particles using hierarchical block timesteps.
//=============================================================================
template <int ndim>
void SphSimulation<ndim>::ComputeBlockTimesteps(void)
{
  int i;                                // Particle counter
  int istep;                            // ??
  int level;                            // Particle timestep level
  int last_level;                       // Previous timestep level
  int level_max_old;                    // Old level_max
  int level_max_sph = 0;                // level_max for SPH particles only
  int level_min_sph = 9999999;          // level_min for SPH particles
  int level_max_nbody = 0;              // level_max for star particles only
  int nstep;                            // ??
  int nfactor;                          // ??
  DOUBLE dt;                            // Aux. timestep variable
  DOUBLE dt_min_sph = big_number_dp;    // Maximum SPH particle timestep
  DOUBLE dt_min_nbody = big_number_dp;  // Maximum N-body particle timestep

  int *ninlevel;

  debug2("[SphSimulation::ComputeBlockTimesteps]");

  // Synchronise all timesteps and reconstruct block timestep structure.
  // ==========================================================================
  if (n == nresync) {

    n = 0;
    timestep = big_number_dp;
    for (i=0; i<sph->Nsph; i++) sph->sphdata[i].dt = big_number_dp;
    for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->dt = big_number_dp;

    // If integrating energy equation, calculate the explicit energy timestep
    if (sph->gas_eos == "energy_eqn") {
      for (i=0; i<sph->Nsph; i++)
	sph->sphdata[i].dt = uint->Timestep(sph->sphdata[i]);
    }

    // Find minimum timestep from all SPH particles
    for (i=0; i<sph->Nsph; i++) {
      dt = min(sph->sphdata[i].dt,
	       sphint->Timestep(sph->sphdata[i],sph->hydro_forces));
      if (dt < timestep) timestep = dt;
      if (dt < dt_min_sph) dt_min_sph = dt;
      sph->sphdata[i].dt = dt;
    }
    
    // Now compute minimum timestep due to stars/systems
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
    level_max_sph = min((int) (invlogetwo*log(dt_max/dt_min_sph)) + 1, 
			level_max);
    level_max_nbody = min((int) (invlogetwo*log(dt_max/dt_min_nbody)) + 1, 
			  level_max);

    // If enforcing a single SPH timestep, set it here.  Otherwise, populate 
    // the timestep levels with SPH particles.
    if (sph_single_timestep == 1) 
      for (i=0; i<sph->Nsph; i++) {
        sph->sphdata[i].level = level_max_sph;
        sph->sphdata[i].levelneib = level_max_sph;
        sph->sphdata[i].nlast = n;
        sph->sphdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
	level_min_sph = min(level_min_sph,sph->sphdata[i].level);
      }
    else {
      for (i=0; i<sph->Nsph; i++) {
        dt = sph->sphdata[i].dt;
        level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
        level = max(level,0);
        sph->sphdata[i].level = level;
        sph->sphdata[i].levelneib = level;
        sph->sphdata[i].nlast = n;
        sph->sphdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
	level_min_sph = min(level_min_sph,sph->sphdata[i].level);
      }
    }

    // Populate timestep levels with N-body particles
    for (i=0; i<nbody->Nnbody; i++) {
      dt = nbody->nbodydata[i]->dt;
      level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
      level = max(level,0);
      nbody->nbodydata[i]->level = max(level,level_max_sph);
      //nbody->nbodydata[i]->level = max(level,level_min_sph);
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

    // Find all SPH particles at the beginning of a new timestep
    // ------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {

      // Skip particles that are not at end of step
      if (sph->sphdata[i].nlast == n) {
        nstep = sph->sphdata[i].nstep;
        last_level = sph->sphdata[i].level;
	
        // Compute new timestep value and level number
        dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
        if (sph->gas_eos == "energy_eqn")
          dt = min(dt,uint->Timestep(sph->sphdata[i]));
        sph->sphdata[i].dt = dt;
        level = max((int) (invlogetwo*log(dt_max/dt)) + 1, 0);
        level = max(level,sph->sphdata[i].levelneib - level_diff_max);

        // Move up one level (if levels are correctly synchronised) or
        // down several levels if required
        if (level < last_level && last_level > 1 && n%(2*nstep) == 0)
          sph->sphdata[i].level = last_level - 1;
        else if (level > last_level)
          sph->sphdata[i].level = level;
        else
          sph->sphdata[i].level = last_level;

        sph->sphdata[i].nlast = n;
        sph->sphdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
      }

      // Find maximum level of all SPH particles
      level_max_sph = max(level_max_sph,sph->sphdata[i].level);
      level_min_sph = min(level_min_sph,sph->sphdata[i].level);
      level_max = max(level_max,sph->sphdata[i].level);
    }
    // ------------------------------------------------------------------------
      

    // Now find all N-body particles at the beginning of a new timestep
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
        level = max(level,level_max_sph);
        //level = max(level,level_min_sph);

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
      }

      // Find maximum level of all N-body particles
      level_max_nbody = max(level_max_nbody,nbody->nbodydata[i]->level);
      level_max = max(level_max,nbody->nbodydata[i]->level);
    }
    // ------------------------------------------------------------------------
      

    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1) {
      for (i=0; i<sph->Nsph; i++) {
        if (sph->sphdata[i].nlast == n)  {
          sph->sphdata[i].level = level_max_sph;
          sph->sphdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
        }
      }
    }

    // For now, don't allow levels to be removed
    level_max = max(level_max,level_max_old);
    level_step = level_max + integration_step - 1;

    for (i=0; i<sph->Nsph; i++) {
      if (sph->sphdata[i].nlast == n)
        sph->sphdata[i].nstep = pow(2,level_step - sph->sphdata[i].level);
    }
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
        for (i=0; i<sph->Nsph; i++) sph->sphdata[i].nlast *= nfactor;
        for (i=0; i<sph->Nsph; i++) sph->sphdata[i].nstep *= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nlast *= nfactor;
        for (i=0; i<nbody->Nnbody; i++) nbody->nbodydata[i]->nstep *= nfactor;
      }
      else if (level_max < level_max_old) {
        nfactor = pow(2,level_max_old - level_max);
        n /= nfactor;
        for (i=0; i<sph->Nsph; i++) sph->sphdata[i].nlast /= nfactor;
        for (i=0; i<sph->Nsph; i++) sph->sphdata[i].nstep /= nfactor;
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

  cout << "Checking timesteps : " << level_max << "   " << level_max_sph << "    " << level_max_nbody << endl;

  for (int l=0; l<level_max+1; l++) ninlevel[l] = 0;
  for (i=0; i<sph->Nsph; i++) ninlevel[sph->sphdata[i].level]++;
  cout << "SPH level occupancy" << endl;
  for (int l=0; l<level_max+1; l++) 
    cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  for (int l=0; l<level_max+1; l++) ninlevel[l] = 0;
  for (i=0; i<nbody->Nstar; i++) ninlevel[nbody->nbodydata[i]->level]++;
  cout << "N-body level occupancy" << endl;
  for (int l=0; l<level_max+1; l++)
    cout << "level : " << l << "     N : " << ninlevel[l] << endl;

  for (int l=0; l<level_max+1; l++) {
    if (ninlevel[l] > 0 && l < level_max_sph) {
      cout << "Something going wrong with timesteps" << endl;
      exit(0);
    }
  }

  delete[] ninlevel;

  return;
}




