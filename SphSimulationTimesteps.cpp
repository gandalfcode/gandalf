// ============================================================================
// SphSimulationTimesteps.cpp
// Contains all main functions controlling the SPH simulation work-flow.
// ============================================================================


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstring>
#include "Precision.h"
#include "Exception.h"
#include "SphSimulation.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;



// ============================================================================
// SphSimulation::ComputeGlobalTimestep
// Computes global timestep for SPH simulation.
// ============================================================================
void SphSimulation::ComputeGlobalTimestep(void)
{
  int i;                                    // Particle counter
  DOUBLE dt;                                // Aux. timestep variable

  debug2("[SphSimulation::ComputeGlobalTimestep]");

  // --------------------------------------------------------------------------
  if (n == nresync) {

    n = 0;
    level_max = 0;
    level_step = level_max + integration_step - 1;
    nresync = integration_step;
    timestep = big_number_dp;

    // Find minimum timestep from all SPH particles
    for (i=0; i<sph->Nsph; i++) {
      dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
      if (dt < timestep) timestep = dt;
    }
    
    // If integrating energy equation, include energy timestep
    if (simparams.stringparams["gas_eos"] == "energy_eqn") {
      for (i=0; i<sph->Nsph; i++) {
	dt = uint->Timestep(sph->sphdata[i]);
	if (dt < timestep) timestep = dt;
      }
    }
    
    // Set all particles to same timestep
    for (i=0; i<sph->Nsph; i++) sph->sphdata[i].level = 0;

  }
  // --------------------------------------------------------------------------

  return;
}



// ============================================================================
// SphSimulation::ComputeBlockTimesteps
// ..
// ============================================================================
void SphSimulation::ComputeBlockTimesteps(void)
{
  int i;                                    // ..
  int istep;                                // ..
  int level;                                // ..
  int last_level;                           // ..
  int level_diff;                           // ..
  int level_old;                            // ..
  int level_max_sph = 0;                    // ..
  DOUBLE dt;                                // Aux. timestep variable
  DOUBLE dt_max_sph = 0.0;                  // ..

  debug2("[SphSimulation::ComputeBlockTimesteps]");

  timestep = big_number;

  // Synchronise all timesteps
  // ==========================================================================
  if (n == nresync) {

    n = 0;

    // Find minimum timestep from all SPH particles
    for (i=0; i<sph->Nsph; i++) {
      dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
      if (dt < timestep) timestep = dt;
      if (dt > dt_max_sph) dt_max_sph = dt;
      sph->sphdata[i].dt = dt;
    }
    
    // If integrating energy equation, include energy timestep
    if (sph->gas_eos == "energy_eqn") {
      for (int i=0; i<sph->Nsph; i++) {
	dt = uint->Timestep(sph->sphdata[i]);
	if (dt < timestep) timestep = dt;
	sph->sphdata[i].dt = min(sph->sphdata[i].dt,dt);
      }
    }

    // Calculate new block timestep levels
    level_max = Nlevels - 1;
    level_step = level_max + integration_step - 1;
    dt_max = timestep*powf(2.0,level_max);
    
    // Calculate the maximum level occupied by all SPH particles
    level_max_sph = max((int) (invlogetwo*log(dt_max/dt_max_sph)) + 1, 0);

    if (sph_single_timestep == 1)
      for (i=0; i<sph->Nsph; i++) sph->sphdata[i].level = level_max_sph;
    else {
      dt = sph->sphdata[i].dt;
      level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);
      level = max(level,0);
      sph->sphdata[i].level = level;
    }

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

  }

  // If not resynchronising, check if any SPH particles need to move up 
  // or down timestep levels
  // ==========================================================================
  else {

    level_old = level_max;
    level_max = 0;

    // Find all SPH particles at the beginning of a new timestep
    // ------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {
      last_level = sph->sphdata[i].level;
      istep = pow(2,level_step - last_level + 1);

      // Skip particles not at end of step
      if (n%istep == 0) {
	last_level = sph->sphdata[i].level;
	dt = sphint->Timestep(sph->sphdata[i],sph->hydro_forces);
	if (sph->gas_eos == "energy_eqn") 
	  dt = min(dt,uint->Timestep(sph->sphdata[i]));
	level = min((int) (invlogetwo*log(dt_max/dt)) + 1, level_max);

	// Move up one level (if levels are correctly synchronised) or 
	// down several levels if required
	if (level < last_level && last_level > 1) sph->sphdata[i].level--;
	else if (level > last_level) sph->sphdata[i].level = level;
      }

      // Find maximum level of all SPH particles
      level_max_sph = max(level_max_sph,sph->sphdata[i].level);
      level_max = max(level_max,sph->sphdata[i].level);
    }
    // ------------------------------------------------------------------------
      

    // Set fixed SPH timestep level here in case maximum has changed
    if (sph_single_timestep == 1)
      for (i=0; i<sph->Nsph; i++) sph->sphdata[i].level = level_max_sph;


    // Update all timestep variables if we have removed or added any levels
    // ------------------------------------------------------------------------
    if (level_max != level_old) {

      // Increase maximum timestep level if correctly synchronised
      istep = pow(2,level_step - level_old + 1);
      if (level_max <= level_old - 1 && level_old > 1 && n%istep == 0)
	level_max = level_old - 1;
      else if (level_max = level_old)
	level_max = level_old;
      level_step = level_max + integration_step - 1;
      level_diff = level_max - level_old;

      // Adjust integer time if levels added or removed
      if (level_max > level_old)
	n *= pow(2,level_max - level_old);
      else if (level_max < level_old)
	n /= pow(2,level_old - level_max);

    }
    // ------------------------------------------------------------------------

    nresync = pow(2,level_step);
    timestep = dt_max / (DOUBLE) nresync;

  }
  // ==========================================================================

#if defined(VERIFY_ALL)
  VerifyBlockTimesteps();
#endif

  return;
}



// ============================================================================
// SphSimulation::IntegerTimestep
// ..
// ============================================================================
int SphSimulation::IntegerTimestep(int level)
{
  return pow(2,level_step - level);
}



// ============================================================================
// SphSimulation::VerifyBlockTimesteps
// ..
// ============================================================================
void SphSimulation::VerifyBlockTimesteps(void)
{
  debug2("[SphSimulation::VerifyBlockTimesteps]");

  // Check integer timestep variables are valid
  if (n < 0 || n > nresync) {
    cout << "Invalid integer timestep value : " 
	 << n << "   " << nresync << endl;
    exit(0);
  }

  // Check all particles occupy valid timestep levels
  for (int i=0; i<sph->Nsph; i++) {
    if (sph->sphdata[i].level > level_max && sph->sphdata[i].level < 0) {
      cout << "Invalid SPH timestep level : " << i << "   "
	   << sph->sphdata[i].level << "   " << level_max << endl;
      exit(0);
    }
  }

  return;
}
