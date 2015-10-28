//=================================================================================================
//  SimulationIC.hpp
//  Contains all routines for generating initial conditions on the fly.
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


#include <iostream>
#include <string>
#include <cstdio>
#include <cstring>
#include <math.h>
#include "Precision.h"
#include "Exception.h"
#include "Simulation.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "RandomNumber.h"
#include "Debug.h"
#include "Ghosts.h"
#include "IC.h"
#if defined(FFTW_TURBULENCE)
#include "fftw3.h"
#endif
using namespace std;



//=================================================================================================
//  Simulation::GenerateIC
/// Generate initial conditions for SPH simulation chosen in parameters file.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::GenerateIC(void)
{
  ifstream f;                                  // Stream of input file
  string in_file;                              // Restart snapshot filename
  string in_file_form;                         // Restart snapshot file format
  string filename;                             // Simulation '.restart' filename
  string ic = simparams->stringparams["ic"];   // Local copy of initial conditions string

  debug2("[Simulation::GenerateIC]");


  // First, check special case of restarting a simulation, in which case
  // determine the name of the last snapshot file to be re-read
  //-----------------------------------------------------------------------------------------------
  if (restart) {
    filename = run_id + ".restart";
    f.open(filename.c_str());

    // If file opens successfully, read snapshot and return
    if (!f.fail()) {
      f >> in_file_form;
      f >> in_file;
      f.close();
      ReadSnapshotFile(in_file,in_file_form);
      ConvertToCodeUnits();
      return;
    }
    // If unsuccessful (e.g. restart file doesn't exist), then make restart flag false
    else {
      restart = false;
    }
  }


  // If not a restart, generate initial conditions either from external file or created on the fly.
  //-----------------------------------------------------------------------------------------------
  Ic<ndim> icGenerator(this, hydro, invndim);

  if (ic == "file") {
    ReadSnapshotFile(simparams->stringparams["in_file"], simparams->stringparams["in_file_form"]);
    rescale_particle_data = true;
    this->initial_h_provided = false;
  }
  //-----------------------------------------------------------------------------------------------
  else if (ic == "bb") {
    icGenerator.BossBodenheimer();
  }
  else if (ic == "binary") {
    icGenerator.BinaryStar();
  }
  else if (ic == "binaryacc") {
    icGenerator.BinaryAccretion();
  }
  else if (ic == "blastwave") {
    icGenerator.BlastWave();
  }
  else if (ic == "bondi") {
    icGenerator.BondiAccretion();
  }
  else if (ic == "box") {
    icGenerator.UniformBox();
  }
  else if (ic == "cdiscontinuity") {
    icGenerator.ContactDiscontinuity();
  }
  else if (ic == "ewaldsine" || ic == "ewaldsine2" ||
           ic == "ewaldslab" ||  ic == "ewaldcylinder") {
    icGenerator.EwaldDensity();
  }
  else if (ic == "gresho") {
    icGenerator.GreshoVortex();
  }
  else if (ic == "khi") {
    icGenerator.KHI();
  }
  else if (ic == "noh") {
    icGenerator.NohProblem();
  }
  else if (ic == "plummer") {
    icGenerator.PlummerSphere();
  }
  else if (ic == "quadruple") {
    icGenerator.QuadrupleStar();
  }
  else if (ic == "rti") {
    icGenerator.RTI();
  }
  else if (ic == "sedov") {
    icGenerator.SedovBlastWave();
  }
  else if (ic == "shearflow") {
    icGenerator.ShearFlow();
  }
  else if (ic == "shocktube") {
    icGenerator.ShockTube();
  }
  else if (ic == "soundwave") {
    icGenerator.SoundWave();
  }
  else if (ic == "sphere") {
    icGenerator.UniformSphere();
  }
  else if (ic == "spitzer") {
    icGenerator.SpitzerExpansion();
  }
  else if (ic == "triple") {
    icGenerator.TripleStar();
  }
  else if (ic == "turbcore") {
    icGenerator.TurbulentCore();
  }
  else if (ic == "isothermsphere") {
    icGenerator.IsothermSphere();
  }
  else if (ic== "rotisothermsphere") {
    icGenerator.RotIsothermSphere();
  }
  else if (ic== "turbisothermsphere") {
    icGenerator.TurbIsothermSphere();
  }
  else if (ic == "evrard"){
	icGenerator.EvrardCollapse() ;
  }
  else if (ic == "python") {
    return;
  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Unrecognised parameter : ic = " + ic;
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------

  // Scale particle data to dimensionless code units if required
  if (rescale_particle_data) ConvertToCodeUnits();

  // Check that the initial conditions are valid
  icGenerator.CheckInitialConditions();

  cout << "Finsihed creating initial conditions" << endl;

  return;
}
