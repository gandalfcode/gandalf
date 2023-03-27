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
#include "Ic.h"
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

  Ic<ndim> * icGenerator = NULL ;

  // Gnerate initial conditions either from external file
  //-----------------------------------------------------------------------------------------------
  if (ic == "file") {
    ReadSnapshotFile(simparams->stringparams["in_file"], simparams->stringparams["in_file_form"]);
    rescale_particle_data = true;
    this->initial_h_provided = false;
  }
  // Generate initial conditions via python scripts
  //-----------------------------------------------------------------------------------------------
  else if (ic == "python") {
    return;
  }
  // Generate initial conditions on-the-fly via various IC class functions
  //-----------------------------------------------------------------------------------------------
  else if (ic == "basic_sine") {
    icGenerator = new BasicIc<ndim>(this, invndim);
  }
  else if (ic == "binaryacc") {
    icGenerator = new BinaryAccretionIc<ndim>(this, invndim);
  }
  else if (ic == "blob") {
    icGenerator = new BlobIc<ndim>(this, invndim);
  }
  else if (ic == "bondi") {
    icGenerator = new BondiAccretionIc<ndim>(this, invndim);
  }
  else if (ic == "bb" || ic == "bossbodenheimer") {
    icGenerator = new BossBodenheimerIc<ndim>(this, invndim);
  }
  else if (ic == "cdiscontinuity") {
    icGenerator = new ContactDiscontinuityIc<ndim>(this, invndim);
  }
  else if (ic == "disc" && ndim > 1) {
    icGenerator = new DiscIc<ndim>(this, invndim);
  }
  else if (ic == "dustybox") {
    icGenerator = new DustyBoxIc<ndim>(this, invndim);
  }
  else if (ic == "evrard") {
    icGenerator = new EvrardCollapseIc<ndim>(this, invndim);
  }
  else if (ic == "ewaldsine" || ic == "ewaldsine2" ||
           ic == "ewaldslab" ||  ic == "ewaldcylinder" ||
           ic == "jeans") {
    icGenerator = new EwaldIc<ndim>(this, invndim);
  }
  else if (ic == "filament") {
    icGenerator = new FilamentIc<ndim>(this, invndim);
  }
  else if (ic == "gaussianring") {
    icGenerator = new GaussianRingIc<ndim>(this, invndim);
  }
  else if (ic == "gresho") {
    icGenerator = new GreshoVortexIc<ndim>(this, invndim);
  }
  else if (ic == "binary" || ic == "triple" || ic == "quadruple") {
    icGenerator = new HierarchicalSystemIc<ndim>(this, invndim);
  }
  else if (ic == "isothermsphere" || ic == "rotisothermsphere" || ic == "turbisothermsphere") {
    icGenerator = new IsothermalSphereIc<ndim>(this, invndim);
  }
  else if (ic == "khi") {
    icGenerator = new KelvinHelmholtzIc<ndim>(this, invndim);
  }
  else if (ic == "noh") {
    icGenerator = new NohIc<ndim>(this, invndim);
  }
  else if (ic == "plummer") {
    icGenerator = new PlummerSphereIc<ndim>(this, invndim);
  }
  else if (ic == "polytrope") {
    icGenerator = new PolytropeIc<ndim>(this, invndim);
  }
  else if (ic == "rti") {
    icGenerator = new RayleighTaylorIc<ndim>(this, invndim);
  }
  else if (ic == "sedov") {
    icGenerator = new SedovBlastwaveIc<ndim>(this, invndim);
  }
  else if (ic == "shocktube") {
    icGenerator = new ShocktubeIc<ndim>(this, invndim);
  }
  else if (ic == "shock2d") {
    icGenerator = new Shock2DIc<ndim>(this, invndim);
  }
  else if (ic == "shearflow") {
    icGenerator = new ShearflowIc<ndim>(this, invndim);
  }
  else if (ic == "silcc") {
    icGenerator = new SilccIc<ndim>(this, invndim);
  }
  else if (ic == "soundwave") {
    icGenerator = new SoundwaveIc<ndim>(this, invndim);
  }
  else if (ic == "spitzer") {
    icGenerator = new SpitzerExpansionIc<ndim>(this, invndim);
  }
  else if (ic == "turbcore") {
    icGenerator = new TurbulentCoreIc<ndim>(this, invndim);
  }
  else if (ic == "box" || ic == "sphere") {
    icGenerator = new UniformIc<ndim>(this, invndim);
  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Unrecognised parameter : ic = " + ic;
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------

  if (icGenerator != NULL) { // We are not reading an IC file

    // Finally, generate the particles for the chosen initial conditions
    icGenerator->Generate();

    // If selected, call all functions for regularising the particle distribution
    if (simparams->intparams["regularise_particle_ics"] == 1) {
      using Regularization::RegularizerFunction;
      using Regularization::ParticleRegularizer;

      // Regularise the particle distribution using the density function providing in the IC class
      RegularizerFunction<ndim> *reg_func = icGenerator->GetParticleRegularizer();
      ParticleRegularizer<ndim>(simparams, icBox)(hydro, neib, nbody, *reg_func);

      // Once regularisation step has finished, (re)set all particle properties
      icGenerator->SetParticleProperties();

      delete reg_func ;
    }
  }
  // Scale particle data to dimensionless code units if required
  if (rescale_particle_data) ConvertToCodeUnits();

  // Check that the initial conditions are valid
  if (icGenerator != NULL) {
    icGenerator->CheckInitialConditions();
    delete icGenerator ;
  }
  cout << "Finished creating initial conditions" << endl;

  return;
}
