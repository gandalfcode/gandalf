//=================================================================================================
//  SimAnalysis.cpp
//  Contains various analysis routines for Simulation object.
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
#include <iomanip>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstring>
#include "Exception.h"
#include "Simulation.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Debug.h"
#ifdef MPI_PARALLEL
#include <stddef.h>
#include "mpi.h"
#endif
using namespace std;



//=================================================================================================
//  Simulation::CalculateDiagnostics
/// Calculates all diagnostic quantities (e.g. conserved quantities),
/// saves to the diagnostic data structure and outputs to screen.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::CalculateDiagnostics(void)
{
  int i;                            // Particle counter
  int k;                            // Dimensionality counter

  debug2("[Simulation::CalculateDiagnostics]");

  diag.Nhydro = hydro->Nhydro;
  diag.Nstar  = nbody->Nstar;
  diag.Ndead  = 0;

  // Zero all diagnostic summation variables
  diag.mtot   = 0.0;
  diag.Etot   = 0.0;
  diag.utot   = 0.0;
  diag.ketot  = 0.0;
  diag.gpetot = 0.0;
  for (k=0; k<ndim; k++) diag.rcom[k]        = 0.0;
  for (k=0; k<ndim; k++) diag.vcom[k]        = 0.0;
  for (k=0; k<ndim; k++) diag.mom[k]         = 0.0;
  for (k=0; k<ndim; k++) diag.force[k]       = 0.0;
  for (k=0; k<3; k++) diag.angmom[k] = 0.0;

  // Loop over all hydro particles and add contributions to all quantities
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    if (part.flags.is_dead()) {
      diag.Ndead++;
      continue;
    }
    diag.mtot   += part.m;
    diag.ketot  += part.m*DotProduct(part.v,part.v,ndim);
    diag.utot   += part.m*part.u;
    diag.gpetot -= part.m*part.gpot;
    for (k=0; k<ndim; k++) {
      diag.rcom[k]        += part.m*part.r[k];
      diag.vcom[k]        += part.m*part.v[k];
      diag.mom[k]         += part.m*part.v[k];
      diag.force[k]       += part.m*part.a[k];
    }
  }

  // Add contributions to angular momentum depending on dimensionality
  if (ndim == 2) {
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      if (part.flags.is_dead()) continue;
      diag.angmom[2] += part.m*(part.r[0]*part.v[1] - part.r[1]*part.v[0]);
    }
  }
  else if (ndim == 3) {
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      if (part.flags.is_dead()) continue;
      diag.angmom[0] += part.m*(part.r[1]*part.v[2] - part.r[2]*part.v[1]);
      diag.angmom[1] += part.m*(part.r[2]*part.v[0] - part.r[0]*part.v[2]);
      diag.angmom[2] += part.m*(part.r[0]*part.v[1] - part.r[1]*part.v[0]);
    }
  }

  // Loop over all star particles and add contributions to all quantities
#if defined MPI_PARALLEL
   Box<ndim> mydomain =  mpicontrol->MyDomain();
#endif
  for (i=0; i<nbody->Nstar; i++) {
#if defined MPI_PARALLEL
	if (!ParticleInBox(nbody->stardata[i], mydomain)  ) continue;
#endif
    diag.mtot   += nbody->stardata[i].m;
    diag.ketot  += nbody->stardata[i].m*DotProduct(nbody->stardata[i].v,nbody->stardata[i].v,ndim);
    diag.gpetot -= nbody->stardata[i].m*nbody->stardata[i].gpot;
    for (k=0; k<ndim; k++) {
      diag.rcom[k]       += nbody->stardata[i].m*nbody->stardata[i].r[k];
      diag.vcom[k]       += nbody->stardata[i].m*nbody->stardata[i].v[k];
      diag.mom[k]        += nbody->stardata[i].m*nbody->stardata[i].v[k];
      diag.force[k]      += nbody->stardata[i].m*nbody->stardata[i].a[k];
    }

  // Add contributions to angular momentum depending on dimensionality
  if (ndim == 2) {
      diag.angmom[2] += nbody->stardata[i].m*(nbody->stardata[i].r[0]*nbody->stardata[i].v[1] -
                                              nbody->stardata[i].r[1]*nbody->stardata[i].v[0]);
  }
  else if (ndim == 3) {
      diag.angmom[0] += nbody->stardata[i].m*(nbody->stardata[i].r[1]*nbody->stardata[i].v[2] -
                                              nbody->stardata[i].r[2]*nbody->stardata[i].v[1]);
      diag.angmom[1] += nbody->stardata[i].m*(nbody->stardata[i].r[2]*nbody->stardata[i].v[0] -
                                              nbody->stardata[i].r[0]*nbody->stardata[i].v[2]);
      diag.angmom[2] += nbody->stardata[i].m*(nbody->stardata[i].r[0]*nbody->stardata[i].v[1] -
                                              nbody->stardata[i].r[1]*nbody->stardata[i].v[0]);
    }
  }

  // Add internal angular momentum (due to sink accretion) and subtract
  // accreted hydro momentum to maintain conservation of individual impulses.
  for (i=0; i<sinks->Nsink; i++) {
#if defined MPI_PARALLEL
	if (!ParticleInBox(nbody->stardata[sinks->sink[i].istar], mydomain )  ) continue;
#endif
    for (k=0; k<3; k++) diag.angmom[k]         += sinks->sink[i].angmom[k];
  }

  // Normalise all quantities and sum all contributions to total energy
  diag.ketot  *= 0.5;
  diag.gpetot *= 0.5;
  diag.Etot   = diag.ketot;
  if (diag.mtot > 0) {
	  for (k=0; k<ndim; k++) diag.rcom[k] /= diag.mtot;
	  for (k=0; k<ndim; k++) diag.vcom[k] /= diag.mtot;
  }
  if (hydro->hydro_forces == 1) diag.Etot += diag.utot;
  if (hydro->self_gravity == 1 || nbody->Nstar > 0) diag.Etot += diag.gpetot;


  // Calculate binary statistics if required
  if (simparams->intparams["binary_stats"] == 1) {
    nbodytree.CreateNbodySystemTree(nbody);
    nbodytree.FindBinarySystems(nbody);
  }


  // For MPI, collect all diagnostic information on root node
#ifdef MPI_PARALLEL
  mpicontrol->CollateDiagnosticsData(diag);
#endif

  RecordDiagnostics();

  return;
}



//=================================================================================================
//  Simulation::OutputDiagnostics
/// Output all diagnostic quantities that are computed in CalculateDiagnostics to screen.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::OutputDiagnostics(void)
{
  debug2("[Simulation::OutputDiagnostics]");

  if (rank != 0) return;

  cout << "Nhydro        : " << diag.Nhydro << endl;
  cout << "Nstar       : " << diag.Nstar << endl;
  if (diag.Ndead > 0)   cout << "Ndead       : " << diag.Ndead << endl;
  cout << "mtot        : " << diag.mtot*simunits.m.outscale << endl;
  cout << "Etot        : " << diag.Etot*simunits.E.outscale << endl;
  cout << "ketot       : " << diag.ketot*simunits.E.outscale << endl;
  if (hydro->hydro_forces == 1)  cout << "utot        : " << diag.utot*simunits.E.outscale << endl;
  if (hydro->self_gravity == 1 || nbody->Nstar > 0) {
    cout << "gpetot      : " << diag.gpetot*simunits.E.outscale << endl;
  }
  if (ndim == 1) {
    cout << "rcom        : " << diag.rcom[0] << endl;
    cout << "vcom        : " << diag.vcom[0] << endl;
    cout << "mom         : " << diag.mom[0] << endl;
    cout << "force       : " << diag.force[0] << endl;
  }
  else if (ndim == 2) {
    cout << "rcom        : " << diag.rcom[0] << "   " << diag.rcom[1] << endl;
    cout << "vcom        : " << diag.vcom[0] << "   " << diag.vcom[1] << endl;
    cout << "ang mom     : " << diag.angmom[2] << endl;
    cout << "mom         : " << diag.mom[0] << "   " << diag.mom[1] << endl;
    cout << "force       : " << diag.force[0] << "   " << diag.force[1] << endl;
  }
  else if (ndim == 3) {
    cout << "rcom        : " << diag.rcom[0] << "   "
         << diag.rcom[1] << "   " << diag.rcom[2] << endl;
    cout << "vcom        : " << diag.vcom[0] << "   "
         << diag.vcom[1] << "   " << diag.vcom[2] << endl;
    cout << "ang mom     : " << diag.angmom[0] << "   "
         << diag.angmom[1] << "   " << diag.angmom[2] << endl;
    cout << "mom         : " << diag.mom[0] << "   "
         << diag.mom[1] << "   " << diag.mom[2] << endl;
    cout << "force       : " << diag.force[0] << "   "
         << diag.force[1] << "   " << diag.force[2] << endl;
  }

  // Calculate binary statistics if required
  if (simparams->intparams["binary_stats"] == 1) nbodytree.OutputBinaryProperties(nbody);

  return;
}



//=================================================================================================
//  Simulation::RecordDiagnostics
/// Output all diagnostic quantities that are computed in CalculateDiagnostics to screen.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::RecordDiagnostics(void)
{
  debug2("[Simulation::RecordDiagnostics]");

  if (rank != 0) return;

  int k;
  ofstream outfile;
  string filename = run_id + ".diag";

  outfile.open(filename.c_str(),std::ofstream::app);

  outfile << t*simunits.t.outscale << "     ";
  outfile << Nsteps << "      ";
  outfile << timestep << "      ";
  outfile << dt_min_hydro << "      ";
  outfile << dt_min_nbody << "      ";
  outfile << level_max << "      ";
  outfile << diag.Nhydro << "     ";
  outfile << diag.Nstar << "     ";
  outfile << diag.Ndead << "     ";
  outfile << diag.mtot*simunits.m.outscale << "     ";
  outfile << diag.Etot*simunits.E.outscale << "     ";
  outfile << diag.ketot*simunits.E.outscale << "     ";
  outfile << diag.utot*simunits.E.outscale << "     ";
  outfile << diag.gpetot*simunits.E.outscale << "     ";
  for (k=0; k<3; k++) outfile << diag.angmom[k]*simunits.angmom.outscale << "     ";
  for (k=0; k<ndim; k++) outfile << diag.rcom[k]*simunits.r.outscale << "     ";
  for (k=0; k<ndim; k++) outfile << diag.vcom[k]*simunits.v.outscale << "     ";
  for (k=0; k<ndim; k++) outfile << diag.mom[k]*simunits.mom.outscale << "     ";
  for (k=0; k<ndim; k++) outfile << diag.force[k] << "     "; //*simunits.f.outscale << "     ";
  outfile << endl;

  outfile.close();

  // Now calculate and output any test specific diagnostics
  OutputTestDiagnostics();

  return;
}



//=================================================================================================
//  Simulation::OutputTestDiagnostics
/// Output all diagnostic quantities related to specific tests.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::OutputTestDiagnostics(void)
{
  int i;
  ofstream outfile;
  ofstream solfile;
  string ic = simparams->stringparams["ic"];

  debug2("[Simulation::OutputTestDiagnostics]");

  if (rank != 0) return;

  //-----------------------------------------------------------------------------------------------
  if (ic == "spitzer") {

    string filename  = run_id + ".ionfront";
    string solname   = run_id + ".spitzer";
    FLOAT temp_ion   = simparams->floatparams["temp_ion"]/simunits.temp.outscale;
    FLOAT mu_ion     = simparams->floatparams["mu_ion"];
    FLOAT cion       = sqrtf(temp_ion/mu_ion);
    FLOAT radius_ion = (FLOAT) 0.0;
    FLOAT m_ion      = (FLOAT) 0.0;
    FLOAT m_if       = (FLOAT) 0.0;
    FLOAT arecomb    = simparams->floatparams["arecomb"]*
      simunits.t.outscale*simunits.t.outcgs/pow(simunits.r.outscale*simunits.r.outcgs,3);
    FLOAT mcloud     = simparams->floatparams["mcloud"]/simunits.m.outscale;
    FLOAT radius     = simparams->floatparams["radius"]/simunits.r.outscale;
    FLOAT NLyC       = simparams->floatparams["NLyC"]*simunits.t.outscale*simunits.t.outSI;
    FLOAT rhoneutral = 3.0*mcloud/(4.0*pi*powf(radius,3.0));
    FLOAT m_h        = m_hydrogen/simunits.m.outscale/simunits.m.outSI;
    FLOAT Rstromgren = powf(0.75*m_h*m_h*NLyC/(pi*arecomb*rhoneutral*rhoneutral),onethird);

    cout << "RSTROMGREN : " << Rstromgren*simunits.r.outscale << endl;
    cout << "Sound speed of ionised gas : " << cion*simunits.v.outscale*simunits.v.outSI << endl;
    cout << "Pressure 1 : " << rhoneutral*temp_ion/mu_ion  << endl;
    cout << "Pressure 2 : " << cion*cion*rhoneutral << endl;


    // Compute location of ionisation front
    //---------------------------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      FLOAT temp = hydro->eos->Temperature(part);
      if (temp > 0.5*temp_ion) m_ion += part.m;
      if (temp > 0.05*temp_ion && temp < 0.95*temp_ion) {
        m_if += part.m;
        radius_ion += part.m*sqrt(DotProduct(part.r, part.r, ndim));
      }
    }
    //---------------------------------------------------------------------------------------------

    // Normalise ionisation front radial position
    radius_ion /= m_if;

    // Open new file if at beginning of simulation.  Otherwise append new data
    if (Nsteps == 0) {
      outfile.open(filename.c_str());
      solfile.open(solname.c_str());
    }
    else {
      outfile.open(filename.c_str(), std::ofstream::app);
      solfile.open(solname.c_str(), std::ofstream::app);
    }
    outfile << t*simunits.t.outscale << "    " << radius_ion*simunits.r.outscale << "    "
            << m_ion*simunits.m.outscale << endl;
    solfile << t*simunits.t.outscale << "   "
            << Rstromgren*powf(1.0 + 7.0*cion*t/(4.0*Rstromgren),(4.0/7.0))*simunits.r.outscale
            << "    " << m_h*m_h*NLyC*powf(1.0 + 7.0*cion*t/(4.0*Rstromgren),(6.0/7.0))*
                         simunits.m.outscale/(arecomb*rhoneutral)
            << endl;
    solfile.close();
    outfile.close();

  }
  //-----------------------------------------------------------------------------------------------

  return;
}
