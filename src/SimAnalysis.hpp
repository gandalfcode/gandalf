//=============================================================================
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
//=============================================================================


#include <iostream>
#include <iomanip>
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



//=============================================================================
//  Simulation::CalculateDiagnostics
/// Calculates all diagnostic quantities (e.g. conserved quantities),
/// saves to the diagnostic data structure and outputs to screen.
//=============================================================================
template <int ndim>
void Simulation<ndim>::CalculateDiagnostics(void)
{
  int i;                            // Particle counter
  int k;                            // Dimensionality counter

  debug2("[SphSimulation::CalculateDiagnostics]");

  diag.Nsph = sph->Nsph;
  diag.Nstar = nbody->Nstar;
  diag.Ndead = 0;

  // Zero all diagnostic summation variables
  diag.mtot = 0.0;
  diag.Etot = 0.0;
  diag.utot = 0.0;
  diag.ketot = 0.0;
  diag.gpetot = 0.0;
  for (k=0; k<ndim; k++) diag.rcom[k] = 0.0;
  for (k=0; k<ndim; k++) diag.vcom[k] = 0.0;
  for (k=0; k<ndim; k++) diag.mom[k] = 0.0;
  for (k=0; k<ndim; k++) diag.force[k] = 0.0;
  for (k=0; k<ndim; k++) diag.force_hydro[k] = 0.0;
  for (k=0; k<ndim; k++) diag.force_grav[k] = 0.0;
  for (k=0; k<3; k++) diag.angmom[k] = 0.0;

  // Loop over all SPH particles and add contributions to all quantities
  for (i=0; i<sph->Nsph; i++) {
    SphParticle<ndim>& part = sph->GetParticleIPointer(i);
    if (part.itype == dead) {
      diag.Ndead++;
      continue;
    }
    diag.mtot += part.m;
    diag.ketot += part.m*
      DotProduct(part.v,part.v,ndim);
    diag.utot += part.m*part.u;
    diag.gpetot -= part.m*part.gpot;
    for (k=0; k<ndim; k++) {
      diag.rcom[k] += part.m*part.r[k];
      diag.vcom[k] += part.m*part.v[k];
      diag.mom[k] += part.m*part.v[k];
      diag.force[k] += part.m*part.a[k];
      diag.force_hydro[k] += part.m*(part.a[k] - part.agrav[k]);
      diag.force_grav[k] += part.m*part.agrav[k];
    }
  }

  // Add contributions to angular momentum depending on dimensionality
  if (ndim == 2) {
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      if (part.itype == dead) continue;
      diag.angmom[2] += part.m*
        (part.r[0]*part.v[1] - part.r[1]*part.v[0]);
    }
  }
  else if (ndim == 3) {
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>& part = sph->GetParticleIPointer(i);
      if (part.itype == dead) continue;
      diag.angmom[0] += part.m*
        (part.r[1]*part.v[2] - part.r[2]*part.v[1]);
      diag.angmom[1] += part.m*
        (part.r[2]*part.v[0] - part.r[0]*part.v[2]);
      diag.angmom[2] += part.m*
        (part.r[0]*part.v[1] - part.r[1]*part.v[0]);
    }
  }

  // Loop over all star particles and add contributions to all quantities
  for (i=0; i<nbody->Nstar; i++) {
    diag.mtot += nbody->stardata[i].m;
    diag.ketot += nbody->stardata[i].m*
      DotProduct(nbody->stardata[i].v,nbody->stardata[i].v,ndim);
    diag.gpetot -= nbody->stardata[i].m*nbody->stardata[i].gpot;
    for (k=0; k<ndim; k++) {
      diag.rcom[k] += nbody->stardata[i].m*nbody->stardata[i].r[k];
      diag.vcom[k] += nbody->stardata[i].m*nbody->stardata[i].v[k];
      diag.mom[k] += nbody->stardata[i].m*nbody->stardata[i].v[k];
      diag.force[k] += nbody->stardata[i].m*nbody->stardata[i].a[k];
      diag.force_grav[k] += nbody->stardata[i].m*nbody->stardata[i].a[k];
    }
  }

  // Add contributions to angular momentum depending on dimensionality
  if (ndim == 2) {
    for (i=0; i<nbody->Nstar; i++)
      diag.angmom[2] += nbody->stardata[i].m*
        (nbody->stardata[i].r[0]*nbody->stardata[i].v[1] -
         nbody->stardata[i].r[1]*nbody->stardata[i].v[0]);
  }
  else if (ndim == 3) {
    for (i=0; i<nbody->Nstar; i++) {
      diag.angmom[0] += nbody->stardata[i].m*
        (nbody->stardata[i].r[1]*nbody->stardata[i].v[2] -
         nbody->stardata[i].r[2]*nbody->stardata[i].v[1]);
      diag.angmom[1] += nbody->stardata[i].m*
        (nbody->stardata[i].r[2]*nbody->stardata[i].v[0] -
         nbody->stardata[i].r[0]*nbody->stardata[i].v[2]);
      diag.angmom[2] += nbody->stardata[i].m*
        (nbody->stardata[i].r[0]*nbody->stardata[i].v[1] -
         nbody->stardata[i].r[1]*nbody->stardata[i].v[0]);
    }
  }

  // Add internal angular momentum (due to sink accretion) and subtract
  // accreted hydro momentum to maintain conservation of individual impulses.
  for (i=0; i<sinks.Nsink; i++) {
    for (k=0; k<3; k++) diag.angmom[k] += sinks.sink[i].angmom[k];
    for (k=0; k<3; k++) diag.force_grav[k] -= sinks.sink[i].fhydro[k];
    for (k=0; k<3; k++) diag.force_hydro[k] += sinks.sink[i].fhydro[k];
  }

  // Normalise all quantities and sum all contributions to total energy
  diag.ketot *= 0.5;
  diag.gpetot *= 0.5;
  diag.Etot = diag.ketot;
  for (k=0; k<ndim; k++) diag.rcom[k] /= diag.mtot;
  for (k=0; k<ndim; k++) diag.vcom[k] /= diag.mtot;
  if (sph->hydro_forces == 1) diag.Etot += diag.utot;
  if (sph->self_gravity == 1 || nbody->Nstar > 0) diag.Etot += diag.gpetot;


  // Calculate binary statistics if required
  if (simparams->intparams["binary_stats"] == 1) {
    nbodytree.CreateNbodySystemTree(nbody);
    nbodytree.FindBinarySystems(nbody);
  }


  // For MPI, collect all diagnostic information on root node
#ifdef MPI_PARALLEL
  mpicontrol->CollateDiagnosticsData(diag);
#endif


  return;
}



//=============================================================================
//  Simulation::OutputDiagnostics
/// Output all diagnostic quantities that are computed in
/// CalculateDiagnostics to screen.
//=============================================================================
template <int ndim>
void Simulation<ndim>::OutputDiagnostics(void)
{
  debug2("[SphSimulation::OutputDiagnostics]");

  if (rank != 0) return;

  cout << "Nsph        : " << diag.Nsph << endl;
  cout << "Nstar       : " << diag.Nstar << endl;
  if (diag.Ndead > 0)   cout << "Ndead       : " << diag.Ndead << endl;
  cout << "mtot        : " << diag.mtot*simunits.m.outscale << endl;
  cout << "Etot        : " << diag.Etot*simunits.E.outscale << endl;
  cout << "ketot       : " << diag.ketot*simunits.E.outscale << endl;
  if (sph->hydro_forces == 1)
    cout << "utot        : " << diag.utot*simunits.E.outscale << endl;
  if (sph->self_gravity == 1 || nbody->Nstar > 0)
    cout << "gpetot      : " << diag.gpetot*simunits.E.outscale << endl;
  if (ndim == 1) {
    cout << "rcom        : " << diag.rcom[0] << endl;
    cout << "vcom        : " << diag.vcom[0] << endl;
    cout << "mom         : " << diag.mom[0] << endl;
    cout << "force       : " << diag.force[0] << endl;
    if (sph->self_gravity == 1 || nbody->Nstar > 0)
      cout << "force_grav  : " << diag.force_grav[0] << endl;
    if (sph->hydro_forces == 1)
      cout << "force_hydro : " << diag.force_hydro[0] << endl;
  }
  else if (ndim == 2) {
    cout << "rcom        : " << diag.rcom[0] << "   " << diag.rcom[1] << endl;
    cout << "vcom        : " << diag.vcom[0] << "   " << diag.vcom[1] << endl;
    cout << "ang mom     : " << diag.angmom[2] << endl;
    cout << "mom         : " << diag.mom[0] << "   " << diag.mom[1] << endl;
    cout << "force       : " << diag.force[0] << "   " << diag.force[1] << endl;
    if (sph->self_gravity == 1 || nbody->Nstar > 0)
      cout << "force_grav  : " << diag.force_grav[0] << "   "
	   << diag.force_grav[1] << endl;
    if (sph->hydro_forces == 1)
      cout << "force_hydro : " << diag.force_hydro[0] << "   "
           << diag.force_hydro[1] << endl;
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
    if (sph->self_gravity == 1 || nbody->Nstar > 0)
      cout << "force_grav  : " << diag.force_grav[0] << "   "
	   << diag.force_grav[1] << "   " << diag.force_grav[2] << endl;
    if (sph->hydro_forces == 1)
      cout << "force_hydro : " << diag.force_hydro[0] << "   "
           << diag.force_hydro[1] << "   " << diag.force_hydro[2] << endl;
  }

  // Calculate binary statistics if required
  if (simparams->intparams["binary_stats"] == 1) {
    nbodytree.OutputBinaryProperties(nbody);
  }

  return;
}
