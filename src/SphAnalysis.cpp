//=============================================================================
//  SphAnalysis.cpp
//  Contains various analysis routines for SphSimulation object.
//=============================================================================


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstring>
#include "Exception.h"
#include "SphSimulation.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Debug.h"
using namespace std;



//=============================================================================
//  SphSimulation::CalculateDiagnostics
/// Calculates all diagnostic quantities (e.g. conserved quantities), 
/// saves to the diagnostic data structure and outputs to screen.
//=============================================================================
template <int ndim>
void SimulationDim<ndim>::CalculateDiagnostics(void)
{
  int i;                            // Particle counter
  int k;                            // Dimensionality counter

  debug2("[SphSimulation::CalculateDiagnostics]");

  // Zero all diagnostic summation variables
  diag.Etot = 0.0;
  diag.utot = 0.0;
  diag.ketot = 0.0;
  diag.gpetot = 0.0;
  for (k=0; k<ndim; k++) diag.mom[k] = 0.0;
  for (k=0; k<ndim; k++) diag.force[k] = 0.0;
  for (k=0; k<ndim; k++) diag.force_grav[k] = 0.0;
  for (k=0; k<3; k++) diag.angmom[k] = 0.0;

  // Loop over all SPH particles and add contributions to all quantities
  for (i=0; i<sph->Nsph; i++) {
    diag.ketot += sph->sphdata[i].m*
      DotProduct(sph->sphdata[i].v,sph->sphdata[i].v,ndim);
    diag.utot += sph->sphdata[i].m*sph->sphdata[i].u;
    diag.gpetot -= sph->sphdata[i].m*sph->sphdata[i].gpot;
    for (k=0; k<ndim; k++) {
      diag.mom[k] += sph->sphdata[i].m*sph->sphdata[i].v[k];
      diag.force[k] += sph->sphdata[i].m*sph->sphdata[i].a[k];
      diag.force_grav[k] += sph->sphdata[i].m*sph->sphdata[i].agrav[k];
    }
  }

  // Add contributions to angular momentum depending on dimensionality
  if (ndim == 2) {
    for (i=0; i<sph->Nsph; i++)
      diag.angmom[2] += sph->sphdata[i].m*
	(sph->sphdata[i].r[0]*sph->sphdata[i].v[1] - 
	 sph->sphdata[i].r[1]*sph->sphdata[i].v[0]);
  }
  else if (ndim == 3) {
    for (i=0; i<sph->Nsph; i++) {
      diag.angmom[0] += sph->sphdata[i].m*
	(sph->sphdata[i].r[1]*sph->sphdata[i].v[2] - 
	 sph->sphdata[i].r[2]*sph->sphdata[i].v[1]);
      diag.angmom[1] += sph->sphdata[i].m*
	(sph->sphdata[i].r[2]*sph->sphdata[i].v[0] - 
	 sph->sphdata[i].r[0]*sph->sphdata[i].v[2]);
      diag.angmom[2] += sph->sphdata[i].m*
	(sph->sphdata[i].r[0]*sph->sphdata[i].v[1] - 
	 sph->sphdata[i].r[1]*sph->sphdata[i].v[0]);
    }
  }

  // Loop over all star particles and add contributions to all quantities
  for (i=0; i<nbody->Nstar; i++) {
    diag.ketot += nbody->stardata[i].m*
      DotProduct(nbody->stardata[i].v,nbody->stardata[i].v,ndim);
    diag.gpetot -= nbody->stardata[i].m*nbody->stardata[i].gpot;
    for (k=0; k<ndim; k++) {
      diag.mom[k] += nbody->stardata[i].m*nbody->stardata[i].v[k];
      diag.force[k] += nbody->stardata[i].m*nbody->stardata[i].a[k];
      diag.force_grav[k] += nbody->stardata[i].m*nbody->stardata[i].a[k];
    }
  }

  // Normalise all quantities and sum all contributions to total energy
  diag.ketot *= 0.5;
  diag.gpetot *= 0.5;
  diag.Etot = diag.ketot;
  if (sph->hydro_forces == 1) diag.Etot += diag.utot;
  if (sph->self_gravity == 1) diag.Etot += diag.gpetot;

  return;
}



//=============================================================================
//  SphSimulation::CalculateDiagnostics
/// Calculates all diagnostic quantities (e.g. conserved quantities),
/// saves to the diagnostic data structure and outputs to screen.
//=============================================================================
template <int ndim>
void SimulationDim<ndim>::OutputDiagnostics(void)
{
  debug2("[SphSimulation::OutputDiagnostics]");

  cout << "Printing out diagnostics" << endl;
  cout << "Nsph       : " << sph->Nsph << endl;
  cout << "Nstar      : " << nbody->Nstar << endl;
  cout << "Etot       : " << diag.Etot << endl;
  cout << "ketot      : " << diag.ketot << endl;
  if (sph->hydro_forces == 1) cout << "utot       : " << diag.utot << endl;
  cout << "gpetot     : " << diag.gpetot << endl;
  if (ndim == 1) {
    cout << "mom        : " << diag.mom[0] << endl;
    cout << "force      : " << diag.force[0] << endl;
    cout << "force_grav : " << diag.force_grav[0] << endl;
  }
  else if (ndim == 2) {
    cout << "mom        : " << diag.mom[0] << "   " << diag.mom[1] << endl;
    cout << "force      : " << diag.force[0] << "   " << diag.force[1] << endl;
    cout << "force_grav : " << diag.force_grav[0] << "   " 
	 << diag.force_grav[1] << endl;
    cout << "ang mom    : " << diag.angmom[2] << endl;
  }
  else if (ndim == 3) {
    cout << "mom        : " << diag.mom[0] << "   " 
	 << diag.mom[1] << "   " << diag.mom[2] << endl;
    cout << "force      : " << diag.force[0] << "   " 
	 << diag.force[1] << "   " << diag.force[2] << endl;
    cout << "force_grav : " << diag.force_grav[0] << "   " 
	 << diag.force_grav[1] << "   " << diag.force_grav[2] << endl;
    cout << "ang mom    : " << diag.angmom[0] << "   "
	 << diag.angmom[1] << "   " << diag.angmom[2] << endl;
  }

  return;
}
