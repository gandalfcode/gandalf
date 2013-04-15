//=============================================================================
// SphAnalysis.cpp
// Contains various analysis routines for SphSimulation object.
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
void SphSimulation<ndim>::CalculateDiagnostics(void)
{
  int k;                            // Dimensionality counter

  debug2("[SphSimulation::CalculateDiagnostics]");

  diag.Etot = 0.0;
  diag.utot = 0.0;
  diag.ketot = 0.0;
  diag.gpetot = 0.0;
  for (k=0; k<ndim; k++) diag.mom[k] = 0.0;
  for (k=0; k<3; k++) diag.angmom[k] = 0.0;
  for (k=0; k<ndim; k++) diag.force[k] = 0.0;
  for (k=0; k<ndim; k++) diag.force_grav[k] = 0.0;

  for (int i=0; i<sph->Nsph; i++) {
    diag.ketot += sph->sphdata[i].m*
      DotProduct(sph->sphdata[i].v,sph->sphdata[i].v,ndim);
    diag.utot += sph->sphdata[i].m*sph->sphdata[i].u;
    diag.gpetot += sph->sphdata[i].m*sph->sphdata[i].gpot;
    for (k=0; k<ndim; k++) {
      diag.mom[k] += sph->sphdata[i].m*sph->sphdata[i].v[k];
      diag.force[k] += sph->sphdata[i].m*sph->sphdata[i].a[k];
      diag.force_grav[k] += sph->sphdata[i].m*sph->sphdata[i].agrav[k];
    }
  }
  diag.ketot *= 0.5;
  diag.gpetot *= 0.5;
  diag.Etot = diag.ketot + diag.utot + diag.gpetot;

  cout << "Printing out diagnostics" << endl;
  cout << "Etot       : " << diag.Etot << endl;
  cout << "utot       : " << diag.utot << endl;
  cout << "ketot      : " << diag.ketot << endl;
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
  }
  else if (ndim == 3) {
    cout << "mom        : " << diag.mom[0] << "   " 
	 << diag.mom[1] << "   " << diag.mom[2] << endl;
    cout << "force      : " << diag.force[0] << "   " 
	 << diag.force[1] << "   " << diag.force[2] << endl;
    cout << "force_grav : " << diag.force_grav[0] << "   " 
	 << diag.force_grav[1] << "   " << diag.force_grav[2] << endl;
  }

  return;
}
