// ============================================================================
// Rendergrid.cpp
// ============================================================================


#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstring>
#include "SphParticle.h"
#include "Sph.h"
#include "SphSnapshot.h"
#include "SphKernel.h"
#include "Exception.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Debug.h"
#include "Render.h"
using namespace std;



// ============================================================================
// Render::Render
// ============================================================================
Render::Render()
{
}



// ============================================================================
// Render::~Render
// ============================================================================
Render::~Render()
{
}




// ============================================================================
// Render::CreateRenderingGrid
// ============================================================================
int Render::CreateRenderingGrid(int ixgrid, int iygrid, string xstring, 
				string ystring, string renderstring, 
				float dx_grid, float xmin, float xmax, 
				float ymin, float ymax, float *rgrid, 
				float *values, SphSnapshot &snap, 
				SphKernel *kern)
{
  int arraycheck = 1;
  int c;
  int i;
  int j;
  int k;
  int idummy;
  int Ngrid = ixgrid*iygrid;
  float dr[2];
  float drsqd;
  float drmag;
  float *xvalues;
  float *yvalues;
  float *rendervalues;
  float *mvalues;
  float *rhovalues;
  float *hvalues;
  float *rendernorm;
  float wnorm;
  float invh;
  float wkern;
  float hrangesqd;

  int ndim = snap.ndim;

  // First, verify x, y and render strings are valid
  snap.ExtractArray(xstring,&xvalues,&idummy); arraycheck *= idummy;
  snap.ExtractArray(ystring,&yvalues,&idummy); arraycheck *= idummy;
  snap.ExtractArray(renderstring,&rendervalues,&idummy); arraycheck *= idummy;
  snap.ExtractArray("m",&mvalues,&idummy); arraycheck *= idummy;
  snap.ExtractArray("rho",&rhovalues,&idummy); arraycheck *= idummy;
  snap.ExtractArray("h",&hvalues,&idummy); arraycheck *= idummy;

  // If any are invalid, exit here with failure code
  if (arraycheck == 0) return -1;

  rendernorm = new float[snap.Nsph];

  for (c=0; c<Ngrid; c++) values[c] = 0.0f;
  for (c=0; c<Ngrid; c++) rendernorm[c] = 0.0f;


  // Loop over all particles in snapshot
  // -----------------------------------------------------------------------------
  for (i=0; i<snap.Nsph; i++) {

    invh = 1.0f/hvalues[i];
    wnorm = mvalues[i]/rhovalues[i]*pow(invh,ndim);
    hrangesqd = kern->kernrangesqd*hvalues[i]*hvalues[i];

    // Now loop over all pixels and add current particles
    // ---------------------------------------------------------------------------
    for (c=0; c<Ngrid; c++) {

      dr[0] = rgrid[2*c] - xvalues[i];
      dr[1] = rgrid[2*c + 1] - yvalues[i];
      drsqd = dr[0]*dr[0] + dr[1]*dr[1];
      
      if (drsqd > hrangesqd) continue;

      drmag = sqrt(hrangesqd);
      wkern = float(kern->w0((FLOAT) drmag*invh));

      values[c] += wnorm*rendervalues[i]*wkern;
      rendernorm[c] += wnorm*wkern;
    }
    // ---------------------------------------------------------------------------

  }
  // -----------------------------------------------------------------------------

  // Normalise all grid cells
  for (c=0; c<Ngrid; c++)
    if (rendernorm[c] > 1.e-10) values[c] /= rendernorm[c];


  return 1;
}
