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
// Render::CreateColumnRenderingGrid
// ============================================================================
int Render::CreateRenderingGrid(int ixgrid, int iygrid, string xstring,
				string ystring, string renderstring,
				string renderunit, float xmin, float xmax,
				float ymin, float ymax, float* values,
				int Ngrid, SphSnapshot &snap,
				Sph *sph, float &scaling_factor)
{
  int arraycheck = 1;                   // ..
  int c;                                // ..
  int i;                                // ..
  int j;                                // ..
  int k;                                // ..
  int idummy;                           // ..
  int ndim = snap.ndim;                 // Local copy of snapshot ndim
  float dr[2];                          // ..
  float drsqd;                          // ..
  float drmag;                          // ..
  float wnorm;                          // ..
  float invh;                           // ..
  float wkern;                          // ..
  float hrangesqd;                      // ..
  float dummyfloat = 0.0;               // ..
  float *xvalues;                       // ..
  float *yvalues;                       // ..
  float *rendervalues;                  // ..
  float *mvalues;                       // ..
  float *rhovalues;                     // ..
  float *hvalues;                       // ..
  float *rendernorm;                    // ..
  float *rgrid;                         // ..
  string dummystring = "";              // ..

  // Check x and y strings are actual co-ordinate strings
  if ((xstring != "x" && xstring != "y" && xstring != "z") ||
	  (ystring != "x" && ystring != "y" && ystring != "z")) return -1;

  // First, verify x, y and render strings are valid
  snap.ExtractArray(xstring,&xvalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;
  snap.ExtractArray(ystring,&yvalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;
  snap.ExtractArray(renderstring,&rendervalues,&idummy,scaling_factor,renderunit); arraycheck *= idummy;
  snap.ExtractArray("m",&mvalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;
  snap.ExtractArray("rho",&rhovalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;
  snap.ExtractArray("h",&hvalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;

  // If any are invalid, exit here with failure code
  if (arraycheck == 0) return -1;

  rendernorm = new float[snap.Nsph];

  // Create grid positions here
  c = 0;
  rgrid = new float[2*Ngrid];
  for (j=iygrid-1; j>=0; j--) {
	  for (i=0; i<ixgrid; i++) {
		  rgrid[2*c] = xmin + ((float) i + 0.5f)*(xmax - xmin)/(float)ixgrid;
		  rgrid[2*c + 1] = ymin + ((float) j + 0.5f)*(ymax - ymin)/(float)iygrid;
		  c++;
	  }
  }

  for (c=0; c<Ngrid; c++) values[c] = 0.0f;
  for (c=0; c<Ngrid; c++) rendernorm[c] = 0.0f;


  // Loop over all particles in snapshot
  // -----------------------------------------------------------------------------
  for (i=0; i<snap.Nsph; i++) {

    invh = 1.0f/hvalues[i];
    wnorm = mvalues[i]/rhovalues[i]*pow(invh,ndim);
    hrangesqd = sph->kernp->kernrangesqd*hvalues[i]*hvalues[i];

    // Now loop over all pixels and add current particles
    // ---------------------------------------------------------------------------
    for (c=0; c<Ngrid; c++) {

      dr[0] = rgrid[2*c] - xvalues[i];
      dr[1] = rgrid[2*c + 1] - yvalues[i];
      drsqd = dr[0]*dr[0] + dr[1]*dr[1];

      if (drsqd > hrangesqd) continue;

      drmag = sqrt(drsqd);
      wkern = float(sph->kernp->w0((FLOAT) (drmag*invh)));

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




// ============================================================================
// Render::CreateSliceRenderingGrid
// ============================================================================
int Render::CreateSliceRenderingGrid(int ixgrid, int iygrid, string xstring,
				string ystring, string zstring, string renderstring,
				string renderunit, float xmin, float xmax,
				float ymin, float ymax, float* values,
				int Ngrid, SphSnapshot &snap, 
				Sph *sph, float &scaling_factor)
{
  int arraycheck = 1;                   // ..
  int c;                                // ..
  int i;                                // ..
  int j;                                // ..
  int k;                                // ..
  int idummy;                           // ..
  int ndim = snap.ndim;                 // Local copy of snapshot ndim
  float dr[2];                          // ..
  float drsqd;                          // ..
  float drmag;                          // ..
  float wnorm;                          // ..
  float invh;                           // ..
  float wkern;                          // ..
  float hrangesqd;                      // ..
  float dummyfloat = 0.0;               // ..
  float *xvalues;                       // ..
  float *yvalues;                       // ..
  float *zvalues;                       // ..
  float *rendervalues;                  // ..
  float *mvalues;                       // ..
  float *rhovalues;                     // ..
  float *hvalues;                       // ..
  float *rendernorm;                    // ..
  float *rgrid;                         // ..
  float zslice = 0.0;
  string dummystring = "";              // ..

  // Check x and y strings are actual co-ordinate strings
  if ((xstring != "x" && xstring != "y" && xstring != "z") ||
	  (ystring != "x" && ystring != "y" && ystring != "z")) return -1;

  // First, verify x, y and render strings are valid
  snap.ExtractArray(xstring,&xvalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;
  snap.ExtractArray(ystring,&yvalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;
  snap.ExtractArray(zstring,&zvalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;
  snap.ExtractArray(renderstring,&rendervalues,&idummy,scaling_factor,renderunit); arraycheck *= idummy;
  snap.ExtractArray("m",&mvalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;
  snap.ExtractArray("rho",&rhovalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;
  snap.ExtractArray("h",&hvalues,&idummy,dummyfloat,dummystring); arraycheck *= idummy;

  // If any are invalid, exit here with failure code
  if (arraycheck == 0) return -1;

  rendernorm = new float[snap.Nsph];

  // Create grid positions here
  c = 0;
  rgrid = new float[2*Ngrid];
  for (j=iygrid-1; j>=0; j--) {
	  for (i=0; i<ixgrid; i++) {
		  rgrid[2*c] = xmin + ((float) i + 0.5f)*(xmax - xmin)/(float)ixgrid;
		  rgrid[2*c + 1] = ymin + ((float) j + 0.5f)*(ymax - ymin)/(float)iygrid;
		  c++;
	  }
  }

  for (c=0; c<Ngrid; c++) values[c] = 0.0f;
  for (c=0; c<Ngrid; c++) rendernorm[c] = 0.0f;


  // Loop over all particles in snapshot
  // -----------------------------------------------------------------------------
  for (i=0; i<snap.Nsph; i++) {

    invh = 1.0f/hvalues[i];
    wnorm = mvalues[i]/rhovalues[i]*pow(invh,ndim);
    hrangesqd = sph->kernp->kernrangesqd*hvalues[i]*hvalues[i];

    // Now loop over all pixels and add current particles
    // ---------------------------------------------------------------------------
    for (c=0; c<Ngrid; c++) {

      dr[0] = rgrid[2*c] - xvalues[i];
      dr[1] = rgrid[2*c + 1] - yvalues[i];
      dr[2] = zslice - zvalues[i];
      drsqd = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      
      if (drsqd > hrangesqd) continue;

      drmag = sqrt(drsqd);
      wkern = float(sph->kernp->w0((FLOAT) (drmag*invh)));

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


