// ============================================================================
// Render.cpp
// Contains all functions for generating rendered images from SPH data.
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
// Calculate column integrated SPH averaged quantities on a grid for use in 
// generated rendered images in python code.
// ============================================================================
int Render::CreateColumnRenderingGrid(int ixgrid, int iygrid, string xstring,
				      string ystring, string renderstring,
				      string renderunit, float xmin, 
				      float xmax,
				      float ymin, float ymax, float* values,
				      int Ngrid, SphSnapshot &snap,
				      Sph *sph, float &scaling_factor)
{
  int arraycheck = 1;                   // Verification flag
  int c;                                // Rendering grid cell counter
  int i;                                // Particle counter
  int j;                                // Aux. counter
  int k;                                // Dimension counter
  int idummy;                           // ..
  int ndim = snap.ndim;                 // Local copy of snapshot ndim
  float dr[2];                          // Rel. position vector on grid plane
  float drsqd;                          // Distance squared on grid plane
  float drmag;                          // Distance
  float wnorm;                          // Kernel normalisation value
  float invh;                           // 1/h
  float wkern;                          // Kernel value
  float hrangesqd;                      // Kernel range squared
  float dummyfloat = 0.0;               // ..
  float *xvalues;                       // Pointer to 'x' array
  float *yvalues;                       // Pointer to 'y' array
  float *rendervalues;                  // Pointer to rendered quantity array
  float *mvalues;                       // Pointer to mass array
  float *rhovalues;                     // Pointer to density array
  float *hvalues;                       // Pointer to smoothing length array
  float *rendernorm;                    // Normalisation array
  float *rgrid;                         // Grid positions
  string dummystring = "";              // ..

  // Check x and y strings are actual co-ordinate strings
  if ((xstring != "x" && xstring != "y" && xstring != "z") ||
      (ystring != "x" && ystring != "y" && ystring != "z")) return -1;

  // First, verify x, y and render strings are valid
  snap.ExtractArray(xstring,&xvalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(ystring,&yvalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(renderstring,&rendervalues,&idummy,scaling_factor,renderunit); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("m",&mvalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("rho",&rhovalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("h",&hvalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);

  // If any are invalid, exit here with failure code
  if (arraycheck == 0) return -1;

  // Allocate temporary memory for creating render grid
  rendernorm = new float[Ngrid];
  rgrid = new float[2*Ngrid];

  // Create grid positions here
  c = 0;
  for (j=iygrid-1; j>=0; j--) {
    for (i=0; i<ixgrid; i++) {
      rgrid[2*c] = xmin + ((float) i + 0.5f)*(xmax - xmin)/(float)ixgrid;
      rgrid[2*c + 1] = ymin + ((float) j + 0.5f)*(ymax - ymin)/(float)iygrid;
      c++;
    }
  }

  // Zero arrays before computing rendering
  for (c=0; c<Ngrid; c++) values[c] = (float) 0.0;
  for (c=0; c<Ngrid; c++) rendernorm[c] = (float) 0.0;


  // Create rendered grid depending on dimensionality
  // ==========================================================================
  if (ndim == 2) {

    // Loop over all particles in snapshot
    // -----------------------------------------------------------------------
#pragma omp parallel for default(shared) private(c,dr,drmag,drsqd,hrangesqd,invh,wkern,wnorm)
    for (i=0; i<snap.Nsph; i++) {
      invh = 1.0f/hvalues[i];
      wnorm = mvalues[i]/rhovalues[i]*pow(invh,ndim);
      hrangesqd = sph->kerntab.kernrangesqd*hvalues[i]*hvalues[i];
      
      // Now loop over all pixels and add current particles
      // ----------------------------------------------------------------------
      for (c=0; c<Ngrid; c++) {
	
    	dr[0] = rgrid[2*c] - xvalues[i];
	dr[1] = rgrid[2*c + 1] - yvalues[i];
	drsqd = dr[0]*dr[0] + dr[1]*dr[1];
	
	if (drsqd > hrangesqd) continue;
	
	drmag = sqrt(drsqd);
	wkern = float(sph->kerntab.w0((FLOAT) (drmag*invh)));
	
#pragma omp atomic
	values[c] += wnorm*rendervalues[i]*wkern;
#pragma omp atomic
	rendernorm[c] += wnorm*wkern;
      }
      // ----------------------------------------------------------------------
      
    }
    // ------------------------------------------------------------------------

    // Normalise all grid cells
    for (c=0; c<Ngrid; c++) {
      if (rendernorm[c] > 1.e-10) values[c] /= rendernorm[c];
    }


  }
  // ==========================================================================
  else if (ndim == 3) {

    // Loop over all particles in snapshot
    // ------------------------------------------------------------------------
    for (i=0; i<snap.Nsph; i++) {
      invh = 1.0f/hvalues[i];
      wnorm = mvalues[i]/rhovalues[i]*pow(invh,(ndim - 1));
      hrangesqd = sph->kerntab.kernrangesqd*hvalues[i]*hvalues[i];
      
      // Now loop over all pixels and add current particles
      // ----------------------------------------------------------------------
      for (c=0; c<Ngrid; c++) {
	
    	dr[0] = rgrid[2*c] - xvalues[i];
	dr[1] = rgrid[2*c + 1] - yvalues[i];
	drsqd = dr[0]*dr[0] + dr[1]*dr[1];
	
	if (drsqd > hrangesqd) continue;
	
	drmag = sqrt(drsqd);
	wkern = float(sph->kerntab.wLOS((FLOAT) (drmag*invh)));
	
	values[c] += wnorm*rendervalues[i]*wkern;
	rendernorm[c] += wnorm*wkern;
      }
      // ----------------------------------------------------------------------

    }
    // ------------------------------------------------------------------------

  }
  // ==========================================================================


  // Free all locally allocated memory
  delete[] rgrid;
  delete[] rendernorm;

  return 1;
}



// ============================================================================
// Render::CreateSliceRenderingGrid
// ============================================================================
int Render::CreateSliceRenderingGrid(int ixgrid, int iygrid, string xstring,
				     string ystring, string zstring, 
				     string renderstring,
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
  snap.ExtractArray(xstring,&xvalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(ystring,&yvalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(zstring,&zvalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(renderstring,&rendervalues,&idummy,scaling_factor,renderunit); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("m",&mvalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("rho",&rhovalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("h",&hvalues,&idummy,dummyfloat,dummystring); arraycheck = min(idummy,arraycheck);

  // If any are invalid, exit here with failure code
  if (arraycheck == 0) return -1;

  rendernorm = new float[Ngrid];

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

  // Zero arrays before computing rendering
  for (c=0; c<Ngrid; c++) values[c] = 0.0f;
  for (c=0; c<Ngrid; c++) rendernorm[c] = 0.0f;


  // Loop over all particles in snapshot
  // --------------------------------------------------------------------------
  for (i=0; i<snap.Nsph; i++) {

    invh = 1.0f/hvalues[i];
    wnorm = mvalues[i]/rhovalues[i]*pow(invh,ndim);
    hrangesqd = sph->kerntab.kernrangesqd*hvalues[i]*hvalues[i];

    // Now loop over all pixels and add current particles
    // ------------------------------------------------------------------------
    for (c=0; c<Ngrid; c++) {

      dr[0] = rgrid[2*c] - xvalues[i];
      dr[1] = rgrid[2*c + 1] - yvalues[i];
      dr[2] = zslice - zvalues[i];
      drsqd = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
      
      if (drsqd > hrangesqd) continue;

      drmag = sqrt(drsqd);
      wkern = float(sph->kerntab.w0((FLOAT) (drmag*invh)));

      values[c] += wnorm*rendervalues[i]*wkern;
      rendernorm[c] += wnorm*wkern;
    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------

  // Normalise all grid cells
  for (c=0; c<Ngrid; c++)
    if (rendernorm[c] > 1.e-10) values[c] /= rendernorm[c];

  return 1;
}


