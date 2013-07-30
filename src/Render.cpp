//=============================================================================
//  Render.cpp
//  Contains all functions for generating rendered images from SPH data.
//=============================================================================


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


//=============================================================================
//  RenderBase::RenderFactory
/// Create new render object for simulation object depending on dimensionality.
//=============================================================================
RenderBase* RenderBase::RenderFactory(int ndim, SimulationBase* sim) {
  RenderBase* render;
  if (ndim==1) {
    render = new Render<1> (sim);
  }
  else if (ndim==2) {
    render = new Render<2> (sim);
  }
  else if (ndim==3) {
    render = new Render<3> (sim);
  }
  else {
    render = NULL;
  }
  return render;
}



//=============================================================================
//  Render::Render
/// Render class constructor
//=============================================================================
template <int ndim>
Render<ndim>::Render(SimulationBase* sim):
sph(static_cast<Sph<ndim>* > (static_cast<Simulation<ndim>* > (sim)->sph))
{
}



//=============================================================================
//  Render::~Render
/// Render class destructor
//=============================================================================
template <int ndim>
Render<ndim>::~Render()
{
}



//=============================================================================
//  Render::CreateColumnRenderingGrid
/// Calculate column integrated SPH averaged quantities on a grid for use in
/// generated rendered images in python code.
//=============================================================================
template <int ndim>
int Render<ndim>::CreateColumnRenderingGrid
(int ixgrid,                       ///< [in] No. of x-grid spacings
 int iygrid,                       ///< [in] No. of y-grid spacings
 string xstring,                   ///< [in] x-axis quantity
 string ystring,                   ///< [in] y-axis quantity
 string renderstring,              ///< [in] Rendered quantity
 string renderunit,                ///< [in] Required unit of rendered quantity
 float xmin,                       ///< [in] Minimum x-extent
 float xmax,                       ///< [in] Maximum x-extent
 float ymin,                       ///< [in] Minimum y-extent
 float ymax,                       ///< [in] Maximum y-extent
 float* values,                    ///< [out] Rendered values for plotting
 int Ngrid,                        ///< [in] No. of grid points (ixgrid*iygrid)
 SphSnapshotBase &snap,            ///< [inout] Snapshot object reference
 float &scaling_factor)            ///< [in] Rendered quantity scaling factor
{
  int arraycheck = 1;              // Verification flag
  int c;                           // Rendering grid cell counter
  int i;                           // Particle counter
  int j;                           // Aux. counter
  int k;                           // Dimension counter
  int idummy;                      // ..
  int Nsph = snap.Nsph;            // ..
  float dr[2];                     // Rel. position vector on grid plane
  float drsqd;                     // Distance squared on grid plane
  float drmag;                     // Distance
  float dummyfloat = 0.0;          // ..
  float hrangesqd;                 // Kernel range squared
  float invh;                      // 1/h
  float wkern;                     // Kernel value
  float wnorm;                     // Kernel normalisation value
  float *xvalues;                  // Pointer to 'x' array
  float *yvalues;                  // Pointer to 'y' array
  float *rendervalues;             // Pointer to rendered quantity array
  float *mvalues;                  // Pointer to mass array
  float *rhovalues;                // Pointer to density array
  float *hvalues;                  // Pointer to smoothing length array
  float *rendernorm;               // Normalisation array
  float *rgrid;                    // Grid positions
  string dummystring = "";         // ..

  // Check x and y strings are actual co-ordinate strings
  if ((xstring != "x" && xstring != "y" && xstring != "z") ||
      (ystring != "x" && ystring != "y" && ystring != "z")) return -1;

  // First, verify x, y and render strings are valid
  snap.ExtractArray(xstring,"sph",&xvalues,&idummy,dummyfloat,dummystring);
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(ystring,"sph",&yvalues,&idummy,dummyfloat,dummystring);
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(renderstring,"sph",&rendervalues,&idummy,
                    scaling_factor,renderunit);
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("m","sph",&mvalues,&idummy,dummyfloat,dummystring);
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("rho","sph",&rhovalues,&idummy,dummyfloat,dummystring);
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("h","sph",&hvalues,&idummy,dummyfloat,dummystring);
  arraycheck = min(idummy,arraycheck);

  // If any are invalid, exit here with failure code
  if (arraycheck == 0) return -1;

  // Allocate temporary memory for creating render grid
  rendernorm = new float[Ngrid];
  rgrid = new float[2*Ngrid];

  // Create grid positions here (need to improve in the future)
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
#pragma omp parallel for default(none) private(c,dr,drmag,drsqd,hrangesqd,invh,i)\
  private(wkern,wnorm) \
  shared(hvalues,mvalues,Nsph,rhovalues,rgrid,xvalues,yvalues,rendervalues,values) \
  shared(rendernorm,Ngrid)
    for (i=0; i<Nsph; i++) {
      invh = 1.0/hvalues[i];
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



//=============================================================================
//  Render::CreateSliceRenderingGrid
/// Calculate gridded SPH properties on a slice for slice-rendering.
//=============================================================================
template <int ndim>
int Render<ndim>::CreateSliceRenderingGrid
(int ixgrid,                       ///< [in] No. of x-grid spacings
 int iygrid,                       ///< [in] No. of y-grid spacings
 string xstring,                   ///< [in] x-axis quantity
 string ystring,                   ///< [in] y-axis quantity
 string zstring,                   ///< [in] z-axis quantity
 string renderstring,              ///< [in] Rendered quantity
 string renderunit,                ///< [in] Required unit of rendered quantity
 float xmin,                       ///< [in] Minimum x-extent
 float xmax,                       ///< [in] Maximum x-extent
 float ymin,                       ///< [in] Minimum y-extent
 float ymax,                       ///< [in] Maximum y-extent
 float zslice,                     ///< [in] z-position of slice
 float* values,                    ///< [out] Rendered values for plotting
 int Ngrid,                        ///< [in] No. of grid points (ixgrid*iygrid)
 SphSnapshotBase &snap,            ///< [inout] Snapshot object reference
 float &scaling_factor)            ///< [in] Rendered quantity scaling factor
{
  int arraycheck = 1;              // Verification flag
  int c;                           // Rendering grid cell counter
  int i;                           // Particle counter
  int j;                           // Aux. counter
  int k;                           // Dimension counter
  int idummy;                      // ..
  float dr[3];                     // Rel. position vector on grid plane
  float drsqd;                     // Distance squared on grid plane
  float drmag;                     // Distance
  float dummyfloat = 0.0;          // ..
  float hrangesqd;                 // Kernel range squared
  float invh;                      // 1/h
  float skern;                     // r/h
  float wkern;                     // Kernel value
  float wnorm;                     // Kernel normalisation value
  float *xvalues;                  // Pointer to 'x' array
  float *yvalues;                  // Pointer to 'y' array
  float *zvalues;                  // Pointer to 'z' array
  float *rendervalues;             // Pointer to rendered quantity array
  float *mvalues;                  // Pointer to mass array
  float *rhovalues;                // Pointer to density array
  float *hvalues;                  // Pointer to smoothing length array
  float *rendernorm;               // Normalisation array
  float *rgrid;                    // Grid positions
  string dummystring = "";         // ..


  // Check x and y strings are actual co-ordinate strings
  if ((xstring != "x" && xstring != "y" && xstring != "z") ||
	  (ystring != "x" && ystring != "y" && ystring != "z")) return -1;

  // First, verify x, y and render strings are valid
  snap.ExtractArray(xstring,"sph",&xvalues,&idummy,dummyfloat,dummystring); 
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(ystring,"sph",&yvalues,&idummy,dummyfloat,dummystring); 
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(zstring,"sph",&zvalues,&idummy,dummyfloat,dummystring); 
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray(renderstring,"sph",&rendervalues,&idummy,
                    scaling_factor,renderunit); 
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("m","sph",&mvalues,&idummy,dummyfloat,dummystring); 
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("rho","sph",&rhovalues,&idummy,dummyfloat,dummystring); 
  arraycheck = min(idummy,arraycheck);
  snap.ExtractArray("h","sph",&hvalues,&idummy,dummyfloat,dummystring); 
  arraycheck = min(idummy,arraycheck);

  // If any are invalid, exit here with failure code
  if (arraycheck == 0) return -1;

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
  for (c=0; c<Ngrid; c++) values[c] = 0.0f;
  for (c=0; c<Ngrid; c++) rendernorm[c] = 0.0f;


  // Loop over all particles in snapshot
  // --------------------------------------------------------------------------
  for (i=0; i<snap.Nsph; i++) {
    invh = 1.0/hvalues[i];
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
      skern = (FLOAT) (drmag*invh);
      wkern = float(sph->kerntab.w0(skern));

#pragma omp atomic
      values[c] += wnorm*rendervalues[i]*wkern;
#pragma omp atomic
      rendernorm[c] += wnorm*wkern;

    }
    // ------------------------------------------------------------------------

  }
  // --------------------------------------------------------------------------

  // Normalise all grid cells
  for (c=0; c<Ngrid; c++)
    if (rendernorm[c] > 1.e-10) values[c] /= rendernorm[c];

  // Free all locally allocated memory
  delete[] rgrid;
  delete[] rendernorm;

  return 1;
}



template class Render<1>;
template class Render<2>;
template class Render<3>;
