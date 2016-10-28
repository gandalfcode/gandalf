//=================================================================================================
//  Render.cpp
//  Contains all functions for generating rendered images from SPH snapshots.
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
#include <sstream>
#include <string>
#include <math.h>
#include <cstdio>
#include <cstring>
#include "Particle.h"
#include "Sph.h"
#include "SphSnapshot.h"
#include "SmoothingKernel.h"
#include "Exception.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Debug.h"
#include "Render.h"
using namespace std;


//=================================================================================================
//  RenderBase::RenderFactory
/// Create new render object for simulation object depending on dimensionality.
//=================================================================================================
RenderBase* RenderBase::RenderFactory
 (int ndim,                            ///< Simulation dimensionality
  SimulationBase* sim)                 ///< Simulation object pointer
{
  RenderBase* render;                  // Pointer to new render object
  if (ndim == 1) {
    render = new Render<1> (sim);
  }
  else if (ndim == 2) {
    render = new Render<2> (sim);
  }
  else if (ndim == 3) {
    render = new Render<3> (sim);
  }
  else {
    render = NULL;
  }
  return render;
}



//=================================================================================================
//  Render::Render
/// Render class constructor.
//=================================================================================================
template <int ndim>
Render<ndim>::Render(SimulationBase* sim):
  hydro(static_cast<Hydrodynamics<ndim>* > (static_cast<Simulation<ndim>* > (sim)->hydro))
{
}



//=================================================================================================
//  Render::~Render
/// Render class destructor.
//=================================================================================================
template <int ndim>
Render<ndim>::~Render()
{
}



//=================================================================================================
//  Render::CreateColumnRenderingGrid
/// Calculate column integrated SPH averaged quantities on a grid for use in
/// generated rendered images in python code.
//=================================================================================================
template <int ndim>
int Render<ndim>::CreateColumnRenderingGrid
 (const int ixgrid,                    ///< [in] No. of x-grid spacings
  const int iygrid,                    ///< [in] No. of y-grid spacings
  const string xstring,                ///< [in] x-axis quantity
  const string ystring,                ///< [in] y-axis quantity
  const string renderstring,           ///< [in] Rendered quantity
  const string renderunit,             ///< [in] Required unit of rendered quantity
  const float xmin,                    ///< [in] Minimum x-extent
  const float xmax,                    ///< [in] Maximum x-extent
  const float ymin,                    ///< [in] Minimum y-extent
  const float ymax,                    ///< [in] Maximum y-extent
  float* values,                       ///< [out] Rendered values for plotting
  const int Ngrid,                     ///< [in] No. of grid points (ixgrid*iygrid)
  SphSnapshotBase &snap,               ///< [inout] Snapshot object reference,
  const string typepart,                         ///< [in] Type of particles to render
  float &scaling_factor)               ///< [in] Rendered quantity scaling factor
{
  int arraycheck = 1;                  // Verification flag
  int c;                               // Rendering grid cell counter
  int i;                               // Particle counter
  int j;                               // Aux. counter
  int ii,jj;                           // ..
  int idummy;                          // Dummy integer to verify valid arrays
  int imin,imax,jmin,jmax;             // Grid limits for column density interpolation
  int Nhydro = snap.GetNparticlesType(typepart);                // No. of SPH particles in snap
  float dx,dy,invdx,invdy;             // ..
  float dr[2];                         // Rel. position vector on grid plane
  float drsqd;                         // Distance squared on grid plane
  float drmag;                         // Distance
  float dummyfloat = 0.0f;             // Dummy variable for function argument
  float invh;                          // 1/h
  float wkern;                         // Kernel value
  float wnorm;                         // Kernel normalisation value
  float *xvalues;                      // Pointer to 'x' array
  float *yvalues;                      // Pointer to 'y' array
  float *rendervalues;                 // Pointer to rendered quantity array
  float *mvalues;                      // Pointer to mass array
  float *rhovalues;                    // Pointer to density array
  float *hvalues;                      // Pointer to smoothing length array
  float *rendernorm;                   // Normalisation array
  float *rgrid;                        // Grid positions
  string dummystring = "";             // Dummy string for function argument

  // Check x and y strings are actual co-ordinate strings
  if ((xstring != "x" && xstring != "y" && xstring != "z") ||
      (ystring != "x" && ystring != "y" && ystring != "z")) return -1;

  cout << "Generating rendered image!" << endl;

  // First, verify x, y, m, rho, h and render strings are valid
  snap.ExtractArray(xstring, typepart, &xvalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray(ystring, typepart, &yvalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray(renderstring, typepart, &rendervalues, &idummy, scaling_factor, renderunit);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray("m", typepart, &mvalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray("rho", typepart, &rhovalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray("h", typepart, &hvalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);

  // If any are invalid, exit here with failure code
  if (arraycheck == 0) return -1;

  // Allocate temporary memory for creating render grid
  rendernorm = new float[Ngrid];
  rgrid = new float[2*Ngrid];

  // Create grid positions here (need to improve in the future)
  dx = (xmax - xmin) / (float) ixgrid;
  dy = (ymax - ymin) / (float) iygrid;
  invdx = 1.0f/dx;
  invdy = 1.0f/dy;
  c = 0;
  for (j=iygrid-1; j>=0; j--) {
    for (i=0; i<ixgrid; i++) {
      rgrid[2*c] = xmin + ((float) i + 0.5f)*dx;
      rgrid[2*c + 1] = ymin + ((float) j + 0.5f)*dy;
      c++;
    }
  }

  // Zero arrays before computing rendered values
  for (c=0; c<Ngrid; c++) values[c] = 0.0f;
  for (c=0; c<Ngrid; c++) rendernorm[c] = 0.0f;



// Loop over all particles in snapshot
//=================================================================================================
#pragma omp parallel for default(none) private(c,dr,drmag,drsqd,i,ii,imax,imin)\
 private(invh,jj,jmax,jmin,wkern,wnorm) shared(cout,dx,dy,invdx,invdy,hvalues,mvalues)\
 shared(Nhydro,rendernorm,rendervalues,rhovalues,rgrid,values,xvalues,yvalues)
  for (i=0; i<Nhydro; i++) {
    const float hrange = hydro->kerntab.kernrange*hvalues[i];
    //const float hrangesqd = hydro->kerntab.kernrangesqd*hvalues[i]*hvalues[i];

    if (xvalues[i] + hrange < xmin || xvalues[i] - hrange > xmax ||
        yvalues[i] + hrange < ymin || yvalues[i] - hrange > ymax) continue;

    // Compute grid coordinate limits for computing kernel convolution
    imin = (xvalues[i] - hrange - xmin)*invdx;  imin = max(0,imin);  imin = min(ixgrid-1,imin);
    imax = (xvalues[i] + hrange - xmin)*invdx;  imax = max(0,imax);  imax = min(ixgrid-1,imax);
    jmin = (yvalues[i] - hrange - ymin)*invdy;  jmin = max(0,jmin);  jmin = min(iygrid-1,jmin);
    jmax = (yvalues[i] + hrange - ymin)*invdy;  jmax = max(0,jmax);  jmax = min(iygrid-1,jmax);

    // If kernel does not overlap rendered image, skip particle entirely
    //if (imin == imax || jmin == jmax) continue;
    invh = 1.0/hvalues[i];
    wnorm = mvalues[i]/rhovalues[i]*pow(invh,ndim);

    // Now loop over all pixels for particle
    //---------------------------------------------------------------------------------------------
    for (jj=jmin; jj<=jmax; jj++) {
      for (ii=imin; ii<=imax; ii++) {
        c = ii + (iygrid - jj - 1)*ixgrid;
        dr[0] = xmin + dx*(float) ii - xvalues[i];
        dr[1] = ymin + dy*(float) jj - yvalues[i];
        drsqd = dr[0]*dr[0] + dr[1]*dr[1];
        //if (drsqd > hrangesqd) continue;
        drmag = sqrt(drsqd);
        if (ndim == 3) wkern = float(hydro->kerntab.wLOS((FLOAT) (drmag*invh)));
        else if (ndim == 2) wkern = float(hydro->kerntab.w0((FLOAT) (drmag*invh)));
        else if (ndim == 1) wkern = 0.0f;
#pragma omp atomic
        values[c] += wnorm*rendervalues[i]*wkern;
#pragma omp atomic
        rendernorm[c] += wnorm*wkern;
      }
    }
   //----------------------------------------------------------------------------------------------

  }
  //===============================================================================================

  // Normalise all grid cells
  for (c=0; c<Ngrid; c++) {
    if (rendernorm[c] > 1.e-10) values[c] /= rendernorm[c];
  }

  // Free all locally allocated memory
  delete[] rgrid;
  delete[] rendernorm;

  return 1;
}



//=================================================================================================
//  Render::CreateSliceRenderingGrid
/// Calculate gridded SPH properties on a slice for slice-rendering.
//=================================================================================================
template <int ndim>
int Render<ndim>::CreateSliceRenderingGrid
 (const int ixgrid,                    ///< [in] No. of x-grid spacings
  const int iygrid,                    ///< [in] No. of y-grid spacings
  const string xstring,                ///< [in] x-axis quantity
  const string ystring,                ///< [in] y-axis quantity
  const string zstring,                ///< [in] z-axis quantity
  const string renderstring,           ///< [in] Rendered quantity
  const string renderunit,             ///< [in] Required unit of rendered quantity
  const float xmin,                    ///< [in] Minimum x-extent
  const float xmax,                    ///< [in] Maximum x-extent
  const float ymin,                    ///< [in] Minimum y-extent
  const float ymax,                    ///< [in] Maximum y-extent
  const float zslice,                  ///< [in] z-position of slice
  float* values,                       ///< [out] Rendered values for plotting
  const int Ngrid,                     ///< [in] No. of grid points (ixgrid*iygrid)
  SphSnapshotBase &snap,               ///< [inout] Snapshot object reference,
  const string typepart,               ///< [in] Type of particles to render
  float &scaling_factor)               ///< [in] Rendered quantity scaling factor
{
  int arraycheck = 1;                  // Verification flag
  int c;                               // Rendering grid cell counter
  int i;                               // Particle counter
  int j;                               // Aux. counter
  int idummy;                          // Dummy integer to verify correct array
  int ii,jj;                           // ..
  int imin,imax,jmin,jmax;             // Grid limits for column density interpolation
  int Nhydro = snap.GetNparticlesType(typepart);                // No. of SPH particles in snap
  float dx,dy,invdx,invdy;             // ..
  float dr[3];                         // Rel. position vector on grid plane
  float drsqd;                         // Distance squared on grid plane
  float drmag;                         // Distance
  float dummyfloat = 0.0f;             // Dummy float for function arguments
  float invh;                          // 1/h
  float wkern;                         // Kernel value
  float wnorm;                         // Kernel normalisation value
  float *xvalues;                      // Pointer to 'x' array
  float *yvalues;                      // Pointer to 'y' array
  float *zvalues;                      // Pointer to 'z' array
  float *rendervalues;                 // Pointer to rendered quantity array
  float *mvalues;                      // Pointer to mass array
  float *rhovalues;                    // Pointer to density array
  float *hvalues;                      // Pointer to smoothing length array
  float *rendernorm;                   // Normalisation array
  float *rgrid;                        // Grid positions
  string dummystring = "";             // Dummy string for function arguments


  // Check x and y strings are actual co-ordinate strings
  if ((xstring != "x" && xstring != "y" && xstring != "z") ||
      (ystring != "x" && ystring != "y" && ystring != "z")) return -1;

  // First, verify x, y, z, m, rho, h and render strings are valid
  snap.ExtractArray(xstring, typepart, &xvalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray(ystring, typepart, &yvalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray(zstring, typepart, &zvalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray(renderstring, typepart, &rendervalues, &idummy, scaling_factor, renderunit);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray("m", typepart, &mvalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray("rho", typepart, &rhovalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);
  snap.ExtractArray("h", typepart, &hvalues, &idummy, dummyfloat, dummystring);
  arraycheck = min(idummy, arraycheck);

  // If any are invalid, exit here with failure code
  if (arraycheck == 0) return -1;

  rendernorm = new float[Ngrid];
  rgrid = new float[2*Ngrid];

  // Create grid positions here (need to improve in the future)
  dx = (xmax - xmin) / (float) ixgrid;
  dy = (ymax - ymin) / (float) iygrid;
  invdx = 1.0f/dx;
  invdy = 1.0f/dy;
  c = 0;
  for (j=iygrid-1; j>=0; j--) {
    for (i=0; i<ixgrid; i++) {
      rgrid[2*c] = xmin + ((float) i + 0.5f)*dx;
      rgrid[2*c + 1] = ymin + ((float) j + 0.5f)*dy;
      c++;
    }
  }


  // Zero arrays before computing rendering
  for (c=0; c<Ngrid; c++) values[c] = 0.0f;
  for (c=0; c<Ngrid; c++) rendernorm[c] = 0.0f;


  // Loop over all particles in snapshot
  //=================================================================================================
#pragma omp parallel for default(none) private(c,dr,drmag,drsqd,i,ii,imax,imin)\
   private(invh,jj,jmax,jmin,wkern,wnorm) shared(cout,dx,dy,invdx,invdy,hvalues,mvalues)\
   shared(Nhydro,rendernorm,rendervalues,rhovalues,rgrid,values,xvalues,yvalues,zvalues)
  for (i=0; i<Nhydro; i++) {
    const float hrange = hydro->kerntab.kernrange*hvalues[i];
    //const float hrangesqd = hydro->kerntab.kernrangesqd*hvalues[i]*hvalues[i];

    if (xvalues[i] + hrange < xmin || xvalues[i] - hrange > xmax ||
        yvalues[i] + hrange < ymin || yvalues[i] - hrange > ymax ||
        zvalues[i] + hrange < zslice || zvalues[i] - hrange > zslice) continue;

    // Compute grid coordinate limits for computing kernel convolution
    imin = (xvalues[i] - hrange - xmin)*invdx;  imin = max(0,imin);  imin = min(ixgrid-1,imin);
    imax = (xvalues[i] + hrange - xmin)*invdx;  imax = max(0,imax);  imax = min(ixgrid-1,imax);
    jmin = (yvalues[i] - hrange - ymin)*invdy;  jmin = max(0,jmin);  jmin = min(iygrid-1,jmin);
    jmax = (yvalues[i] + hrange - ymin)*invdy;  jmax = max(0,jmax);  jmax = min(iygrid-1,jmax);

    // If kernel does not overlap rendered image, skip particle entirely
    //if (imin == imax || jmin == jmax) continue;
    invh = 1.0f/hvalues[i];
    wnorm = mvalues[i]/rhovalues[i]*pow(invh,ndim);

    // Now loop over all pixels for particle
    //---------------------------------------------------------------------------------------------
    for (jj=jmin; jj<=jmax; jj++) {
      for (ii=imin; ii<=imax; ii++) {
        c = ii + (iygrid - jj - 1)*ixgrid;
        dr[0] = xmin + dx*(float) ii - xvalues[i];
        dr[1] = ymin + dy*(float) jj - yvalues[i];
        dr[2] = zslice - zvalues[i];
        drsqd = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        //if (drsqd > hrangesqd) continue;
        drmag = sqrt(drsqd);
        wkern = float(hydro->kerntab.w0((FLOAT) (drmag*invh)));
#pragma omp atomic
        values[c] += wnorm*rendervalues[i]*wkern;
#pragma omp atomic
        rendernorm[c] += wnorm*wkern;
      }
    }
   //----------------------------------------------------------------------------------------------

  }
  //===============================================================================================


  // Normalise all grid cells
  for (c=0; c<Ngrid; c++) {
    if (rendernorm[c] > 1.e-10) values[c] /= rendernorm[c];
  }

  // Free all locally allocated memory
  delete[] rgrid;
  delete[] rendernorm;

  return 1;
}



template class Render<1>;
template class Render<2>;
template class Render<3>;
