//=================================================================================================
//  PlaneParallelRadiation.cpp
//  ...
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


#include "Radiation.h"
#include "chealpix.h"


//=================================================================================================
//  PlaneParallelRadiation::PlaneParallelRadiation
/// Constructor for main PlaneParallelRadiation object
//=================================================================================================
template <int ndim, template<int> class ParticleType>
PlaneParallelRadiation<ndim,ParticleType>::PlaneParallelRadiation
 (Parameters *params, SmoothingKernel<ndim> *_kern, SimUnits *_units,
  NeighbourSearch<ndim> *_neib, CodeTiming *_timing):
  Radiation<ndim>()
{
  debug2("[PlaneParallelRadiation::PlaneParallelRadiation]");

  kern               = _kern;
  neib               = _neib;
  timing             = _timing;
  units              = _units;
  minRayDivisions    = params->intparams["min_ray_divisions"];
  radiationDirection = params->intparams["radiation_direction"];
  NLyC               = params->floatparams["NLyC"];
  arecomb            = params->floatparams["arecomb"];
  FLOAT gammam1      = params->floatparams["gamma_eos"] - 1.0;
  FLOAT mu_ion       = params->floatparams["mu_ion"];
  FLOAT temp_ion     = params->floatparams["temp_ion"];
  uion               = temp_ion/gammam1/mu_ion;
  maxIntegral        = NLyC / arecomb;

  tree = static_cast<OctTree<ndim,ParticleType,OctTreeCell>* > (neib->GetTree());
}



//=================================================================================================
//  PlaneParallelRadiation::~PlaneParallelRadiation
/// Destructor for PlaneParallelRadiation object
//=================================================================================================
template <int ndim, template<int> class ParticleType>
PlaneParallelRadiation<ndim,ParticleType>::~PlaneParallelRadiation()
{
}



//=================================================================================================
//  PlaneParallelRadiation::UpdateRadiationField
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void PlaneParallelRadiation<ndim,ParticleType>::UpdateRadiationField
 (int Nhydro,                                    ///< [in] No. of hydro particles
  int Nnbody,                                    ///< [in] No. of N-body particles
  int Nsink,                                     ///< [in] No. of sink particles
  Particle<ndim> *part_gen,                      ///< [in] Generic hydro particle data array
  NbodyParticle<ndim> **nbodydata,               ///< [in] N-body data array
  SinkParticle<ndim> *sinkdata)                  ///< [in] Sink data array
{
  int Nraysmax = 2*Nhydro;
  int Nrays = 0;
  FLOAT dir[ndim];
  PlanarRay<ndim> *rayStack = new PlanarRay<ndim>[Nraysmax];
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* >(part_gen);

  debug2("[PlaneParallelRadiation::UpdateRadiationField]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("RADIATION_PLANAR");


  // Set all particles to neutral before computing radiation properties
  for (int i=0; i<Nhydro; i++) partdata[i].ionstate = 0;
  numIonised = 0;

  // Start by creating rays from the root cell
  const OctTreeCell<ndim> &rootCell = tree->celldata[0];
  const FLOAT rootCellSize = tree->rootCellSize;
  const int ltotTree = tree->ltot;
  const FLOAT leafCellSize = rootCellSize / pow(2, ltotTree);
  PlanarRay<ndim> &rootRay = rayStack[Nrays++];
  CreateRootRay(rootCell, rootCellSize, rootRay);

  // Add small offset to ensure ray in inside computational domain
  for (int k=0; k<ndim; k++) rootRay.r[k] += 0.00001*leafCellSize*dir[k];

  //std::cout << "RAY : " << rootRay.dir[0] << "  " << radiationDirection << "   r : " << rootRay.r[0] << std::endl;

  // Evaluate all rays on the stack until all are terminated (e.g. by reaching the maximim
  // ionisation integral value) or leave the computational domain.
  //-----------------------------------------------------------------------------------------------
  while (Nrays > 0) {
    PlanarRay<ndim> ray = rayStack[--Nrays];
    bool terminateRay = false;

    // Integrate current ray while (i) the integral is less than the limit which defines the
    // location of the ionisation front, (ii) the ray is still inside the computational domain,
    // and (iii) the ray resolution is sufficient (i.e. no larger than the cell size).
    //---------------------------------------------------------------------------------------------
    do {

      const int cc = tree->FindLeafCell(ray.r);
      const OctTreeCell<ndim> &cell = tree->celldata[cc];

      //std::cout << "Investigating cell : " << cc << "  " << cell.rcell[0] << "  " <<  cell.cellBox.min[0] << "  " << cell.cellBox.max[0] << "  N : " << cell.N << std::endl;
      //std::cout << "level : " << cell.level << std::endl;

      // If leaf cell level is deeper than ray (i.e. so cell/particle size is much smaller than
      // the ray beam width), then split ray into 2^(ndim - 1) child rays
      if (cell.level > ray.level) {
        const FLOAT rayWidth = rootCellSize / pow(2, ray.level);
        SplitRay(ray, rayWidth, Nrays, rayStack);
        terminateRay = true;
      }
      else {
        const FLOAT cellSize = rootCellSize / pow(2, cell.level);
        terminateRay = CellRayIntegration(cell, cellSize, ray, partdata);
        //std::cout << "ray;  level : " << ray.level << "   rayIntegral : " << ray.rayIntegral << "  " << maxIntegral << std::endl;
      }

    } while (!terminateRay && ray.rayIntegral < maxIntegral && PointInBox(ray.r, rootCell.cellBox));
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  //int k;
  //cin >> k;

  // Find initial position of ray based on the size of the tree root cell
  /*for (int k=0; k<ndim; k++) ray.r[k] = (FLOAT) 0.5;
  ray.r[0] = xmin;

  // Initial ray points in positive x-direction
  for (int k=0; k<ndim; k++) dir[k] = (FLOAT) 0.0;
  dir[0] = (FLOAT) 1.0;

  // For now, use minimum smoothing length of particles as constant step-size
  FLOAT hmin = big_number;
  FLOAT hmax = big_number;
  for (int i=0; i<Nhydro; i++) {
    hmin = min(hmin, partdata[i].h);
    hmax = max(hmax, partdata[i].h);
  }
  FLOAT step = (FLOAT) 0.2*hmin;
  FLOAT invhmaxsqd = (FLOAT) 1.0 / hmax / hmax;
  FLOAT hfactor = pow((FLOAT) 1.0 / hmax, ndim);

  // Compute density at initial place
  FLOAT rho0 = (FLOAT) 0.0;
  for (int i=0; i<Nhydro; i++) {
    FLOAT dr[ndim];
    for (int k=0; k<ndim; k++) dr[k] = ray.r[k] - partdata[i].r[k];
    FLOAT drsqd = DotProduct(dr, dr, ndim);
    FLOAT ssqd  = drsqd*invhmaxsqd;
    rho0       += partdata[i].m*kern->w0_s2(ssqd);
  }
  rho0 *= hfactor;

  // Ray-trace through the computational domain computing the ionisation integral
  //-----------------------------------------------------------------------------------------------
  do {

    for (int k=0; k<ndim; k++) ray.r[k] += step*dir[k];

    // Compute density at new point
    FLOAT rho1 = (FLOAT) 0.0;
    for (int i=0; i<Nhydro; i++) {
      FLOAT dr[ndim];
      for (int k=0; k<ndim; k++) dr[k] = ray.r[k] - partdata[i].r[k];
      FLOAT drsqd = DotProduct(dr, dr, ndim);
      FLOAT ssqd  = drsqd*invhmaxsqd;
      rho1       += partdata[i].m*kern->w0_s2(ssqd);
    }
    rho1 *= hfactor;

    FLOAT dIntegral = (FLOAT) 0.5*(rho0*rho0 + rho1*rho1)*step;
    ray.rayIntegral += dIntegral;
    numSteps++;
    rho0 = rho1;

  } while (ray.rayIntegral < maxIntegral && ray.r[0] < xmax);
  //-----------------------------------------------------------------------------------------------

  // Set ionisation fractions of particles based on position relative to ionisation front
  for (int i=0; i<Nhydro; i++) {
    if (partdata[i].r[0] > ray.r[0]) {
      partdata[i].ionstate = 0;
    }
    else {
      partdata[i].ionstate = 1;
      partdata[i].u = uion;
      numIonised++;
      //std::cout << "FOUND IONISED PARTICLE : " << i << "  " << partdata[i].r[0] << std::endl;
    }
  }

  std::cout << "NO. OF STEPS : " << numSteps << "  " << ray.rayIntegral << "   " << maxIntegral << std::endl;
  std::cout << "NO. IONISED  : " << numIonised << "    r : " << ray.r[0] << std::endl;
  */

  delete[] rayStack;

  return;
}



//=================================================================================================
//  PlaneParallelRadiation::CellRayIntegration
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void PlaneParallelRadiation<ndim,ParticleType>::CreateRootRay
 (const OctTreeCell<ndim> &rootCell,             ///< ..
  const FLOAT rootCellSize,                      ///< ..
  PlanarRay<ndim> &rootRay)                      ///< ..
{
  for (int k=0; k<ndim; k++) rootRay.r[k] = rootCell.rcell[k];
  rootRay.level = 0;

  // +ve x-direction
  if (radiationDirection == 0) {
    rootRay.r[0] = rootCell.cellBox.min[0];
    rootRay.dir[0] = (FLOAT) 1.0;
    kx = 0;
    ky = 1;
    kz = 2;
  }
  // -ve x-direction
  else if (radiationDirection == 1) {
    rootRay.r[0] = rootCell.cellBox.max[0];
    rootRay.dir[0] = -(FLOAT) 1.0;
    kx = 0;
    ky = 1;
    kz = 2;
  }
  // +ve y-direction
  else if (ndim > 1 && radiationDirection == 2) {
    rootRay.r[1] = rootCell.cellBox.min[1];
    rootRay.dir[1] = (FLOAT) 1.0;
    kx = 1;
    ky = 0;
    kz = 2;
  }
  // -ve y-direction
  else if (ndim > 1 && radiationDirection == 3) {
    rootRay.r[1] = rootCell.cellBox.max[1];
    rootRay.dir[1] = -(FLOAT) 1.0;
    kx = 1;
    ky = 0;
    kz = 2;
  }
  // +ve z-direction
  else if (ndim == 3 && radiationDirection == 4) {
    rootRay.r[2] = rootCell.cellBox.min[2];
    rootRay.dir[2] = (FLOAT) 1.0;
    kx = 2;
    ky = 0;
    kz = 1;
  }
  // -ve x-direction
  else if (ndim == 3 && radiationDirection == 5) {
    rootRay.r[2] = rootCell.cellBox.max[2];
    rootRay.dir[2] = -(FLOAT) 1.0;
    kx = 2;
    ky = 0;
    kz = 1;
  }
  return;
}



//=================================================================================================
//  PlaneParallelRadiation::CellRayIntegration
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
bool PlaneParallelRadiation<ndim,ParticleType>::CellRayIntegration
 (const OctTreeCell<ndim> &cell,                 ///< ..
  const FLOAT cellSize,                          ///< ..
  PlanarRay<ndim> &ray,                          ///< ..
  ParticleType<ndim> *partdata)                  ///< ..
{
  // If leaf-cell is not empty, then compute mean density and compute integral contribution
  if (cell.N > 0) {
    const FLOAT rho = cell.m / pow(cellSize, ndim);
    const FLOAT dIntegral = (FLOAT) 0.5*rho*cellSize;
    if (ray.rayIntegral + dIntegral > maxIntegral) {
      const FLOAT frac = (maxIntegral - ray.rayIntegral) / dIntegral;
      for (int k=0; k<ndim; k++) ray.r[k] += frac*cellSize*ray.dir[k];
      int i = cell.ifirst;
      while (i != -1) {
        FLOAT dr[ndim];
        for (int k=0; k<ndim; k++) dr[k] = ray.r[k] - partdata[i].r[k];
        const FLOAT dot = DotProduct(dr, ray.dir, ndim);
        if (dot > (FLOAT) 0.0) {
          partdata[i].ionstate = 1;
          partdata[i].u = uion;
          //std::cout << "IONISING : " << i << "   r : " << partdata[i].r[0] << std::endl;
          numIonised++;
        }
        if (i == cell.ilast) break;
        i = tree->inext[i];
      };
      //std::cout << "FRONT : " << ray.r[0] << "    numIonised : " << numIonised << std::endl;
      return true;
    }
    else {
      ray.rayIntegral += dIntegral;
      for (int k=0; k<ndim; k++) ray.r[k] += cellSize*ray.dir[k];
      //std::cout << "MOVING RAY : " << ray.r[0] << "   dir : " << ray.dir[0] << std::endl;
      int i = cell.ifirst;
      while (i != -1) {
        partdata[i].ionstate = 1;
        partdata[i].u = uion;
        numIonised++;
        //std::cout << "IONISING : " << i << "   r : " << partdata[i].r[0] << std::endl;
        if (i == cell.ilast) break;
        i = tree->inext[i];
      };
      return false;
    }
  }
  else {
    for (int k=0; k<ndim; k++) ray.r[k] += cellSize*ray.dir[k];
    return false;
  }
}



//=================================================================================================
//  PlaneParallelRadiation::SplitRay
/// ...
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void PlaneParallelRadiation<ndim,ParticleType>::SplitRay
 (const PlanarRay<ndim> &ray,                    ///< [in] ..
  const FLOAT rayWidth,                          ///< [in] ..
  int &Nrays,                                    ///< [inout] ..
  PlanarRay<ndim> *rayStack)                     ///< [inout] ..
{
  // In 1d, trivial replacement with identical ray at different level
  if (ndim == 1) {
    PlanarRay<ndim> &childRay = rayStack[Nrays++];
    childRay = ray;
    childRay.level++;
  }
  // In 2d, split ray into 2 child rays separated by cell size
  else if (ndim == 2) {
    //std::cout << "Ray pos : " << ray.r[0] << "  " << ray.r[1] << std::endl;
    PlanarRay<ndim> &childRay1 = rayStack[Nrays++];
    childRay1 = ray;
    childRay1.level++;
    //std::cout << "Child1 pos : " << childRay1.r[0] << "  " << childRay1.r[1] << std::endl;
    childRay1.r[ky] += (FLOAT) 0.25*rayWidth;
    //std::cout << "Ray pos (again) : " << ray.r[0] << "  " << ray.r[1] << std::endl;
    PlanarRay<ndim> &childRay2 = rayStack[Nrays++];
    childRay2 = ray;
    childRay2.level++;
    //std::cout << "Child2 pos : " << childRay2.r[0] << "  " << childRay2.r[1] << std::endl;
    childRay2.r[ky] -= (FLOAT) 0.25*rayWidth;
    //std::cout << "Splitting ray; level : " << ray.level << "  rayWidth : " << rayWidth << "   ry : " << ray.r[ky] << "   child r : " << childRay1.r[ky] << "  " << childRay2.r[ky] << std::endl;
    //int k;
    //cin >> k;
  }
  else if (ndim == 3) {
    PlanarRay<ndim> &childRay1 = rayStack[Nrays++];
    childRay1 = ray;
    childRay1.level++;
    childRay1.r[ky] += (FLOAT) 0.25*rayWidth;
    childRay1.r[kz] += (FLOAT) 0.25*rayWidth;

    PlanarRay<ndim> &childRay2 = rayStack[Nrays++];
    childRay2 = ray;
    childRay2.level++;
    childRay2.r[ky] -= (FLOAT) 0.25*rayWidth;
    childRay2.r[kz] += (FLOAT) 0.25*rayWidth;

    PlanarRay<ndim> &childRay3 = rayStack[Nrays++];
    childRay3 = ray;
    childRay3.level++;
    childRay3.r[ky] += (FLOAT) 0.25*rayWidth;
    childRay3.r[kz] -= (FLOAT) 0.25*rayWidth;

    PlanarRay<ndim> &childRay4 = rayStack[Nrays++];
    childRay4 = ray;
    childRay4.level++;
    childRay4.r[ky] -= (FLOAT) 0.25*rayWidth;
    childRay4.r[kz] -= (FLOAT) 0.25*rayWidth;
  }

  return;
}



template class PlaneParallelRadiation<3, GradhSphParticle>;
template class PlaneParallelRadiation<2, GradhSphParticle>;
template class PlaneParallelRadiation<1, GradhSphParticle>;

template class PlaneParallelRadiation<3, MeshlessFVParticle>;
template class PlaneParallelRadiation<2, MeshlessFVParticle>;
template class PlaneParallelRadiation<1, MeshlessFVParticle>;
