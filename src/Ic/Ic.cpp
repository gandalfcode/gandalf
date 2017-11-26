//=================================================================================================
//  IC.cpp
//  Implementation of the initial conditions generation
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


#include <fstream>
#include <sstream>
#include "Precision.h"
#include "Debug.h"
#include "Ic.h"
using namespace std;



//=================================================================================================
//  Ic::CalculateMassInBox
/// Calculates total mass throughout volume of given box.
//=================================================================================================
template <int ndim>
FLOAT Ic<ndim>::CalculateMassInBox
 (const Box<ndim> &box,
  const int ptype)
{
  int maxGridSize = 512;
  int grid[ndim];
  FLOAT dr[ndim];
  FLOAT r[ndim];
  FLOAT mtot      = (FLOAT) 0.0;
  FLOAT dV        = (FLOAT) 1.0;
  FLOAT lengthmax = (FLOAT) 0.0;

  // Calculate the maximum box length and then the appropriate grid spacing in each direction
  for (int k=0; k<ndim; k++) {
    lengthmax = max(lengthmax, (box.max[k] - box.min[k]));
  }
  for (int k=0; k<ndim; k++) {
    grid[k] = (int) ((FLOAT) maxGridSize*(box.max[k] - box.min[k])/lengthmax);
    dr[k]   = (box.max[k] - box.min[k]) / (FLOAT) grid[k];
    dV      *= dr[k];
  }

  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    for (int i=0; i<grid[0]; i++) {
      r[0] = box.min[0] + ((FLOAT) i + (FLOAT) 0.5)*dr[0];
      mtot += this->GetDensity(r, ptype);
    }
  }
  //-----------------------------------------------------------------------------------------------
  if (ndim == 2) {
    for (int i=0; i<grid[0]; i++) {
      r[0] = box.min[0] + ((FLOAT) i + (FLOAT) 0.5)*dr[0];
      for (int j=0; j<grid[1]; j++) {
        r[1] = box.min[1] + ((FLOAT) j + (FLOAT) 0.5)*dr[1];
        mtot += this->GetDensity(r, ptype);
      }
    }
  }
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {
    for (int i=0; i<grid[0]; i++) {
      r[0] = box.min[0] + ((FLOAT) i + (FLOAT) 0.5)*dr[0];
      for (int j=0; j<grid[1]; j++) {
        r[1] = box.min[1] + ((FLOAT) j + (FLOAT) 0.5)*dr[1];
        for (int k=0; k<grid[2]; k++) {
          r[2] = box.min[2] + ((FLOAT) k + (FLOAT) 0.5)*dr[2];
          mtot += this->GetDensity(r, ptype);
        }
      }
    }
  }
  //-----------------------------------------------------------------------------------------------

  mtot *= dV;

  return mtot;
}



//=================================================================================================
//  Ic::GetSmoothedDensity
/// Calculates the smoothed value of the required quantity 'var' at the position 'rsmooth' for
/// a smoothing kernel 'kern' and smoothing length 'h'.  Numerically calculates the smoothed
/// integral by discretising over the smoothing kernel extent with a grid of size 'gridsize' in
/// each dimension.  Returns the final smoothed quantity.
//=================================================================================================
template <int ndim>
FLOAT Ic<ndim>::GetSmoothedDensity
 (const FLOAT rsmooth[ndim],
  const int ptype,
  const FLOAT h,
  SmoothingKernel<ndim> *kern)
{
  const int gridSize  = 3;
  const FLOAT hrange  = kern->kernrange*h;
  const FLOAT invh    = (FLOAT) 1.0/h;
  const FLOAT invhsqd = invh*invh;
  const FLOAT hfactor = pow(invh, ndim);
  const FLOAT dV      = pow(2.0*hrange/(FLOAT) gridSize, ndim);
  FLOAT sumValue      = (FLOAT) 0.0;
  FLOAT dr[ndim];
  FLOAT drsqd;
  FLOAT r[ndim];

  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    for (int i=0; i<gridSize; i++) {
      dr[0] = ((FLOAT) i + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
      drsqd = dr[0]*dr[0];
      r[0]  = rsmooth[0] + dr[0];
      sumValue += hfactor*kern->w0_s2(drsqd*invhsqd)*this->GetDensity(r, ptype);
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 2) {
    for (int i=0; i<gridSize; i++) {
      for (int j=0; j<gridSize; j++) {
        dr[0] = ((FLOAT) i + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
        dr[1] = ((FLOAT) j + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
        r[0]  = rsmooth[0] + dr[0];
        r[1]  = rsmooth[1] + dr[1];
        drsqd = dr[0]*dr[0] + dr[1]*dr[1];
        sumValue += hfactor*kern->w0_s2(drsqd*invhsqd)*this->GetDensity(r, ptype);
      }
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 3) {
    for (int i=0; i<gridSize; i++) {
      for (int j=0; j<gridSize; j++) {
        for (int k=0; k<gridSize; k++) {
          dr[0] = ((FLOAT) i + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
          dr[1] = ((FLOAT) j + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
          dr[2] = ((FLOAT) k + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
          r[0]  = rsmooth[0] + dr[0];
          r[1]  = rsmooth[1] + dr[1];
          r[2]  = rsmooth[2] + dr[2];
          drsqd = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
          sumValue += hfactor*kern->w0_s2(drsqd*invhsqd)*this->GetDensity(r, ptype);
        }
      }
    }
  }
  //-----------------------------------------------------------------------------------------------

  // Normalise summation/integral with the volume of each individual summation point
  sumValue *= dV;

  return sumValue;
}



//=================================================================================================
//  Ic::CheckInitialConditions
/// Performs some simple sanity checks on all initial conditions
//=================================================================================================
template <int ndim>
void Ic<ndim>::CheckInitialConditions(void)
{
  bool valid_ic = true;                      // Valid initial conditions flag
  int i,k;                                   // Particle and dimension counter
  DomainBox<ndim>& icBox = sim->icBox;
  const bool cutbox       = simparams->intparams["cut_box"];

  debug2("[Ic::CheckInitialConditions]");

  // Loop through all particles performing various checks
  //-----------------------------------------------------------------------------------------------
  if (hydro != NULL) {
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      bool boxflag = true;
      bool okflag  = true;

      // Check that main particle properties are valid floating point numbers
      if (!isnormal(part.m)) okflag = false;
      if (!isfinite(part.u)) okflag = false;
      for (k=0; k<ndim; k++) {
        if (!isfinite(part.r[k])) okflag = false;
        if (!isfinite(part.v[k])) okflag = false;
      }

      // Check that all particles reside inside any defined boundaries
      for (k=0; k<ndim; k++) {
        if (part.r[k] < icBox.min[k]) {
          if (icBox.boundary_lhs[k] == periodicBoundary || icBox.boundary_lhs[k] == mirrorBoundary) {
            if (cutbox) {
              part.flags.set(dead);
            }
            else {
              boxflag = false;
            }
          }
        }
        if (part.r[k] > icBox.max[k]) {
          if (icBox.boundary_rhs[k] == periodicBoundary || icBox.boundary_lhs[k] == mirrorBoundary) {
            if (cutbox) {
              part.flags.set(dead);
            }
            else {
              boxflag = false;
            }
          }
        }
      }

      // If flag indicates a problem, print error and quit
      if (!okflag) {
        cout << "Invalid particle : " << i << "   m : " << part.m << "   u : " << part.u << "    r : ";
        for (k=0; k<ndim; k++) cout << part.r[k] << "  ";
        cout << "    v : ";
        for (k=0; k<ndim; k++) cout << part.v[k] << "  ";
        cout << std::endl;
        ExceptionHandler::getIstance().raise("Error : Invalid floating point values for particles");
      }
      if (!boxflag) {
        cout << "Particle " << i << " not inside periodic box" << endl;
        for (k=0; k<ndim; k++)
          cout << "r[" << k << "] : " << part.r[k] << "    "
          << icBox.min[k] << "    " << icBox.max[k] << endl;
      }

      valid_ic = valid_ic && okflag;

    }

    // Check particles are sorted in type order
    int ptype =  -1 ;
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      if (part.ptype < ptype) {
        ExceptionHandler::getIstance().raise("Error Particles must be ordered by ptype");
      }
      ptype = part.ptype;
    }
  }
  //-----------------------------------------------------------------------------------------------

  if (!valid_ic) {
    string message = "Invalid initial conditions for SPH particles";
    ExceptionHandler::getIstance().raise(message);
  }

  return;
}



//=================================================================================================
//  Ic::AddBinaryStar
/// Add a binary star of given mass, eccentricity and separation.
/// (Code provided courtesy of S. P. Goodwin; 29/09/2013)
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddBinaryStar
 (FLOAT sma,                          ///< Semi-major axis
  FLOAT eccent,                       ///< Orbital eccentricity
  FLOAT m1,                           ///< Mass of star 1
  FLOAT m2,                           ///< Mass of star 2
  FLOAT h1,                           ///< Smoothing length of star 1
  FLOAT h2,                           ///< Smoothing length of star 2
  FLOAT phirot,                       ///< 'phi' Euler rotation angle
  FLOAT thetarot,                     ///< 'theta' Euler rotation angle
  FLOAT phase,                        ///< Phase angle
  FLOAT psirot,                       ///< 'tpsi' rotation angle
  FLOAT *rbinary,                     ///< Position of COM of binary
  FLOAT *vbinary,                     ///< Velocity of COM of binary
  NbodyParticle<ndim> &s1,             ///< Star 1
  NbodyParticle<ndim> &s2)             ///< Star 2
{

  //-----------------------------------------------------------------------------------------------
  if (ndim > 1) {
    debug2("[Ic::AddBinaryStar]");

    int k;                               // Dimension counter
    FLOAT mbin = m1 + m2;                // Total binary mass

    // randomly sample M to get theta
    FLOAT M = 2.0*pi*sim->randnumb->floatrand();

    // from this solve to get eccentric anomoly E - e sin(theta) = M
    // N-R method x_1 = x_0 - f(x_0)/f'(x_0)
    FLOAT Ee = M;
    FLOAT old = M;
    do {
      old = Ee;
      Ee = Ee - (Ee - eccent*sin(Ee) - M)/(1.0 - eccent*cos(Ee));
    } while (fabs(old - Ee) > (FLOAT) 1.0e-6);

    // Next get theta tan(t/2) = sqrt((1+e)/(1-e))*tan(E/2)
    FLOAT theta = sqrt(((FLOAT) 1.0 + eccent)/((FLOAT) 1.0 - eccent))*tan((FLOAT) 0.5*Ee);
    theta = (FLOAT) 2.0*atan(theta);

    // Total separation
    FLOAT sep = sma*((FLOAT) 1.0 - eccent*eccent)/((FLOAT) 1.0 + eccent*cos(theta));

    // Get velocity
    FLOAT vel = (m1 + m2)*((FLOAT) 2.0/sep - (FLOAT) 1.0/sma);
    vel = sqrt(vel);

    // ..
    FLOAT hc = sqrt(((FLOAT) 1.0 + eccent*cos(theta))/((FLOAT) 2.0 - sep/sma));
    FLOAT phi = acos(hc);

    // Set properties of star 1
    // put on x-y plane to start
    for (k=0; k<ndim; k++) s1.r[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) s1.v[k] = (FLOAT) 0.0;
    s1.m = m1;
    s1.h = h1;
    s1.invh = (FLOAT) 1.0 / s1.h;
    s1.r[0] += sep*cos(theta)*m2/mbin;
    s1.r[1] += sep*sin(theta)*m2/mbin;
    s1.v[0] += -vel*cos((FLOAT) 0.5*pi - theta + phi)*m2/mbin;
    s1.v[1] += vel*sin((FLOAT) 0.5*pi - theta + phi)*m2/mbin;

    // Set properties of star 2
    for (k=0; k<ndim; k++) s2.r[k] = (FLOAT) 0.0;
    for (k=0; k<ndim; k++) s2.v[k] = (FLOAT) 0.0;
    s2.m = m2;
    s2.h = h2;
    s2.invh = (FLOAT) 1.0 / s2.h;
    s2.r[0] -= sep*cos(theta)*m1/mbin;
    s2.r[1] -= sep*sin(theta)*m1/mbin;
    s2.v[0] -= -vel*cos((FLOAT) 0.5*pi - theta + phi)*m1/mbin;
    s2.v[1] -= vel*sin((FLOAT) 0.5*pi - theta + phi)*m1/mbin;

    // Rotate binary to given orientation using Euler angles
    EulerAngleRotation(phirot,thetarot,psirot,s1.r);
    EulerAngleRotation(phirot,thetarot,psirot,s1.v);
    EulerAngleRotation(phirot,thetarot,psirot,s2.r);
    EulerAngleRotation(phirot,thetarot,psirot,s2.v);

    // Now move binary to given centre of position and velocity
    for (k=0; k<ndim; k++) s1.r[k] += rbinary[k];
    for (k=0; k<ndim; k++) s1.v[k] += vbinary[k];
    for (k=0; k<ndim; k++) s2.r[k] += rbinary[k];
    for (k=0; k<ndim; k++) s2.v[k] += vbinary[k];

  }
  //-----------------------------------------------------------------------------------------------
  else {
    string message = "Binary test not available in 1D";
    ExceptionHandler::getIstance().raise(message);
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::GetMaximumDensityField
/// Populate given bounding box with random particles.
//=================================================================================================
template <int ndim>
FLOAT Ic<ndim>::GetMaximumDensity
 (const int ptype,                     ///< [in] Particle type for max. density field
  const Box<ndim> &box,                ///< [in] Bounding box containing particles
  RandomNumber *randnumb)              ///< [inout] Pointer to random number generator
{
  int numSamples = 1000000;
  FLOAT rhoMax = (FLOAT) 0.0;
  FLOAT rrand[ndim];

  for (int i=0; i<numSamples; i++) {
    for (int k=0; k<ndim; k++) {
      rrand[k] = box.min[k] + (box.max[k] - box.min[k])*randnumb->floatrand();
    }
    rhoMax = max(rhoMax, GetDensity(rrand, ptype));
  }

  return rhoMax;
}



//=================================================================================================
//  Ic::AddMonteCarloDensityField
/// Populate given bounding box with random particles.
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddMonteCarloDensityField
 (const int Npart,                     ///< [in] No. of particles
  const int ptype,                     ///< [in] Type of particle
  const Box<ndim> &box,                ///< [in] Bounding box containing particles
  FLOAT *r,                            ///< [out] Positions of particles
  RandomNumber *randnumb)              ///< [inout] Pointer to random number generator
{
  FLOAT rho;
  FLOAT rhoMax = GetMaximumDensity(ptype, box, randnumb);
  FLOAT rrand[ndim];

  debug2("[Ic::AddMonteCarloDensityField]");
  assert(r);

  //-----------------------------------------------------------------------------------------------
  for (int i=0; i<Npart; i++) {

    // Iterate with random Monte-Carlo sampling of the density field
    do {
      for (int k=0; k<ndim; k++) {
        rrand[k] = box.min[k] + (box.max[k] - box.min[k])*randnumb->floatrand();
      }
      rho = rhoMax*randnumb->floatrand();
    } while (GetDensity(rrand, ptype) < rho);

    // Copy "correct" particle position to array for IC generation
    for (int k=0; k<ndim; k++) r[ndim*i + k] = rrand[k];
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::AddRandomBox
/// Populate given bounding box with random particles.
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddRandomBox
 (const int Npart,                     ///< [in] No. of particles
  const DomainBox<ndim> &box,          ///< [in] Bounding box containing particles
  FLOAT *r,                            ///< [out] Positions of particles
  RandomNumber *randnumb)              ///< [inout] Pointer to random number generator
{
  debug2("[Ic::AddRandomBox]");
  assert(r);

  for (int i=0; i<Npart; i++) {
    for (int k=0; k<ndim; k++) {
      r[ndim*i + k] = box.min[k] + (box.max[k] - box.min[k])*randnumb->floatrand();
    }
  }

  return;
}



//=================================================================================================
//  Ic::AddRandomSphere
/// Add random sphere of particles
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddRandomSphere
 (const int Npart,                     ///< [in] No. of particles in sphere
  const FLOAT rcentre[ndim],           ///< [in] Position of sphere centre
  const FLOAT radius,                  ///< [in] Radius of sphere
  FLOAT *r,                            ///< [out] Positions of particles in sphere
  RandomNumber *randnumb)              ///< [inout] Pointer to random number generator
{
  int i,k;                             // Particle and dimension counters
  FLOAT rad;                           // Radius of particle
  FLOAT rpos[ndim];                    // Random position of new particle

  debug2("[Ic::AddRandomSphere]");
  assert(r);

  // Loop over all required particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {

    // Continously loop until random particle lies inside sphere
    do {
      for (k=0; k<ndim; k++)
      rpos[k] = (FLOAT) 1.0 - (FLOAT) 2.0*randnumb->floatrand();
      rad = DotProduct(rpos,rpos,ndim);
    } while (rad > 1.0);

    for (k=0; k<ndim; k++) r[ndim*i + k] = rcentre[k] + radius*rpos[k];
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::AddLatticeSphere
/// Add sphere of particles cut-out of regular lattice
//=================================================================================================
template <int ndim>
int Ic<ndim>::AddLatticeSphere
 (const int Npart,                     ///< [in] Desired/maximum no. of particles in sphere
  const FLOAT rcentre[ndim],           ///< [in] Position of sphere centre
  const FLOAT radius,                  ///< [in] Radius of sphere
  const string particle_dist,          ///< [in] String of lattice type
  FLOAT *r,                            ///< [out] Positions of particles in sphere
  RandomNumber *randnumb)              ///< [inout] Pointer to random number generator
{
  int i,k;                             // Particle and dimension counters
  int Naux;                            // Aux. particle number
  int Nlattice[3];                     // Lattice size
  FLOAT theta;                         // Euler angle for random rotation of lattice
  FLOAT phi;                           //        "                      "
  FLOAT psi;                           //        "                      "
  FLOAT *raux;                         // Temp. array to hold particle positions
  DomainBox<ndim> box1;                // Bounding box

  debug2("[Ic::AddLatticeSphere]");
  assert(r);


  // Set parameters for box and lattice to ensure it contains enough particles
  for (k=0; k<3; k++) Nlattice[k] = 1;
  for (k=0; k<ndim; k++) Nlattice[k] = (int) (3.0*powf((FLOAT) Npart, (FLOAT)1/ndim));
  for (k=0; k<ndim; k++) box1.min[k] = -2.0;
  for (k=0; k<ndim; k++) box1.max[k] = 2.0;
  Naux = Nlattice[0]*Nlattice[1]*Nlattice[2];
  raux = new FLOAT[ndim*Naux];

  // Create a bounding box to hold lattice sphere
  if (particle_dist == "cubic_lattice") {
    AddCubicLattice(Naux, Nlattice, box1, true, raux);
  }
  else if (particle_dist == "hexagonal_lattice") {
    AddHexagonalLattice(Naux, Nlattice, box1, true, raux);
  }
  else {
    string message = "Invalid particle distribution option";
    ExceptionHandler::getIstance().raise(message);
  }

  // Now cut-out sphere from lattice containing exact number of particles
  // (unless lattice structure prevents this).
  Naux = CutSphere(Npart, Naux, box1, false, raux);
  assert(Naux <= Npart);

  // Rotate sphere through random Euler angles (to prevent alignment problems
  // during tree construction)
  if (ndim == 2) {
    FLOAT rtemp[ndim];
    theta = twopi*randnumb->floatrand();
    for (i=0; i<Naux; i++) {
      for (k=0; k<ndim; k++) rtemp[k] = raux[ndim*i + k];
      raux[ndim*i] = rtemp[0]*cos(theta) - rtemp[1]*sin(theta);
      raux[ndim*i + 1] = rtemp[0]*sin(theta) + rtemp[1]*cos(theta);
    }
  }
  else if (ndim == 3) {
    theta = acos(sqrtf(randnumb->floatrand()));
    phi   = twopi*randnumb->floatrand();
    psi   = twopi*randnumb->floatrand();
    EulerAngleArrayRotation(Naux, phi, theta, psi, raux);
  }

  // Copy particle positions to main position array to be returned
  for (i=0; i<Naux; i++) {
    for (k=0; k<ndim; k++) r[ndim*i + k] = rcentre[k] + radius*raux[ndim*i + k];
  }

  // Free allocated memory
  delete[] raux;

  return Naux;
}



//=================================================================================================
//  Ic::Addr2Sphere
/// Add r^-2 sphere of particles
//=================================================================================================
template <int ndim>
void Ic<ndim>::Addr2Sphere
 (int Npart,                           ///< [in] No. of particles in sphere
  FLOAT *r,                            ///< [out] Positions of particles in sphere
  FLOAT rcentre[ndim],                 ///< [in] Position of sphere centre
  FLOAT radius)                        ///< [in] Radius of sphere
{
  int i;                               // Particle counter
  FLOAT phi;                           // ..
  FLOAT theta;                         // ..
  FLOAT sintheta;                      // ..
  FLOAT costheta;                      // ..
  FLOAT rpart;                         // ..

  debug2("[Ic::Addr2Sphere]");
  assert(r);

  // Loop over all required particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Npart; i++) {

    // Continously loop until random particle lies inside sphere
    phi      = (FLOAT) 2.0*pi*randnumb->floatrand();
    costheta = (FLOAT) 2.0*randnumb->floatrand() - (FLOAT) 1.0;
    theta    = acos(costheta);
    sintheta = sin(theta);
    rpart    = radius*randnumb->floatrand();
    r[ndim*i + 0] = rpart*sintheta*cos(phi);
    r[ndim*i + 1] = rpart*sintheta*sin(phi);
    r[ndim*i + 2] = rpart*costheta;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::AddCubicLattice
/// Add regular (cubic) lattice of particles
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddCubicLattice
 (const int Npart,                     ///< [in] No. of particles in lattice
  const int Nlattice[ndim],            ///< [in] Ptcls per dimension in lattice
  const DomainBox<ndim> &box,          ///< [in] Bounding box of particles
  const bool normalise,                ///< [in] Normalise lattice shape and size
  FLOAT *r)                            ///< [out] Positions of particles
{
  int i,k;                             // Particle and dimension counters
  int ii,jj,kk;                        // Aux. lattice counters
  FLOAT spacing[ndim];                 // Lattice spacing in each direction

  debug2("[Ic::AddCubicLattice]");
  assert(r);

  // If normalised, ensure equal spacing between all lattice layers.
  // Otherwise set spacing to fit bounding box
  if (normalise) {
    for (k=0; k<ndim; k++) {
      spacing[k] = (box.max[0] - box.min[0])/(FLOAT) Nlattice[0];
    }
  }
  else {
    for (k=0; k<ndim; k++) {
      spacing[k] = (box.max[k] - box.min[k])/(FLOAT) Nlattice[k];
    }
  }


  // Create lattice depending on dimensionality
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.min[0] + ((FLOAT)ii + (FLOAT) 0.5)*spacing[0];
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
        i = jj*Nlattice[0] + ii;
        r[ndim*i] = box.min[0] + ((FLOAT)ii + (FLOAT) 0.5)*spacing[0];
        r[ndim*i + 1] = box.min[1] + ((FLOAT)jj + (FLOAT) 0.5)*spacing[1];
      }
    }
  }
  //-----------------------------------------------------------------------------------------------
  else if (ndim == 3) {
#pragma omp parallel for default(none) shared(cout,box,Nlattice,r,spacing) private(i,ii,jj,kk)
    for (kk=0; kk<Nlattice[2]; kk++) {
      for (jj=0; jj<Nlattice[1]; jj++) {
        for (ii=0; ii<Nlattice[0]; ii++) {
          i = kk*Nlattice[0]*Nlattice[1] + jj*Nlattice[0] + ii;
          r[ndim*i] = box.min[0] + ((FLOAT)ii + (FLOAT) 0.5)*spacing[0];
          r[ndim*i + 1] = box.min[1] + ((FLOAT)jj + (FLOAT) 0.5)*spacing[1];
          r[ndim*i + 2] = box.min[2] + ((FLOAT)kk + (FLOAT) 0.5)*spacing[2];
        }
      }
    }
  }

  return;
}



//=================================================================================================
//  Ic::AddHexagonalLattice
/// Create simple hexagonal-packed lattice using A-B-A-B pattern in z-direction
/// N.B. the box is scaled to fit to the x-boxsize.
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddHexagonalLattice
 (const int Npart,                     ///< [in] No. of particles in lattice
  const int Nlattice[3],               ///< [in] Ptcls per dimension in lattice
  const DomainBox<ndim> &box,          ///< [in] Bounding box of particles
  const bool normalise,                ///< [in] Normalise lattice shape and size
  FLOAT *r)                            ///< [out] Positions of particles
{
  int i,k;                             // Particle and dimension counters
  int ii,jj,kk;                        // Aux. lattice counters
  FLOAT rad[ndim];                     // 'Radius' of particle in lattice

  debug2("[Ic::AddHexagonalLattice]");
  assert(r);

  // If normalised, ensure equal spacing between all particles.
  // Otherwise set spacing to fit bounding box.
  if (normalise) {
    for (k=0; k<ndim; k++) rad[k] = (FLOAT) 0.5*(box.max[0] - box.min[0])/(FLOAT) Nlattice[0];
  }
  else {
    for (k=0; k<ndim; k++) rad[k] = (FLOAT) 0.5*(box.max[k] - box.min[k])/(FLOAT) Nlattice[k];
  }


  // Create lattice depending on dimensionality
  //-----------------------------------------------------------------------------------------------
  if (ndim == 1) {
    for (ii=0; ii<Nlattice[0]; ii++) {
      i = ii;
      r[i] = box.min[0] + (FLOAT) 0.5*rad[0] + (FLOAT) 2.0*(FLOAT)ii*rad[0];
    }
  }

  //-----------------------------------------------------------------------------------------------
  else if (ndim == 2) {
    for (jj=0; jj<Nlattice[1]; jj++) {
      for (ii=0; ii<Nlattice[0]; ii++) {
        i = jj*Nlattice[0] + ii;
        r[ndim*i] = box.min[0] +
          (FLOAT) 0.5*rad[0] + ((FLOAT) 2.0*(FLOAT)ii + (FLOAT)(jj%2))*rad[0];
        r[ndim*i + 1] = box.min[1] +
          (FLOAT) 0.5*sqrt((FLOAT) 3.0)*rad[1] + (FLOAT)jj*sqrt(3.0)*rad[1];
      }
    }
  }

  //-----------------------------------------------------------------------------------------------
  else if (ndim == 3) {
#pragma omp parallel for default(none) shared(box,Nlattice,r,rad) private(i,ii,jj,kk)
    for (kk=0; kk<Nlattice[2]; kk++) {
      for (jj=0; jj<Nlattice[1]; jj++) {
        for (ii=0; ii<Nlattice[0]; ii++) {
          i = kk*Nlattice[0]*Nlattice[1] + jj*Nlattice[0] + ii;
          r[ndim*i] = box.min[0] + (FLOAT) 0.5*rad[0] +
            ((FLOAT) 2.0*(FLOAT) ii + (FLOAT) (jj%2) + (FLOAT) ((kk+1)%2))*rad[0];
          r[ndim*i + 1] = box.min[1] + (FLOAT) 0.5*sqrt((FLOAT) 3.0)*rad[1] +
            (FLOAT) jj*sqrt((FLOAT) 3.0)*rad[1] + (FLOAT) (kk%2)*rad[1]/sqrt((FLOAT) 3.0);
          r[ndim*i + 2] = box.min[2] + sqrt((FLOAT) 6.0)*rad[2]/(FLOAT) 3.0 +
            (FLOAT) kk*(FLOAT) 2.0*sqrt((FLOAT) 6.0)*rad[2]/(FLOAT) 3.0;
        }
      }
    }
  }

  return;
}



//=================================================================================================
//  Ic::CutSphere
/// Cut-out a sphere containing exactly 'Nsphere' particles from a uniform box of particles.
//=================================================================================================
template <int ndim>
int Ic<ndim>::CutSphere
 (const int Nsphere,                   ///< [in] Desired no. of particles in sphere
  const int Naux,                      ///< [in] No. of input particles
  const DomainBox<ndim> &box,          ///< [in] Bounding box of particles
  const bool exact,                    ///< [in] Flag if we require exact numbers (not used)
  FLOAT *r)                            ///< [inout] Positions of particles
{
  int i,k;                             // Particle and dimension counters
  int Ninterior = 0;                   // No. of particle
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT r_low = 0.0;                   // Lower-bound for bisection iteration
  FLOAT r_high;                        // Upper-bound for bisection iteration
  FLOAT radius;                        // Current radius containing Nsphere ptcls
  FLOAT rcentre[ndim];                 // Centre of sphere

  debug2("[Ic::CutSphere]");

  assert(r);

  // Find centre and shortest edge-length of bounding box
  r_high = (FLOAT) big_number;
  for (k=0; k<ndim; k++) {
    rcentre[k] = (FLOAT) 0.5*(box.min[k] + box.max[k]);
    r_high = min(r_high, (FLOAT) 0.5*(box.max[k] - box.min[k]));
  }

  // Bisection iteration to determine the radius containing the desired
  // number of particles
  //-----------------------------------------------------------------------------------------------
  do {
    radius = (FLOAT) 0.5*(r_low + r_high);
    Ninterior = 0;

    // Count how many particles lie inside current radius
    for (i=0; i<Naux; i++) {
      for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= radius*radius) Ninterior++;
    }

    // If it's impossible to converge on the desired number of particles, due
    // to lattice configurations, then exit iteration with approximate number
    // of particles (must be less than Nsphere due to memory).
    if (Ninterior < Nsphere && fabs(r_high - r_low)/radius < 1.e-8) break;

    // Otherwise, continue bisection iteration to find radius
    if (Ninterior > Nsphere) r_high = radius;
    if (Ninterior < Nsphere) r_low = radius;

  } while (Ninterior != Nsphere);


  // Now that the radius containing require number has been identified,
  // record only the particles inside the sphere.
  //-----------------------------------------------------------------------------------------------
  Ninterior = 0;
  for (i=0; i<Naux; i++) {
    for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
    drsqd = DotProduct(dr,dr,ndim);
    if (drsqd <= radius*radius) {
      for (k=0; k<ndim; k++) r[ndim*Ninterior + k] = r[ndim*i + k]/radius;
      Ninterior++;
    }
  }

  return Ninterior;
}



//=================================================================================================
//  Ic::AddAzimuthalDensityPerturbation
/// Add an azimuthal density perturbation for implementing Boss-Bodenheimer-type initial conditions
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddAzimuthalDensityPerturbation
 (const int Npart,                       ///< [in] No. of particles in sphere
  const int mpert,                       ///< [in] Perturbation mode
  const FLOAT amp,                       ///< [in] Amplitude of perturbation
  const FLOAT rcentre[ndim],             ///< [in] Position of sphere centre
  FLOAT *r)                              ///< [inout] Positions of particles
{
  if (ndim > 1) {
    int i,k;                             // Particle and dimension counters
    int j;                               // Aux. counter
    FLOAT phi,phi1,phi2,phiprime;        // Aux. azimuthal angle variables
    FLOAT Rsqd;                          // Radial distance (from z-axis) squared
    FLOAT Rmag;                          // Radial distance (from z-axis)
    FLOAT rpos[2];                       // Random position of new particle
    const int tabtot = 2048;             // No of elements in tables
    const FLOAT invmpert = (FLOAT) 1.0/(FLOAT) mpert;  // 1 / mpert
    const FLOAT spacing = twopi/(FLOAT)(tabtot - 1);   // Table angular spacing

    debug2("[Ic::AddAzimuthalDensityPerturbation]");
    assert(r);


    // Loop over all required particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(r,rcentre)\
  private(i,j,k,phi,phiprime,phi1,phi2,rpos,Rmag,Rsqd)
    for (i=0; i<Npart; i++) {
      for (k=0; k<2; k++) rpos[k] = r[ndim*i + k] - rcentre[k];

      // Calculate distance from z-axis and
      Rsqd = rpos[0]*rpos[0] + rpos[1]*rpos[1];
      Rmag = sqrt(Rsqd);

      // Find azimuthal angle around z-axis correcting for which quadrant
      if (Rmag > small_number) phi = asin(fabs(rpos[1])/Rmag);
      else phi = (FLOAT) 0.0;
      phiprime = (FLOAT) 0.0;

      if (rpos[0] < (FLOAT) 0.0 && rpos[1] > (FLOAT) 0.0) phi = pi - phi;
      else if (rpos[0] < (FLOAT) 0.0 && rpos[1] < (FLOAT) 0.0) phi = pi + phi;
      else if (rpos[0] > (FLOAT) 0.0 && rpos[1] < (FLOAT) 0.0) phi = twopi - phi;

      // Wrap angle to fit between 0 and two*pi
      if (phi < amp*invmpert) phi = phi + twopi;

      // Numerically find new phi angle for perturbation.  Search through grid of values,
      // find upper and lower bounds, then use linear interpolation to find new value of phi.
      for (j=1; j<tabtot; j++) {
        phi1 = spacing*(FLOAT) (j - 1);
        phi2 = spacing*(FLOAT) j;
        phi1 = phi1 + amp*cos((FLOAT) mpert*phi1)*invmpert;
        phi2 = phi2 + amp*cos((FLOAT) mpert*phi2)*invmpert;

        if (phi2 >= phi && phi1 < phi) {
          phiprime = spacing*(FLOAT)(j - 1) + spacing*(phi - phi1) / (phi2 - phi1);
          r[ndim*i] = rcentre[0] + Rmag*cos(phiprime);
          r[ndim*i + 1] = rcentre[1] + Rmag*sin(phiprime);
          break;
        }
      }

    }
    //---------------------------------------------------------------------------------------------

  }

  return;
}



//=================================================================================================
//  Ic::AddSinusoidalDensityPerturbation
/// Add a 1D sinusoidal density perturbation (in x-direction) to given uniform density field.
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddSinusoidalDensityPerturbation
 (int Npart,                     ///< [in] No. of particles in sphere
  FLOAT amp,                     ///< [in] Amplitude of perturbation
  FLOAT lambda,                  ///< [in] Wave number of perturbation
  FLOAT *r)                            ///< [inout] Positions of particles
{
  int i;                               // Particle counter
  FLOAT diff;                          // Convergence error/difference
  FLOAT xold;                          // Old particle x-position in iteration
  FLOAT xnew;                          // New particle position
  FLOAT kwave = twopi/lambda;          // Sine wave-number

  debug2("[Ic::AddSinusoidalDensityPerturbation]");
  assert(r);


  // Loop over all required particles
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(lambda,Npart,amp,kwave,r) private(diff,i,xnew,xold)
  for (i=0; i<Npart; i++) {
    xnew = r[ndim*i];

    // Solve iterative procedure for particle positions in sinusoid
    do {
      xold = xnew;
      xnew = r[ndim*i] - amp*((FLOAT) 1.0 - cos(kwave*xnew))/kwave;
      diff = fabs((xnew - xold)/lambda);
    } while (diff > (FLOAT) 1.0e-12);

    if (xnew > icBox.max[0]) xnew -= icBox.size[0];
    if (xnew < icBox.min[0]) xnew += icBox.size[0];

    r[ndim*i] = xnew;
  }
  //-----------------------------------------------------------------------------------------------


  return;
}



//=================================================================================================
//  Ic::AddRotationalVelocityField
/// Add a solid-body rotational velocity field
//=================================================================================================
template <int ndim>
void Ic<ndim>::AddRotationalVelocityField
 (const int Npart,                     ///< [in] No. of particles in sphere
  const FLOAT angvelaux,               ///< [in] Angular velocity of cloud
  const FLOAT rcentre[ndim],           ///< [in] Position of sphere centre
  const FLOAT *r,                      ///< [in] Positions of particles
  FLOAT *v)                            ///< [out] Velocities of particles
{
  // Only compile for 2 or 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim > 1) {

    int i,k;                           // Particle and dimension counters
    FLOAT Rmag;                        // Distance from z-axis
    FLOAT Rsqd;                        // Distance squared from z-axis
    FLOAT dr[3];                       // Relative position vector

    debug2("[Ic::AddRotationalVelocityField]");
    assert(r);
    assert(v);


    // Loop over all required particles
    //---------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) shared(r,rcentre,v) private(dr,i,k,Rmag,Rsqd)
    for (i=0; i<Npart; i++) {
      for (k=0; k<ndim; k++) dr[k] = r[ndim*i + k] - rcentre[k];
      for (k=0; k<ndim; k++) v[ndim*i + k] = (FLOAT) 0.0;

      // Calculate distance from z-axis and
      Rsqd = dr[0]*dr[0] + dr[1]*dr[1] + small_number;
      Rmag = sqrt(Rsqd);

      // Find azimuthal angle around z-axis correcting for which quadrant
      if (Rmag > small_number) {
        dr[0] = dr[0]/Rmag;
        dr[1] = dr[1]/Rmag;

        v[ndim*i] = -angvelaux*Rmag*dr[1];
        v[ndim*i + 1] = angvelaux*Rmag*dr[0];
      }
    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  Ic::GenerateTurbulentVelocityField
/// Generates turbulent velocity field using FFTW library.
/// Based on original code by A. McLeod.
//=================================================================================================
template <int ndim>
void Ic<ndim>::GenerateTurbulentVelocityField
 (const int field_type,                ///< [in] Type of turbulent velocity field
  const int gridsize,                  ///< [in] Size of velocity grid
  const DOUBLE power_turb,             ///< [in] Power spectrum index
  DOUBLE *vfield)                      ///< [out] Array containing velocity field grid
{
#if defined(FFTW_TURBULENCE)

  // Only valid for 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    bool divfree;                      // Select div-free turbulence
    bool curlfree;                     // Select curl-free turbulence
    int kmax;                          // Max. extent of k (in 3D)
    int kmin;                          // Min. extent of k (in 3D)
    int krange;                        // Range of k values (kmax - kmin + 1)
    int i,j,k;                         // Grid counters
    int ii,jj,kk;                      // Aux. grid counters
    int kmod;                          // ..
    int d;                             // Dimension counter
    DOUBLE F[ndim];                    // Fourier vector component
    DOUBLE unitk[3];                   // Unit k-vector
    DOUBLE **power,**phase;            // ..
    DOUBLE Rnd[3];                     // Random numbers, variable in Gaussian calculation
    //DOUBLE k_rot[ndim];                // bulk rotation modes
    //DOUBLE k_com[ndim];                // bulk compression modes
    fftw_plan plan;                    // ??
    fftw_complex *incomplexfield;      // ..
    fftw_complex *outcomplexfield;     // ..

    debug2("[Ic::GenerateTurbulentVelocityField]");


    divfree = false;
    curlfree = false;
    if (field_type == 1) curlfree = true;
    if (field_type == 2) divfree = true;

    kmin = -(gridsize/2 - 1);
    kmax = gridsize/2;
    krange = kmax - kmin + 1;
//    cout << "kmin : " << kmin << "    kmax : " << kmax << "    krange : "
//         << krange << "    gridsize : " << gridsize << endl;
    if (krange != gridsize) {
      string msg="Error : krange != gridsize in Ic::GenerateTurbulentVelocityField";
      ExceptionHandler::getIstance().raise(msg);
    }
    power = new DOUBLE*[3];  phase = new DOUBLE*[3];
    for (d=0; d<3; d++) power[d] = new DOUBLE[krange*krange*krange];
    for (d=0; d<3; d++) phase[d] = new DOUBLE[krange*krange*krange];
    incomplexfield = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gridsize*gridsize*gridsize);
    outcomplexfield = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*gridsize*gridsize*gridsize);

    for (d=0; d<3; d++) {
      for (i=0; i<krange*krange*krange; i++) power[d][i] = 0.0;
      for (i=0; i<krange*krange*krange; i++) phase[d][i] = 0.0;
    }
    for (i=0; i<3*gridsize*gridsize*gridsize; i++) vfield[i] = 0.0;


    // Define wave vectors in Fourier space.  Each wave vector has coordinates in Fourier space,
    // random phases in three dimensions, and power in three dimensions, giving a power and a
    // phase vector field in Fourier space.  With a 'natural' field type there is no coordination
    // between x,y,z components.  With a div-free or curl-free field, there is a relation between
    // coordinates in Fourier space and power vector. For a div-free field, power vectors are
    // perpendicular to the Fourier coordinate vector, and for a curl-free field power vectors
    // are (anti-)parallel to the Fourier coordinate vector
    // (i,j,k) is the Fourier coordinate vector
    // power(1:3,i,j,k) is the power vector at Fourier coordinates i,j,k
    // phase(1:3,i,j,k) is the phase vector at Fourier coordinates i,j,k
    //---------------------------------------------------------------------------------------------
    for (kmod=0; kmod<=kmax; kmod++) {

      for (i=kmin; i<=kmax; i++) {
        for (j=kmin; j<=kmax; j++) {
          for (k=kmin; k<=kmax; k++) {

            // Skip any k-vectors that have already been calculated
            if (abs(i) != kmod && abs(j) != kmod && abs(k) != kmod) continue;
            if (abs(i) > kmod || abs(j) > kmod || abs(k) > kmod) continue;

            //continue;
            ii = (i + krange)%krange;
            jj = (j + krange)%krange;
            kk = (k + krange)%krange;

            // cycle antiparallel k-vectors
            //if (k < 0) continue;
            //if (k == 0) {
            //  if (j < 0) continue;
            //  if (j == 0 && i < 0) continue;
            //}

            // Central power = 0
            if (i == 0 && j == 0 && k == 0) continue;
            if (i*i + j*j + k*k >= kmax*kmax) continue;

            // Power value, to be multipled by random power chosen from a Gaussian
            // This is what gives the slope of the eventual power spectrum
            for (d=0; d<3; d++) F[d] = sqrt(pow(sqrt((DOUBLE)(i*i + j*j + k*k)),power_turb));

            for (d=0; d<3; d++) {
              Rnd[0] = randnumb->floatrand();

              // Random phase between 0 and 2*pi (actually -pi and +pi).
              phase[d][ii + krange*jj + krange*krange*kk] = (2.0*Rnd[0] - 1.0)*pi;

              // Create Gaussian distributed random numbers
              Rnd[1] = randnumb->gaussrand(0.0,1.0);
              Rnd[2] = randnumb->gaussrand(0.0,1.0);
              F[d] = Rnd[1]*F[d];
            }

            // Helmholtz decomposition!
            unitk[0] = (DOUBLE) i;
            unitk[1] = (DOUBLE) j;
            unitk[2] = (DOUBLE) k;
            DOUBLE ksqd = DotProduct(unitk,unitk,3);
            for (d=0; d<3; d++) unitk[d] /= sqrt(ksqd);

            // For curl free turbulence, vector F should be
            // parallel/anti-parallel to vector k
            if (curlfree) {
              for (d=0; d<3; d++) {
                power[d][ii + krange*jj + krange*krange*kk] = unitk[d]*DotProduct(F,unitk,3);
              }
            }

            // For divergence free turbulence, vector F should be perpendicular to vector k
            else if (divfree) {
              for (d=0; d<3; d++) {
                power[d][ii + krange*jj + krange*krange*kk] = F[d] - unitk[d]*DotProduct(F,unitk,3);
              }
            }
            else {
              for (d=0; d<3; d++) power[d][ii + krange*jj + krange*krange*kk] = F[d];
            }

          }
        }
      }

    }
    //---------------------------------------------------------------------------------------------


    DOUBLE power_spectrum[kmax+1][3];
    for (i=0; i<kmax+1; i++) {
      for (d=0; d<3; d++) power_spectrum[i][d] = 0.0;
    }

    for (i=kmin; i<=kmax; i++) {
      for (j=kmin; j<=kmax; j++) {
        for (k=kmin; k<=kmax; k++) {
          ii = (i + krange)%krange;
          jj = (j + krange)%krange;
          kk = (k + krange)%krange;
          DOUBLE kmag = sqrt((DOUBLE)(i*i + j*j + k*k));
          if (kmag >= kmax) kmag = kmax;
          int ibin = (int) (1.0001*kmag);
          power_spectrum[ibin][0] += pow(power[0][ii + krange*jj + krange*krange*kk],2);
          power_spectrum[ibin][1] += pow(power[1][ii + krange*jj + krange*krange*kk],2);
          power_spectrum[ibin][2] += pow(power[2][ii + krange*jj + krange*krange*kk],2);
        }
      }
    }

    plan = fftw_plan_dft_3d(gridsize, gridsize, gridsize, incomplexfield,
                            outcomplexfield, FFTW_BACKWARD, FFTW_ESTIMATE);


    // reorder array: positive wavenumbers are placed in ascending order along
    // first half of dimension, i.e. 0 to k_max, negative wavenumbers are placed
    // along second half of dimension, i.e. -k_min to 1.
    for (d=0; d<3; d++) {

      for (i=kmin; i<=kmax; i++) {
        for (j=kmin; j<=kmax; j++) {
          for (k=kmin; k<=kmax; k++) {
            ii = (i + krange)%krange;
            jj = (j + krange)%krange;
            kk = (k + krange)%krange;
            //ii = i - kmin;
            //jj = j - kmin;
            //kk = k - kmin;
            incomplexfield[ii +krange*jj + krange*krange*kk][0] =
              power[d][ii + krange*jj + krange*krange*kk]*
              cos(phase[d][ii + krange*jj + krange*krange*kk]);
            incomplexfield[ii +krange*jj + krange*krange*kk][1] =
              power[d][ii + krange*jj + krange*krange*kk]*
              sin(phase[d][ii + krange*jj + krange*krange*kk]);

          }
        }
      }

      //fftw_execute_dft(plan, complexfield, complexfield);
      fftw_execute_dft(plan, incomplexfield, outcomplexfield);


      for (i=0; i<krange; i++) {
        for (j=0; j<krange; j++) {
          for (k=0; k<krange; k++) {
            vfield[d + 3*i + 3*krange*j + 3*krange*krange*k] =
              outcomplexfield[i + krange*j + krange*krange*k][0];

          }
        }
      }

    }


    fftw_destroy_plan(plan);

    for (d=0; d<3; d++) delete[] phase[d];
    for (d=0; d<3; d++) delete[] power[d];
    delete[] phase;
    delete[] power;
    fftw_free(outcomplexfield);
    fftw_free(incomplexfield);


  }
  //-----------------------------------------------------------------------------------------------

#endif

  return;
}



//=================================================================================================
//  Ic::InterpolateVelocityField
/// Calculate Interpolated velocity from uniform grid onto particle positions.
//=================================================================================================
template <int ndim>
void Ic<ndim>::InterpolateVelocityField
 (const int Npart,                       ///< [in] No of particles
  const int Ngrid,                       ///< [in] Size (per dim) of velocity grid
  const FLOAT xmin,                      ///< [in] Minimum position
  const FLOAT dxgrid,                    ///< [in] Grid size
  const FLOAT *r,                        ///< [in] Positions of particles
  const DOUBLE *vfield,                   ///< [in] Tabulated velocity field
  FLOAT *v)                              ///< [out] Interpolated particle velocity
{

  // Only compile routine for 3 dimensions
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3) {

    int i,j,k;                           // Grid coordinates
    int kk;                              // Dimension counter
    int p;                               // Particle counter
    FLOAT dx[3];                         // Position relative to grid point
    FLOAT vint[8];                       // Interpolated velocity


    // Now interpolate velocity field onto particle positions
    //---------------------------------------------------------------------------------------------
#pragma omp parallel for default(none) private(dx,i,j,k,kk,p,vint) \
    shared(cout,r,v,vfield)
    for (p=0; p<Npart; p++) {
      for (kk=0; kk<ndim; kk++) dx[kk] = (r[ndim*p + kk] - xmin)/dxgrid;

      i = (int) dx[0];
      j = (int) dx[1];
      k = (int) dx[2];

      if (i > Ngrid || j > Ngrid || k > Ngrid || i < 0 || j < 0 || k < 0)  {
        cout << "Problem with velocity interpolation grid!! : "
             << i << "    " << j << "    " << k << "   " << Ngrid << endl;
        ExceptionHandler::getIstance().raise("Problem with velocity interpolation grid");
      }

      for (kk=0; kk<ndim; kk++) dx[kk] -= (int) dx[kk];

      // Interpolate to get more accurate velocities
      vint[0] = ((FLOAT) 1.0 - dx[0])*((FLOAT) 1.0 - dx[1])*((FLOAT) 1.0 - dx[2]);
      vint[1] = ((FLOAT) 1.0 - dx[0])*((FLOAT) 1.0 - dx[1])*dx[2];
      vint[2] = ((FLOAT) 1.0 - dx[0])*dx[1]*((FLOAT) 1.0 - dx[2]);
      vint[3] = ((FLOAT) 1.0 - dx[0])*dx[1]*dx[2];
      vint[4] = dx[0]*((FLOAT) 1.0 - dx[1])*((FLOAT) 1.0 - dx[2]);
      vint[5] = dx[0]*((FLOAT) 1.0 - dx[1])*dx[2];
      vint[6] = dx[0]*dx[1]*((FLOAT) 1.0 - dx[2]);
      vint[7] = dx[0]*dx[1]*dx[2];

      v[ndim*p]   = vint[0]*vfield[3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[1]*vfield[3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[2]*vfield[3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[3]*vfield[3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)] +
                    vint[4]*vfield[3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[5]*vfield[3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[6]*vfield[3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[7]*vfield[3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)];

      v[ndim*p+1] = vint[0]*vfield[1 + 3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[1]*vfield[1 + 3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[2]*vfield[1 + 3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[3]*vfield[1 + 3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)] +
                    vint[4]*vfield[1 + 3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[5]*vfield[1 + 3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[6]*vfield[1 + 3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[7]*vfield[1 + 3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)];

      v[ndim*p+2] = vint[0]*vfield[2 + 3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[1]*vfield[2 + 3*i + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[2]*vfield[2 + 3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[3]*vfield[2 + 3*i + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)] +
                    vint[4]*vfield[2 + 3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*k] +
                    vint[5]*vfield[2 + 3*(i+1) + 3*Ngrid*j + 3*Ngrid*Ngrid*(k+1)] +
                    vint[6]*vfield[2 + 3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*k] +
                    vint[7]*vfield[2 + 3*(i+1) + 3*Ngrid*(j+1) + 3*Ngrid*Ngrid*(k+1)];

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class Ic<1>;
template class Ic<2>;
template class Ic<3>;
