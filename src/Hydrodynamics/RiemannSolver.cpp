//=================================================================================================
//  RiemannSolver.cpp
//  Contains all available Riemann solver functions.
//=================================================================================================


#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <map>
#include <string>
#include "Constants.h"
#include "RiemannSolver.h"
using namespace std;


static const int Niterationmax = 100;            // Max. no. of Riemann solver iterations
static const FLOAT tolerance = 1.0e-6;           // Iteration tolerance


//=================================================================================================
//  RiemannSolver::RiemannSolver
/// Constructor for RiemannSolver class.
//=================================================================================================
template <int ndim>
RiemannSolver<ndim>::RiemannSolver(FLOAT gamma_aux, bool _zeroMassFlux):
  gamma(gamma_aux),
  invgamma(1.0/gamma_aux),
  g1(0.5*(gamma_aux - 1.0)/gamma_aux),
  g2(0.5*(gamma_aux + 1.0)/gamma_aux),
  g3(2.0*gamma_aux/(gamma_aux - 1.0)),
  g4(2.0/(gamma_aux - 1.0)),
  g5(2.0/(gamma_aux + 1.0)),
  g6((gamma_aux - 1.0)/(gamma_aux + 1.0)),
  g7(0.5*(gamma_aux - 1.0)),
  g8(gamma_aux - 1.0),
  g9(1.0/(gamma_aux - 1.0)),
  zeroMassFlux(_zeroMassFlux)
{
}



//=================================================================================================
//  RiemannSolver::ComputeRotationMatrices
/// Compute all rotation matrices from unit vector that are required for transforming from original
/// coordinate frame to x-axis frame and back again once flux terms have been computed.
//=================================================================================================
template <int ndim>
void RiemannSolver<ndim>::ComputeRotationMatrices
 (const FLOAT runit[ndim],             ///< [in] Directional unit vector
  FLOAT rotMat[ndim][ndim],            ///< [out] Rotational matrix
  FLOAT invRotMat[ndim][ndim])         ///< [out] Inverse rotational matrix
{
  if (ndim == 1) {
    rotMat[0][0]    = runit[0];
    invRotMat[0][0] = runit[0];
  }
  else if (ndim == 2) {
    FLOAT theta = atan2(runit[1], runit[0]);
    rotMat[0][0] = cos(theta);
    rotMat[0][1] = -sin(theta);
    rotMat[1][0] = sin(theta);
    rotMat[1][1] = cos(theta);
    invRotMat[0][0] = cos(theta);
    invRotMat[0][1] = sin(theta);
    invRotMat[1][0] = -sin(theta);
    invRotMat[1][1] = cos(theta);
  }
  else if (ndim == 3) {
    const FLOAT phi = atan2(runit[1], runit[0]);
    const FLOAT delta = atan2(runit[2], sqrt(runit[0]*runit[0] + runit[1]*runit[1]));
    invRotMat[0][0] = cos(delta)*cos(phi);
    invRotMat[0][1] = cos(delta)*sin(phi);
    invRotMat[0][2] = sin(delta);
    invRotMat[1][0] = -sin(phi);
    invRotMat[1][1] = cos(phi);
    invRotMat[1][2] = 0.0;
    invRotMat[2][0] = -sin(delta)*cos(phi);
    invRotMat[2][1] = -sin(delta)*sin(phi);
    invRotMat[2][2] = cos(delta);
    rotMat[0][0] = cos(delta)*cos(phi);
    rotMat[0][1] = -sin(phi);
    rotMat[0][2] = -sin(delta)*cos(phi);
    rotMat[1][0] = cos(delta)*sin(phi);
    rotMat[1][1] = cos(phi);
    rotMat[1][2] = -sin(delta)*sin(phi);
    rotMat[2][0] = sin(delta);
    rotMat[2][1] = 0.0;
    rotMat[2][2] = cos(delta);
  }

  return;
}



//=================================================================================================
//  RiemannSolver::RotateVector
/// Rotate vector using the given rotation matrix
//=================================================================================================
template <int ndim>
void RiemannSolver<ndim>::RotateVector
 (const FLOAT rotMat[ndim][ndim],      ///< [in] Rotation mtrix
  FLOAT vec[ndim])                     ///< [inout] Vector to be rotated
{
  if (ndim == 1) {
    vec[0] = rotMat[0][0]*vec[0];
  }
  else if (ndim == 2) {
    FLOAT oldVec[ndim];
    for (int k=0; k<ndim; k++) oldVec[k] = vec[k];
    vec[0] = rotMat[0][0]*oldVec[0] + rotMat[0][1]*oldVec[1];
    vec[1] = rotMat[1][0]*oldVec[0] + rotMat[1][1]*oldVec[1];
  }
  else if (ndim == 3) {
    FLOAT oldVec[ndim];
    for (int k=0; k<ndim; k++) oldVec[k] = vec[k];
    vec[0] = rotMat[0][0]*oldVec[0] + rotMat[0][1]*oldVec[1] + rotMat[0][2]*oldVec[2];
    vec[1] = rotMat[1][0]*oldVec[0] + rotMat[1][1]*oldVec[1] + rotMat[1][2]*oldVec[2];
    vec[2] = rotMat[2][0]*oldVec[0] + rotMat[2][1]*oldVec[1] + rotMat[2][2]*oldVec[2];
  }

  return;
}



//=================================================================================================
//  ExactRiemannSolver::SolveRiemannProblem
/// Exact Riemann solver, based on approach outlined by Toro (1999).
//=================================================================================================
template <int ndim>
void ExactRiemannSolver<ndim>::ComputeStarRegion
 (const FLOAT pl,                      ///< [in] LHS pressure
  const FLOAT pr,                      ///< [in] RHS pressure
  const FLOAT dl,                      ///< [in] LHS density
  const FLOAT dr,                      ///< [in] RHS density
  const FLOAT cl,                      ///< [in] LHS sound speed
  const FLOAT cr,                      ///< [in] RHS sound speed
  const FLOAT ul,                      ///< [in] LHS velocity
  const FLOAT ur,                      ///< [in] RHS velocity
  FLOAT &pstar,                        ///< [out] Intermediate pressure state
  FLOAT &ustar)                        ///< [out] Velocity of intermediate state
{
  const FLOAT Al = g5/dl;              // ..
  const FLOAT Ar = g5/dr;              // ..
  const FLOAT Bl = pl*g6;              // ..
  const FLOAT Br = pr*g6;              // ..
  int iteration = 0;                   // No. of iterations (should be configured in parameter file)
  FLOAT cup;                           // Primitive variable ???
  FLOAT fl;                            // Left-state variable
  FLOAT flprime;                       // Left-state variable gradient
  FLOAT fr;                            // Right-state variable
  FLOAT frprime;                       // Right-state variable gradient
  FLOAT gel;                           // Two-shock variable
  FLOAT ger;                           // Two-shock variable
  FLOAT quser = 2.0;                   // First guess of solution
  FLOAT pmin;                          // Max. pressure
  FLOAT pmax;                          // Min. pressure
  FLOAT pold;                          // Old intermediate pressure
  FLOAT ppv;                           // Primitive Variable pressure
  // Two-rarefaction RS variables
  FLOAT pq;                            // ..
  FLOAT um;                            // ..
  FLOAT ptl;                           // ..
  FLOAT ptr;                           // ..

  // Compute guess from Primitive Variable RS
  cup  = (FLOAT) 0.25*(dl + dr)*(cl + cr);
  ppv  = (FLOAT) 0.5*(pl + pr) + 0.5*(ul - ur)*cup;
  ppv  = max((FLOAT) 0.0, ppv);
  pmin = min(pl, pr);
  pmax = max(pl, pr);

  // Check for vacuum condition
  if (cl + cr - g7*(ur - ul) <= (FLOAT) 0.0) {
    pstar = (FLOAT) 0.0;
    ustar = (FLOAT) 0.0;
    return;
  }

  // Select guess from PVRS
  if (pmax/pmin <= quser && pmin <= ppv && ppv <= pmax) {
    pstar = ppv;
  }
  // Select two-rarefaction RS
  else if (ppv < pmin) {
    pq    = pow(pl/pr,g1);
    um    = (pq*ul/cl + ur/cr + g4*(pq - 1.0))/(pq/cl + 1.0/cr);
    ptl   = 1.0 + g7*(ul - um)/cl;
    ptr   = 1.0 + g7*(um - ur)/cr;
    pstar = 0.5*(pl*pow(ptl,g3) + pr*pow(ptr,g3));
    //FLOAT pstar2 = pow((cl + cr - g7*(ur - ul))/(cl/pow(pl,g1) + cr/pow(pr,g1)), g3);
  }
  // Select two-shock RS with PVRS as estimate
  else {
    gel   = sqrt((g5/dl)/(g6*pl + ppv));
    ger   = sqrt((g5/dr)/(g6*pr + ppv));
    pstar = (gel*pl + ger*pr - (ur - ul))/(gel + ger);
  }

  pstar = max(pstar, small_number);
  assert(pstar > 0.0);

  // Main iteration loop
  //-----------------------------------------------------------------------------------------------
  do {
    iteration++;

    // Calculate contribution to f and fprime for LHS
    if (pstar > pl) {
      fl = (pstar - pl)*sqrt(Al/(pstar + Bl));
      flprime = sqrt(Al/(pstar + Bl))*(1.0 - 0.5*(pstar - pl)/(pstar + Bl));
    }
    else {
      fl = g4*cl*(pow(pstar/pl,g1) - 1.0);
      flprime = pow(pstar/pl,-g2)/(dl*cl);
    }

    // Calculate contribution to f and fprime for RHS
    if (pstar > pr) {
      fr = (pstar - pr)*sqrt(Ar/(pstar + Br));
      frprime = sqrt(Ar/(pstar + Br))*(1.0 - 0.5*(pstar - pr)/(pstar + Br));
    }
    else {
      fr = g4*cr*(pow(pstar/pr,g1) - 1.0);
      frprime = pow(pstar/pr,-g2)/(dr*cr);
    }

    // Perform Newton-Raphson iteration
    pold = pstar;
    pstar = pstar - (fl + fr + ur - ul)/(flprime + frprime);

    // Check if convergence has been achieved
    if (pstar < small_number) pstar = small_number; //this could also be smaller than zero
    else if (2.0*fabs(pstar - pold)/(pstar + pold) < tolerance) break;

    // Check the star variables have not become NaNs
    if (!isfinite(pstar)) {
      std::cout << "Checking Riemann values : " << pstar
                << "   iteration : " << iteration << std::endl;
      std::cout << "rho : " << dl << "   " << dr << "     vel : " << ul << "    "
                << ur << "    press : " << pl << "    " << pr << std::endl;
      std::cout << "f/fprime : " << fl << "   " << fr
                << "   " << flprime << "   " << frprime << endl;
      ExceptionHandler::getIstance().raise("Error : Invalid star values in Riemann solver");
    }


  } while (iteration < Niterationmax);
  //-----------------------------------------------------------------------------------------------

  // Compute velocity of star region
  ustar = 0.5*(ul + ur) + 0.5*(fr - fl);
  assert(pstar > 0.0);

  return;
}



//=================================================================================================
//  ExactRiemannSolver::SampleExactSolution
/// Sample the solution of the exact Riemann solver at the dimensionless position 's'.
//=================================================================================================
template <int ndim>
void ExactRiemannSolver<ndim>::SampleExactSolution
 (const FLOAT pstar,                   ///< [in] Pressure of central 'star' region
  const FLOAT ustar,                   ///< [in] Velocity in central 'star' region
  const FLOAT s,                       ///< [in] Dimensionless position to sample solution at
  const FLOAT pl,                      ///< [in] LHS pressure
  const FLOAT pr,                      ///< [in] RHS pressure
  const FLOAT dl,                      ///< [in] LHS density
  const FLOAT dr,                      ///< [in] RHS density
  const FLOAT cl,                      ///< [in] LHS sound speed
  const FLOAT cr,                      ///< [in] RHS sound speed
  const FLOAT ul,                      ///< [in] LHS velocity
  const FLOAT ur,                      ///< [in] RHS velocity
  FLOAT &p,                            ///< [out] Sampled pressure state
  FLOAT &d,                            ///< [out] Sampled density state
  FLOAT &u)                            ///< [out] Sampled velocity state
{

  FLOAT cm;                            // ..
  FLOAT sl,sr;                         // ..
  FLOAT sh, st;                        // Rarefaction head and tail speed
  FLOAT c;                             // Parameter to determine fan position


  // If sampling point lies to the left of the contact discontinuity
  //-----------------------------------------------------------------------------------------------
  if (s <= ustar) {

    // Left rarefaction
    //---------------------------------------------------------------------------------------------
    if (pstar <= pl) {
      sh = ul - cl;

      // Sampled point is left data state
      if (s <= sh) {
        d = dl;
        u = ul;
        p = pl;
      }
      else {
        cm = cl*pow(pstar/pl,g1);
        st = ustar - cm;

        // Sample point is star left state
        if (s > st) {
          d = dl*pow(pstar/pl,invgamma);
          u = ustar;
          p = pstar;
        }
        // Sampled point is inside left fan
        else {
          u = g5*(cl + g7*ul + s);
          c = g5*(cl + g7*(ul - s));
          d = dl*pow(c/cl,g4);
          p = pl*pow(c/cl,g3);
        }

      }

    }

    // left shock
    //---------------------------------------------------------------------------------------------
    else {
      FLOAT pml = pstar/pl;
      sl = ul - cl*sqrt(g2*pml + g1);

      // Sampled point is left data state
      if (s <= sl) {
        d = dl;
        u = ul;
        p = pl;
      }
      // Sampled point is star left state
      else {
        d = dl*(pml + g6)/(pml*g6 + 1.0);
        u = ustar;
        p = pstar;
      }
    }
    //---------------------------------------------------------------------------------------------

  }
  // Sampling point lies to the right of the contact discontinuity
  //-----------------------------------------------------------------------------------------------
  else {

    // Right shock
    //---------------------------------------------------------------------------------------------
    if (pstar >= pr) {

      FLOAT pmr = pstar/pr;
      sr = ur + cr*sqrt(g2*pmr + g1);

      // Sampled point is right data state
      if (s >= sr) {
        d = dr;
        u = ur;
        p = pr;
      }
      // Sampled point is star right state
      else {
        d = dr*(pmr + g6)/(pmr*g6 + 1.0);
        u = ustar;
        p = pstar;
      }

    }
    // Right rarefaction
    //---------------------------------------------------------------------------------------------
    else {

      sh = ur + cr;

      // Sampled point is left data state
      if (s >= sh) {
        d = dr;
        u = ur;
        p = pr;
      }
      else {
        cm = cr*pow(pstar/pr,g1);
        st = ustar + cm;

        // Sample point is star right state
        if (s <= st) {
          d = dr*pow(pstar/pr,invgamma);
          u = ustar;
          p = pstar;
        }
        // Sampled point is inside left fan
        else {
          u = g5*(-cr + g7*ur + s);
          c = g5*(cr - g7*(ur - s));
          d = dr*pow(c/cr,g4);
          p = pr*pow(c/cr,g3);
        }

      }

    }

  }
  //-----------------------------------------------------------------------------------------------

  // Check pressure is positive.  If not, print to screen and quit.
  if (p < 0.0) {
    cout << "s : " << s << "   " << ustar << "   " << pstar << "   " << pl << "    " << pr << endl;
    cout << "p : " << p << "    d : " << d << "    c : " << c << "    u : " << u << endl;
    ExceptionHandler::getIstance().raise("Error : Negative pstar in Riemann solver");
  }

  return;
}



//=================================================================================================
//  ExactRiemannSolver::ComputeFluxes
/// Exact Riemann solver, based on approach outlined by Toro (1999).
//=================================================================================================
template <int ndim>
void ExactRiemannSolver<ndim>::ComputeFluxes
 (const FLOAT Wleft[nvar],             ///< [in] LHS primitive state
  const FLOAT Wright[nvar],            ///< [in] RHS primitive state
  const FLOAT runit[ndim],             ///< [in] Unit vector pointing from left to right state
  FLOAT vface[ndim],                   ///< [in] Velocity of the working face
  FLOAT flux[nvar][ndim])              ///< [out] Flux vector
{
  const FLOAT cl = sqrt(gamma*Wleft[ipress]/Wleft[irho]);      // LHS sound speed
  const FLOAT cr = sqrt(gamma*Wright[ipress]/Wright[irho]);    // RHS sound speed
  int k,kv;                            // Dimension counters
  FLOAT p,d,u;                         // Primitive variables at s=0 from Riemann solver
  FLOAT rotMat[ndim][ndim];            // Rotation matrix
  FLOAT invRotMat[ndim][ndim];         // Inverse rotation matrix
  FLOAT uleft[ndim];                   // Left velocity state
  FLOAT uright[ndim];                  // Right velocity state
  FLOAT Wface[nvar];                   // Primitive vector at working face
  FLOAT pstar = small_number;          // Pressure in star region
  FLOAT ustar = (FLOAT) 0.0;           // Velocity in star region

  assert(Wleft[ipress] > 0.0);
  assert(Wleft[irho] > 0.0);
  assert(Wright[ipress] > 0.0);
  assert(Wright[irho] > 0.0);

  for (k=0; k<ndim; k++) uleft[k] = Wleft[k];
  for (k=0; k<ndim; k++) uright[k] = Wright[k];

  // Compute rotation matrices and rotate both left and right velocity states
  this->ComputeRotationMatrices(runit, rotMat, invRotMat);
  this->RotateVector(invRotMat, uleft);
  this->RotateVector(invRotMat, uright);

  // Compute p and u values at interface (in star region)
  ComputeStarRegion(Wleft[ipress], Wright[ipress], Wleft[irho], Wright[irho],
                    cl, cr, uleft[0], uright[0], pstar, ustar);


  // If not a vacuum state (i.e. pstar > 0), then sample solution
  //-----------------------------------------------------------------------------------------------
  if (pstar > 0.0) {
    SampleExactSolution(pstar, ustar, 0.0, Wleft[ipress], Wright[ipress], Wleft[irho],
                        Wright[irho], cl, cr, uleft[0], uright[0], p, d, u);
    assert(p >= 0.0);
    assert(d >= 0.0);

    for (kv=0; kv<ndim; kv++) Wface[kv] = (FLOAT) 0.0;
    Wface[irho]   = d;
    Wface[ivx]    = u;
    Wface[ipress] = p;

    // If using face that moves with star velocity, then mass flux should be zero.
    // (Used for Meshless finite-mass scheme of Hopkins 2015)
    if (zeroMassFlux) {

      Wface[ivx] = 0.0;
      for (k=0; k<ndim; k++) vface[k] += u*runit[k];

    }
    else {

      // Calculate transverse velocity fluxes from upwind direction
      if (u > 0.0) {
        for (k=1; k<ndim; k++) Wface[k] = uleft[k];
      }
      else {
        for (k=1; k<ndim; k++) Wface[k] = uright[k];
      }
    }


    // Rotate velocity to original frame
    this->RotateVector(rotMat, Wface);


    for (k=0; k<ndim; k++) Wface[k] += vface[k];
    FLOAT etot = (FLOAT) 0.0;
    for (kv=0; kv<ndim; kv++) etot += Wface[kv]*Wface[kv];
    etot = (FLOAT) 0.5*etot + Wface[ipress]/(gamma - (FLOAT) 1.0)/Wface[irho];
    for (k=0; k<ndim; k++) {
      for (kv=0; kv<ndim; kv++) flux[kv][k] = Wface[irho]*(Wface[k] - vface[k])*Wface[kv];
      flux[k][k]     = Wface[irho]*Wface[k]*(Wface[k] - vface[k]) + Wface[ipress];
      flux[irho][k]  = Wface[irho]*(Wface[k] - vface[k]);
      flux[ietot][k] = Wface[irho]*etot*(Wface[k] - vface[k]) + Wface[ipress]*Wface[k];
      //(Wface[ipress]/(gamma - 1.0) + 0.5*Wface[irho]*ekin)*(Wface[k])
    }


    // Compute fluxes in moving frame
    /*FLOAT ekin = (FLOAT) 0.0;
    for (kv=0; kv<ndim; kv++) ekin += Wface[kv]*Wface[kv];
    for (k=0; k<ndim; k++) {
      for (kv=0; kv<ndim; kv++) flux[kv][k] = Wface[irho]*Wface[k]*Wface[kv];
      flux[k][k]     = Wface[irho]*Wface[k]*Wface[k] + Wface[ipress];
      flux[irho][k]  = Wface[irho]*Wface[k];
      flux[ietot][k] = (Wface[ipress]/(gamma - (FLOAT) 1.0) + (FLOAT) 0.5*Wface[irho]*ekin)*Wface[k]
        + Wface[ipress]*Wface[k];
    }

    // Add corrections for transforming back to original lab frame
    for (k=0; k<ndim; k++) {
      flux[ietot][k] += (FLOAT) 0.5*DotProduct(vface, vface, ndim)*flux[irho][k] +
        DotProduct(vface, flux[k], ndim);
    }
    for (k=0; k<ndim; k++) {
      for (kv=0; kv<ndim; kv++) flux[kv][k] += vface[kv]*flux[irho][k];
    }*/

  }
  // Otherwise assume vacuum state conditions
  //-----------------------------------------------------------------------------------------------
  else {
    d = 0.0;
    u = 0.0;
    p = 0.0;

    for (kv=0; kv<ndim; kv++) Wface[kv] = (FLOAT) 0.0;
    for (int var=0; var<ndim+2; var++) {
      for (k=0; k<ndim; k++) flux[var][k] = (FLOAT) 0.0;
    }
  }
  //-----------------------------------------------------------------------------------------------


  return;
}




//=================================================================================================
//  ShocktubeSolution::ShocktubeSolution
/// Constructor for RiemannSolver class.
//=================================================================================================
ShocktubeSolution::ShocktubeSolution
 (FLOAT _rhol, FLOAT _rhor, FLOAT _vl, FLOAT _vr, FLOAT _pl, FLOAT _pr,
  FLOAT _xl, FLOAT _x0, FLOAT _xr, FLOAT _t, int _nvalues, FLOAT _gamma) :
  rhol(_rhol), rhor(_rhor), vl(_vl), vr(_vr), pl(_pl), pr(_pr),
  xl(_xl), x0(_x0), xr(_xr), t(_t), nvalues(_nvalues), gamma(_gamma)
{
  const FLOAT cl = sqrt(gamma*pl/rhol);        // LHS sound speed
  const FLOAT cr = sqrt(gamma*pr/rhor);        // RHS sound speed

  riemann = new ExactRiemannSolver<1>(_gamma, false);
  riemann->ComputeStarRegion(pl, pr, rhol, rhor, cl, cr, vl, vr, pstar, ustar);
}



//=================================================================================================
//  ShocktubeSolution::ShocktubeSolution
/// Destructor for RiemannSolver class.
//=================================================================================================
ShocktubeSolution::~ShocktubeSolution()
{
  delete riemann;
}



//=================================================================================================
//  ShocktubeSolution::ComputeShocktubeSolution
/// Destructor for RiemannSolver class.
//=================================================================================================
void ShocktubeSolution::ComputeShocktubeSolution
  (const std::string quantity,                  ///< [in] ..
//  (char *quantity,                              ///< [in]
#ifdef GANDALF_SNAPSHOT_SINGLE_PRECISION
   float* vals,                                ///< [out] ..
   int N)                                      ///< [in] ..
#else
   double* vals,
   int N)
#endif
{
  const FLOAT cl = sqrt(gamma*pl/rhol);        // LHS sound speed
  const FLOAT cr = sqrt(gamma*pr/rhor);        // RHS sound speed
  FLOAT d,p,s,u,x;

  //-----------------------------------------------------------------------------------------------
  for (int i=0; i<N; i++) {
    x = xl + (xr - xl)*((FLOAT) i) / (FLOAT) (N - 1);
    s = (x - x0)/t;
    riemann->SampleExactSolution(pstar, ustar, s, pl, pr, rhol, rhor, cl, cr, vl, vr, p, d, u);

    if (quantity == "x") vals[i] = x;
    else if (quantity == "rho") vals[i] = d;
    else if (quantity == "vx") vals[i] = u;
    else if (quantity == "press") vals[i] = p;
    else if (quantity == "u") vals[i] = (p/d/(gamma - 1.0));
    else vals[i] = 0.0f;
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class RiemannSolver<1>;
template class RiemannSolver<2>;
template class RiemannSolver<3>;
template class ExactRiemannSolver<1>;
template class ExactRiemannSolver<2>;
template class ExactRiemannSolver<3>;
template class HllcRiemannSolver<1>;
template class HllcRiemannSolver<2>;
template class HllcRiemannSolver<3>;
