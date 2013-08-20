//=============================================================================
//  RiemannSolver.cpp
//  Contains all available Riemann solver functions.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G Rosotti
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
//=============================================================================


#include <math.h>
#include <map>
#include <string>
#include "Constants.h"
#include "Exception.h"
#include "RiemannSolver.h"
using namespace std;



//=============================================================================
//  RiemannSolver::RiemannSolver
/// Constructor for RiemannSolver class.
//=============================================================================
RiemannSolver::RiemannSolver(FLOAT gamma_aux):
  gamma(gamma_aux),
  g1(gamma_aux - 1.0),
  g2((gamma_aux - 1.0)/(gamma_aux + 1.0)),
  g3(0.5*(gamma_aux + 1.0)/gamma_aux),
  g4(0.5*(gamma_aux - 1.0)/gamma_aux),
  g5(1.0/g4)
{
}



//=============================================================================
//  ExactRiemannSolver::SolveRiemannProblem
/// Exact Riemann solver, based on approach outlined by Toro (19??).
//=============================================================================
void ExactRiemannSolver::SolveRiemannProblem
(FLOAT pl,                            ///< LHS pressure
 FLOAT pr,                            ///< RHS pressure
 FLOAT rhol,                          ///< LHS density
 FLOAT rhor,                          ///< RHS density
 FLOAT soundl,                        ///< LHS sound speed
 FLOAT soundr,                        ///< RHS sound speed
 FLOAT vl,                            ///< LHS velocity
 FLOAT vr,                            ///< RHS velocity
 FLOAT &pstar,                        ///< Intermediate pressure state
 FLOAT &vstar)                        ///< Velocity of intermediate state
{
  int niteration = 0;                 // No. of iterations
  FLOAT fl;                           // Left-state variable
  FLOAT flprime;                      // Left-state variable gradient
  FLOAT fr;                           // Right-state variable
  FLOAT frprime;                      // Right-state variable gradient
  FLOAT pold;                         // Old intermediate pressure
  FLOAT tolerance = 1.0e-6;           // Iteration tolerance
  FLOAT Al = 2.0/(gamma + 1.0)/rhol;  // Aux. variable for efficiency
  FLOAT Ar = 2.0/(gamma + 1.0)/rhor;  // ""
  FLOAT Bl = pl*g2;                   // ""
  FLOAT Br = pr*g2;                   // ""


  // Make guess of solution for first iteration
  // For now, use two-rarefaction waves approximation
  pstar = (soundl + soundr - 0.5*g1*(vr - vl))/
    (soundl/pow(pl,g4) + soundr/pow(pr,g4));
  pstar = pow(pstar,1.0/g4);

  // Main iteration loop
  // --------------------------------------------------------------------------
  do {
    niteration++;

    // Calculate contribution to f and fprime for LHS
    if (pstar > pl) {
      fl = (pstar - pl)*sqrt(Al/(pstar + Bl));
      flprime = sqrt(Al/(pstar + Bl))*(1.0 - 0.5*(pstar - pl)/(pstar + Bl));
    }
    else {
      fl = 2.0*soundl*(pow(pstar/pl,g4) - 1.0)/g1;
      flprime = pow(pstar/pl,-g3)/(rhol*soundl);
    }

    // Calculate contribution to f and fprime for RHS
    if (pstar > pr) {
      fr = (pstar - pr)*sqrt(Ar/(pstar + Br));
      frprime = sqrt(Ar/(pstar + Br))*(1.0 - 0.5*(pstar - pr)/(pstar + Br));
    }
    else {
      fr = 2.0*soundr*(pow(pstar/pr,g4) - 1.0)/g1;
      frprime = pow(pstar/pr,-g3)/(rhor*soundr);
    }

    // Perform Newton-Raphson iteration
    pold = pstar;
    pstar = pstar - (fl + fr + vr - vl)/(flprime + frprime);

    // Check if convergence has been achieved
    if (pstar < small_number) pstar = small_number;
    else if (2.0*fabs(pstar - pold)/(pstar + pold) < tolerance) break;

  } while(2 > 1);
  // --------------------------------------------------------------------------

  // Compute velocity of star region
  vstar = 0.5*(vl + vr) + 0.5*(fr - fl);

  return;
}



//=============================================================================
//  HllcRiemannSolver::SolveRiemannProblem
/// HLLC Riemann solver for pstar and vstar.
//=============================================================================
void HllcRiemannSolver::SolveRiemannProblem
(FLOAT pl,                          ///< LHS pressure
 FLOAT pr,                          ///< RHS pressure
 FLOAT rhol,                        ///< LHS density
 FLOAT rhor,                        ///< RHS density
 FLOAT soundl,                      ///< LHS sound speed
 FLOAT soundr,                      ///< RHS sound speed
 FLOAT vl,                          ///< LHS velocity
 FLOAT vr,                          ///< RHS velocity
 FLOAT &pstar,                      ///< Intermediate pressure state
 FLOAT &vstar)                      ///< Velocity of intermediate state
{
  FLOAT Sl = vl - soundl;           // LHS wave speed estimator
  FLOAT Sr = vr + soundr;           // RHS wave speed estimator

  // Compute intermediate ('star') velocity and pressure
  vstar = (pr - pl + rhol*vl*(Sl - vl) - rhor*vr*(Sr - vr))/
    (rhol*(Sl - vl) - rhor*(Sr - vr));
  pstar = ((Sr - vr)*rhor*pl - (Sl - vl)*rhol*pr + rhol*rhor*(Sr - vr)*
	   (Sl - vl)*(vr - vl))/((Sr - vr)*rhor - (Sl - vl)*rhol);

  return;
}

