//=============================================================================
//  RiemannSolver.h
//  ..
//=============================================================================


#ifndef _RIEMANN_SOLVER_H_
#define _RIEMANN_SOLVER_H_


#include "Precision.h"
#include "Constants.h"
#include "Parameters.h"



//=============================================================================
//  Class RiemannSolver
/// \brief   Main parent class for Riemann solvers
/// \details ..
/// \author  D. A. Hubber
/// \date    26/04/2013
//=============================================================================
class RiemannSolver
{
 public:

  RiemannSolver(FLOAT);
  ~RiemannSolver();

  virtual void SolveRiemannProblem(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT,  
				   FLOAT, FLOAT, FLOAT &, FLOAT &) = 0;

  const FLOAT gamma;
  const FLOAT g1;
  const FLOAT g2;
  const FLOAT g3;
  const FLOAT g4;
  const FLOAT g5;

};



//=============================================================================
//  Class ExactRiemannSolver
/// \brief   Exact Riemann solver as described by Toro (19??)
/// \details ..
/// \author  D. A. Hubber
/// \date    26/04/2013
//=============================================================================
class ExactRiemannSolver: public RiemannSolver
{
  using RiemannSolver::gamma;
  using RiemannSolver::g1;
  using RiemannSolver::g2;
  using RiemannSolver::g3;
  using RiemannSolver::g4;
  using RiemannSolver::g5;

 public:

  ExactRiemannSolver(FLOAT gamma_aux): RiemannSolver(gamma_aux) {};
    //~ExactRiemannSolver() {};

  void SolveRiemannProblem(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, 
			   FLOAT, FLOAT, FLOAT, FLOAT &, FLOAT &);

};



//=============================================================================
//  Class HllcRiemannSolver
/// \brief   HLLC Riemann Solver
/// \details ..
/// \author  D. A. Hubber
/// \date    26/04/2013
//=============================================================================
class HllcRiemannSolver: public RiemannSolver
{
  using RiemannSolver::gamma;
  using RiemannSolver::g1;
  using RiemannSolver::g2;
  using RiemannSolver::g3;
  using RiemannSolver::g4;
  using RiemannSolver::g5;

 public:

  HllcRiemannSolver(FLOAT gamma_aux): RiemannSolver(gamma_aux) {};
    //~HllcRiemannSolver() {};

  void SolveRiemannProblem(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, 
			   FLOAT, FLOAT, FLOAT, FLOAT &, FLOAT &);

};
#endif
