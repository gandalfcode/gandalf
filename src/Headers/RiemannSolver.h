//=============================================================================
//  RiemannSolver.h
//  Contains all class definitions for Riemann solvers used in Godunov scheme.
//=============================================================================


#ifndef _RIEMANN_SOLVER_H_
#define _RIEMANN_SOLVER_H_


#include "Precision.h"
#include "Constants.h"
#include "InlineFuncs.h"


//=============================================================================
//  Class RiemannSolver
//=============================================================================
class RiemannSolver
{
 public:


  RiemannSolver(FLOAT);
  ~RiemannSolver();

  virtual void ComputeStarRegion(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT,
                                 FLOAT, FLOAT, FLOAT &, FLOAT &) {};
  virtual void ComputeFluxes(const FLOAT *, const FLOAT *, const FLOAT *, FLOAT *) = 0;
  //virtual void SolveRiemannProblem(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT,
  //                                 FLOAT, FLOAT, FLOAT &, FLOAT &) = 0;

  const FLOAT gamma;
  const FLOAT invgamma;
  const FLOAT g1;
  const FLOAT g2;
  const FLOAT g3;
  const FLOAT g4;
  const FLOAT g5;
  const FLOAT g6;
  const FLOAT g7;
  const FLOAT g8;
  const FLOAT g9;

};



//=============================================================================
//  Class ExactRiemannSolver
//=============================================================================
template <int ndim>
class ExactRiemannSolver: public RiemannSolver
{
  using RiemannSolver::gamma;
  using RiemannSolver::g1;
  using RiemannSolver::g2;
  using RiemannSolver::g3;
  using RiemannSolver::g4;
  using RiemannSolver::g5;
  using RiemannSolver::g6;
  using RiemannSolver::g7;
  using RiemannSolver::g8;
  using RiemannSolver::invgamma;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;

 public:

  ExactRiemannSolver(FLOAT gamma_aux): RiemannSolver(gamma_aux) {};
    //~ExactRiemannSolver() {};

  virtual void ComputeStarRegion(const FLOAT, const FLOAT, const FLOAT,
                                 const FLOAT, const FLOAT, const FLOAT,
                                 const FLOAT, const FLOAT, FLOAT &, FLOAT &);
  virtual void ComputeFluxes(const FLOAT *, const FLOAT *, const FLOAT *, FLOAT *);
  void SampleExactSolution(const FLOAT, const FLOAT, const FLOAT, const FLOAT,
                           const FLOAT, const FLOAT, const FLOAT, const FLOAT,
                           const FLOAT, const FLOAT, const FLOAT,
                           FLOAT &, FLOAT &, FLOAT &);

};



//=============================================================================
//  Class HllcRiemannSolver
//=============================================================================
template <int ndim>
class HllcRiemannSolver: public RiemannSolver
{
  using RiemannSolver::gamma;
  using RiemannSolver::g1;
  using RiemannSolver::g2;
  using RiemannSolver::g3;
  using RiemannSolver::g4;
  using RiemannSolver::g5;
  using RiemannSolver::g6;
  using RiemannSolver::g7;
  using RiemannSolver::g8;
  using RiemannSolver::invgamma;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;

 public:

  HllcRiemannSolver(FLOAT gamma_aux): RiemannSolver(gamma_aux) {};
    //~HllcRiemannSolver() {};

  virtual void ComputeFluxes(const FLOAT *, const FLOAT *, const FLOAT *, FLOAT *);

};
#endif
