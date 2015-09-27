//=================================================================================================
//  RiemannSolver.h
//  Contains all class definitions for Riemann solvers used in Godunov scheme.
//=================================================================================================


#ifndef _RIEMANN_SOLVER_H_
#define _RIEMANN_SOLVER_H_


#include "Precision.h"
#include "Constants.h"
#include "InlineFuncs.h"


//=================================================================================================
//  Class RiemannSolver
/// \brief   Virtual parent class for all Riemann solver classes.
/// \details Virtual parent class for all Riemann solver classes.
/// \author  D. A. Hubber, S. Heigl, J. Ngoumou
/// \date    01/10/2014
//=================================================================================================
template <int ndim>
class RiemannSolver
{
 public:

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
  const bool zeroMassFlux;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;


  RiemannSolver(FLOAT, bool);
  virtual ~RiemannSolver() {};

  virtual void ComputeStarRegion(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT,
                                 FLOAT, FLOAT, FLOAT &, FLOAT &) {};
  virtual void ComputeFluxes(const FLOAT [ndim+2], const FLOAT [ndim+2],
                             const FLOAT [ndim], FLOAT [ndim], FLOAT [ndim+2][ndim]) = 0;

  void ComputeRotationMatrices(const FLOAT dr[ndim], FLOAT rotMat[ndim][ndim], FLOAT invRotMat[ndim][ndim]);
  void RotateVector(const FLOAT rotMat[ndim][ndim], FLOAT vec[ndim]);

};



//=================================================================================================
//  Class ExactRiemannSolver
/// \brief   Exact Riemann solver solution; based on algorithm presented in Toro (1999).
/// \details Exact Riemann solver solution; based on algorithm presented in Toro (1999).
/// \author  D. A. Hubber, S. Heigl, J. Ngoumou
/// \date    01/10/2014
//=================================================================================================
template <int ndim>
class ExactRiemannSolver: public RiemannSolver<ndim>
{
  using RiemannSolver<ndim>::gamma;
  using RiemannSolver<ndim>::g1;
  using RiemannSolver<ndim>::g2;
  using RiemannSolver<ndim>::g3;
  using RiemannSolver<ndim>::g4;
  using RiemannSolver<ndim>::g5;
  using RiemannSolver<ndim>::g6;
  using RiemannSolver<ndim>::g7;
  using RiemannSolver<ndim>::g8;
  using RiemannSolver<ndim>::g9;
  using RiemannSolver<ndim>::invgamma;
  using RiemannSolver<ndim>::zeroMassFlux;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;

 public:

  ExactRiemannSolver(FLOAT gamma_aux, bool _zeroMassFlux) : RiemannSolver<ndim>(gamma_aux, _zeroMassFlux) {};
  virtual ~ExactRiemannSolver() {};

  virtual void ComputeStarRegion(const FLOAT, const FLOAT, const FLOAT, const FLOAT, const FLOAT,
                                 const FLOAT, const FLOAT, const FLOAT, FLOAT &, FLOAT &);
  virtual void ComputeFluxes(const FLOAT [nvar], const FLOAT [nvar],
                             const FLOAT [ndim], FLOAT [ndim], FLOAT [nvar][ndim]);
  void SampleExactSolution(const FLOAT, const FLOAT, const FLOAT, const FLOAT, const FLOAT,
                           const FLOAT, const FLOAT, const FLOAT, const FLOAT, const FLOAT,
                           const FLOAT, FLOAT &, FLOAT &, FLOAT &);
  FLOAT Prefun(const FLOAT, const FLOAT, const FLOAT, const FLOAT, FLOAT &);

};



//=================================================================================================
//  Class HllcRiemannSolver
/// \brief   HLLC approximate Riemann solver.
/// \details HLLC approximate Riemann solver.
/// \author  S. Heigl
/// \date    01/10/2014
//=================================================================================================
template <int ndim>
class HllcRiemannSolver: public RiemannSolver<ndim>
{
  using RiemannSolver<ndim>::gamma;
  using RiemannSolver<ndim>::g1;
  using RiemannSolver<ndim>::g2;
  using RiemannSolver<ndim>::g3;
  using RiemannSolver<ndim>::g4;
  using RiemannSolver<ndim>::g5;
  using RiemannSolver<ndim>::g6;
  using RiemannSolver<ndim>::g7;
  using RiemannSolver<ndim>::g8;
  using RiemannSolver<ndim>::g9;
  using RiemannSolver<ndim>::invgamma;
  using RiemannSolver<ndim>::zeroMassFlux;

  static const int nvar = ndim + 2;
  static const int ivx = 0;
  static const int ivy = 1;
  static const int ivz = 2;
  static const int irho = ndim;
  static const int ietot = ndim + 1;
  static const int ipress = ndim + 1;

 public:

  HllcRiemannSolver(FLOAT gamma_aux, bool _zeroMassFlux) : RiemannSolver<ndim>(gamma_aux, zeroMassFlux) {};
  virtual ~HllcRiemannSolver() {};

  virtual void ComputeFluxes(const FLOAT [nvar], const FLOAT [nvar],
                             const FLOAT [ndim], FLOAT [ndim], FLOAT [nvar][ndim]);

};



//=================================================================================================
//  Class ShocktubeSolution
/// \brief   Simple wrapper class to communicate Shocktube problem solution to python module.
/// \details Simple wrapper class to communicate Shocktube problem solution to python module.
/// \author  D. A. Hubber, S. Heigl, J. Ngoumou
/// \date    25/09/2015
//=================================================================================================
class ShocktubeSolution
{
 public:
  const int nvalues;
  const FLOAT gamma,pl,pr,rhol,rhor,t,vl,vr,x0,xl,xr;
  FLOAT pstar,ustar;
  ExactRiemannSolver<1>* riemann;

  ShocktubeSolution(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT,
                    FLOAT, FLOAT, FLOAT, FLOAT, int, FLOAT);
  ~ShocktubeSolution();

  //void ComputeShocktubeSolution(const std::string, int, float *);
  void ComputeShocktubeSolution(const std::string, float* vals, int N);

};
#endif
