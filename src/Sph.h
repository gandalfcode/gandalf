//=============================================================================
// Sph.h
// Contains main parent virtual class plus child classes for various SPH 
// algorithms that are implemented.
//=============================================================================


#ifndef _SPH_H_
#define _SPH_H_


#include <string>
#include "Precision.h"
#include "Constants.h"
#include "SphParticle.h"
#include "SphKernel.h"
#include "NbodyParticle.h"
#include "Nbody.h"
#include "Parameters.h"
#include "EOS.h"
#include "RiemannSolver.h"
using namespace std;


enum aviscenum{noneav, mon97};
enum acondenum{noneac, wadsley2008, price2008};



//=============================================================================
//  Class Sph
/// \brief   Main parent Sph class.
/// \details Different SPH implementations (e.g. grad-h SPH, Saitoh &
///          Makino 2012) are derived from this class.  Each implementation
///          requires defining its own version of each function (e.g. ComputeH
///          for its own method of computing smoothing lengths).
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim>
class Sph
{
 public:

  const aviscenum avisc;
  const acondenum acond;

  // Constructor
  // --------------------------------------------------------------------------
  Sph(int hydro_forces_aux, int self_gravity_aux, FLOAT alpha_visc_aux, 
      FLOAT beta_visc_aux, FLOAT h_fac_aux, FLOAT h_converge_aux, 
      aviscenum avisc_aux, acondenum acond_aux, string gas_eos_aux, 
      string KernelName);


  // SPH functions for computing SPH sums with neighbouring particles 
  // (fully coded in each separate SPH implementation, and not in Sph.cpp)
  // --------------------------------------------------------------------------
  virtual int ComputeH(int, int, FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
                       SphParticle<ndim> &, Nbody<ndim> *) = 0;
  virtual void ComputeSphHydroForces(int, int, int *, FLOAT *, FLOAT *, 
                                     FLOAT *, SphParticle<ndim> &,
                                     SphParticle<ndim> *) = 0;
  virtual void ComputeSphHydroGravForces(int, int, int *, SphParticle<ndim> &, 
					 SphParticle<ndim> *) = 0;
  virtual void ComputeSphGravForces(int, int, int *, SphParticle<ndim> &,
				    SphParticle<ndim> *) = 0;
  virtual void ComputeDirectGravForces(int, int, int *, FLOAT *, FLOAT *, 
				       SphParticle<ndim> &, 
                                       SphParticle<ndim> *) = 0;
  virtual void ComputeSphNeibDudt(int, int, int *, FLOAT *, FLOAT *,
				  FLOAT *, SphParticle<ndim> &, 
                                  SphParticle<ndim> *) = 0;
  virtual void ComputeSphDerivatives(int, int, int *, FLOAT *, FLOAT *,
				     FLOAT *, SphParticle<ndim> &, 
                                     SphParticle<ndim> *) = 0;
  virtual void ComputePostHydroQuantities(SphParticle<ndim> &) = 0;
  virtual void ComputeStarGravForces(int, NbodyParticle<ndim> **, 
				     SphParticle<ndim> &) = 0;


  // SPH array memory allocation functions
  // --------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);
  void SphBoundingBox(FLOAT *, FLOAT *, int);
  void InitialSmoothingLengthGuess(void);


  // SPH particle counters and main particle data array
  // --------------------------------------------------------------------------
  bool allocated;                     ///< Is SPH memory allocated?
  int Ngather;                        ///< Average no. of gather neighbours
  int Nsph;                           ///< No. of SPH particles in simulation
  int Nghost;                         ///< No. of ghost SPH particles
  int Ntot;                           ///< No. of real + ghost particles
  int Nsphmax;                        ///< Max. no. of SPH particles in array
  int Nghostmax;                      ///< Max. allowed no. of ghost particles

  const FLOAT alpha_visc;             ///< alpha artificial viscosity parameter
  const FLOAT beta_visc;              ///< beta artificial viscosity parameter
  const FLOAT h_fac;                  ///< Smoothing length-density factor
  const FLOAT h_converge;             ///< h-rho iteration tolerance
  const string gas_eos;               ///< Gas EOS option
  const int hydro_forces;             ///< Compute hydro forces?
  const int self_gravity;             ///< Compute gravitational forces?
  static const FLOAT invndim=1./ndim; ///< Copy of 1/ndim
  int create_sinks;                   ///< ..
    FLOAT mmean;                        ///< ..

  string riemann_solver;              ///< Selected Riemann solver
  string slope_limiter;               ///< Selected slope limiter
  int riemann_order;                  ///< Order of Riemann solver
  FLOAT kernfac;                      ///< Kernel range neighbour fraction
  FLOAT kernfacsqd;                   ///< Kernel range neib. fraction squared

  FLOAT *rsph;                           ///< Position array (for efficiency)
  struct SphParticle<ndim> *sphdata;  ///< Main SPH particle data array
  SphKernel<ndim> *kernp;             ///< Pointer to chosen kernel object
  TabulatedKernel<ndim> kerntab;      ///< Tabulated version of chosen kernel
  EOS<ndim> *eos;                     ///< Equation-of-state
  RiemannSolver *riemann;             ///< Riemann solver

};



//=============================================================================
//  Class GradhSph
/// \brief   Class definition for conservative 'grad-h' SPH simulations.
/// \details Class definition for conservative 'grad-h' SPH simulations 
///          (as derived from the parent Sph class).  Full code for each of 
///          these class functions written in 'GradhSph.cpp'.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
#if !defined(SWIG)
template <int ndim, template<int> class kernelclass>
class GradhSph: public Sph<ndim>
{
  using Sph<ndim>::allocated;
  using Sph<ndim>::Nsph;
  using Sph<ndim>::eos;
  using Sph<ndim>::h_fac;
  using Sph<ndim>::invndim;
  using Sph<ndim>::h_converge;
  using Sph<ndim>::hydro_forces;
  using Sph<ndim>::self_gravity;
  using Sph<ndim>::avisc;
  using Sph<ndim>::beta_visc;
  using Sph<ndim>::alpha_visc;
  using Sph<ndim>::acond;
  using Sph<ndim>::create_sinks;

 public:

  GradhSph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
           aviscenum, acondenum, string, string);
  ~GradhSph();

  int ComputeH(int, int, FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               SphParticle<ndim> &, Nbody<ndim> *);
  void ComputeSphGravForces(int, int, int *, SphParticle<ndim> &,
			    SphParticle<ndim> *);
  void ComputeSphHydroGravForces(int, int, int *, SphParticle<ndim> &, 
				 SphParticle<ndim> *);
  void ComputeSphHydroForces(int, int, int *, FLOAT *, FLOAT *, FLOAT *,
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphNeibDudt(int, int, int *, FLOAT *, FLOAT *,
			  FLOAT *, SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphDerivatives(int, int, int *, FLOAT *, FLOAT *, FLOAT *, 
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeDirectGravForces(int, int, int *, FLOAT *, FLOAT *,
			       SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputePostHydroQuantities(SphParticle<ndim> &);
  void ComputeStarGravForces(int, NbodyParticle<ndim> **, SphParticle<ndim> &);

  kernelclass<ndim> kern;                  ///< SPH kernel

};



//=============================================================================
//  Class SM2012Sph
/// \brief   Class definition for Saitoh & Makino (2012) SPH simulations
/// \details Class definition for Saitoh & Makino (2012) SPH simulations 
///          (as derived from the parent Sph class).  Full code for each of 
///          these class functions written in 'SM2012Sph.cpp'.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=============================================================================
template <int ndim, template<int> class kernelclass>
class SM2012Sph: public Sph<ndim>
{
  using Sph<ndim>::allocated;
  using Sph<ndim>::Nsph;
  using Sph<ndim>::eos;
  using Sph<ndim>::h_fac;
  using Sph<ndim>::invndim;
  using Sph<ndim>::h_converge;
  using Sph<ndim>::hydro_forces;
  using Sph<ndim>::avisc;
  using Sph<ndim>::beta_visc;
  using Sph<ndim>::alpha_visc;
  using Sph<ndim>::acond;

 public:

  SM2012Sph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
            aviscenum, acondenum, string, string);
  ~SM2012Sph();

  int ComputeH(int, int, FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               SphParticle<ndim> &, Nbody<ndim> *);
  void ComputeSphHydroForces(int, int, int *, FLOAT *, FLOAT *, FLOAT *,
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphHydroGravForces(int, int, int *, SphParticle<ndim> &, 
				 SphParticle<ndim> *);
  void ComputeSphGravForces(int, int, int *,
			    SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphNeibDudt(int, int, int *, FLOAT *, FLOAT *,
			  FLOAT *, SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphDerivatives(int, int, int *, FLOAT *, FLOAT *, FLOAT *, 
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeDirectGravForces(int, int, int *, FLOAT *, FLOAT *, 
                               SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputePostHydroQuantities(SphParticle<ndim> &);
  void ComputeStarGravForces(int, NbodyParticle<ndim> **, SphParticle<ndim> &);

  kernelclass<ndim> kern;                  ///< SPH kernel

};



//=============================================================================
//  Class GodunovSph
/// Class definition for Godunov SPH (Inutsuka 2002) algorithm.
/// Full code for each of these class functions
/// written in 'GodunovSph.cpp'.
//=============================================================================
template <int ndim, template<int> class kernelclass>
class GodunovSph: public Sph<ndim>
{
  using Sph<ndim>::allocated;
  using Sph<ndim>::Nsph;
  using Sph<ndim>::eos;
  using Sph<ndim>::h_fac;
  using Sph<ndim>::invndim;
  using Sph<ndim>::h_converge;
  using Sph<ndim>::hydro_forces;
  using Sph<ndim>::avisc;
  using Sph<ndim>::beta_visc;
  using Sph<ndim>::alpha_visc;
  using Sph<ndim>::acond;
  using Sph<ndim>::riemann;
  using Sph<ndim>::riemann_solver;
  using Sph<ndim>::riemann_order;
  using Sph<ndim>::slope_limiter;

 public:

  GodunovSph(int, int, FLOAT, FLOAT, FLOAT, FLOAT,
             aviscenum, acondenum, string, string);
  ~GodunovSph();

  int ComputeH(int, int, FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *,
               SphParticle<ndim> &, Nbody<ndim> *);
  void ComputeSphHydroForces(int, int, int *, FLOAT *, FLOAT *, FLOAT *, 
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphHydroGravForces(int, int, int *, SphParticle<ndim> &, 
				 SphParticle<ndim> *);
  void ComputeSphGravForces(int, int, int *,
			    SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphNeibDudt(int, int, int *, FLOAT *, FLOAT *,
  			  FLOAT *, SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeSphDerivatives(int, int, int *, FLOAT *, FLOAT *, FLOAT *, 
			     SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputeDirectGravForces(int, int, int *, FLOAT *, FLOAT *, 
                               SphParticle<ndim> &, SphParticle<ndim> *);
  void ComputePostHydroQuantities(SphParticle<ndim> &);
  void InitialiseRiemannProblem(SphParticle<ndim>, SphParticle<ndim>, FLOAT *,
                                FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT &, 
                                FLOAT &, FLOAT &, FLOAT &, FLOAT &, FLOAT &);
  void ComputeStarGravForces(int, NbodyParticle<ndim> **, SphParticle<ndim> &);

  kernelclass<ndim> kern;                 ///< SPH kernel

};
#endif


#endif
