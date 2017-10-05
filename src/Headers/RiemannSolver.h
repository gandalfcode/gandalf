//=================================================================================================
//  RiemannSolver.h
//  Contains all class definitions for Riemann solvers used in Godunov scheme.
//=================================================================================================


#ifndef _RIEMANN_SOLVER_H_
#define _RIEMANN_SOLVER_H_


#include "Precision.h"
#include "Constants.h"
#include "InlineFuncs.h"

enum RiemannSolverEnum {
  exact = 1,
  hllc  = 2  
} ;


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

};



//=================================================================================================
//  Class HllcRiemannSolver
/// \brief   HLLC approximate Riemann solver.
/// \details HLLC approximate Riemann solver.
/// \author  S. Heigl, R. Booth
/// \date    01/10/2014
//=================================================================================================
template<int ndim>
class HllcRiemannSolver
{
 public:
	HllcRiemannSolver(double gamma,
	                  bool zero_mass_flux=false,
	                  bool isothermal=false)
   : _gamma(gamma), _zmf(zero_mass_flux), _isothermal(isothermal)
  { } ;


  void ComputeFluxes(const FLOAT Wl[ndim+2], const FLOAT Wr[ndim+2],
                     const FLOAT n[ndim], FLOAT vface[ndim],
                     FLOAT flux[ndim+2][ndim]) {

    FLOAT flux_tmp[ndim+2] ;
    solve(Wl, Wr, n, vface, flux_tmp) ;

    // return vector flux:
    for (int i(0); i < ndim+2; ++i)
      for (int j(0); j < ndim; ++j)
        flux[i][j] = flux_tmp[i]*n[j] ;
  }



 private:
  void solve(const FLOAT Wl[ndim+2], const FLOAT Wr[ndim+2],
             const FLOAT n[ndim],
             FLOAT vface[ndim], FLOAT flux[ndim+2]) const  {


    // Compute speed of the 3 waves
    HLLC_State Sl(Wl, n, _gamma, _isothermal) ;
    HLLC_State Sr(Wr, n, _gamma, _isothermal) ;

    double Smin, Smax, vm ;
    HLL_Speeds(Sl, Sr, Smin, Smax) ;
    vm = compute_central_wave_speed(Sl, Sr, Smin, Smax) ;

    // Move to frame of contact discontinuity
    if (_zmf) {
      Smin -= vm ;
      Smax -= vm ;
      Sl.vline -= vm ;
      Sr.vline -= vm ;

      for (int i(0); i < ndim; ++i) {
        Sl.v[i]  -= vm*n[i] ;
        Sr.v[i]  -= vm*n[i] ;
        vface[i] += vm*n[i] ;
      }
      vm = 0 ;
    }

    // Compute the fluxes
    if (Smax <= 0)
      Hydro_Flux(Sr, n, flux) ;
    else if (Smin >= 0)
      Hydro_Flux(Sl, n, flux) ;
    else {
      // Compute the flux on the correct side of the contact wave{
      if (vm > 0) {
          Hydro_Flux(Sl, n, flux) ;
          add_RH_flux(Sl, n, Smin, vm, flux) ;
        }
        else {
          Hydro_Flux(Sr, n, flux) ;
          add_RH_flux(Sr, n, Smax, vm, flux) ;
        }
    }

    if (_zmf) {
      flux[irho] = 0 ;
    }

    // Convert back to the lab frame
    for (int i(0); i < ndim; ++i) {
      flux[iE] += flux[i]    * vface[i] ;
      flux[i]  += flux[irho] * vface[i] ;
    }
    flux[iE] += flux[irho] * 0.5*DotProduct(vface, vface, ndim) ;
  }

 private:

  class HLLC_State {
  public:
    HLLC_State(const FLOAT W[ndim+2], const FLOAT n[ndim], double gamma,
               bool isothermal)
      : rho(W[irho]),
        press(W[ipress]),
        vline(DotProduct(W, n, ndim))
     {
      if (not isothermal)
        cs = sqrt(gamma * press / rho) ;
      else
        cs = sqrt(        press / rho) ;

      e = 0 ;
      for (int i(0); i < ndim; ++i) {
        v[i] = W[i] ;
        e += v[i]*v[i] ;
      }
      if (not isothermal)
        e = 0.5 * rho * e + press / (gamma - 1) ;
      else
        e = 0.5 * rho * e + press ;
     }
    FLOAT v[ndim] ;
    FLOAT rho, press, e, cs ;
    FLOAT vline ;
  };

  void add_RH_flux(const HLLC_State& S, const FLOAT n[ndim],
                   double vwave, double vm,
		               FLOAT flux[ndim+2]) const {

    // Conserved quantities of edge state:
    FLOAT Q[nvar] ;
    for (int i(0); i < ndim; ++i)
      Q[i] = S.rho * S.v[i] ;

    Q[irho] = S.rho ;
    Q[iE]   = S.e ;


    // Compute conserved quantities of the starred state:
    double vs  = S.vline ;
    double dms = S.rho*(vs - vwave) ;

    FLOAT Qs[nvar] ;
    Qs[irho] = Q [irho]*(vwave - vs) / (vwave - vm) ;
    Qs[iE]   = Qs[irho]*(Q[iE] / Q[irho] +
			 (vm - vs)*(vm - S.press / dms));

    for (int i(0); i < ndim; ++i)
      Qs[i] = Qs[irho] * (S.v[i] + (vm - vs)*n[i]) ;

    // Add extra terms from RH conditions
    for (int i(0); i < nvar; ++i)
      flux[i] += vwave * (Qs[i] - Q[i]) ;
  }



  /* Roe_average_HLL_speeds
   *
   * Compute the Roe-averaged HLL wave speeds. Max/min wave-speed estimates
   * are from Einfeldt et al (1991), Batten et al (1997) or Toro (1999).
   */
  void HLL_Speeds(const HLLC_State& Sl, const HLLC_State& Sr,
                  double& Smin, double& Smax) const {

    double R = sqrt(Sr.rho / Sl.rho) ;
    double fl = 1. / (1. + R) ;
    double fr = 1. - fl ;

    // Average velocity
    double
      vl = Sl.vline ,
      vr = Sr.vline ;

    double v_av = fl*vl + fr*vr ;

    double cs_l = Sl.cs ;
    double cs_r = Sr.cs ;


    // Compute Average sound-speed
    double cs_av ;
    if (!_isothermal) {
      // Compute the sound speed from the Roe-averaged enthalpy in the standard
      // way.
      double dv2 = 0 ;
      for (int i(0); i < ndim; ++i) {
        double dvi = Sl.v[i] - Sr.v[i] ;
        dv2 += dvi * dvi ;
      }
      cs_av =
      sqrt(fl*cs_l*cs_l + fr*cs_r*cs_r + 0.5*fl*fr*(_gamma-1)*dv2) ;
    }
    else {
      // The sound speeds should just be the same, but they might differ
      // slightly (if e.g. a locally isothermal equations of state are used). We
      // use a Roe averaged thermal energy as a sensible guess in this case.
      cs_av = sqrt(fl*cs_l*cs_l + fr*cs_r*cs_r) ;
    }

    // Final Wave-speeds
    Smin = min(vl - cs_l, v_av - cs_av);
    Smax = max(vr + cs_r, v_av + cs_av);
  }

  // Middle wave-speed : Contact discontinuity
  double compute_central_wave_speed(const HLLC_State& Sl,
                                    const HLLC_State& Sr,
                                    double Smin, double Smax) const {

    double vl = Sl.vline, vr = Sr.vline ;

    double
      dml = Sl.rho * (vl - Smin),
      dmr = Sr.rho * (vr - Smax);

    if (not _isothermal) {
      double
        Pl = vl*dml + Sl.press,
        Pr = vr*dmr + Sr.press;

      return (Pr - Pl) / (dmr - dml) ;
    }
    else {
      return (Smax*dml - Smin*dmr) / (dml - dmr) ;
    }
  }

  /* Evaluate the hydrodynamic flux from state s in direction [0] */
  void Hydro_Flux(const HLLC_State& S, const FLOAT n[ndim],
                  FLOAT flux[ndim+2]) const
  {
    double rho = S.rho ;
    double P   = S.press ;

    double vj = S.vline ;
    double E = S.e ;

    for (int i(0); i < ndim; ++i)
      flux[i] = rho*vj*S.v[i] + P*n[i] ;

    flux[irho] = rho * vj ;
    flux[ipress]  = (P + E)*vj ;
  }

 private:
  double _gamma ;
  bool _zmf, _isothermal ;

  static const int nvar = ndim + 2;
  static const int irho   = ndim ;
  static const int ipress = ndim + 1 ;
  static const int iE     = ipress ;

} ;
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
#ifdef GANDALF_SNAPSHOT_SINGLE_PRECISION
  void ComputeShocktubeSolution(const std::string, float* vals, int N);
  static const bool single=true;
#else
  void ComputeShocktubeSolution(const std::string, double* vals, int N);
  static const bool single=false;
#endif

};

//=================================================================================================
//  Class ViscousFlux
/// \brief   Computes the viscous flux at an interface
/// \details Simple arithmetic average for the interface flux in a diffusion problem. Assumes
///          a constant kinematic viscosity (shear and bulk). Provides an interface for more
///          complex viscosities as well
/// \author  R. A. Booth
/// \date    17/5/2017
//=================================================================================================
template<int ndim>
class ViscousFlux
{
 public:

  ViscousFlux(FLOAT nu_shear, FLOAT nu_bulk=0)
   : _nu_shear0(nu_shear), _nu_bulk0(nu_bulk)
 { }

  virtual ~ViscousFlux(){} ;

  void ComputeViscousFlux(const FLOAT Wl[ndim+2], const FLOAT Wr[ndim+2],
                          const FLOAT gradv_l[][ndim], const FLOAT gradv_r[][ndim],
                          FLOAT flux[ndim+2][ndim]) const {

    // Compute average face state
    FLOAT W[ndim+2], gradv[ndim][ndim], div_v(0) ;
    for (int i=0; i < ndim+2; i++) W[i] = (Wl[i] + Wr[i])/2 ;
    for (int i=0; i < ndim; i++) {
      for (int j=0; j < ndim; j++) {
        gradv[i][j] = (gradv_l[i][j] + gradv_r[i][j])/2 ;
      }
      div_v += gradv[i][i];
    }

    // Face viscosity
    FLOAT eta_s = eta_shear(W);
    FLOAT eta_b = eta_bulk(W) ;


    FLOAT stress[ndim][ndim];
    for (int i=0; i < ndim; i++) {
      for (int j=0; j < ndim; j++) {
        stress[i][j] = eta_s * (gradv[i][j] + gradv[j][i]);
      }
      stress[i][i] += (eta_b - 2*eta_s/3)*div_v;
    }

    // Add the fluxes
    for (int i=0; i < ndim; i++) {
      for (int j=0; j < ndim; j++) {
        flux[i][j]  -= stress[i][j];
        flux[iE][j] -= stress[i][j]*W[i];
      }
    }
  }

protected:
  virtual FLOAT eta_shear(const FLOAT *W) const {
    return _nu_shear0 * W[irho];
  }
  virtual FLOAT eta_bulk(const FLOAT *W) const  {
    return _nu_bulk0 * W[irho];
  }

private:
  FLOAT _nu_shear0, _nu_bulk0 ;

  static const int nvar = ndim + 2;
  static const int irho   = ndim ;
  static const int ipress = ndim + 1 ;
  static const int iE     = ipress ;
};




#endif
