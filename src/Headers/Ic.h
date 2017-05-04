//=================================================================================================
//  Ic.h
//  Contains class for generating the initial conditions
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


#ifndef _IC_H_
#define _IC_H_


#include "Exception.h"
#include "Precision.h"
#include "Hydrodynamics.h"
#include "InlineFuncs.h"
#include "Simulation.h"
#include "SimUnits.h"
#include "SmoothingKernel.h"
#include "RandomNumber.h"
#include "Parameters.h"
#include "Particle.h"
#if defined(FFTW_TURBULENCE)
#include "fftw3.h"
#endif


// Fwd Declarations
template <int ndim> class NeighbourSearch ;
template <int ndim> class Nbody ;
template <int ndim> struct Particle ;
namespace Regularization {
  template <int ndim> class RegularizerFunction ;
}



//=================================================================================================
//  Class Ic
/// \brief   Class containing functions for generating initial conditions.
/// \details Class containing functions for generating initial conditions.
/// \author  D. A. Hubber, G. Rosotti
/// \date    02/02/2015
//=================================================================================================
template <int ndim>
class Ic
{
protected:

  Simulation<ndim>* const sim;               ///< Simulation class pointer
  Hydrodynamics<ndim>* const hydro;          ///< Hydrodynamics algorithm pointer
  Nbody<ndim>* const nbody;                  ///< Nbody algorithm pointer
  NeighbourSearch<ndim>* const neib;         ///< Neighbour search algorithm pointer
  const FLOAT invndim;                       ///< 1/ndim
  const SimUnits& simunits;                  ///< Reference to main simunits object
  const DomainBox<ndim>& icBox;              ///< Reference to IC box
  const DomainBox<ndim>& simbox;             ///< Reference to simulation bounding box object

  Parameters* simparams;                     ///< Pointer to parameters object
  RandomNumber *randnumb;                    ///< Random number object pointer


  // Helper routines
  //-----------------------------------------------------------------------------------------------
  void AddAzimuthalDensityPerturbation(const int, const int, const FLOAT, const FLOAT *, FLOAT *);
  void AddBinaryStar(FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT,
                     FLOAT *, FLOAT *, NbodyParticle<ndim> &, NbodyParticle<ndim> &);
  void AddRotationalVelocityField(const int, const FLOAT, const FLOAT *, const FLOAT *, FLOAT *);
  void Addr2Sphere(int, FLOAT *, FLOAT *, FLOAT);
  void AddSinusoidalDensityPerturbation(int, FLOAT, FLOAT, FLOAT *);
  void GenerateTurbulentVelocityField(const int, const int, const DOUBLE, DOUBLE *);
  void InterpolateVelocityField(const int, const int, const FLOAT, const FLOAT,
                                const FLOAT *, const DOUBLE *, FLOAT *);

  void AddMonteCarloDensityField(const int, const int, const Box<ndim> &, FLOAT *, RandomNumber *);
  virtual FLOAT GetMaximumDensity(const int, const Box<ndim> &, RandomNumber *);


public:

  // Constructor for setting important variables and pointers
  //-----------------------------------------------------------------------------------------------
  Ic(Simulation<ndim>* _sim, FLOAT _invndim) :
    sim(_sim), hydro(_sim->hydro), nbody(_sim->nbody), neib(_sim->neib), icBox(_sim->icBox),
    invndim(_invndim), simunits(_sim->simunits), simbox(_sim->simbox),
    simparams(_sim->simparams), randnumb(_sim->randnumb)
  {
  };
  virtual ~Ic() {};


  // Virtual functions
  //-----------------------------------------------------------------------------------------------
  virtual void Generate(void) {
    ExceptionHandler::getIstance().raise("Ic::Generate function not implemented");
  }
  virtual FLOAT GetDensity(const FLOAT *, const int) const {
    ExceptionHandler::getIstance().raise("Ic::GetDensity function not implemented");
    return (FLOAT) 0.0;
  }
  virtual void SetParticleProperties() {
    ExceptionHandler::getIstance().raise("Ic::SetParticleProperties function not implemented");
  }
  virtual Regularization::RegularizerFunction<ndim>* GetParticleRegularizer() const {
    ExceptionHandler::getIstance().raise("GetParticleRegularizer not implemented");
    return NULL;
  }


  // Other common functions
  //-----------------------------------------------------------------------------------------------
  void CheckInitialConditions(void);
  FLOAT CalculateMassInBox(const Box<ndim>&, const int);
  FLOAT GetSmoothedDensity(const FLOAT *, const int, const FLOAT, SmoothingKernel<ndim> *);


  // Static functions which can be used outside of Ic class
  // (e.g. generating new particles on the fly in simulations)
  //-----------------------------------------------------------------------------------------------
  static void AddCubicLattice(const int, const int *, const DomainBox<ndim> &, const bool, FLOAT *);
  static void AddHexagonalLattice(const int, const int *, const DomainBox<ndim> &, const bool, FLOAT *);
  static int AddLatticeSphere(const int, const FLOAT *, const FLOAT, const string, FLOAT *, RandomNumber *);
  static void AddRandomBox(const int, const DomainBox<ndim> &, FLOAT *, RandomNumber *);
  static void AddRandomSphere(const int, const FLOAT *, const FLOAT, FLOAT *, RandomNumber *);
  static int CutSphere(const int, const int, const DomainBox<ndim> &, const bool, FLOAT *);


};



//=================================================================================================
//  Class NullIc
/// \brief   ..
/// \details ..
/// \author  D. A. Hubber
/// \date    20/03/2016
//=================================================================================================
template <int ndim>
class NullIc : public Ic<ndim>
{
public:

  NullIc(Simulation<ndim>* _sim, FLOAT _invndim) :
    Ic<ndim>(_sim, _invndim) {};

};



//=================================================================================================
//  Class BinaryAccretionIc
/// \brief   Binary accretion simulation IC class.
/// \details Binary accretion simulation IC class.
/// \author  D. A. Hubber
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class BinaryAccretionIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  BinaryAccretionIc(Simulation<ndim>* _sim, FLOAT _invndim);
  ~BinaryAccretionIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class BlobIc
/// \brief   ...
/// \details ...
/// \author  ...
/// \date    09/12/2016
//=================================================================================================
template <int ndim>
class BlobIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  BlobIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~BlobIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class BondiAccretionIc
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    17/11/2016
//=================================================================================================
template <int ndim>
class BondiAccretionIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  BondiAccretionIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~BondiAccretionIc() {};

  virtual void Generate(void);
  void ComputeBondiSolution(int, FLOAT *, FLOAT *, FLOAT *, FLOAT *);

};



//=================================================================================================
//  Class BossBodenheimerIc
/// \brief   Boss-Bodenheimer simulation IC class.
/// \details Boss-Bodenheimer simulation IC class.
/// \author  D. A. Hubber
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class BossBodenheimerIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;

  FLOAT amp;                               ///< Azimuthal density perturbation amplitude
  FLOAT angvel;                            ///< Angular velocity of rotating cloud
  FLOAT mcloud;                            ///< Mass of cloud
  FLOAT radius;                            ///< Radius of cloud
  FLOAT rho0;                              ///< Average density of cloud
  FLOAT temp0;                             ///< Temperature of cloud
  FLOAT u0;                                ///< Specific internal energy of gas


public:

  BossBodenheimerIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~BossBodenheimerIc() {};

  virtual void Generate();
  virtual FLOAT GetDensity(const FLOAT *, const int) const;
  virtual void SetParticleProperties();
  virtual Regularization::RegularizerFunction<ndim>* GetParticleRegularizer() const;

};



//=================================================================================================
//  Class ContactDiscontinuityIc
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    19/11/2016
//=================================================================================================
template <int ndim>
class ContactDiscontinuityIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  ContactDiscontinuityIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~ContactDiscontinuityIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class DustyBoxIc
/// \brief   ...
/// \details ...
/// \author  ...
/// \date    16/11/2016
//=================================================================================================
template <int ndim>
class DustyBoxIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  DustyBoxIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~DustyBoxIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class EvrardCollapseIc
/// \brief   ...
/// \details ...
/// \author  ...
/// \date    16/11/2016
//=================================================================================================
template <int ndim>
class EvrardCollapseIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  EvrardCollapseIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~EvrardCollapseIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class EwaldIc
/// \brief   ...
/// \details ...
/// \author  ...
/// \date    16/11/2016
//=================================================================================================
template <int ndim>
class EwaldIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  EwaldIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~EwaldIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class FilamentIc
/// \brief   Class to generate a simple filament for ICs.
/// \details Class to generate a simple filament for ICs.
/// \author  D. A. Hubber
/// \date    12/10/2016
//=================================================================================================
template <int ndim>
class FilamentIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;

  FLOAT aconst;
  FLOAT n0;
  FLOAT r0;
  FLOAT rho0;
  FLOAT Rfilament;
  FLOAT Lfilament;
  FLOAT temp0;
  FLOAT v_cyl_infall;
  FLOAT v_rad_infall;

  FLOAT mtot;
  FLOAT u0;
  FLOAT mp;


public:

  FilamentIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~FilamentIc() {};

  virtual void Generate();
  virtual FLOAT GetDensity(const FLOAT *, const int) const;
  virtual void SetParticleProperties();
  virtual Regularization::RegularizerFunction<ndim>* GetParticleRegularizer() const ;

};



//=================================================================================================
//  Class GaussianRingIc
/// \brief   ...
/// \details ...
/// \author  ...
/// \date    09/12/2016
//=================================================================================================
template <int ndim>
class GaussianRingIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  GaussianRingIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~GaussianRingIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class GreshoVortexIc
/// \brief   Class to generate initial conditions for 2d Gresho vortex-type simulations.
/// \details Class to generate initial conditions for 2d Gresho vortex-type simulations.
/// \author  D. A. Hubber
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class GreshoVortexIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  GreshoVortexIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~GreshoVortexIc() {};

  virtual void Generate(void);
  virtual FLOAT GetValue(const std::string, const FLOAT *);

};



//=================================================================================================
//  Class HierarchicalSystemIc
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    18/11/2016
//=================================================================================================
template <int ndim>
class HierarchicalSystemIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  HierarchicalSystemIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~HierarchicalSystemIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class IsothermalSphereIc
/// \brief   Class to generate ...
/// \details Class to generate ...
/// \author  ...
/// \date    17/11/2016
//=================================================================================================
template <int ndim>
class IsothermalSphereIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  IsothermalSphereIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~IsothermalSphereIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class KelvinHelmholtzIc
/// \brief   Class to generate initial conditions for 2d Kelvin-Helmholtz instabilty.
/// \details Class to generate initial conditions for a 2d Kelvin-Helmholtz instabilty.
/// \author  D. A. Hubber
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class KelvinHelmholtzIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  KelvinHelmholtzIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~KelvinHelmholtzIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class NohIc
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    17/11/2016
//=================================================================================================
template <int ndim>
class NohIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  NohIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~NohIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class PlummerSphereIc
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    17/11/2016
//=================================================================================================
template <int ndim>
class PlummerSphereIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  PlummerSphereIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~PlummerSphereIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class PolytropeIc
/// \brief   Class to generate initial conditions with a Polytrope density profile
/// \details Class to generate initial conditions with a Polytrope density profile
/// \author  D. A. Hubber
/// \date    22/06/2016
//=================================================================================================
template <int ndim>
class PolytropeIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;

  int Ntablemax;
  FLOAT *xiArray;
  FLOAT *psiArray;
  FLOAT *phiArray;
  FLOAT *muArray;
  FLOAT *pressArray;
  FLOAT *rhoArray;
  FLOAT *thetaArray;
  FLOAT *massArray;


public:

  PolytropeIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~PolytropeIc();

  virtual void Generate(void);
  virtual FLOAT GetValue(const std::string, const FLOAT *);
  static void ComputeIsothermalLaneEmdenSolution(const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *);
  static void ComputeLaneEmdenSolution(const int, const FLOAT, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *);


};



//=================================================================================================
//  Class RayleighTaylorIc
/// \brief   Class to generate Rayleigh-Taylor instability initial conditions.
/// \details Class to generate Rayleigh-Taylor instability initial conditions.
/// \author  D. A. Hubber
/// \date    18/11/2016
//=================================================================================================
template <int ndim>
class RayleighTaylorIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  RayleighTaylorIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~RayleighTaylorIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class SedovBlastwaveIc
/// \brief   Class to generate Sedov-Taylor blastwave initial conditions.
/// \details Class to generate Sedov-Taylor blastwave initial conditions.
/// \author  D. A. Hubber
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class SedovBlastwaveIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  SedovBlastwaveIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~SedovBlastwaveIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class ShearflowIc
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    19/11/2016
//=================================================================================================
template <int ndim>
class ShearflowIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  ShearflowIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~ShearflowIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class ShocktubeIc
/// \brief   Class to generate shocktube problem initial conditions.
/// \details Class to generate shocktube problem initial conditions.
/// \author  D. A. Hubber
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class ShocktubeIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  ShocktubeIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~ShocktubeIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class SilccIc
/// \brief   TODO: NEEDS A MORE GENERAL DESCRIPTION
/// \details TODO: NEEDS A MORE GENERAL DESCRIPTION
/// \author  D. A. Hubber, S. Walch
/// \date    20/03/2016
//=================================================================================================
template <int ndim>
class SilccIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;

  FLOAT a_midplane;                    // Location of midplane
  FLOAT box_area;                      // Area of x-y plane of simulation box
  FLOAT h_midplane;                    // Scale-height of midplane density region
  FLOAT m_box;                         // Total gas mass in box
  FLOAT m_exp;                         // Total gas mass in exponential profile region
  FLOAT m_uniform;                     // Total gas mass in uniform density region
  FLOAT rho_a;                         // Density at edge of exponential midplane profile
  FLOAT rho_midplane;                  // Density at the midplane
  FLOAT rho_star;                      // Stellar density at the midplane
  FLOAT sigma_star;                    // Surface density of stars
  FLOAT temp0;                         // Temperature
  FLOAT u0;                            // ??
  FLOAT z_d;                           // ??


public:

  SilccIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~SilccIc() {};

  virtual void Generate();
  virtual FLOAT GetDensity(const FLOAT *, const int) const;
  virtual void SetParticleProperties();
  virtual Regularization::RegularizerFunction<ndim>* GetParticleRegularizer() const ;

};



//=================================================================================================
//  Class SoundwaveIc
/// \brief   Class to generate simple 1d soundwave initial conditions.
/// \details Class to generate simple 1d soundwave initial conditions.
/// \author  D. A. Hubber, S. Walch
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class SoundwaveIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  SoundwaveIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~SoundwaveIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class SpitzerExpansionIc
/// \brief   ...
/// \details ...
/// \author  D. A. Hubber
/// \date    19/11/2016
//=================================================================================================
template <int ndim>
class SpitzerExpansionIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  SpitzerExpansionIc(Simulation<ndim>* _sim, FLOAT _invndim);
  ~SpitzerExpansionIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class TurbulentCoreIc
/// \brief   Class to generate simple turbulent core ICs.
/// \details Class to generate simple turbulent core ICs.
/// \author  D. A. Hubber
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class TurbulentCoreIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  TurbulentCoreIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~TurbulentCoreIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class UniformIc
/// \brief
/// \details
/// \author  D. A. Hubber
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class UniformIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  UniformIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~UniformIc() {};

  virtual void Generate(void);

};


//=================================================================================================
//  Class CMZIc
/// \brief
/// \details
/// \author  J. E. Dale, G. P. Rosotti
/// \date    04/05/2017
//=================================================================================================
template <int ndim>
class CMZIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::icBox;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  CMZIc(Simulation<ndim>* _sim, FLOAT _invndim);
  virtual ~CMZIc() {};

  virtual void Generate(void);

};



namespace Regularization {

// Base class for regularization funciton interface
template<int ndim>
class RegularizerFunction {
public:
  virtual ~RegularizerFunction() { } ;

  virtual FLOAT operator()(const Particle<ndim>&, const Particle<ndim>&) const  = 0 ;
  //virtual FLOAT GetSmoothedDensity(const FLOAT *, const FLOAT) = 0;
};


// Standard regularizer based upon density and particle disorder
template<int ndim, class DensityFunc>
class DefaultRegularizerFunction : public RegularizerFunction<ndim>
{
  const bool regularise_smooth_density;
  const FLOAT alphaReg, rhoReg;
  const DensityFunc* _rho_func;
  SmoothingKernel<ndim>* _kern;

public:
  DefaultRegularizerFunction(SmoothingKernel<ndim> *kern, Parameters* simparams, const DensityFunc* rho)
  : regularise_smooth_density(simparams->intparams["regularise_smooth_density"]),
    alphaReg(simparams->floatparams["alpha_reg"]),
    rhoReg(simparams->floatparams["rho_reg"]),
    _rho_func(rho),
    _kern(kern)
  {} ;

  virtual ~DefaultRegularizerFunction() {};


  virtual FLOAT operator()(const Particle<ndim>& part,
                           const Particle<ndim>& neibpart) const
  {
    FLOAT rhotrue;

    if (regularise_smooth_density) {
      const FLOAT h       = neibpart.h;
      const int ptype     = neibpart.ptype;
      const int gridSize  = 5;
      const FLOAT hrange  = _kern->kernrange*h;
      const FLOAT invh    = (FLOAT) 1.0/h;
      const FLOAT invhsqd = invh*invh;
      const FLOAT hfactor = pow(invh, ndim);
      const FLOAT dV      = pow(2.0*hrange/(FLOAT) gridSize, ndim);
      FLOAT sumValue      = (FLOAT) 0.0;
      FLOAT dr[ndim];
      FLOAT drsqd;
      FLOAT r[ndim];

      //-------------------------------------------------------------------------------------------
      if (ndim == 1) {
        for (int i=0; i<gridSize; i++) {
          dr[0] = ((FLOAT) i + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
          drsqd = dr[0]*dr[0];
          r[0]  = neibpart.r[0] + dr[0];
          sumValue += hfactor*_kern->w0_s2(drsqd*invhsqd)*_rho_func->GetDensity(r, ptype);
        }
      }
      //-------------------------------------------------------------------------------------------
      else if (ndim == 2) {
        for (int i=0; i<gridSize; i++) {
          for (int j=0; j<gridSize; j++) {
            dr[0] = ((FLOAT) i + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
            dr[1] = ((FLOAT) j + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
            r[0]  = neibpart.r[0] + dr[0];
            r[1]  = neibpart.r[1] + dr[1];
            drsqd = dr[0]*dr[0] + dr[1]*dr[1];
            sumValue += hfactor*_kern->w0_s2(drsqd*invhsqd)*_rho_func->GetDensity(r, ptype);
          }
        }
      }
      //-------------------------------------------------------------------------------------------
      else if (ndim == 3) {
        for (int i=0; i<gridSize; i++) {
          for (int j=0; j<gridSize; j++) {
            for (int k=0; k<gridSize; k++) {
              dr[0] = ((FLOAT) i + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
              dr[1] = ((FLOAT) j + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
              dr[2] = ((FLOAT) k + (FLOAT) 0.5)*2.0*hrange/(FLOAT) gridSize - hrange;
              r[0]  = neibpart.r[0] + dr[0];
              r[1]  = neibpart.r[1] + dr[1];
              r[2]  = neibpart.r[2] + dr[2];
              drsqd = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
              sumValue += hfactor*_kern->w0_s2(drsqd*invhsqd)*_rho_func->GetDensity(r, ptype);
            }
          }
        }
      }
      //-------------------------------------------------------------------------------------------

      // Normalise summation/integral with the volume of each individual summation point
      sumValue *= dV;
      rhotrue = sumValue;
    }
    else {
      rhotrue = _rho_func->GetDensity(neibpart.r, neibpart.ptype);
    }

    FLOAT rhofrac = (neibpart.rho - rhotrue)/(rhotrue + small_number);
    rhofrac = std::min(std::max(rhofrac, -(FLOAT) 0.1), (FLOAT) 10.0);

    return rhoReg*rhofrac + alphaReg;
  }

};



//=================================================================================================
//  Class ParticleRegularizer
/// \brief Regularizes the particle distribution based upon the provided regularization criteria.
//=================================================================================================
template<int ndim>
class ParticleRegularizer
{
  DomainBox<ndim> localBox;
  int Nreg;


public:
  ParticleRegularizer(Parameters* simparams, const DomainBox<ndim>& icbox);

  void operator()(Hydrodynamics<ndim>* hydro, NeighbourSearch<ndim> *neib, Nbody<ndim>* nbody,
                  const RegularizerFunction<ndim>& regularizer) const;
};

} // namespace Regularization


#endif
