//=================================================================================================
//  IC.h
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


#include "Precision.h"
#include "Hydrodynamics.h"
#include "InlineFuncs.h"
#include "Simulation.h"
#include "RandomNumber.h"
#if defined(FFTW_TURBULENCE)
#include "fftw3.h"
#endif



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

  Simulation<ndim>* const sim;              ///< Simulation class pointer
  Hydrodynamics<ndim>* const hydro;         ///< Hydrodynamics algorithm pointer
  const FLOAT invndim;                      ///< 1/ndim
  const SimUnits& simunits;                 ///< Reference to main simunits object
  const DomainBox<ndim>& simbox;            ///< Reference to simulation bounding box object

  int Ntable;                               ///< No. of table elements
  FLOAT *xTable;                            ///< Tabulated position values
  FLOAT *mTable;                            ///< Tabulated integrated mass values
  FLOAT *mFracTable;                        ///< Tabulated fractional integrated mass values
  std::string posQuantity;                  ///< Position quantity string (e.g. x, y, r)
  Parameters* simparams;                    ///< Pointer to parameters object
  RandomNumber *randnumb;                   ///< Random number object pointer


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


public:

  // Constructor for setting important variables and pointers
  //-----------------------------------------------------------------------------------------------
  Ic(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim) :
    sim(_sim), hydro(_hydro), invndim(_invndim),
    simunits(_sim->simunits), simbox(_sim->simbox),
    simparams(_sim->simparams), randnumb(_sim->randnumb)
  {
  };
  virtual ~Ic() {};


  // Virtual functions
  //-----------------------------------------------------------------------------------------------
  virtual void Generate(void) {};
  virtual FLOAT GetValue(const std::string, const FLOAT *) {return (FLOAT) 0.0;}
  virtual FLOAT GetDensity(const FLOAT rad) {return (FLOAT) 0.0;}


  // Other common functions
  //-----------------------------------------------------------------------------------------------
  void CheckInitialConditions(void);
  void CalculateMassTable(std::string, FLOAT, FLOAT);
  FLOAT CalculateMassInBox(const int*, const Box<ndim>&);
  FLOAT FindMassIntegratedPosition(FLOAT);



  // Static functions which can be used outside of Ic class
  // (e.g. generating new particles on the fly in simulations)
  //-----------------------------------------------------------------------------------------------
  static void AddCubicLattice(const int, const int *, const DomainBox<ndim> &, const bool, FLOAT *);
  static void AddHexagonalLattice(const int, const int *, const DomainBox<ndim> &, const bool, FLOAT *);
  static int AddLatticeSphere(const int, const FLOAT *, const FLOAT, const string, FLOAT *, RandomNumber *);
  static void AddRandomBox(const int, const DomainBox<ndim> &, FLOAT *, RandomNumber *);
  static void AddRandomSphere(const int, const FLOAT *, const FLOAT, FLOAT *, RandomNumber *);
  static void ComputeIsothermalLaneEmdenSolution(const int, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *);
  static void ComputeLaneEmdenSolution(const int, const FLOAT, const FLOAT, FLOAT *, FLOAT *, FLOAT *, FLOAT *);
  static int CutSphere(const int, const int, const DomainBox<ndim> &, const bool, FLOAT *);



  // Initial conditions routines
  //-----------------------------------------------------------------------------------------------
  void ContactDiscontinuity(void);
  void GaussianRing(void);
  void SedovBlastWave(void);
  void BlobTest(void);
  void UniformBox(void);
  void UniformSphere(void);

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

  NullIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim) :
    Ic<ndim>(_sim, _hydro, _invndim) {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  BinaryAccretionIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~BinaryAccretionIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  BondiAccretionIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~BondiAccretionIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  BossBodenheimerIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~BossBodenheimerIc() {};

  virtual void Generate(void);

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  ContactDiscontinuityIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~ContactDiscontinuityIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  DustyBoxIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~DustyBoxIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  EvrardCollapseIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~EvrardCollapseIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  EwaldIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~EwaldIc() {};

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

  FLOAT mtot;
  FLOAT u0;
  FLOAT mp;


public:

  FilamentIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~FilamentIc() {};

  virtual void Generate(void);
  virtual FLOAT GetValue(const std::string, const FLOAT *);

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  GreshoVortexIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~GreshoVortexIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  HierarchicalSystemIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~HierarchicalSystemIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  IsothermalSphereIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~IsothermalSphereIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class KhiIc
/// \brief   Class to generate initial conditions for 2d Kelvin-Helmholtz instabilty.
/// \details Class to generate initial conditions for a 2d Kelvin-Helmholtz instabilty.
/// \author  D. A. Hubber
/// \date    15/11/2016
//=================================================================================================
template <int ndim>
class KhiIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  KhiIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~KhiIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  NohIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~NohIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  PlummerSphereIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~PlummerSphereIc() {};

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

  PolytropeIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~PolytropeIc();

  virtual void Generate(void);
  virtual FLOAT GetValue(const std::string, const FLOAT *);

};



//=================================================================================================
//  Class RtiIc
/// \brief   Class to generate Rayleigh-Taylor instability initial conditions.
/// \details Class to generate Rayleigh-Taylor instability initial conditions.
/// \author  D. A. Hubber
/// \date    18/11/2016
//=================================================================================================
template <int ndim>
class RtiIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  RtiIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~RtiIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  SedovBlastwaveIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~SedovBlastwaveIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  ShearflowIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~ShearflowIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  ShocktubeIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~ShocktubeIc() {};

  virtual void Generate(void);

};



//=================================================================================================
//  Class SilccIc
/// \brief   Class to generate SILCC-like initial conditions
/// \details Class to generate SILCC-like initial conditions
/// \author  D. A. Hubber, S. Walch
/// \date    20/03/2016
//=================================================================================================
template <int ndim>
class SilccIc : public Ic<ndim>
{
protected:

  using Ic<ndim>::hydro;
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

  SilccIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~SilccIc() {};

  virtual void Generate(void);
  virtual FLOAT GetValue(const std::string, const FLOAT *);

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  SoundwaveIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~SoundwaveIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  SpitzerExpansionIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
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
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  TurbulentCoreIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~TurbulentCoreIc() {};

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
  using Ic<ndim>::invndim;
  using Ic<ndim>::randnumb;
  using Ic<ndim>::sim;
  using Ic<ndim>::simbox;
  using Ic<ndim>::simparams;
  using Ic<ndim>::simunits;


public:

  UniformIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  ~UniformIc() {};

  virtual void Generate(void);

};
#endif
