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
  void ComputeBondiSolution(int, FLOAT *, FLOAT *, FLOAT *, FLOAT *);
  void GenerateTurbulentVelocityField(const int, const int, const DOUBLE, DOUBLE *);
  void InterpolateVelocityField(const int, const int, const FLOAT, const FLOAT,
                                const FLOAT *, const DOUBLE *, FLOAT *);


public:

  Ic(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim) :
    sim(_sim), hydro(_hydro), invndim(_invndim),
    simunits(_sim->simunits), simbox(_sim->simbox),
    simparams(_sim->simparams), randnumb(_sim->randnumb)
  {
  };

  virtual void Generate(void) {};
  virtual FLOAT GetValue(std::string, FLOAT *) {return (FLOAT) 0.0;}


  // Initial conditions routines
  //-----------------------------------------------------------------------------------------------
  void BinaryAccretion(void);
  void BinaryStar(void);
  void BlastWave(void);
  void BondiAccretion(void);
  void BossBodenheimer(void);
  void CheckInitialConditions(void);
  void ContactDiscontinuity(void);
  void EwaldDensity(void);
  void GaussianRing(void);
  void GreshoVortex(void);
  void KHI(void);
  void NohProblem(void);
  void PlummerSphere(void);
  void QuadrupleStar(void);
  void RTI(void);
  void SedovBlastWave(void);
  void ShearFlow(void);
  void ShockTube(void);
  void Silcc(void);
  void SoundWave(void);
  void SpitzerExpansion(void);
  void TripleStar(void);
  void BlobTest(void);
  void TurbulentCore(void);
  void UniformBox(void);
  void UniformSphere(void);
  void IsothermSphere(void);
  void RotIsothermSphere(void);
  void TurbIsothermSphere(void);
  void EvrardCollapse(void);
  void DustyBox(void);

  // Static functions which can be used outside of Ic class
  // (e.g. generating new particles on the fly in simulations)
  //-----------------------------------------------------------------------------------------------
  static void AddCubicLattice(const int, const int *, const DomainBox<ndim> &, const bool, FLOAT *);
  static void AddHexagonalLattice(const int, const int *, const DomainBox<ndim> &, const bool, FLOAT *);
  static int AddLatticeSphere(const int, const FLOAT *, const FLOAT, const string, FLOAT *, RandomNumber *);
  static void AddRandomBox(const int, const DomainBox<ndim>, FLOAT *, RandomNumber *);
  static void AddRandomSphere(const int, const FLOAT *, const FLOAT, FLOAT *, RandomNumber *);
  static int CutSphere(const int, const int, const DomainBox<ndim>, const bool, FLOAT *);

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

  virtual void Generate(void) {};
  virtual FLOAT GetValue(std::string, FLOAT *) {return (FLOAT) 0.0;}

};


static FLOAT sech(FLOAT arg) {
  return (FLOAT) 2.0 / (exp(arg) + exp(-arg));
}

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

  FLOAT a_midplane;                    // ..
  FLOAT box_area;                      // Area of x-y plane of simulation box
  FLOAT h_midplane;                    // ..
  FLOAT m_box;                         // Total gas mass in box
  FLOAT m_exp;                         // Total gas mass in exponential profile region
  FLOAT m_uniform;                     // Total gas mass in uniform density region
  FLOAT rho_a;                         // Density at edge of exponential midplane profile
  FLOAT rho_midplane;                  // ..
  FLOAT rho_star;                      // Stellar density at the midplane
  FLOAT sigma_star;                    // ..
  FLOAT z_d;                           // ..


public:


  SilccIc(Simulation<ndim>* _sim, Hydrodynamics<ndim>* _hydro, FLOAT _invndim);
  virtual ~SilccIc() {};

  virtual void Generate(void);
  virtual FLOAT GetValue(std::string, FLOAT *);

};
#endif
