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
/// \brief   Class to contain all generic initial conditions functions.
/// \details Class to contain all generic initial conditions functions.
/// \author  D. A. Hubber & G. Rosotti
/// \date    15/04/2013
//=================================================================================================
template <int ndim>
class Ic
{
private:

  Simulation<ndim>* const sim;              ///< Simulation class pointer
  Hydrodynamics<ndim>* const hydro;         ///< Hydrodynamics algorithm pointer
  const FLOAT invndim;                      ///< 1/ndim
  const SimUnits& simunits;                 ///< Reference to main simunits object
  const DomainBox<ndim>& simbox;            ///< Reference to simulation bounding box object

  Parameters* simparams;                    ///< Pointer to parameters object
  RandomNumber *randnumb;                   ///< Random number object pointer


  // Helper routines
  //-----------------------------------------------------------------------------------------------
  void AddAzimuthalDensityPerturbation(int, int, FLOAT, FLOAT *, FLOAT *);
  void AddBinaryStar(DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE,
                     DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE *,DOUBLE *,
                     NbodyParticle<ndim> &, NbodyParticle<ndim> &);
  void AddCubicLattice(int, int *, FLOAT *, DomainBox<ndim>, bool);
  void AddHexagonalLattice(int, int *, FLOAT *, DomainBox<ndim>, bool);
  int AddLatticeSphere(int, FLOAT *, FLOAT *, FLOAT, string);
  void AddRotationalVelocityField(int, FLOAT, FLOAT *, FLOAT *, FLOAT *);
  void AddRandomBox(int, FLOAT *, DomainBox<ndim>);
  void AddRandomSphere(int, FLOAT *, FLOAT *, FLOAT);
  void Addr2Sphere(int, FLOAT *, FLOAT *, FLOAT);
  void AddSinusoidalDensityPerturbation(int, FLOAT, FLOAT, FLOAT *);
  int CutSphere(int, int, FLOAT *, DomainBox<ndim>, bool);
  void ComputeBondiSolution(int, FLOAT *, FLOAT *, FLOAT *, FLOAT *);
  void GenerateTurbulentVelocityField(int, int, DOUBLE, DOUBLE *);
  void InterpolateVelocityField(int, int, FLOAT, FLOAT, FLOAT *, FLOAT *,FLOAT *);


public:

  Ic(Simulation<ndim>* sim_aux, Hydrodynamics<ndim>* hydro_aux, FLOAT invndim_aux) :
    sim(sim_aux), hydro(hydro_aux), invndim(invndim_aux),
    simunits(sim_aux->simunits), simbox(sim_aux->simbox),
    simparams(sim_aux->simparams), randnumb(sim_aux->randnumb)
  {
  };


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
  void GreshoVortex(void);
  void KHI(void);
  void NohProblem(void);
  void PlummerSphere(void);
  void QuadrupleStar(void);
  void RTI(void);
  void ShockTube(void);
  void SedovBlastWave(void);
  void ShearFlow(void);
  void SoundWave(void);
  void SpitzerExpansion(void);
  void TripleStar(void);
  void TurbulentCore(void);
  void UniformBox(void);
  void UniformSphere(void);
  void IsothermSphere(void);
  void RotIsothermSphere(void);
  void TurbIsothermSphere(void);

};
#endif
