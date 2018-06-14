//=================================================================================================
//  Radiation.h
//  Contains definitions for all classes that control the transport of
//  radiation through the computational domain.
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


#ifndef _RADIATION_H_
#define _RADIATION_H_


#include <map>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <ostream>
#include "CodeTiming.h"
#include "EnergyEquation.h"
#include "EOS.h"
#include "DomainBox.h"
#include "KDRadiationTree.h"
#include "Nbody.h"
#include "Precision.h"
#include "Parameters.h"
#include "RandomNumber.h"
#include "SimUnits.h"
#include "Sinks.h"
#include "SmoothingKernel.h"
#include "SphNeighbourSearch.h"
#include "Particle.h"
#include "NbodyParticle.h"
using namespace std;



//=================================================================================================
//  enum radiationSourceType
/// ...
//=================================================================================================
enum radiationSourceType {
  isotropicSource,
  pointSource,
  planarSource,
  numRadiationSourceTypes
};



//=================================================================================================
//  Struct PhotonPacket
/// Radiation photon packet data structure
//=================================================================================================
template <int ndim>
struct PhotonPacket {
  int c;                               ///< Current cell of photon packet
  int cnext;                           ///< Next cell for photon
  FLOAT energy;                        ///< Total energy carried by packet
  FLOAT r[ndim];                       ///< Position of ray
  FLOAT eray[ndim];                    ///< Unit vector direction of ray
  FLOAT inveray[ndim];                 ///< 1/eray
};



//=================================================================================================
//  Struct RadiationSource
/// Radiation photon packet data structure
//=================================================================================================
template <int ndim>
struct RadiationSource {
  string sourcetype;                   ///< Type of radiation source
  int c;                               ///< i.d. of cell containing source
  FLOAT luminosity;                    ///< Source luminosity
  FLOAT r[ndim];                       ///< Position of radiation source
  FLOAT esource[ndim];                 ///< Unit vector of radiation from source
                                       ///< (for uni-directional sources)
};



//=================================================================================================
//  Struct ionpar
/// ..
//=================================================================================================
struct ionpar
{
  int sink;                         // Is particle sink
  int fionised;
  int neighstorcont;
  double x;                         // Particle x,y,z co-ordinates
  double y;
  double z;
  double rho;                       // Density
  double t;                         // Temperature
  double h;                         // Smoothing length
  double u;                         // Specific internal energy
  int *checked;
  int *ionised;                     // Is particle ionised by source?
  int *neigh;                       // Part. neib array (Neibs closest to sources)
  int *neighstor;
  double *angle;
  double *prob;                     // Prob. of transmition from each source
  double *photons;                  // No. of photons lost up to this point
  double *rad_pre_acc;
};



//=================================================================================================
//  Class Radiation
/// \brief   Main base radiation class
/// \details Main base radiation class from which child classes containing
///          implementations are inherited.
/// \author  D. A. Hubber
/// \date    21/04/2014
//=================================================================================================
template <int ndim>
class Radiation
{
 public:

  Radiation() {};
  ~Radiation() {};

  virtual void UpdateRadiationField(int, int, int, SphParticle<ndim> *,
                                    NbodyParticle<ndim> **, SinkParticle<ndim> *) = 0;

  CodeTiming *timing;                  ///< Pointer to code timing object
  SimUnits *units;                     ///< Pointer to code units object

};



//=================================================================================================
//  Class MultipleSourceIonisiation
/// \brief   Radiation Scheme to treat ionising radiation for multiple sources
/// \details ..
/// \author  S. K. Balfour
/// \date    24/04/2014
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class MultipleSourceIonisation : public Radiation<ndim>
{
 public:

  MultipleSourceIonisation(SphNeighbourSearch<ndim> *, float, float,
                           float, float, double, float, float, float, double);
  //MultipleSourceIonisation(SphNeighbourSearch<ndim> *,float,float,
  //                         float,float,float,float,float,float);
  ~MultipleSourceIonisation();

  virtual void UpdateRadiationField(int, int, int, SphParticle<ndim> *,
                                    NbodyParticle<ndim> **, SinkParticle<ndim> *);

  void ionisation_intergration(int, int, NbodyParticle<ndim> **,
                               SphParticle<ndim> *, double, double,
                               SphNeighbourSearch<ndim> *, double, double,
                               double, double, double, double);
  void photoncount(ionpar *,int *,double *, int &,int &,int &,int &);
  double lost(ionpar *,int *,double *,int &, int &,int &,int &,int &);
  void probs(int &, ionpar *, int *, int &, double *);


  SphNeighbourSearch<ndim> *sphneib;
  float mu_bar,temp0,mu_ion,temp_ion,gamma_eos,scale,tempscale; //cmscott
  double rad_cont,Ndotmin; //cmscott
  //float mu_bar,temp0,mu_ion,temp_ion,Ndotmin,gamma_eos,scale,tempscale;
  vector< vector<int> > ionisation_fraction;
};



//=================================================================================================
//  Class NullRadiation
/// \brief   Empty radiation class when no radiation object is selected
/// \details Empty radiation class when no radiation object is selected
/// \author  D. A. Hubber
/// \date    21/04/2014
//=================================================================================================
template <int ndim>
class NullRadiation : public Radiation<ndim>
{
  using Radiation<ndim>::timing;

 public:

  NullRadiation(): Radiation<ndim>() {};

  virtual void UpdateRadiationField(int, int, int, SphParticle<ndim> *,
                                    NbodyParticle<ndim> **, SinkParticle<ndim> *) {};
};
#endif
