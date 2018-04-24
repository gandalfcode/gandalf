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
#include "OctTree.h"
#include "Precision.h"
#include "Parameters.h"
#include "RandomNumber.h"
#include "SimUnits.h"
#include "Sinks.h"
#include "SmoothingKernel.h"
#include "SphNeighbourSearch.h"
#include "Tree.h"
#include "Particle.h"
#include "NbodyParticle.h"
using namespace std;


// Compiler flag for optimised version (in terms of memory allocation) for Scott's ionisation algorithm
// (Comment out for now due to OpenMP bug, but will fix asap).
//#define RAD_OPTIMISE


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
#ifdef RAD_OPTIMISE
static const int maxSources=16;
#endif
struct ionpar
{
  int sink;                         // Is particle sink
  int fionised;
  int neighstorcont;
  DOUBLE x;                         // Particle x,y,z co-ordinates
  DOUBLE y;
  DOUBLE z;
  DOUBLE rho;                       // Density
  DOUBLE t;                         // Temperature
  DOUBLE h;                         // Smoothing length
  DOUBLE u;                         // Specific internal energy
#ifdef RAD_OPTIMISE
  int checked[maxSources];
  int ionised[maxSources];                     // Is particle ionised by source?
  int neigh[maxSources];                       // Part. neib array (Neibs closest to sources)
  int neighstor[200];
  DOUBLE angle[maxSources];
  DOUBLE prob[maxSources];                     // Prob. of transmition from each source
  DOUBLE photons[maxSources];                  // No. of photons lost up to this point
  DOUBLE rad_pre_acc[3];
#else
  int *checked;
  int *ionised;                     // Is particle ionised by source?
  int *neigh;                       // Part. neib array (Neibs closest to sources)
  int *neighstor;
  DOUBLE *angle;
  DOUBLE *prob;                     // Prob. of transmition from each source
  DOUBLE *photons;                  // No. of photons lost up to this point
  DOUBLE *rad_pre_acc;
#endif
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

  virtual void UpdateRadiationField(int, int, int, Particle<ndim> *,
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

  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  MultipleSourceIonisation(NeighbourSearch<ndim> *, FLOAT, FLOAT, FLOAT,
                           FLOAT, FLOAT, DOUBLE, FLOAT, FLOAT, FLOAT, FLOAT, DOUBLE);
  ~MultipleSourceIonisation();


  // Function prototypes
  //-----------------------------------------------------------------------------------------------
  virtual void UpdateRadiationField(int, int, int, Particle<ndim> *,
                                    NbodyParticle<ndim> **, SinkParticle<ndim> *);

  void ionisation_intergration(int, int, NbodyParticle<ndim> **, Particle<ndim> *,
                               DOUBLE, DOUBLE, NeighbourSearch<ndim> *,
                               DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE, DOUBLE);
  void photoncount(ionpar *, int *, DOUBLE *, int &, int &, int &, int &);
  DOUBLE lost(ionpar *, int *, DOUBLE *, int &, int &, int &, int &, int &);
  void probs(int &, ionpar *, int *, int &, DOUBLE *);


  // Variables
  //-----------------------------------------------------------------------------------------------
  FLOAT mu_bar;
  FLOAT X_comp;
  FLOAT temp0;
  FLOAT mu_ion;
  FLOAT temp_ion;
  FLOAT gamma_eos;
  FLOAT arecomb;
  FLOAT scale;
  FLOAT tempscale;
  DOUBLE rad_cont;
  DOUBLE Ndotmin;
  vector< vector<int> > ionisation_fraction;
  NeighbourSearch<ndim> *neib;
};



//=================================================================================================
//  Class TreeMonteCarlo
/// \brief   Class that propagates radiation through KD-tree
/// \details Class that propagates radiation through KD-tree
/// \author  D. A. Hubber, A. P. Whitworth
/// \date    25/04/2014
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
class TreeMonteCarlo : public Radiation<ndim>
{
public:
  using Radiation<ndim>::timing;


  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  TreeMonteCarlo(int, int, RandomNumber *);
  ~TreeMonteCarlo();


  // Function prototypes
  //-----------------------------------------------------------------------------------------------
  virtual void UpdateRadiationField(int, int, int, Particle<ndim> *,
                                    NbodyParticle<ndim> **, SinkParticle<ndim> *) ;
  void IterateRadiationField(int, int, int, int, Particle<ndim> *,
                             NbodyParticle<ndim> **, SinkParticle<ndim> *) ;
  PhotonPacket<ndim> GenerateNewPhotonPacket(RadiationSource<ndim> &);
  void ScatterPhotonPacket(PhotonPacket<ndim> &);


  // Variables
  //-----------------------------------------------------------------------------------------------
  int Nphoton;                                  // No. of photon packets
  FLOAT boundaryradius;                         // Radius from which isotropic
                                                // photons are emitted.
  FLOAT packetenergy;                           // Energy in photon packet
  RandomNumber *randnumb;                       // Random number object pointer

  KDRadiationTree<ndim,nfreq,ParticleType,CellType> *radtree;  // Rad. tree

};



//=================================================================================================
//  Class MonochromaticIonisationMonteCarlo
/// \brief   Class that propagates radiation through KD-tree
/// \details Class that propagates radiation through KD-tree
/// \author  D. A. Hubber, A. P. Whitworth
/// \date    25/04/2014
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
class MonochromaticIonisationMonteCarlo : public Radiation<ndim>
{
  using Radiation<ndim>::timing;

 public:


  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  MonochromaticIonisationMonteCarlo(int, int, int, FLOAT, FLOAT, FLOAT, DOUBLE,
                                    string, SimUnits *, EOS<ndim> *);
  ~MonochromaticIonisationMonteCarlo();


  // Function prototypes
  //-----------------------------------------------------------------------------------------------
  virtual void UpdateRadiationField(int, int, int, Particle<ndim> *,
                                    NbodyParticle<ndim> **, SinkParticle<ndim> *);
  void InterpolateParticleProperties(const int, const int, Particle<ndim> *);
  void IterateRadiationField(const int, const int, const int, const int, const int,
                             Particle<ndim> *, NbodyParticle<ndim> **, SinkParticle<ndim> *);
  PhotonPacket<ndim> GenerateNewPhotonPacket(const RadiationSource<ndim> &, RandomNumber *);
  void ScatterPhotonPacket(PhotonPacket<ndim> &, RandomNumber *);
  bool UpdateIonisationFraction(const int, const int);
  void UpdateCellOpacity(CellType<ndim,nfreq> &, ParticleType<ndim> *);


  // Variables
  //-----------------------------------------------------------------------------------------------
  FLOAT Nphotonratio;                  // Ratio of photons to radiation cells
  int Nraditerations;                  // No. of iterations of radiation field
  int Nradlevels;                      // No. of tree levels to converge over
  int Nthreads;                        // No. of OpenMP threads
  FLOAT boundaryradius;                // Radius from which isotropic photons are emitted.
  FLOAT across;                        // Photoionisation cross-section
  FLOAT arecomb;                       // Recombination coefficient
  FLOAT Eion;                          // Ionisation energy
  FLOAT invEion;                       // 1 / Eion
  //FLOAT packetenergy;                  // Energy in photon packet
  FLOAT temp_ion;                      // ..
  DOUBLE invmh;                        // ..
  DOUBLE ionconst;                     // ..
  DOUBLE NLyC;                         // No. of ionising photons per second
  RandomNumber **randNumbArray;        // Random number object array (for OpenMP threads)
  SimUnits *units;                     // ..
  EOS<ndim> *eos;                      // Pointer to main EOS object

  KDRadiationTree<ndim,nfreq,ParticleType,CellType> *radtree;  // Rad. tree

};



//===============================================================================================
//  Class PlaneParallelRadiation
/// \brief   Radiation photon packet data structure
/// \author  D. A. Hubber
/// \date    11/12/2017
//===============================================================================================
template <int ndim>
class PlanarRay
{
public:

  int level;                                   ///< Level of ray (with respect to OctTree level)
  FLOAT rayIntegral;                           ///< Ionisation integral
  FLOAT r[ndim];                               ///< Current position of ray
  FLOAT dir[ndim];                             ///< Radiation propogation direction

  PlanarRay()
  {
    level = 0;
    rayIntegral = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) r[k] = (FLOAT) 0.0;
    for (int k=0; k<ndim; k++) dir[k] = (FLOAT) 0.0;
  }

  inline void operator= (const PlanarRay<ndim> &other)
  {
    level = other.level;
    rayIntegral = other.rayIntegral;
    for (int k=0; k<ndim; k++) r[k] = other.r[k];
    for (int k=0; k<ndim; k++) dir[k] = other.dir[k];
  }

};



//=================================================================================================
//  Class PlaneParallelRadiation
/// \brief   ...
/// \author  D. A. Hubber
/// \date    11/12/2017
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class PlaneParallelRadiation : public Radiation<ndim>
{
private:

  int kx, ky, kz;
  int minRayDivisions;
  int numIonised;
  int radiationDirection;
  FLOAT arecomb;
  FLOAT FLyC;
  FLOAT maxIntegral;
  FLOAT rayStepMult;
  FLOAT uion;
  FLOAT h_fac;
  FLOAT invndim;

  CodeTiming *timing;
  SimUnits *units;
  SmoothingKernel<ndim> *kern;
  NeighbourSearch<ndim> *neib;
  OctTree<ndim,ParticleType,OctTreeCell> *tree;

  int DivideRay(PlanarRay<ndim> &, PlanarRay<ndim> *);


public:

  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  PlaneParallelRadiation(Parameters *, SmoothingKernel<ndim> *, SimUnits*,
                         NeighbourSearch<ndim> *, CodeTiming *);
  ~PlaneParallelRadiation();


  // Function prototypes
  //-----------------------------------------------------------------------------------------------
  virtual void UpdateRadiationField(int, int, int, Particle<ndim> *,
                                    NbodyParticle<ndim> **, SinkParticle<ndim> *);

  FLOAT CalculateDensity(const int, ParticleType<ndim> *, FLOAT *, FLOAT &, FLOAT &);
  bool CellRayIntegration(const OctTreeCell<ndim> &, const FLOAT, const int, FLOAT &, PlanarRay<ndim> &, ParticleType<ndim> *);
  bool CheapCellRayIntegration(const OctTreeCell<ndim> &, const FLOAT, const int, FLOAT &, PlanarRay<ndim> &, ParticleType<ndim> *);
  void CreateRootRay(const OctTreeCell<ndim> &, const FLOAT, PlanarRay<ndim> &);
  bool ExpensiveCellRayIntegration(const OctTreeCell<ndim> &, const FLOAT, const int, FLOAT &, PlanarRay<ndim> &, ParticleType<ndim> *);
  void SplitRay(const PlanarRay<ndim> &, const FLOAT, int &, PlanarRay<ndim> *);

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

  virtual void UpdateRadiationField(int, int, int, Particle<ndim> *,
                                    NbodyParticle<ndim> **, SinkParticle<ndim> *) {};
};
#endif
