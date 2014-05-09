//=============================================================================
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
//=============================================================================


#ifndef _RADIATION_H_
#define _RADIATION_H_


#include <map>
#include <string>
#include <list>
#include <iostream>
#include <ostream>
#include "CodeTiming.h"
#include "EnergyEquation.h"
#include "DomainBox.h"
#include "KDRadiationTree.h"
#include "Nbody.h"
#include "Precision.h"
#include "Parameters.h"
#include "SimUnits.h"
#include "Sinks.h"
#include "SphKernel.h"
#include "SphNeighbourSearch.h"
#include "SphParticle.h"
using namespace std;



//=============================================================================
//  Struct PhotonPacket
/// Radiation photon packet data structure
//=============================================================================
template <int ndim>
struct PhotonPacket {
  int c;                            ///< Current cell of photon packet
  int cnext;                        ///< Next cell for photon
  FLOAT energy;                     ///< Total energy carried by packet
  FLOAT r[ndim];                    ///< Position of ray
  FLOAT eray[ndim];                 ///< Unit vector direction of ray
  FLOAT inveray[ndim];              ///< 1/eray
};



//=============================================================================
//  Struct RadiationSource
/// Radiation photon packet data structure
//=============================================================================
template <int ndim>
struct RadiationSource {
  string sourcetype;                ///< Type of radiation source
  int c;                            ///< i.d. of cell containing source
  FLOAT luminosity;                 ///< Source luminosity
  FLOAT r[ndim];                    ///< Position of radiation source
  FLOAT esource[ndim];              ///< Unit vector of radiation from source 
                                    ///< (for uni-directional sources)
};



//=============================================================================
//  Class Radiation
/// \brief   Main base radiation class
/// \details Main base radiation class from which child classes containing 
///          implementations are inherited.
/// \author  D. A. Hubber
/// \date    21/04/2014
//=============================================================================
template <int ndim>
class Radiation
{
 public:

  Radiation() {};
  ~Radiation() {};

  virtual void UpdateRadiationField(int, int, int, SphParticle<ndim> *, 
                                    NbodyParticle<ndim> **, 
                                    SinkParticle<ndim> *) = 0;

  CodeTiming *timing;               ///< Pointer to code timing object

};



//=============================================================================
//  Class TreeMonteCarlo
/// \brief   Class that propagates radiation through KD-tree
/// \details Class that propagates radiation through KD-tree
/// \author  D. A. Hubber, A. P. Whitworth
/// \date    25/04/2014
//=============================================================================
template <int ndim, template<int> class ParticleType>
class TreeMonteCarlo : public Radiation<ndim>
{
  using Radiation<ndim>::timing;

 public:


  TreeMonteCarlo(int, int);
  ~TreeMonteCarlo();


  virtual void UpdateRadiationField(int, int, int, SphParticle<ndim> *, 
                                    NbodyParticle<ndim> **, 
                                    SinkParticle<ndim> *) ;
  PhotonPacket<ndim> GenerateNewPhotonPacket(RadiationSource<ndim> &);
  int FindRayExitFace(KDRadTreeCell<ndim> &, FLOAT *, 
                      FLOAT *, FLOAT *, FLOAT &);
  int FindAdjacentCell(int, FLOAT *);
  void ScatterPhotonPacket(PhotonPacket<ndim> &);


  int Nphoton;                                  // No. of photon packets
  FLOAT packetenergy;                           // Energy in photon packet
  KDRadiationTree<ndim,ParticleType> *radtree;  // Radiation tree

};




//=============================================================================
//  Class NullRadiation
/// \brief   Empty radiation class when no radiation object is selected
/// \details Empty radiation class when no radiation object is selected
/// \author  D. A. Hubber
/// \date    21/04/2014
//=============================================================================
template <int ndim>
class NullRadiation : public Radiation<ndim>
{
  using Radiation<ndim>::timing;

 public:

  NullRadiation():Radiation<ndim>() {};


  virtual void UpdateRadiationField(int, int, int, SphParticle<ndim> *, 
                                    NbodyParticle<ndim> **, 
                                    SinkParticle<ndim> *) {};
};
#endif
