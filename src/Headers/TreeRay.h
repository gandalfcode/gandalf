//=================================================================================================
//  TreeRay.h
//  ...
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


#ifndef _TREE_RAY_H_
#define _TREE_RAY_H_


#include <map>
#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <ostream>
#include "CodeTiming.h"
#include "Constants.h"
#include "Debug.h"
#include "EnergyEquation.h"
#include "EOS.h"
#include "DomainBox.h"
#include "KDRadiationTree.h"
#include "Nbody.h"
#include "Particle.h"
#include "Precision.h"
#include "Parameters.h"
#include "Radiation.h"
#include "RandomNumber.h"
#include "SimUnits.h"
#include "Sinks.h"
#include "SmoothingKernel.h"
#include "SphNeighbourSearch.h"
#include "Particle.h"
#include "NbodyParticle.h"
#include "OctTree.h"
using namespace std;


/*static const int nfreq=1;


//=================================================================================================
//  Struct TreeRayCell
/// TreeRay cell data structure
//=================================================================================================
template <int ndim>
struct TreeRayCell : public OctTreeCell<ndim>
{
  FLOAT volume;
};



//=================================================================================================
//  Struct ..
/// TreeRay cell data structure
//=================================================================================================
//template <int ndim, int nfreq>
template <int ndim>
struct OsTreeRayCell : public TreeRayCell<ndim>
{
  FLOAT srcEUV[nfreq];                        ///< Source function of EUV radiation
  FLOAT erdEUVold[nfreq];                     ///< ..
  FLOAT erdEUV[nfreq];                        ///< Radiation energy density
};
*/


//=================================================================================================
//  Struct ..
/// ...
//=================================================================================================
//template <int nfreq>
struct Rays {
  FLOAT mass;                          ///< ..
  FLOAT volume;                        ///< ..
  FLOAT rho;                           ///< ..
  FLOAT srcF[nfreq];                   ///< ..
  FLOAT erad[nfreq];                   ///< Radiation energy DENSITY
  FLOAT Erad[nfreq];                   ///< Radiation energy contained in volume
};



//=================================================================================================
//  Class TreeRay
/// \brief
/// \details
/// \author  R. Wunsch & D. A. Hubber
/// \date    16/03/2015
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int> class TreeCell>
class TreeRay : public Radiation<ndim>
{
 public:

#define IIL(i,iNS,iPhi,iTheta) (i + iNS*ilNI + iPhi*ilNI*ilNNS + iTheta*ilNI*ilNNS*(ilNPhi + 1))
#define IRNM(iR,iFR,iL) (iR + iFR*bhNR + iL*bhNR*bhNR*nFineR)
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))

  static const int nFineR = 10;        // Richard's MAGIC number!!
  static const int NL = 40;            // No. of levels in tree (replace later)

  // Constant parameters
  const bool onTheSpot;                ///< 'On the spot' approximation
  const int nSide;                     ///< HEALPix nside
  const int ilNR;                      ///< Intersection list no. of points along ray
  const int ilNTheta;                  ///< No. of points along theta direction
  const int ilNPhi;                    ///< ""                  phi
  const int ilNNS;                     ///< ""                  angular node size
  const long int ilFinePix;            ///< Extra HEALPix levels for fine ??
  const FLOAT maxDist;                 ///< Maximum distance
  const FLOAT rayRadRes;               ///< Radial ray resolution (relative to node size)
  const FLOAT relErr;                  ///< Iteration error tolerance
  const string errControl;             ///< Iteration control method/mode

  int bhNR;                            ///< ..
  int ilNI;                            ///< ..
  int NTBLevels;                       ///< ..
  long int nPix;                       ///< ..
  FLOAT bhLocRelErr;                   ///< ..
  FLOAT bhMaxRelEradErr;               ///< ..
  FLOAT bhLocEradTot;                  ///< ..
  FLOAT bhLocMionTot;                  ///< ..
  FLOAT ilNSSampFac;                   ///< ???
  FLOAT ilNSSampFacI;                  ///< ..
  FLOAT max_ray_length;                ///< ..
  FLOAT minCellSize;                   ///< ..

  int *radNodeMapIndex;                ///< ..
  FLOAT *intersectList;                ///< ..
  FLOAT *radNodeMapValue;              ///< ..
  FLOAT *rayR;                         ///< ..
  FLOAT *rayR2;                        ///< ..
  FLOAT *rayRi;                        ///< ..
  FLOAT *rayR2i;                       ///< ..

  OctTree<ndim,ParticleType,TreeCell> *tree;   ///< ..
  Rays **rays;                         ///< ..


  //-----------------------------------------------------------------------------------------------
  TreeRay(bool, int, int, int, int, int, int, FLOAT, FLOAT, FLOAT, string,
          DomainBox<ndim> &, SimUnits *, Parameters *);
  ~TreeRay();

  virtual void UpdateRadiationField(int, int, int, SphParticle<ndim> *,
                                    NbodyParticle<ndim> **, SinkParticle<ndim> *);

  void AddRadiationSource(SinkParticle<ndim> &);
  void FixRay(Rays *ray);
  void FinaliseCell(TreeCell<ndim> &, FLOAT **eflux, FLOAT **cdMaps);
  void FinaliseIteration(void);
  void GenerateIntersectList(void);
  void GenerateRadNodeMapping(void);
  void IntegrateRay(Rays *ray, FLOAT eflux[nfreq]);
  FLOAT NodeKernel(const FLOAT, const FLOAT, const FLOAT);
  void NodeContribution(TreeCell<ndim> &targetNode, TreeCell<ndim> &contributingNode,
                        FLOAT dr[ndim], FLOAT drsqd, FLOAT cellSize);
  void RadToGas(TreeCell<ndim> &, ParticleType<ndim> *, FLOAT timestep);
  void TreeWalk(TreeCell<ndim> &);
  bool TreeWalkEnd(void);

};
#endif
