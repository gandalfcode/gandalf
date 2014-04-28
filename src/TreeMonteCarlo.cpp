//=============================================================================
//  TreeMonteCarlo.cpp
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
//=============================================================================


#include "Precision.h"
#include "Radiation.h"
#include "Debug.h"
using namespace std;



//=============================================================================
//  TreeMonteCarlo::TreeMonteCarlo()
/// Constructor for KD-tree radiation class
//=============================================================================
template <int ndim, template<int> class ParticleType>
TreeMonteCarlo<ndim,ParticleType>::TreeMonteCarlo(int Nphotonaux)
{
  Nphoton = Nphotonaux;
}



//=============================================================================
//  TreeMonteCarlo::~TreeMonteCarlo()
/// Destructor for KD-tree radiation class
//=============================================================================
template <int ndim, template<int> class ParticleType>
TreeMonteCarlo<ndim,ParticleType>::~TreeMonteCarlo()
{
}




//=============================================================================
//  TreeMonteCarlo::UpdateRadiationField
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
void TreeMonteCarlo<ndim,ParticleType>::UpdateRadiationField
(int Nsph,                          ///< No. of SPH particle
 int Nnbody,                        ///< No. of N-body particles
 int Nsink,                         ///< No. of sink particles
 SphParticle<ndim> *sph_gen,        ///< Generic SPH particle data array
 NbodyParticle<ndim> **nbodydata,   ///< N-body data array
 SinkParticle<ndim> *sinkdata)      ///< Sink data array
{
  int c;                            // Cell counter
  int iphoton;                      // Photon counter
  int isource;                      // Radiation source i.d.
  int k;                            // Dimension counter
  FLOAT dpath;                      // Path to next cell boundary
  FLOAT packetenergy;               // Energy carried by one photon packet
  PhotonPacket<ndim> photon;        // Current photon packet
  RadiationSource<ndim> source;     // Current radiation source
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* >(sph_gen);


  debug2("[TreeMonteCarlo::UpdateRadiationField]");
  
 
  // Re-build radiation tree from scratch
  radtree->BuildTree(Nsph,Nsph,sphdata);

  // Emit photon packets from single source (for now)
  for (k=0; k<ndim; k++) source.r[k] = 0.0;
  source.c = radtree->FindCell(0,source.r);
  source.luminosity = 1.0;
  packetenergy = source.luminosity/(FLOAT) Nphoton;


  // Now emit all photons from radiation sources, updating the radiation field
  //---------------------------------------------------------------------------
  for (iphoton=0; iphoton<Nphoton; iphoton++) {

    // Initialise new photon packet from chosen source 
    // (e.g. pick random direction, frequency, etc..)
    photon.c = source.c;
    photon.energy = packetenergy;
    for (k=0; k<ndim; k++) photon.r[k] = source.r[k];


    // Main photon transmission/scattering/absorption-reemission iteration loop
    //-------------------------------------------------------------------------
    do {

      // Find i.d. of next (parent) cell the path length in current cell
      photon.cnext = FindRayExitFace(radtree->radcell[photon.c],
                                     photon.r,photon.eray,dpath);

      // Propagate photon packet to edge of cell and add contribution to 
      // radiation field of cell
      for (k=0; k<ndim; k++) photon.r[k] += dpath*photon.eray[k];
      radtree->radcell[photon.c].uphoton += packetenergy*dpath;

      // Exit loop if we've reached the edge of the computational domain
      if (photon.cnext == -1) break;

      // Find i.d. of next cell from the parent cell
      photon.cnext = FindAdjacentCell(photon.cnext,photon.r);


    } while (photon.c != -1);
    //-------------------------------------------------------------------------
    

  }
  //---------------------------------------------------------------------------


  // Normalise photon energy density for all cells
  for (c=0; c<radtree->Ncell; c++)
    radtree->radcell[c].uphoton /= radtree->radcell[c].volume;


  return;
} 



//=============================================================================
//  TreeMonteCarlo::FindRayExitFace
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
int TreeMonteCarlo<ndim,ParticleType>::FindRayExitFace
(KDRadTreeCell<ndim> &cell,         ///< [in] Reference to cell
 FLOAT rp[ndim],                    ///< [in] Position of point/ray
 FLOAT eray[ndim],                  ///< [in] Unit vector direction of ray
 FLOAT &dpath)                      ///< [out] Length of ray path across cell
{
  int cexit;                        // i.d. of cell that ray is travelling to
  int k;                            // Dimension counter
  FLOAT daux;                       // Aux. value of face-intersection distance
  
  // Initialise variables before finding face
  dpath = big_number;

  // Check each cell boundary
  //---------------------------------------------------------------------------
  for (k=0; k<ndim; k++) {

    // If radiation if travelling in +ve k-direction
    if (eray[0] > 0.0) {
      daux = (cell.bbmax[k] - rp[k])/eray[k];
      if (daux < dpath) {
        dpath = daux;
        cexit = cell.cexit[1][k];
      }
    }

    // If radiation if travelling in -ve k-direction
    else {
      daux = (cell.bbmin[k] - rp[k])/eray[k];
      if (daux < dpath) {
        dpath = daux;
        cexit = cell.cexit[0][k];
      }
    }

  }
  //---------------------------------------------------------------------------

  return cexit;
}



//=============================================================================
//  TreeMonteCarlo::FindAdjacentCell
/// Find i.d. of cell adjacent to current cell that the radiation packet is 
/// travelling into.
//=============================================================================
template <int ndim, template<int> class ParticleType>
int TreeMonteCarlo<ndim,ParticleType>::FindAdjacentCell
(int cparent,                       ///< [in] i.d. of larger parent cell
 FLOAT rp[ndim])                    ///< [in] Position of point/ray
{
  int c = cparent;                  // Cell i.d.
  int c1;                           // i.d. of 1st cell child
  int k_divide;                     // Dimension of cell division

  // Walk back down through tree to bottom level
  //---------------------------------------------------------------------------
  while (radtree->radcell[c].level > radtree->ltot) {
    c1 = c + 1;
    k_divide = radtree->radcell[c].k_divide;

    // If point is left of divide, pick 1st child cell.  Else pick 2nd child.
    if (rp[k_divide] < radtree->radcell[c1].bbmax[k_divide])
      c = c1;
    else
      c = radtree->radcell[c].c2;

  };
  //---------------------------------------------------------------------------

  return c;
}



template class TreeMonteCarlo<1,GradhSphParticle>;
template class TreeMonteCarlo<2,GradhSphParticle>;
template class TreeMonteCarlo<3,GradhSphParticle>;
template class TreeMonteCarlo<1,SM2012SphParticle>;
template class TreeMonteCarlo<2,SM2012SphParticle>;
template class TreeMonteCarlo<3,SM2012SphParticle>;
template class TreeMonteCarlo<1,GodunovSphParticle>;
template class TreeMonteCarlo<2,GodunovSphParticle>;
template class TreeMonteCarlo<3,GodunovSphParticle>;
