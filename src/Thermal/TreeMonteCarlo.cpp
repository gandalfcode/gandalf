//=============================================================================
//  TreeMonteCarlo.cpp
//  Class for controlling Monte-Carlo radiation transport algorithms walking
//  a KD-tree constructed from the gas positions.
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


#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <math.h>
#include "Precision.h"
#include "Radiation.h"
#include "Debug.h"
using namespace std;



//=============================================================================
//  TreeMonteCarlo::TreeMonteCarlo()
/// Constructor for KD-tree radiation class
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
TreeMonteCarlo<ndim,nfreq,ParticleType,CellType>::TreeMonteCarlo
(int Nphotonaux, int Nleafmaxaux, RandomNumber *randaux)
{
  Nphoton = Nphotonaux;
  randnumb = randaux;
  radtree = new KDRadiationTree<ndim,nfreq,ParticleType,CellType>(Nleafmaxaux);
}



//=============================================================================
//  TreeMonteCarlo::~TreeMonteCarlo()
/// Destructor for KD-tree radiation class
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
TreeMonteCarlo<ndim,nfreq,ParticleType,CellType>::~TreeMonteCarlo()
{
}



//=============================================================================
//  TreeMonteCarlo::UpdateRadiationField
/// Update the radiation field by iterating the MCRT method on the tree
/// several times until acceptable convergence has been achieved.
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void TreeMonteCarlo<ndim,nfreq,ParticleType,CellType>::UpdateRadiationField
(int Nhydro,                          ///< No. of SPH particle
 int Nnbody,                        ///< No. of N-body particles
 int Nsink,                         ///< No. of sink particles
 SphParticle<ndim> *sph_gen,        ///< Generic SPH particle data array
 NbodyParticle<ndim> **nbodydata,   ///< N-body data array
 SinkParticle<ndim> *sinkdata)      ///< Sink data array
{
  int c;                            // Cell counter
  int k;                            // Dimension counter
  int level;                        // Level to walk tree on
  RadiationSource<ndim> source;     // Current radiation source
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* >(sph_gen);


  debug2("[TreeMonteCarlo::UpdateRadiationField]");


  // Re-build radiation tree from scratch
  timing->StartTimingSection("RADTREE_BUILD");
  radtree->BuildTree(Nhydro,Nhydro,sphdata);
  timing->EndTimingSection("RADTREE_BUILD");

  // Compute maximum radius of cell extent from co-ordinate centre, in order
  // to know where to emit external photons from
  boundaryradius = 0.0;
  for (k=0; k<ndim; k++) boundaryradius +=
    max(pow(radtree->radcell[0].bbmax[k] - radtree->radcell[0].rcell[k],2),
	pow(radtree->radcell[0].rcell[k] - radtree->radcell[0].bbmin[k],2));
  boundaryradius = 1.1*sqrt(boundaryradius);


  // Perform single iteration step on computing radiation field
  timing->StartTimingSection("RADTREE_MONTE_CARLO");
  level = radtree->ltot;
  IterateRadiationField(level,Nhydro,Nnbody,Nsink,sph_gen,nbodydata,sinkdata);
  timing->EndTimingSection("RADTREE_MONTE_CARLO");


  // Update the radiation field on all higher cells
  radtree->SumRadiationField(level,radtree->radcell[0]);


  // Normalise photon energy density for all cells
  //for (c=0; c<radtree->Ncell; c++)
  //radtree->radcell[c].uphoton =
  //  radtree->radcell[c].lsum/radtree->radcell[c].volume;


  // Output info to file for plotting
  for (int l=0; l<level; l++) {

    ofstream outfile;
    string filename;
    string nostring;
    stringstream ss;
    nostring = "";
    ss << setfill('0') << setw(2) << l;
    nostring = ss.str();

    filename = "uphoton-l" + nostring + ".dat";
    outfile.open(filename.c_str());
    for (c=0; c<radtree->Ncell; c++) {
      if (radtree->radcell[c].level != l ||
	  radtree->radcell[c].N == 0) continue;
      outfile << sqrt(DotProduct(radtree->radcell[c].rcell,
				 radtree->radcell[c].rcell,ndim)) << "   "
	      << radtree->radcell[c].lsum[0]/radtree->radcell[c].volume
	      << "    " << radtree->radcell[c].Nphoton << "    "
	      << radtree->radcell[c].N << "    "
	      << radtree->radcell[c].volume << endl;
    }
    outfile.close();

  }

  return;
}



//=============================================================================
//  TreeMonteCarlo::IterateRadiationField
/// Single iteration step for updating the radiation field using the
/// tree-MCRT method.
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void TreeMonteCarlo<ndim,nfreq,ParticleType,CellType>::IterateRadiationField
(int level,                         ///< Level to walk tree on
 int Nhydro,                          ///< No. of SPH particle
 int Nnbody,                        ///< No. of N-body particles
 int Nsink,                         ///< No. of sink particles
 SphParticle<ndim> *sph_gen,        ///< Generic SPH particle data array
 NbodyParticle<ndim> **nbodydata,   ///< N-body data array
 SinkParticle<ndim> *sinkdata)      ///< Sink data array
{
  int iphoton;                      // Photon counter
  int k;                            // Dimension counter
  int kfreq;                        // Frequency bin counter
  long int Ncellcount = 0;          // Count no. of cells passed through
  long int Nscattercount = 0;       // No. of scattering events
  FLOAT dpath;                      // Path to next cell boundary
  FLOAT rand1;                      // Random number
  FLOAT taumax;                     // Optical depth travelled by photon
  FLOAT tau;                        // Current value of photon optical depth
  PhotonPacket<ndim> photon;        // Current photon packet
  RadiationSource<ndim> source;     // Current radiation source

  if (level < 0 && level > radtree->ltot) {
    cout << "Walking invalid tree level" << endl;
    exit(0);
  }

  // Emit photon packets from single source (for now)
  source.sourcetype = "pointsource";
  for (k=0; k<ndim; k++) source.r[k] = 0.0;
  source.c = radtree->FindCell(0,level,source.r);
  source.luminosity = 1.0;
  packetenergy = source.luminosity/(FLOAT) Nphoton;
  kfreq = 0;


  // Now emit all photons from radiation sources, updating the radiation field
  //===========================================================================
#pragma omp parallel for default(none) reduction(+:Ncellcount,Nscattercount) \
  shared(kfreq,level,source) private(dpath,iphoton,k,photon,rand1,tau,taumax)
  for (iphoton=0; iphoton<Nphoton; iphoton++) {

    // Initialise new photon packet from single source (modify later)
    photon = GenerateNewPhotonPacket(source);

    // Calculate optical depth to be travelled by photon
    rand1 = randnumb->floatrand();
    taumax = -log(rand1);
    tau = 0.0;
    Nscattercount++;


    // Main photon transmission/scattering/absorption-reemission iteration loop
    //-------------------------------------------------------------------------
    do {

      // Increase cell counter
      Ncellcount++;

      // Find i.d. of next (parent) cell the path length in current cell
      photon.cnext = radtree->FindRayExitFace(radtree->radcell[photon.c],
                                              photon.r,photon.eray,photon.inveray,dpath);

      // Check if maximum optical depth has been reached in order to
      // scatter or absorb/re-emit photon.
      //-----------------------------------------------------------------------
      if (tau + dpath*radtree->radcell[photon.c].opacity[kfreq] > taumax) {

        // Propagate photon packet until absorption/scattering event
        dpath = (taumax - tau)/radtree->radcell[photon.c].opacity[kfreq];
        for (k=0; k<ndim; k++) photon.r[k] += dpath*photon.eray[k];
#pragma omp atomic
        radtree->radcell[photon.c].lsum[kfreq] += packetenergy*dpath;
#pragma omp atomic
        radtree->radcell[photon.c].Nphoton++;

        // Scatter photon (isotropic scattering for now)
        ScatterPhotonPacket(photon);

        // Calculate new optical depth to be travelled by scattered photon
        rand1 = randnumb->floatrand();
        taumax = -log(rand1);
        tau = 0.0;
        Nscattercount++;

      }

      // Otherwise, photon continues through cell and exits to adjacent cell
      //-----------------------------------------------------------------------
      else {

        // Propagate photon packet to edge of cell and add contribution to
        // radiation field of cell
        for (k=0; k<ndim; k++) photon.r[k] += dpath*photon.eray[k];
#pragma omp atomic
        radtree->radcell[photon.c].lsum[kfreq] += packetenergy*dpath;
#pragma omp atomic
        radtree->radcell[photon.c].Nphoton++;
        //radtree->radcell[photon.c].uphoton += packetenergy*dpath;
        tau += dpath*radtree->radcell[photon.c].opacity[kfreq];


#ifdef OUTPUT_ALL
        cout << "Found path length : " << dpath << "     cnext : " << photon.cnext << endl;
        cout << "Photon exitting cell at : " << photon.r[0] << "   "
             << photon.r[1] << "   " << photon.r[2] << endl;
        cout << "Checking distance : " << photon.r[0]/photon.eray[0] << "   "
             << photon.r[1]/photon.eray[1] << "   " << photon.r[2]/photon.eray[2] << endl;
#endif

        // Exit loop if we've reached the edge of the computational domain
        if (photon.cnext == -1) break;

        // Find i.d. of next cell from the parent cell
        photon.c = radtree->FindAdjacentCell(photon.cnext,level,photon.r);

      }
      //-----------------------------------------------------------------------


    } while (photon.c != -1);
    //-------------------------------------------------------------------------


  }
  //===========================================================================


  // Normalise photon energy density for all cells
  //for (c=0; c<radtree->Ncell; c++)
  //radtree->radcell[c].uphoton =
  //  radtree->radcell[c].lsum/radtree->radcell[c].volume;


#ifdef OUTPUT_ALL
  cout << "Radiation field : " << radtree->Ntot << "   " << radtree->Ntotmax
       << "   " << radtree->Ncell << "   " << radtree->ltot
       << "   " << radtree->radcell[0].volume << endl;
#endif
  cout << "No. of photons propagated       : " << Nphoton << endl;
  cout << "Total no. of cells crossed      : " << Ncellcount << endl;
  cout << "Average no. of cells per photon : " << Ncellcount/Nphoton
       << "    " << pow(radtree->gtot,0.33333333) << endl;
  cout << "Average no. of scatter events   : " << Nscattercount/Nphoton
       << "    " << pow(radtree->radcell[0].opacity[kfreq],2) << endl;

  return;
}



//=================================================================================================
//  TreeMonteCarlo::GenerateNewPhotonPacket
/// Create a new photon packet for a given radiation source.
/// Currently only implements point sources.
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
PhotonPacket<ndim> TreeMonteCarlo<ndim,nfreq,ParticleType,CellType>::GenerateNewPhotonPacket
 (RadiationSource<ndim> &source)       ///< [in] Source for generating new packet
{
  int k;                               // Dimension counter
  FLOAT theta;                         // Random angle for photon direction
  PhotonPacket<ndim> photon;           // Photon packet type


  // Initialise new photon packet from chosen source
  // (e.g. pick random direction, frequency, etc..)
  //-----------------------------------------------------------------------------------------------
  if (source.sourcetype == "pointsource") {
    photon.c = source.c;
    photon.cnext = source.c;
    photon.energy = packetenergy;
    for (k=0; k<ndim; k++) photon.r[k] = source.r[k];

    // Generate random direction for photon
    theta = pi*(2.0*randnumb->floatrand() - 1.0);
    if (ndim == 3) {
      photon.eray[2] = 2.0*randnumb->floatrand() - 1.0;
      photon.eray[0] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*cos(theta);
      photon.eray[1] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*sin(theta);
    }
    for (k=0; k<ndim; k++) photon.inveray[k] = 1.0/(photon.eray[k] + small_number);

  }

  // Isotropic source
  //-----------------------------------------------------------------------------------------------
  else if (source.sourcetype == "isotropic") {
    cout << "Isotropic radiation field not yet implemented" << endl;
    exit(0);

  }

  // Isotropic source
  //-----------------------------------------------------------------------------------------------
  else if (source.sourcetype == "planar") {
    cout << "Planar radiation field not yet implemented" << endl;
    exit(0);

  }
  //-----------------------------------------------------------------------------------------------


#ifdef OUTPUT_ALL
  cout << "Emitting photon with direction " << photon.eray[0] << "   "
       << photon.eray[1] << "   " << photon.eray[2] << endl;
#endif


  return photon;
}



//=============================================================================
//  TreeMonteCarlo::ScatterPhotonPacket
/// Scatter photon packet inot random, isotropic direction.
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void TreeMonteCarlo<ndim,nfreq,ParticleType,CellType>::ScatterPhotonPacket
(PhotonPacket<ndim> &photon)        ///< [inout] Reference to photon packet
{
  int k;                            // Dimension counter
  FLOAT theta;                      // Random angle for photon direction

  // Generate random direction for photon
  //theta = pi*(2.0*((FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX) - 1.0);
  //photon.eray[2] = 2.0*((FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX) - 1.0;
  theta = pi*(2.0*randnumb->floatrand() - 1.0);
  if (ndim == 3) {
    photon.eray[2] = 2.0*randnumb->floatrand() - 1.0;
    photon.eray[0] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*cos(theta);
    photon.eray[1] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*sin(theta);
  }
  for (k=0; k<ndim; k++) photon.inveray[k] = 1.0/(photon.eray[k] + small_number);

  return;
}



template class TreeMonteCarlo<1,1,GradhSphParticle,KDRadTreeCell>;
template class TreeMonteCarlo<2,1,GradhSphParticle,KDRadTreeCell>;
template class TreeMonteCarlo<3,1,GradhSphParticle,KDRadTreeCell>;
//template class TreeMonteCarlo<1,SM2012SphParticle>;
//template class TreeMonteCarlo<2,SM2012SphParticle>;
//template class TreeMonteCarlo<3,SM2012SphParticle>;
