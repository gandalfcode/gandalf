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


#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
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
template <int ndim, template<int> class ParticleType>
TreeMonteCarlo<ndim,ParticleType>::TreeMonteCarlo
(int Nphotonaux, int Nleafmaxaux)
{
  Nphoton = Nphotonaux;
  radtree = new KDRadiationTree<ndim,ParticleType>(Nleafmaxaux);
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
  long int Ncellcount = 0;          // Count no. of cells passed through
  long int Nscattercount = 0;       // No. of scattering events
  FLOAT dpath;                      // Path to next cell boundary
  FLOAT randnumb;                   // Random number
  FLOAT taumax;                     // Optical depth travelled by photon
  FLOAT tau;                        // Current value of photon optical depth
  FLOAT theta;                      // Random angle for photon direction
  PhotonPacket<ndim> photon;        // Current photon packet
  RadiationSource<ndim> source;     // Current radiation source
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* >(sph_gen);

  ofstream outfile;
  string filename = "uphoton.dat";


  debug2("[TreeMonteCarlo::UpdateRadiationField]");
 

  // Re-build radiation tree from scratch
  timing->StartTimingSection("RADTREE_BUILD",2);
  radtree->BuildTree(Nsph,Nsph,sphdata);
  timing->EndTimingSection("RADTREE_BUILD");

  timing->StartTimingSection("TREE_MONTE_CARLO",2);


  // Emit photon packets from single source (for now)
  source.sourcetype = "pointsource";
  for (k=0; k<ndim; k++) source.r[k] = 0.0;
  source.c = radtree->FindCell(0,source.r);
  source.luminosity = 1.0;
  packetenergy = source.luminosity/(FLOAT) Nphoton;

  cout << "Source located at : " << source.c << "   " 
       << radtree->radcell[source.c].r[0] << "   " 
       << radtree->radcell[source.c].r[1] << "   " 
       << radtree->radcell[source.c].r[2] << endl;
 

  // Now emit all photons from radiation sources, updating the radiation field
  //===========================================================================
#pragma omp parallel for default(none) reduction(+:Ncellcount,Nscattercount) \
  shared(source) private(dpath,iphoton,k,photon,randnumb,tau,taumax) 
  for (iphoton=0; iphoton<Nphoton; iphoton++) {

    // Initialise new photon packet from single source (modify later)
    photon = GenerateNewPhotonPacket(source);
    
    // Calculate optical depth to be travelled by photon
    randnumb = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
    taumax = -log(randnumb);
    tau = 0.0;
    Nscattercount++;


    // Main photon transmission/scattering/absorption-reemission iteration loop
    //-------------------------------------------------------------------------
    do {

      // Increase cell counter
      Ncellcount++;

      // Find i.d. of next (parent) cell the path length in current cell
      photon.cnext = FindRayExitFace(radtree->radcell[photon.c],photon.r,
                                     photon.eray,photon.inveray,dpath);

      // Check if maximum optical depth has been reached in order to 
      // scatter or absorb/re-emit photon.
      //-----------------------------------------------------------------------
      if (tau + dpath*radtree->radcell[photon.c].opacity > taumax) {

        // Propagate photon packet until absorption/scattering event
        dpath = (taumax - tau)/radtree->radcell[photon.c].opacity;
        for (k=0; k<ndim; k++) photon.r[k] += dpath*photon.eray[k];
#pragma omp atomic
        radtree->radcell[photon.c].uphoton += packetenergy*dpath;

        // Scatter photon (isotropic scattering for now)
        ScatterPhotonPacket(photon);

	// Calculate new optical depth to be travelled by scattered photon
	randnumb = (FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX;
	taumax = -log(randnumb);
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
        radtree->radcell[photon.c].uphoton += packetenergy*dpath;
	tau += dpath*radtree->radcell[photon.c].opacity;


#ifdef OUTPUT_ALL
	cout << "Found path length : " << dpath << "     cnext : " 
	     << photon.cnext << endl;
	cout << "Photon exitting cell at : " << photon.r[0] << "   " 
	     << photon.r[1] << "   " << photon.r[2] << endl;
	cout << "Checking distance : " << photon.r[0]/photon.eray[0] 
	     << "   " << photon.r[1]/photon.eray[1] << "   " 
	     << photon.r[2]/photon.eray[2] << endl;
#endif

	// Exit loop if we've reached the edge of the computational domain
	if (photon.cnext == -1) break;
	
	// Find i.d. of next cell from the parent cell
	photon.c = FindAdjacentCell(photon.cnext,photon.r);

      }
      //-----------------------------------------------------------------------


    } while (photon.c != -1);
    //-------------------------------------------------------------------------
    

  }
  //===========================================================================


  // Normalise photon energy density for all cells
  for (c=0; c<radtree->Ncell; c++)
    radtree->radcell[c].uphoton /= radtree->radcell[c].volume;

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
       << "    " << pow(radtree->radcell[0].opacity,2) << endl;

  // Output info to file for plotting
  outfile.open(filename.c_str());
  for (c=0; c<radtree->Ncell; c++) {
    if (radtree->radcell[c].level != radtree->ltot || 
	radtree->radcell[c].N == 0) continue;
    outfile << sqrt(DotProduct(radtree->radcell[c].rcell,
			       radtree->radcell[c].rcell,ndim))
	    << "   " << radtree->radcell[c].uphoton << "    " 
	    << radtree->radcell[c].volume << "    " 
	    << radtree->radcell[c].N << endl;
  }
  outfile.close();

  timing->EndTimingSection("TREE_MONTE_CARLO");

  return;
} 



//=============================================================================
//  TreeMonteCarlo::GenerateNewPhotonPacket
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
PhotonPacket<ndim> TreeMonteCarlo<ndim,ParticleType>::GenerateNewPhotonPacket
(RadiationSource<ndim> &source)
{
  int k;                            // Dimension counter
  FLOAT theta;                      // Random angle for photon direction
  PhotonPacket<ndim> photon;        // Photon packet type
  
  
  // Initialise new photon packet from chosen source 
  // (e.g. pick random direction, frequency, etc..)
  //---------------------------------------------------------------------------
  if (source.sourcetype == "pointsource") {
    photon.c = source.c;
    photon.cnext = source.c;
    photon.energy = packetenergy;
    for (k=0; k<ndim; k++) photon.r[k] = source.r[k];
    
    // Generate random direction for photon
    theta = pi*(2.0*((FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX) - 1.0);
    photon.eray[2] = 2.0*((FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX) - 1.0;
    photon.eray[0] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*cos(theta);
    photon.eray[1] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*sin(theta);
    for (k=0; k<ndim; k++) 
      photon.inveray[k] = 1.0/(photon.eray[k] + small_number);

  }

  // Isotropic source
  //---------------------------------------------------------------------------
  else if (source.sourcetype == "isotropic") {
    cout << "Isotropic radiation field not yet implemented" << endl;
    exit(0);

  }

  // Isotropic source
  //---------------------------------------------------------------------------
  else if (source.sourcetype == "planar") {
    cout << "Planar radiation field not yet implemented" << endl;
    exit(0);

  }
  //---------------------------------------------------------------------------


#ifdef OUTPUT_ALL
  cout << "Emitting photon " << iphoton << " with direction " 
       << photon.eray[0] << "   " << photon.eray[1] << "   " 
       << photon.eray[2] << endl;
#endif
  

  return photon;
}



//=============================================================================
//  TreeMonteCarlo::ScatterPhotonPacket
/// Scatter photon packet inot random, isotropic direction.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void TreeMonteCarlo<ndim,ParticleType>::ScatterPhotonPacket
(PhotonPacket<ndim> &photon)        ///< [inout] Reference to photon packet
{
  int k;                            // Dimension counter
  FLOAT theta;                      // Random angle for photon direction

  // Generate random direction for photon
  theta = pi*(2.0*((FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX) - 1.0);
  photon.eray[2] = 2.0*((FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX) - 1.0;
  photon.eray[0] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*cos(theta);
  photon.eray[1] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*sin(theta);
  for (k=0; k<ndim; k++) 
    photon.inveray[k] = 1.0/(photon.eray[k] + small_number);

  return;
}



//=============================================================================
//  TreeMonteCarlo::FindRayExitFace
/// Find face in current cell that photon packet will intercept first.
/// Also computes the path length through the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
int TreeMonteCarlo<ndim,ParticleType>::FindRayExitFace
(KDRadTreeCell<ndim> &cell,         ///< [in] Reference to cell
 FLOAT rp[ndim],                    ///< [in] Position of point/ray
 FLOAT eray[ndim],                  ///< [in] Unit vector direction of ray
 FLOAT inveray[ndim],               ///< [in] 1/eray
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
    if (eray[k] > 0.0) {
      daux = (cell.bbmax[k] - rp[k])*inveray[k];
      if (daux < dpath) {
        dpath = daux;
        cexit = cell.cexit[1][k];
      }
    }

    // If radiation if travelling in -ve k-direction
    else {
      daux = (cell.bbmin[k] - rp[k])*inveray[k];
      if (daux < dpath) {
        dpath = daux;
        cexit = cell.cexit[0][k];
      }
    }

#ifdef OUTPUT_ALL
    if (daux < 0.0) {
      cout << "Problem with ray path length : " << daux << "   " << k 
	   << "   " << eray[k] << "   " << cell.bbmin[k] << "   " 
	   << cell.bbmax[k] << "   " <<  rp[k] << endl;
      cout << "LH ray : " << (cell.bbmax[k] - rp[k])/eray[k] << endl;
      cout << "RH ray : " << (cell.bbmin[k] - rp[k])/eray[k] << endl;
      exit(0);
    }
#endif

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
  while (radtree->radcell[c].level < radtree->ltot) {
    c1 = c + 1;
    k_divide = radtree->radcell[c].k_divide;

    // If point is left of divide, pick 1st child cell.  Else pick 2nd child.
    if (rp[k_divide] < radtree->radcell[c1].bbmax[k_divide])
      c = c1;
    else
      c = radtree->radcell[c].c2;

  };
  //---------------------------------------------------------------------------


#ifdef OUTPUT_ALL
  cout << "Looking for cell containing : " 
       << rp[0] << "  " << rp[1] << "  " << rp[2] << endl;
  cout << "Cell x-range : " << radtree->radcell[c].bbmin[0] 
       << "   " << radtree->radcell[c].bbmax[0] << endl;
  cout << "Cell y-range : " << radtree->radcell[c].bbmin[1] 
       << "   " << radtree->radcell[c].bbmax[1] << endl;
  cout << "Cell z-range : " << radtree->radcell[c].bbmin[2] 
       << "   " << radtree->radcell[c].bbmax[2] << endl;
  cout << "Cell level : " << radtree->radcell[c].level << endl;
#endif

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
