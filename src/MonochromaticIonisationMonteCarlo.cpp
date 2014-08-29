//=============================================================================
//  MonochromaticIonisationMonteCarlo.cpp
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
#include "Constants.h"
#include "Radiation.h"
#include "Debug.h"
using namespace std;



//=============================================================================
//  MonochromaticIonisationMonteCarlo::MonochromaticIonisationMonteCarlo()
/// Constructor for ..
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::MonochromaticIonisationMonteCarlo
(int Nphotonaux, int Nleafmaxaux, FLOAT tempionaux, DOUBLE NLyCaux,
 RandomNumber *randaux, SimUnits *unitsaux, EOS<ndim> *eosaux)
{
  randnumb = randaux;
  units    = unitsaux;
  Nphoton  = Nphotonaux;
  NLyC     = 1.0e49; //NLyCaux;
  temp_ion = tempionaux/units->temp.outscale;

  ionconst = NLyC*m_hydrogen*1.0e3*7.9e-18*units->r.outscale*units->r.outcgs/
    (units->m.outscale*units->m.outcgs*2.7e-13*(FLOAT) Nphoton);

  cout << "Rstromgren   : " << pow(3.0*m_hydrogen*m_hydrogen*NLyC/
                                   (4.0*pi*2.7e-19*1.6e-17*1.6e-17),onethird)/r_pc << endl;
  cout << "Opacity      : " << 7.9e-22*1.6e-17/m_hydrogen << endl;


  // Set all physical values here (scaled to dimensionless units)
  NLyC         = NLyC*units->t.outscale*units->t.outSI;
  Eion         = 13.6*e_charge/units->E.outscale/units->E.outSI;
  invEion      = 1.0/Eion;
  packetenergy = Eion*NLyC/(FLOAT) Nphoton;
  across       = 7.9e-18/pow(units->r.outscale*units->r.outcgs,2);
  arecomb      = 2.7e-13*units->t.outscale*units->t.outcgs/
    pow(units->r.outscale*units->r.outcgs,3);
  invmh        = units->m.outscale*units->m.outSI/m_hydrogen;

  cout << "Rstromgren2  : " << pow(3.0*NLyC/(4.0*pi*arecomb*invmh*invmh*
                                   pow(1.6e-17/(units->rho.outscale*units->rho.outSI),2)),onethird)*
                                   units->r.outscale << endl;


  cout << "Eion         : " << Eion << endl;
  cout << "invEion      : " << invEion << endl;
  cout << "packetenergy : " << packetenergy << "   " << NLyC
       << "    " << NLyCaux << "   " << Nphoton << endl;
  cout << "across       : " << across << endl;
  cout << "arecomb      : " << arecomb << endl;
  cout << "invmh        : " << invmh << endl;

  radtree = new KDRadiationTree<ndim,nfreq,ParticleType,CellType>(Nleafmaxaux);
}



//=============================================================================
//  MonochromaticIonisationMonteCarlo::~MonochromaticIonisationMonteCarlo()
/// Destructor for KD-tree radiation class
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::~MonochromaticIonisationMonteCarlo()
{
}



//=============================================================================
//  MonochromaticIonisationMonteCarlo::UpdateRadiationField
/// Update the radiation field by iterating the MCRT method on the tree
/// several times until acceptable convergence has been achieved.
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::UpdateRadiationField
(int Nsph,                          ///< No. of SPH particle
 int Nnbody,                        ///< No. of N-body particles
 int Nsink,                         ///< No. of sink particles
 SphParticle<ndim> *sph_gen,        ///< Generic SPH particle data array
 NbodyParticle<ndim> **nbodydata,   ///< N-body data array
 SinkParticle<ndim> *sinkdata)      ///< Sink data array
{
  bool converged;                   // Is radiation field converged?
  int c;                            // Cell counter
  int i;                            // ..
  int isource;                      // Radiation source i.d.
  int it;                           // Iteration counter
  int k;                            // Dimension counter
  int level;                        // Level to walk tree on
  RadiationSource<ndim> source;     // Current radiation source
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* >(sph_gen);


  debug2("[MonochromaticIonisationMonteCarlo::UpdateRadiationField]");


  // Re-build radiation tree from scratch
  timing->StartTimingSection("IONTREE_BUILD",2);
  radtree->BuildTree(Nsph,Nsph,sphdata);
  timing->EndTimingSection("IONTREE_BUILD");

  // Compute maximum radius of cell extent from co-ordinate centre, in order
  // to know where to emit external photons from
  boundaryradius = 0.0;
  for (k=0; k<ndim; k++) boundaryradius +=
    max(pow(radtree->radcell[0].bbmax[k] - radtree->radcell[0].rcell[k],2),
	pow(radtree->radcell[0].rcell[k] - radtree->radcell[0].bbmin[k],2));
  boundaryradius = 1.1*sqrt(boundaryradius);


  // Perform single iteration step on computing radiation field
  timing->StartTimingSection("RADTREE_MONTE_CARLO",2);
  level = radtree->ltot;

  int Nit = 1;

  // Set initial opacities of cells
  for (c=0; c<radtree->Ncell; c++) {
    radtree->radcell[c].opacity[0] = (1.0 - radtree->radcell[c].Xold)*
      across*radtree->radcell[c].rho*invmh;
  }


  // Main Monte Carlo iteration loop
  //---------------------------------------------------------------------------
  for (it=0; it<Nit; it++) {

    // Zero all path length sums before new iteration
    for (c=0; c<radtree->Ncell; c++) radtree->radcell[c].lsum[0] = 0.0;

    // Perform one iteration step of the radiation field
    IterateRadiationField(level,Nsph,Nnbody,Nsink,sph_gen,nbodydata,sinkdata);

    // Update the radiation field on all higher cells
    radtree->SumRadiationField(level,radtree->radcell[0]);

    // Update ionisation fractions and check if converged
    converged = UpdateIonisationFraction();

    // If converged early, exit iteration loop
    if (converged) {
      cout << "Converged radiation field in " << it << " iteration(s)" << endl;
      //break;
    }

  }
  //---------------------------------------------------------------------------


  timing->EndTimingSection("RADTREE_MONTE_CARLO");


  // Set thermal properties of all particles in leaf cells
  for (c=0; c<radtree->Ncell; c++) {
    if (radtree->radcell[c].level == level) {
      i = radtree->radcell[c].ifirst;
      while (i != -1) {
        sphdata[i].ionfrac = radtree->radcell[c].Xion;
        if (i == radtree->radcell[c].ilast) break;
        i = radtree->inext[i];
      };
    }
  }


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
      if (radtree->radcell[c].level != l || radtree->radcell[c].N == 0) continue;
      outfile << sqrt(DotProduct(radtree->radcell[c].rcell,radtree->radcell[c].rcell,ndim)) << "   "
              << radtree->radcell[c].Xion << "   "
              << radtree->radcell[c].Xold << "   "
              << radtree->radcell[c].opacity[0]/(units->r.outscale*units->r.outSI) << "   "
              << pow(units->r.outscale*units->r.outcgs,2)*
                     radtree->radcell[c].opacity[0]/radtree->radcell[c].rho/invmh << "   "
              << radtree->radcell[c].lsum[0] << "   "
              << radtree->radcell[c].temp << "   "
              << radtree->radcell[c].rho << "   "
              << radtree->radcell[c].Nphoton << "    "
              << radtree->radcell[c].N << "    "
              << radtree->radcell[c].volume << endl;
    }
    outfile.close();

  }

  return;
}



//=============================================================================
//  MonochromaticIonisationMonteCarlo::IterateRadiationField
/// Single iteration step for updating the radiation field using the
/// tree-MCRT method.
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::IterateRadiationField
(int level,                         ///< Level to walk tree on
 int Nsph,                          ///< No. of SPH particle
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
  int kfreq;                        // Frequency bin counter
  long int Ncellcount = 0;          // Count no. of cells passed through
  long int Nscattercount = 0;       // No. of scattering events
  FLOAT dpath;                      // Path to next cell boundary
  FLOAT rand1;                      // Random number
  FLOAT taumax;                     // Optical depth travelled by photon
  FLOAT tau;                        // Current value of photon optical depth
  FLOAT theta;                      // Random angle for photon direction
  PhotonPacket<ndim> photon;        // Current photon packet
  RadiationSource<ndim> source;     // Current radiation source
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* >(sph_gen);

  if (level < 0 && level > radtree->ltot) {
    cout << "Walking invalid tree level" << endl;
    exit(0);
  }

  // Emit photon packets from single source (for now)
  source.sourcetype = "pointsource";
  for (k=0; k<ndim; k++) source.r[k] = 0.0;
  source.c = radtree->FindCell(0,level,source.r);
  source.luminosity = 1.0;
  kfreq = 0;

  cout << "Source located at : " << source.c << "   "
       << radtree->radcell[source.c].r[0] << "   "
       << radtree->radcell[source.c].r[1] << "   "
       << radtree->radcell[source.c].r[2] << endl;


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
      photon.cnext = radtree->FindRayExitFace(radtree->radcell[photon.c],photon.r,
                                              photon.eray,photon.inveray,dpath);


      // Check if maximum optical depth has been reached in order to
      // scatter or absorb/re-emit photon.
      //-----------------------------------------------------------------------
      if (tau + dpath*radtree->radcell[photon.c].opacity[kfreq] > taumax) {

        // Propagate photon packet until absorption/scattering event
        dpath = (taumax - tau)/radtree->radcell[photon.c].opacity[kfreq];
        for (k=0; k<ndim; k++) photon.r[k] += dpath*photon.eray[k];
#pragma omp atomic
        radtree->radcell[photon.c].lsum[kfreq] += dpath; //packetenergy*dpath;
#pragma omp atomic
        radtree->radcell[photon.c].Nphoton++;

        break;
      }

      // Otherwise, photon continues through cell and exits to adjacent cell
      //-----------------------------------------------------------------------
      else {

        // Propagate photon packet to edge of cell and add contribution to
        // radiation field of cell
        for (k=0; k<ndim; k++) photon.r[k] += dpath*photon.eray[k];
#pragma omp atomic
        radtree->radcell[photon.c].lsum[kfreq] += dpath; //packetenergy*dpath;
#pragma omp atomic
        radtree->radcell[photon.c].Nphoton++;
        //radtree->radcell[photon.c].uphoton += packetenergy*dpath;
        tau += dpath*radtree->radcell[photon.c].opacity[kfreq];


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



//=============================================================================
//  MonochromaticIonisationMonteCarlo::UpdateIonisationFraction
/// Create a new photon packet for a given radiation source.
/// Currently only implements point sources.
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
bool MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::UpdateIonisationFraction(void)
{
  bool converged = true;
  int c;
  int k;
  int Nconverged = 0;
  FLOAT nfrac;
  MonoIonTreeCell<ndim,nfreq> *cell;


  //---------------------------------------------------------------------------
  for (c=0; c<radtree->Ncell; c++) {
    cell = &(radtree->radcell[c]);

    // Record old ionisation fraction to later test convergence
    for (int j=0; j<20; j++) {
      cell->Xold = cell->Xion + small_number;

      // Compute new ionisation fraction and temperature
      //nfrac = NLyC*across*cell->lsum[0]/(Nphoton*arecomb*invmh*cell->volume*cell->rho*cell->Xold);
      nfrac = ionconst*cell->lsum[0]/(cell->volume*cell->rho*cell->Xold);
      cell->Xion = nfrac/(1.0 + nfrac);
      cell->temp = temp_ion*cell->Xion;
      cell->opacity[0] = ((FLOAT) 1.0 - cell->Xion)*across*cell->rho*invmh;

    }


    // Check if ionisation fraction has converged or not
    if (cell->Xion < small_number && cell->Xold < small_number) {
      Nconverged++;
      continue;
    }
    else if (cell->Xion > small_number &&
	     fabs(cell->Xion - cell->Xold)/cell->Xion > 0.1)
      converged = false;
    else if (cell->Xold > small_number &&
	     fabs(cell->Xion - cell->Xold)/cell->Xold > 0.1)
      converged = false;
    else
      Nconverged++;

  }
  //---------------------------------------------------------------------------

  cout << "Cell ionisation fractions converged? : " << converged << "   "
       << Nconverged << "    " << radtree->Ncell << "   "
       << (FLOAT) Nconverged / (FLOAT) radtree->Ncell << endl;

  return converged;
}



//=============================================================================
//  MonochromaticIonisationMonteCarlo::GenerateNewPhotonPacket
/// Create a new photon packet for a given radiation source.
/// Currently only implements point sources.
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
PhotonPacket<ndim> MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::GenerateNewPhotonPacket
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
    theta = pi*(2.0*randnumb->floatrand() - 1.0);
    photon.eray[2] = 2.0*randnumb->floatrand() - 1.0;
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
//  MonochromaticIonisationMonteCarlo::ScatterPhotonPacket
/// Scatter photon packet inot random, isotropic direction.
//=============================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::ScatterPhotonPacket
(PhotonPacket<ndim> &photon)        ///< [inout] Reference to photon packet
{
  int k;                            // Dimension counter
  FLOAT theta;                      // Random angle for photon direction

  // Generate random direction for photon
  //theta = pi*(2.0*((FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX) - 1.0);
  //photon.eray[2] = 2.0*((FLOAT)(rand()%RAND_MAX)/(FLOAT)RAND_MAX) - 1.0;
  theta = pi*(2.0*randnumb->floatrand() - 1.0);
  photon.eray[2] = 2.0*randnumb->floatrand() - 1.0;
  photon.eray[0] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*cos(theta);
  photon.eray[1] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*sin(theta);
  for (k=0; k<ndim; k++)
    photon.inveray[k] = 1.0/(photon.eray[k] + small_number);

  return;
}



template class MonochromaticIonisationMonteCarlo<1,1,GradhSphParticle,MonoIonTreeCell>;
template class MonochromaticIonisationMonteCarlo<2,1,GradhSphParticle,MonoIonTreeCell>;
template class MonochromaticIonisationMonteCarlo<3,1,GradhSphParticle,MonoIonTreeCell>;
