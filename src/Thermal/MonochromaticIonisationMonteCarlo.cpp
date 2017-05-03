//=================================================================================================
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
//=================================================================================================


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
#include "RandomNumber.h"
#include "Debug.h"
using namespace std;



//=================================================================================================
//  MonochromaticIonisationMonteCarlo::MonochromaticIonisationMonteCarlo()
/// Constructor for ..
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::MonochromaticIonisationMonteCarlo
 (int Nleafmaxaux, int Nraditerationsaux, int Nradlevelsaux,
  FLOAT Nphotonratioaux, FLOAT tempionaux, FLOAT arecombaux, DOUBLE NLyCaux,
  string rand_algorithm, SimUnits *unitsaux, EOS<ndim> *eosaux)
{
  units          = unitsaux;
  //Nphoton        = Nphotonaux;
  Nphotonratio   = Nphotonratioaux;
  Nraditerations = Nraditerationsaux;
  Nradlevels     = Nradlevelsaux;
  NLyC           = NLyCaux;
  arecomb        = arecombaux;
  temp_ion       = tempionaux/units->temp.outscale;

  //ionconst = NLyC*m_hydrogen*1.0e3*7.9e-18*units->r.outscale*units->r.outcgs/
  //  (units->m.outscale*units->m.outcgs*2.7e-13*(FLOAT) Nphoton);
  ionconst = NLyC*m_hydrogen*1.0e3*7.9e-18*units->r.outscale*units->r.outcgs/
    (units->m.outscale*units->m.outcgs*arecomb);

  // Set all physical values here (scaled to dimensionless units)
  NLyC         = NLyC*units->t.outscale*units->t.outSI;
  Eion         = 13.6*e_charge/units->E.outscale/units->E.outSI;
  invEion      = 1.0/Eion;
  //packetenergy = Eion*NLyC/(FLOAT) Nphoton;
  across       = 7.9e-18/pow(units->r.outscale*units->r.outcgs,2);
  arecomb      = arecomb*units->t.outscale*units->t.outcgs/pow(units->r.outscale*units->r.outcgs,3);
  invmh        = units->m.outscale*units->m.outSI/m_hydrogen;

  cout << "Opacity      : " << 7.9e-22*1.6e-17/m_hydrogen << endl;

  /*cout << "Rstromgren2  : " << pow(3.0*NLyC/(4.0*pi*arecomb*invmh*invmh*
                                   pow(1.6e-17/(units->rho.outscale*units->rho.outSI),2)),onethird)*
                                   units->r.outscale << endl;
  cout << "Eion         : " << Eion << endl;
  cout << "invEion      : " << invEion << endl;
  cout << "packetenergy : " << packetenergy << "   " << NLyC
       << "    " << NLyCaux << "   " << Nphoton << endl;
  cout << "across       : " << across << endl;
  cout << "arecomb      : " << arecomb << endl;
  cout << "invmh        : " << invmh << endl;*/

#if defined _OPENMP
  Nthreads = omp_get_max_threads();
#else
  Nthreads = 1;
#endif

  // Generate random number objects based on number of OpenMP threads
  randNumbArray = new RandomNumber*[Nthreads];
  if (rand_algorithm == "xorshift") {
    for (int i=0; i<Nthreads; i++) randNumbArray[i] = new XorshiftRand(i);
  }
  else if (rand_algorithm == "none") {
    for (int i=0; i<Nthreads; i++) randNumbArray[i] = new DefaultSystemRand(i);
  }
  else {
    string message = "Unrecognised parameter : rand_algorithm= " + rand_algorithm;
    ExceptionHandler::getIstance().raise(message);
  }

  radtree = new KDRadiationTree<ndim,nfreq,ParticleType,CellType>(Nleafmaxaux);
}



//=================================================================================================
//  MonochromaticIonisationMonteCarlo::~MonochromaticIonisationMonteCarlo()
/// Destructor for KD-tree radiation class
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::~MonochromaticIonisationMonteCarlo()
{
  delete radtree;
  for (int i=0; i<Nthreads; i++) delete randNumbArray[i];
  delete[] randNumbArray;
}



//=================================================================================================
//  MonochromaticIonisationMonteCarlo::UpdateRadiationField
/// Update the radiation field by iterating the MCRT method on the tree
/// several times until acceptable convergence has been achieved.
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::UpdateRadiationField
 (int Nhydro,                          ///< [in] No. of hydro particle
  int Nnbody,                          ///< [in] No. of N-body particles
  int Nsink,                           ///< [in] No. of sink particles
  Particle<ndim> *part_gen,            ///< [in] Generic hydro particle data array
  NbodyParticle<ndim> **nbodydata,     ///< [in] N-body data array
  SinkParticle<ndim> *sinkdata)        ///< [in] Sink data array
{
  bool converged;                      // Is radiation field converged?
  int c;                               // Cell counter
  //int i;                               // ..
  int it;                              // Iteration counter
  int k;                               // Dimension counter
  int level;                           // Level to walk tree on
  int Nit;                             // No. of iterations of radiation field
  int Nphoton;                         // No. of photon packets
  RadiationSource<ndim> source;        // Current radiation source
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* >(part_gen);


  debug2("[MonochromaticIonisationMonteCarlo::UpdateRadiationField]");


  // Re-build radiation tree from scratch
  CodeTiming::BlockTimer timer1 = timing->StartNewTimer("IONTREE_BUILD");
  radtree->BuildTree(Nhydro, Nhydro, partdata);
  timer1.EndTiming();

  // Compute maximum radius of cell extent from co-ordinate centre, in order
  // to know where to emit external photons from
  boundaryradius = 0.0;
  for (k=0; k<ndim; k++) boundaryradius +=
    max(pow(radtree->radcell[0].bbmax[k] - radtree->radcell[0].rcell[k],2),
        pow(radtree->radcell[0].rcell[k] - radtree->radcell[0].bbmin[k],2));
  boundaryradius = 1.1*sqrt(boundaryradius);


  // Perform single iteration step on computing radiation field
  CodeTiming::BlockTimer timer2 = timing->StartNewTimer("RADTREE_MONTE_CARLO");
  Nit = Nraditerations + Nradlevels - 1;
  level = radtree->ltot - Nradlevels + 1;
  //level = radtree->ltot;

  // Set initial opacities of cells
  UpdateCellOpacity(radtree->radcell[0], partdata);
  //for (c=0; c<radtree->Ncell; c++) {
  //  radtree->radcell[c].opacity[0] = (1.0 - radtree->radcell[c].Xold)*
  //    across*radtree->radcell[c].rho*invmh;
  //}


  // Main Monte Carlo iteration loop
  //-----------------------------------------------------------------------------------------------
  for (it=0; it<Nit; it++) {

    Nphoton = (int) (Nphotonratio*(FLOAT) pow(2,level));

    cout << "Propagating photons on level : " << level << "     ltot : " << radtree->ltot << endl;

    // Zero all path length sums before new iteration
    for (c=0; c<radtree->Ncell; c++) radtree->radcell[c].lsum[0] = 0.0;

    // Perform one iteration step of the radiation field
    IterateRadiationField(level, Nphoton, Nhydro, Nnbody, Nsink, part_gen, nbodydata, sinkdata);

    // Update the radiation field on all higher cells
    radtree->SumRadiationField(level, radtree->radcell[0]);

    // Update ionisation fractions and check if converged
    converged = UpdateIonisationFraction(level, Nphoton);

    // If converged early on lowest tree level, exit iteration loop
    if (converged && level == radtree->ltot) {
      cout << "Converged radiation field in " << it << " iteration(s)" << endl;
      //break;
    }

    // Move to next level (if not already on lowest level)
    if (level != radtree->ltot) level++;

  }
  //-----------------------------------------------------------------------------------------------


  timer2.EndTiming();

  // Set thermal properties of all particles in leaf cells by interpolation from grid values
  InterpolateParticleProperties(level, Nhydro, part_gen);


  // Output info to file for plotting
  //-----------------------------------------------------------------------------------------------
  /*for (int l=0; l<=level; l++) {
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

  static int Nfiles = 0;
  ofstream outfile;
  string filename;
  string nostring;
  stringstream ss;
  nostring = "";
  ss << setfill('0') << setw(5) << Nfiles;
  nostring = ss.str();

  filename = "uphoton." + nostring + ".dat";
  outfile.open(filename.c_str());
  for (c=0; c<radtree->Ncell; c++) {
    if (radtree->radcell[c].level != level || radtree->radcell[c].N == 0) continue;
    outfile << radtree->radcell[c].rcell[0] << "   "
            << radtree->radcell[c].rcell[1] << "   "
            << sqrt(DotProduct(radtree->radcell[c].rcell,radtree->radcell[c].rcell,ndim)) << "   "
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
  Nfiles++;*/


  //cin >> it;

  return;
}



//=================================================================================================
//  MonochromaticIonisationMonteCarlo::IterateRadiationField
/// Single iteration step for updating the radiation field using the tree-MCRT method.
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::IterateRadiationField
 (const int level,                     ///< Level to walk tree on
  const int Nphoton,                   ///< No. of photon packets
  const int Nhydro,                      ///< No. of hydro particle
  const int Nnbody,                    ///< No. of N-body particles
  const int Nsink,                     ///< No. of sink particles
  Particle<ndim> *part_gen,            ///< Generic hydro particle data array
  NbodyParticle<ndim> **nbodydata,     ///< N-body data array
  SinkParticle<ndim> *sinkdata)        ///< Sink data array
{
  int kfreq;                           // Frequency bin counter
  long int Ncellcount = 0;             // Count no. of cells passed through
  long int Nscattercount = 0;          // No. of scattering events
  RadiationSource<ndim> source;        // Current radiation source

  if (level < 0 && level > radtree->ltot) {
    cout << "Walking invalid tree level" << endl;
    string msg = "Error : walking invalid tree level in "
      "MonochromaticIonisationMonteCarlo::IterateRadiationField";
    ExceptionHandler::getIstance().raise(msg);
  }

  // Emit photon packets from single source (for now)
  source.sourcetype = "pointsource";
  for (int k=0; k<ndim; k++) source.r[k] = 0.0;
  source.c = radtree->FindCell(0,level,source.r);
  source.luminosity = 1.0;
  kfreq = 0;


  // Now emit all photons from radiation sources, updating the radiation field
  //===============================================================================================
#pragma omp parallel default(none) reduction(+:Ncellcount,Nscattercount) shared(kfreq,source)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();        // OpenMP thread i.d.
#else
    const int ithread = 0;
#endif
    int iphoton;                                     // Photon counter
    int k;                                           // Dimension counter
    FLOAT dpath;                                     // Path to next cell boundary
    FLOAT taumax;                                    // Optical depth travelled by photon
    FLOAT tau;                                       // Current value of photon optical depth
    RandomNumber *randnumb = randNumbArray[ithread];  // ..


    // Main photon packet loop counter
    //---------------------------------------------------------------------------------------------
#pragma omp for
    for (iphoton=0; iphoton<Nphoton; iphoton++) {

      // Initialise new photon packet from single source (modify later)
      PhotonPacket<ndim> photon = GenerateNewPhotonPacket(source, randnumb);

      // Calculate optical depth to be travelled by photon
      taumax = -log((FLOAT) 1.0 - randnumb->floatrand());
      tau    = (FLOAT) 0.0;
      Nscattercount++;


      // Main photon transmission/scattering/absorption-reemission iteration loop
      //-------------------------------------------------------------------------------------------
      do {

        // Increase cell counter
        Ncellcount++;

        // Find i.d. of next (parent) cell and the maximum path length travelled in current cell
        photon.cnext = radtree->FindRayExitFace(radtree->radcell[photon.c], photon.r,
                                                photon.eray, photon.inveray, dpath);


        // Check if max optical depth has been reached in order to scatter or absorb/re-emit photon.
        //-----------------------------------------------------------------------------------------
        if (tau + dpath*radtree->radcell[photon.c].opacity[kfreq] > taumax) {

          // Propagate photon packet until absorption/scattering event
          dpath = (taumax - tau)/radtree->radcell[photon.c].opacity[kfreq];
          for (k=0; k<ndim; k++) photon.r[k] += dpath*photon.eray[k];
#pragma omp atomic
          radtree->radcell[photon.c].lsum[kfreq] += dpath;
#pragma omp atomic
          radtree->radcell[photon.c].Nphoton++;

          break;
        }

        // Otherwise, photon continues through cell and exits to adjacent cell
        //-----------------------------------------------------------------------------------------
        else {

        // Propagate photon packet to edge of cell and add contribution to radiation field of cell
        for (k=0; k<ndim; k++) photon.r[k] += dpath*photon.eray[k];
#pragma omp atomic
          radtree->radcell[photon.c].lsum[kfreq] += dpath;
#pragma omp atomic
          radtree->radcell[photon.c].Nphoton++;
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
          photon.c = radtree->FindAdjacentCell(photon.cnext, level, photon.r);

        }
        //-----------------------------------------------------------------------------------------


      } while (photon.c != -1);
      //-------------------------------------------------------------------------------------------

    }
    //---------------------------------------------------------------------------------------------

  }
  //===============================================================================================


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
//  MonochromaticIonisationMonteCarlo::UpdateIonisationFraction
/// Create a new photon packet for a given radiation source.
/// Currently only implements point sources.
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
bool MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::UpdateIonisationFraction
 (const int level,
  const int Nphoton)
{
  bool converged = true;               // ..
  int c;                               // ..
  int Nlivecells = 0;                  // ..
  int Nconverged = 0;                  // ..
  FLOAT nfrac;                         // ..
  MonoIonTreeCell<ndim,nfreq> *cell;   // ..

  //const FLOAT invNphoton = 1.0 / Nphoton;


  //-----------------------------------------------------------------------------------------------
  for (c=0; c<radtree->Ncell; c++) {
    cell = &(radtree->radcell[c]);

    // For cells on upper levels of tree (where radiation is prop)
    //---------------------------------------------------------------------------------------------
    if (cell->level <= level) {

      // Skip any empty cells
      if (cell->N == 0) continue;
      Nlivecells++;

      // Record old ionisation fraction to later test convergence
      for (int j=0; j<20; j++) {
        cell->Xold = cell->Xion + small_number;

        // Compute new ionisation fraction and temperature
        //nfrac = NLyC*across*cell->lsum[0]/(Nphoton*arecomb*invmh*cell->volume*cell->rho*cell->Xold);
        nfrac            = ionconst*cell->lsum[0]/(Nphoton*cell->volume*cell->rho*cell->Xold);
        cell->Xion       = nfrac/((FLOAT) 1.0 + nfrac);
        cell->temp       = temp_ion*cell->Xion;
        cell->opacity[0] = ((FLOAT) 1.0 - cell->Xion)*across*cell->rho*invmh;

        /*cout << "Xion : " << cell->Xion << "    " << cell->Xold << endl;
        cout << "nums : " << ionconst << "    " << cell->lsum[0] << "    " << Nphoton*cell->volume*cell->rho*cell->Xold << endl;
        int jj;
        cin >> jj;*/

        if (nfrac != nfrac || cell->Xion != cell->Xion) {
          cout << "NaN detected while computing Xion : " << c << "    " << nfrac << "   "
               << cell->Xion << "   " << cell->Xold << "   " << cell->lsum[0] << "    "
               << cell->rho << "    " << cell->volume << "    " << cell->N << "    "
               << cell->m << endl;
          ExceptionHandler::getIstance().raise("Error : NaN detected while computing Xion");
        }

      }


      // Check if ionisation fraction has converged or not
      if (cell->Xion < small_number && cell->Xold < small_number) {
        Nconverged++;
        continue;
      }
      else if (cell->Xion > small_number && fabs(cell->Xion - cell->Xold)/cell->Xion > 0.02) {
        converged = false;
      }
      else if (cell->Xold > small_number && fabs(cell->Xion - cell->Xold)/cell->Xold > 0.02) {
        converged = false;
      }
      else {
        Nconverged++;
      }

    }


    //---------------------------------------------------------------------------------------------
    if (cell->level >= level && cell->level < radtree->ltot) {

      radtree->radcell[cell->c1].Xion = cell->Xion;
      radtree->radcell[cell->c1].Xold = cell->Xold;
      radtree->radcell[cell->c1].temp = cell->temp;
      radtree->radcell[cell->c1].opacity[0] = cell->opacity[0];

      radtree->radcell[cell->c2].Xion = cell->Xion;
      radtree->radcell[cell->c2].Xold = cell->Xold;
      radtree->radcell[cell->c2].temp = cell->temp;
      radtree->radcell[cell->c2].opacity[0] = cell->opacity[0];

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------


  cout << "Cell ionisation fractions converged? : " << converged << "   " << Nconverged
       << "    " << radtree->Ncell << "   " << (FLOAT) Nconverged / (FLOAT) Nlivecells << endl;

  return converged;
}



//=================================================================================================
//  MonochromaticIonisationMonteCarlo::GenerateNewPhotonPacket
/// Create a new photon packet for a given radiation source.
/// Currently only implements point sources.
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
PhotonPacket<ndim> MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::GenerateNewPhotonPacket
 (const RadiationSource<ndim> &source,     ///< ..
  RandomNumber *randnumb)                  ///< ..
{
  int k;                                   // Dimension counter
  PhotonPacket<ndim> photon;               // Photon packet type


  // Initialise new photon packet from chosen source (e.g. pick random direction, frequency, etc..)
  //-----------------------------------------------------------------------------------------------
  if (source.sourcetype == "pointsource") {
    photon.c = source.c;
    photon.cnext = source.c;
    //photon.energy = packetenergy;
    for (k=0; k<ndim; k++) photon.r[k] = source.r[k];

    // Generate random direction for photon
    if (ndim == 1) {
      if (randnumb->floatrand() < (FLOAT) 0.5) {
        photon.eray[0] = -1.0;
      }
      else {
        photon.eray[0] = 1.0;
      }
    }
    else if (ndim == 2) {
      FLOAT theta = twopi*randnumb->floatrand();
      photon.eray[0] = cos(theta);
      photon.eray[1] = sin(theta);
    }
    else if (ndim == 3) {
      FLOAT theta = pi*(2.0*randnumb->floatrand() - 1.0);
      photon.eray[2] = 2.0*randnumb->floatrand() - 1.0;
      photon.eray[0] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*cos(theta);
      photon.eray[1] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*sin(theta);
    }
    for (k=0; k<ndim; k++) photon.inveray[k] = 1.0/(photon.eray[k] + small_number);

  }

  // Isotropic source
  //-----------------------------------------------------------------------------------------------
  else if (source.sourcetype == "isotropic") {
    ExceptionHandler::getIstance().raise("Error : isotropic radiation field not yet implemented");
  }
  // Planar source
  //-----------------------------------------------------------------------------------------------
  else if (source.sourcetype == "planar") {
    ExceptionHandler::getIstance().raise("Error : planar radiation field not yet implemented");
  }
  //-----------------------------------------------------------------------------------------------


#ifdef OUTPUT_ALL
  cout << "Emitting photon with direction " << photon.eray[0] << "   "
       << photon.eray[1] << "   " << photon.eray[2] << endl;
#endif


  return photon;
}



//=================================================================================================
//  MonochromaticIonisationMonteCarlo::ScatterPhotonPacket
/// Scatter photon packet inot random, isotropic direction.
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::ScatterPhotonPacket
 (PhotonPacket<ndim> &photon,          ///< [inout] Reference to photon packet
  RandomNumber *randnumb)              ///< ..
{
  int k;                               // Dimension counter

  // Generate random direction for photon
  if (ndim == 1) {
    if (randnumb->floatrand() < (FLOAT) 0.5) {
      photon.eray[0] = -1.0;
    }
    else {
      photon.eray[0] = 1.0;
    }
  }
  else if (ndim == 2) {
    FLOAT theta = twopi*randnumb->floatrand();
    photon.eray[0] = cos(theta);
    photon.eray[1] = sin(theta);
  }
  else if (ndim == 3) {
    FLOAT theta = pi*(2.0*randnumb->floatrand() - 1.0);
    photon.eray[2] = 2.0*randnumb->floatrand() - 1.0;
    photon.eray[0] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*cos(theta);
    photon.eray[1] = sqrt(1.0 - photon.eray[2]*photon.eray[2])*sin(theta);
  }
  for (k=0; k<ndim; k++) photon.inveray[k] = 1.0/(photon.eray[k] + small_number);

  return;
}



//=================================================================================================
//  MonochromaticIonisationMonteCarlo::UpdateCellOpacity
/// ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::UpdateCellOpacity
 (CellType<ndim,nfreq> &cell,          ///< Reference to current tree cell
  ParticleType<ndim> *partdata)        ///< Particle data array
{
  int i;                               // Aux. child cell counter

  // If cell is not leaf, stock child cells
  if (cell.level != radtree->ltot) {
    for (i=0; i<2; i++) {
      if (i == 0) UpdateCellOpacity(radtree->radcell[cell.c1], partdata);
      else if (i == 1) UpdateCellOpacity(radtree->radcell[cell.c2], partdata);
    }
  }

  // Stock node once all children are stocked
  //-----------------------------------------------------------------------------------------------
  if (cell.level == radtree->ltot) {
    cell.Xion = 0.0;
    i = cell.ifirst;

    // Loop over all particles in cell summing their contributions
    while (i != -1) {
      cell.Xion += partdata[i].m*partdata[i].Xion;  //ionfrac;
      if (i == cell.ilast) break;
      i = radtree->inext[i];
    };
    if (cell.m > 0.0) cell.Xion /= cell.m;
    cell.opacity[0] = ((FLOAT) 1.0 - cell.Xion)*across*cell.rho*invmh;
  }

  // For non-leaf cells, sum together two children cells
  //-----------------------------------------------------------------------------------------------
  else {
    cell.Xion = 0.0;
    if (cell.m > 0.0) {
      cell.Xion = radtree->radcell[cell.c1].m*radtree->radcell[cell.c1].Xion +
        radtree->radcell[cell.c2].m*radtree->radcell[cell.c2].Xion;
      cell.Xion /= cell.m;
    }
    cell.opacity[0] = ((FLOAT) 1.0 - cell.Xion)*across*cell.rho*invmh;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  MonochromaticIonisationMonteCarlo::InterpolateParticleProperties
/// ...
//=================================================================================================
template <int ndim, int nfreq, template<int> class ParticleType, template<int,int> class CellType>
void MonochromaticIonisationMonteCarlo<ndim,nfreq,ParticleType,CellType>::InterpolateParticleProperties
 (const int level,                     ///< Level to walk tree on
  const int Nhydro,                    ///< No. of hydro particle
  Particle<ndim> *part_gen)            ///< Generic hydro particle data array
{
  int c;
  int i;
  int k;
  int Ngather;
  int Ngathermax = 1024;
  int *celllist = new int[Ngathermax];
  FLOAT xion = 0.0;
  FLOAT xionNorm = 0.0;
  FLOAT dr[ndim];
  FLOAT drmag;
  FLOAT weight;
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* >(part_gen);

  for (c=0; c<radtree->Ncell; c++) {
    MonoIonTreeCell<ndim,nfreq> &cell = radtree->radcell[c];
    //cout << "Xion : " << cell.Xion << endl;
    if (cell.Xion != cell.Xion) {
      cout << "NaN detected : " << c << "   " << cell.Xion << endl;
      string msg = "Nan detected in MonochromaticIonisationMonteCarlo::InterpolateParticleProperties";
      ExceptionHandler::getIstance().raise(msg);
    }

  }

  //-----------------------------------------------------------------------------------------------
  for (c=0; c<radtree->Ncell; c++) {
    if (radtree->radcell[c].level == level && radtree->radcell[c].N != 0) {
      i = radtree->radcell[c].ifirst;

      while (i != -1) {
        ParticleType<ndim> &part = partdata[i];
        do {
          Ngather = radtree->ComputeGatherCellList(part.r, 3.0*part.h, Ngathermax, celllist);
          if (Ngather == 0) {
            ExceptionHandler::getIstance().raise("Ngather = 0 in radtree->ComputeGatherCellList");
          }
          if (Ngather > 0) break;
          delete[] celllist;
          Ngathermax *= 2;
          celllist = new int[Ngathermax];
        } while (Ngather < 0);

        /*xion = 0.0;
        xionNorm = 0.0;
        for (int jj=0; jj<Ngather; jj++) {
          int cc = celllist[jj];
          MonoIonTreeCell<ndim,nfreq> &cell = radtree->radcell[cc];

          for (k=0; k<ndim; k++) dr[k] = cell.rcell[k] - part.r[k];
          drmag    = sqrt(DotProduct(dr,dr,ndim) + small_number);
          if (drmag*part.invh < 2.0) {
            //weight = 1.0 - drmag*part.invh;
            weight = exp(-drmag*drmag*part.invh*part.invh*2.0*2.0);
          }
          else {
            weight = 0.0;
          }
          xion     += cell.Xion*weight;
          xionNorm += weight;
        }

        if (xion == 0.0) {
          part.ionfrac = radtree->radcell[c].Xion;
        }
        else {
          part.ionfrac = xion/xionNorm;
        }
        part.Xion = radtree->radcell[c].Xion;
        part.ionfrac = part.Xion;*/

        // Use Inverse distance weighting (Shepard's method) to interpolate thermal proeprties
        xion = 0.0;
        xionNorm = 0.0;
        FLOAT pshep = 8.0;
        for (int jj=0; jj<Ngather; jj++) {
          int cc = celllist[jj];
          MonoIonTreeCell<ndim,nfreq> &cell = radtree->radcell[cc];

          for (k=0; k<ndim; k++) dr[k] = cell.rcell[k] - part.r[k];
          drmag    = sqrt(DotProduct(dr,dr,ndim) + small_number);
          /*if (drmag*part.invh < 2.0) {
            weight = max(0.0, 2.0*part.h - drmagh;
          }
          else {
            weight = 0.0;
          }*/
          weight = powf(max(0.0, 2.0*part.h - drmag)/(2.0*part.h*(drmag + 0.00001*part.h)), pshep);
          xion     += cell.Xion*weight;
          xionNorm += weight;
        }

        if (xion == 0.0) {
          part.ionfrac = radtree->radcell[c].Xion;
        }
        else {
          part.ionfrac = xion/xionNorm;
        }
        part.Xion = radtree->radcell[c].Xion;
        part.ionfrac = part.Xion;


        if (i == radtree->radcell[c].ilast) break;
        i = radtree->inext[i];
      };
    }
  }
  //-----------------------------------------------------------------------------------------------

  return;
}



template class MonochromaticIonisationMonteCarlo<1,1,GradhSphParticle,MonoIonTreeCell>;
template class MonochromaticIonisationMonteCarlo<2,1,GradhSphParticle,MonoIonTreeCell>;
template class MonochromaticIonisationMonteCarlo<3,1,GradhSphParticle,MonoIonTreeCell>;
