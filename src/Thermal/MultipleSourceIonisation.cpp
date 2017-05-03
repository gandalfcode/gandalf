//=================================================================================================
//  MultipleSourceIonisation.cpp
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



#include "Constants.h"
#include "Precision.h"
#include "Particle.h"
#include "Debug.h"
#include "NeighbourSearch.h"
#include "Radiation.h"
#include "Sinks.h"
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include "Exception.h"
using namespace std;


//=================================================================================================
//	MultipleSourceIonisation::MultipleSourceIonisation
//	MultipleSourceIonisation class constructor
//=================================================================================================

template <int ndim, template<int> class ParticleType>
MultipleSourceIonisation<ndim,ParticleType>::MultipleSourceIonisation(
  NeighbourSearch<ndim> *partneibaux ,
  FLOAT mu_baraux,
  FLOAT X_compaux,
  FLOAT mu_ionaux,
  FLOAT temp0aux,
  FLOAT temp_ionaux,
  DOUBLE Ndotminaux,
  FLOAT gamma_eosaux,
  FLOAT arecombaux,
  FLOAT scaleaux,
  FLOAT tempscaleaux,
  DOUBLE rad_contaux)
{
  neib = partneibaux;
  mu_bar=mu_baraux;
  X_comp=X_compaux;
  mu_ion=mu_ionaux;
  temp0=temp0aux;
  temp_ion=temp_ionaux;
  Ndotmin=Ndotminaux;
  gamma_eos=gamma_eosaux;
  arecomb=arecombaux;
  scale=scaleaux;
  tempscale=tempscaleaux;
  rad_cont=rad_contaux;
#ifdef MPI_PARALLEL
  string message = "Multiple Source Ionisations does not work with MPI";
  ExceptionHandler::getIstance().raise(message);
#endif

}

//=================================================================================================
//  MultipleSourceIonisation::~MultipleSourceIonisation
///	MultipleSourceIonisation class destructor
//=================================================================================================

template <int ndim, template<int> class ParticleType>
MultipleSourceIonisation<ndim,ParticleType>::~MultipleSourceIonisation()
{
}



//=================================================================================================
//  MultipleSourceIonisation::UpdateRadiationFieldMMS
/// Calculates the internal energy of particles due to ionising radiation.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::UpdateRadiationField
 (int N,
  int nos,
  int aux,
  Particle<ndim>  * partgen,
  NbodyParticle<ndim> ** ndata,
  SinkParticle<ndim> * partaux)
{
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (partgen);

  debug2("[MultipleSourceIonisation::UpdateRadiationField]");

  ionisation_intergration(nos,N,ndata,partdata,scale,tempscale,neib,
                          temp0,mu_bar,mu_ion,temp_ion,Ndotmin,1./gamma_eos);

  return;
}



//=================================================================================================
//  MultipleSourceIonisation::probs
/// Works out the fraction of ionisation each star is responsible for.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::probs
 (int &nos,
  ionpar *ionisedpart,
  int *sinkid,
  int &testpart,
  DOUBLE *ndot)
{
  int pp;
  DOUBLE fluxcontrole[nos];                                            //To store temporary radii
  DOUBLE sum=0;	                                                       //Controle to stop devision by zero

  // Loop over sources and add photon flux at current location
  for (pp=0;pp<nos;pp++) {

    // If the particle is ionised by the source we are testing, calculate contribution
    if (ionisedpart[ionisedpart[testpart].neigh[pp]].ionised[pp] == 1) {
      fluxcontrole[pp] = ndot[pp] - ionisedpart[ionisedpart[testpart].neigh[pp]].photons[pp];
    }
    // If not then the contribution will be zero
    else {
      fluxcontrole[pp] = 0;
    }

   // Add to controle parameter to ensure we carry out the scaling
    sum = sum + fluxcontrole[pp];
  }

  //Scale so total fraction of used photons is one
  for (pp=0;pp<nos;pp++) {
    if (sum > 0) {                                               // Do we have to scale
      ionisedpart[testpart].prob[pp]=fluxcontrole[pp]/sum;       // Scale to one
    }
    else {
      ionisedpart[testpart].prob[pp]=fluxcontrole[pp];           // Pass through 0s as no photons are received
    }
  }

  return;
}



//=================================================================================================
//  MultipleSourceIonisation::lost
/// Recursive function that works out number of photons used in ionisation along the path
/// from the source to the particle.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
DOUBLE MultipleSourceIonisation<ndim, ParticleType>::lost
 (ionpar *ionisedpart,
  int *sinkid,
  DOUBLE *ndot,
  int &N,
  int &pp,
  int &testpart,
  int &nos,
  int &change)
{
  DOUBLE absorbed = 0.0;                         // Total amount absorbed
  DOUBLE d1;                                     // Distance between the particle and the source
  DOUBLE d2;                                     // Distance between the neighbour and the source

  // If the particle has not been checked
  //-----------------------------------------------------------------------------------------------
  if (ionisedpart[testpart].checked[pp] == 0) {

    // Works out amount of photons lost along the path
    // Is the particle a sink (this stops recursion when we get to the source)
    if (ionisedpart[testpart].sink == 0) {

      // Call the probs function to work out ionisation fraction of each source
      probs(nos, ionisedpart, sinkid, testpart, ndot);

      // Does the particle have a neighbour for this source?
      if (ionisedpart[testpart].neigh[pp] != N) {

        d1 = sqrt(pow(ionisedpart[testpart].x - ionisedpart[sinkid[pp]].x, 2) +
                  pow(ionisedpart[testpart].y - ionisedpart[sinkid[pp]].y, 2) +
                  pow(ionisedpart[testpart].z - ionisedpart[sinkid[pp]].z, 2));
        d2 = sqrt(pow(ionisedpart[ionisedpart[testpart].neigh[pp]].x - ionisedpart[sinkid[pp]].x, 2) +
                  pow(ionisedpart[ionisedpart[testpart].neigh[pp]].y - ionisedpart[sinkid[pp]].y, 2) +
                  pow(ionisedpart[ionisedpart[testpart].neigh[pp]].z - ionisedpart[sinkid[pp]].z, 2));

        if (ionisedpart[ionisedpart[testpart].neigh[pp]].sink == 0) {
          absorbed = ((pow((ionisedpart[testpart].rho + ionisedpart[ionisedpart[testpart].neigh[pp]].rho)/2.,2.))/3.)*
                     (pow(d1,3.)-pow(d2,3.))*(ionisedpart[testpart].prob[pp]) +
                      lost(ionisedpart, sinkid, ndot, N, pp, ionisedpart[testpart].neigh[pp], nos, change);
        }
        else {
          absorbed = ((pow((ionisedpart[testpart].rho),2.))/3.)*(pow(d1,3.)-pow(d2,3.))*
                     (ionisedpart[testpart].prob[pp]) +
                      lost(ionisedpart, sinkid, ndot, N, pp, ionisedpart[testpart].neigh[pp], nos, change);
        }
          //absorbed=rho**2.*[d1**3-d2**3]/3*transmitonfrac+absorbed(previous particle in chain)
      }

      // All photons used up as we cant have links to the dummy partcile
      else {
        absorbed=ndot[pp];
      }
    }
    ionisedpart[testpart].photons[pp] = absorbed;

    if ((ndot[pp] - absorbed) > 0) {

      // Record if the particle is changing state (for convergence)
      if (ionisedpart[testpart].ionised[pp] == 0) {
        change = change + 1;
      }

      // Set particle as source ionised
      ionisedpart[testpart].ionised[pp] = 1;
    }
    else {

      // Record if the particle is changing state (for convergence)
      if (ionisedpart[testpart].ionised[pp] == 1) {
        change=change+1;
      }

      // Set particle as not source ionised
      ionisedpart[testpart].ionised[pp] = 0;
    }
    ionisedpart[testpart].checked[pp] = 1;
    return absorbed;
  }
  else {
    return ionisedpart[testpart].photons[pp];
  }
}



//=================================================================================================
//  MultipleSourceIonisation::photoncount
/// Works out if a particle is ionised using partner functions.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::photoncount
 (ionpar *ionisedpart,
  int *sinkid,
  DOUBLE *ndot,
  int &N,
  int &nos,
  int &testpart,
  int &change)
{
  int pp;

  // Set as not ionised initially
  ionisedpart[testpart].fionised = 0;

  // Looping over all sources to test if the amount lost smaller than the amount available
  for (pp=0;pp<nos;pp++) {
    lost(ionisedpart, sinkid, ndot, N, pp, testpart, nos, change);
  }

  // Check to see if the particle should be ionised
  ionisedpart[testpart].fionised = 0;
  for (pp=0;pp<nos;pp++) {
    if (ionisedpart[testpart].ionised[pp]==1) ionisedpart[testpart].fionised = 1;
  }

  return;
}



//=================================================================================================
//  MultipleSourceIonisation::ionisation_intergration
/// Main contole routine.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MultipleSourceIonisation<ndim, ParticleType>::ionisation_intergration
 (int newnos,                           ///< Number of ionising sources
  int N,                                ///< Number of SPH particles
  NbodyParticle<ndim> ** ndata,         ///< Source Data
  Particle<ndim> * partgen,             ///< SPH particle data
  DOUBLE scale,                         ///< Scaling
  DOUBLE tempscale,                     ///< Temperature scaling
  NeighbourSearch<ndim> * neib,         ///< Neighbour Search Routine
  DOUBLE tn,                            ///< Neutral gas temperature
  DOUBLE mu_bar,                        ///< Average neutral gas mass
  DOUBLE mu_ion,                        ///< Average ionised gas mass
  DOUBLE ti,                            ///< Ionised gas temperature
  DOUBLE Ndotmin,                       ///< Minimum Ionising output
  DOUBLE gammam1)                       ///< 1/gamma
{
  int ii,jj,pp,tt; //Integer allocation for loops
  int debug=0,smoothing=1; //Debug mode controler and maximum number of neighbours allowed
  FLOAT delta=0;
  const FLOAT m_per_hatom=m_hydrogen*1000.0/X_comp;
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (partgen);

  struct timeval start, end;
  gettimeofday(&start, NULL);

  //Check that the stellar.dat file is present
  bool check;
  struct stat buf;
  stat("stellar.dat", &buf);
  check = S_ISREG(buf.st_mode);
  if (check == 0) {
    string message = "stellar.dat is not present in run directory, ionisation will not be included";
    ExceptionHandler::getIstance().raise(message);
    return;
  }

  // Checks if there are currently any sinks in gandalf and if not returns
  if (newnos == 0) {
    cout << "No stars" << endl;
    return;
  }

  int nos = 0;                                   ///< No. of sources
#ifndef RAD_OPTIMISE
  int *newnosid = new int[newnos];               ///< Which sinks are active sources
#else
  if (nos > maxSources) {
    string message = "Too many sources for optimised mode.";
    ExceptionHandler::getIstance().raise(message);
  }
  int newnosid[maxSources];
#endif

  // Determines which sinks are active sources based on user choices
  for(ii=0; ii<newnos; ii++) {
    if (ndata[ii]->NLyC >= Ndotmin) {
      nos = nos + 1;
      newnosid[nos - 1] = ii;
    }
  }


  // Checks if the sinks are of large enough size
  if (nos == 0) {
    cout << "No stars of suitable mass" << endl;
    return;
  }

  if (debug == 1) {
    cout << "# of sources followed is " << nos << ". ";
  }

  // Increases N to accomidate sinks
  N = N + nos;

  if (debug == 1) {
    gettimeofday(&end, NULL);

    delta = (((end.tv_sec  - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec) / 1.e6)-delta;
    cout<<delta<<"s to ";
    cout<<"Starting"<<endl; //Debug message
  }


  // Fill ionisedpart particle array
  // Create sink id table and ndot table
#ifndef RAD_OPTIMISE
  int *sinkid = new int[nos];
  DOUBLE *ndot = new DOUBLE[nos];
#else
  int sinkid[maxSources];
  DOUBLE ndot[maxSources];
#endif

  // Create the ionisedpart array and resize to the number of particles
  // (Number of simulation particiles, inluding sinks. Also the dummy particle.)
  ionpar *ionisedpart = new ionpar[N+1];

  if (ionisation_fraction.size() == 0) {
    ionisation_fraction.resize(N+1);
    for (ii=0; ii<N+1; ii++) {
      ionisation_fraction[ii].resize(newnos);
      for (pp=0; pp<nos; pp++) {
        ionisation_fraction[ii][pp] = 0;
      }
    }
  }
  ionisation_fraction.resize(N + 1);

  // Add ionisedpart particle data
#pragma omp parallel for private(ii,jj)
  for (ii=0; ii<N-nos; ii++) {
    ionisedpart[ii].x             = partdata[ii].r[0];   // Particle x from gandalf
    ionisedpart[ii].y             = partdata[ii].r[1];   // Particle y from gandalf
    ionisedpart[ii].z             = partdata[ii].r[2];   // Particle z from gandalf
    ionisedpart[ii].rho           = partdata[ii].rho;    // Particle density from gandalf
    ionisedpart[ii].h             = partdata[ii].h;      // Particle h from gandalf
    ionisedpart[ii].sink          = 0;                   // Is the particle a sink
    ionisedpart[ii].t             = tn;                  // Neutral gas temp
    ionisedpart[ii].fionised      = 0;                   // Is the particle ionised at all
    ionisedpart[ii].neighstorcont = 0;                   // Control var. for building neibstore
#ifndef RAD_OPTIMISE
    ionisedpart[ii].neighstor     = new int[200];  //Array containing references to all particles that consider this one a neighbour
    ionisedpart[ii].angle         = new DOUBLE[nos];
    ionisedpart[ii].neigh         = new int[nos];
    ionisedpart[ii].photons       = new DOUBLE[nos];
    ionisedpart[ii].prob          = new DOUBLE[nos];
    ionisedpart[ii].checked       = new int[nos];
    ionisedpart[ii].ionised       = new int[nos];
    ionisedpart[ii].rad_pre_acc   = new DOUBLE[ndim];
#endif

    for(jj=0; jj<ndim; jj++) {
      ionisedpart[ii].rad_pre_acc[jj] = 0;
    }
    ionisation_fraction[ii].resize(newnos); //Resize operation

    // Correctly filling new spaces with 0
    for (jj=0; jj<(newnos - ionisation_fraction[ii].size()); jj++) {
      ionisedpart[ii].ionised[ionisation_fraction[ii].size()+jj-1]=0;
    }

    // Copying in current values for active sources
    for (jj=0; jj<nos; jj++) {
      ionisedpart[ii].ionised[jj] = ionisation_fraction[ii][newnosid[jj]];
    }

    for (jj=0; jj<nos; jj++) {
      ionisedpart[ii].angle[jj]   = 2.*pi;
      ionisedpart[ii].neigh[jj]   = N;
      ionisedpart[ii].checked[jj] = 0;
      ionisedpart[ii].prob[jj]    = 0;
    }
    for(jj=0; jj<200; jj++) {
      ionisedpart[ii].neighstor[jj]=N;
    }
  }


  // Add sink propertys to sink particles
#pragma omp parallel for private(ii,jj)
  for (ii=0; ii<nos; ii++) {
    ionisedpart[N-nos+ii].t             = ti;                        // Set stars to ionised gas temp for smoothing
    ionisedpart[N-nos+ii].sink          = 1;                         // Is the particle a sink
    ionisedpart[N-nos+ii].x             = ndata[newnosid[ii]]->r[0]; // Source x from gandalf
    ionisedpart[N-nos+ii].y             = ndata[newnosid[ii]]->r[1]; // Source y from gandalf
    ionisedpart[N-nos+ii].z             = ndata[newnosid[ii]]->r[2]; // Source z from gandalf
    ionisedpart[N-nos+ii].fionised      = 1;                         // As this is source it is marked as ionised
    ionisedpart[N-nos+ii].neighstorcont = 0;                         // Controle varible for building of neighstore
    sinkid[ii] = N - nos + ii;
    ndot[ii] = pow(m_per_hatom ,2)*ndata[newnosid[ii]]->NLyC/(4.*pi*arecomb)*scale;  //Ndot table s^-1

#ifndef RAD_OPTIMISE
    ionisedpart[N-nos+ii].neighstor   = new int[200];          //Array containing references to all particles that consider this one a neighbour ionisedpart[N-nos+ii].angle=new DOUBLE[nos];
    ionisedpart[N-nos+ii].neigh       = new int[nos];
    ionisedpart[N-nos+ii].photons     = new DOUBLE[nos];
    ionisedpart[N-nos+ii].checked     = new int[nos];
    ionisedpart[N-nos+ii].prob        = new DOUBLE[nos];
    ionisedpart[N-nos+ii].ionised     = new int[nos];
    ionisedpart[N-nos+ii].angle       = new DOUBLE[nos];
    ionisedpart[N-nos+ii].rad_pre_acc = new DOUBLE[3];
#endif
    for(jj=0; jj<3; jj++) {
      ionisedpart[N-nos+ii].rad_pre_acc[jj]=0;
    }
    for(jj=0;jj<nos;jj++) {
      ionisedpart[N-nos+ii].angle[jj]   = 2.*pi;
      ionisedpart[N-nos+ii].neigh[jj]   = N;
      ionisedpart[N-nos+ii].checked[jj] = 0;
      ionisedpart[N-nos+ii].prob[jj]    = 0;
      ionisedpart[N-nos+ii].ionised[jj] = 0;
    }
    ionisedpart[N-nos+ii].ionised[ii] =1;                   //Set the sink location in the ionised array to 1
    for(jj=0; jj<200; jj++) {
      ionisedpart[N-nos+ii].neighstor[jj] = N;
    }
  }


  // Add control ionisedpart to which all particles are initally linked
  // Large values are used so it is never linked in the furture
  ionisedpart[N].x    = 100000.;
  ionisedpart[N].y    = 100000.;
  ionisedpart[N].z    = 100000.;
  ionisedpart[N].rho  = 1e100;
  ionisedpart[N].sink = 0;
#ifndef RAD_OPTIMISE
  ionisedpart[N].angle        = new DOUBLE[nos];
  ionisedpart[N].neigh        = new int[nos];
  ionisedpart[N].photons      = new DOUBLE[nos];
  ionisedpart[N].checked      = new int[nos];
  ionisedpart[N].prob         = new DOUBLE[nos];
  ionisedpart[N].ionised      = new int[nos];
  ionisedpart[N].rad_pre_accv = new DOUBLE[ndim];
#endif
  for (jj=0; jj<ndim; jj++) {
    ionisedpart[N].rad_pre_acc[jj] = 0;
  }
  for (jj=0; jj<nos; jj++) {
    ionisedpart[N].angle[jj]   = 2.*pi;
    ionisedpart[N].neigh[jj]   = N;
    ionisedpart[N].checked[jj] = 0;
    ionisedpart[N].prob[jj]    = 0;
    ionisedpart[N].ionised[jj] = 0;
  }


  // Debug message
  if (debug == 1) {
    gettimeofday(&end, NULL);

    delta = (((end.tv_sec  - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec) / 1.e6)-delta;
    cout << delta << "s to Particle arrays created" << endl;
  }

  // Find the closest source neighbour in chain for each particle
  //-----------------------------------------------------------------------------------------------
#ifdef RAD_OPTIMISE
  int maxneigh = 200;
  int Nneigb;
  int *current_paricle_nn = new int[maxneigh];
#else
  int maxneigh = max(N/8, 2000);
  int current_paricle_nn[maxneigh];              // Working particle neighbour list
  int Nneigb;                                    // No. of neighbours found by neib roughtine
#endif
  DOUBLE dot,mag,angletest; 				//holding variables
  DOUBLE distanceii,distancejj,temp_radius,temp_radius2;	//Distance between the source of the test particle ii and the neighbour particle jj

  // Begin neighbour find.  Loop over all particles
#pragma omp parallel for private(current_paricle_nn,ii,jj,tt,Nneigb,temp_radius,temp_radius2,dot,mag,angletest,pp,distanceii,distancejj) //Initiate openmp
  for (ii=0; ii<N; ii++) {

    // Find NNnumber nearest neighbours
    Nneigb = neib->GetGatherNeighbourList(partdata[ii].r, 2.0*partdata[ii].h,
                                          partgen, N, maxneigh, current_paricle_nn);
#ifdef RAD_OPTIMISE
    while (Nneigb == -1) {
      maxneigh = maxneigh*2;
      delete [] current_paricle_nn;
      current_paricle_nn = new int[maxneigh];
      Nneigb = neib->GetGatherNeighbourList(partdata[ii].r,partdata[ii].h*3.,partgen,N,maxneigh,current_paricle_nn); //Find NNnumber nearest neighbours
    }
#else
    if (!(Nneigb >= 0 && Nneigb <= maxneigh)) {
      cout << "Invalid value of Nneighb : " << Nneigb << "   " << maxneigh << endl;
      exit(0);
    }
#endif

    // Checks which sources lie within the smoothing length
    for (tt=0;tt<nos;tt++){
      if(sqrt(pow((ionisedpart[sinkid[tt]].x - ionisedpart[ii].x), 2) +
              pow((ionisedpart[sinkid[tt]].y - ionisedpart[ii].y), 2) +
              pow((ionisedpart[sinkid[tt]].z - ionisedpart[ii].z), 2)) <= ionisedpart[ii].h*2.) {
      current_paricle_nn[Nneigb]=sinkid[tt];
      Nneigb = Nneigb + 1;
      }
    }

    // For each of these neighbours write my id to there neighstorcont array
    for (jj=0; jj<Nneigb; jj++) {

      // If there is still room in the neighstorecont array.
      // (This should always be true but this ensures no seg fault)
      if (ionisedpart[current_paricle_nn[jj]].neighstorcont < 200) {
        ionisedpart[current_paricle_nn[jj]].neighstor[ionisedpart[current_paricle_nn[jj]].neighstorcont] = ii;
        ionisedpart[current_paricle_nn[jj]].neighstorcont += 1;
      }
      else {
        temp_radius = sqrt(pow(ionisedpart[ii].x - ionisedpart[current_paricle_nn[jj]].x, 2) +
                           pow(ionisedpart[ii].y - ionisedpart[current_paricle_nn[jj]].y, 2) +
                           pow(ionisedpart[ii].z - ionisedpart[current_paricle_nn[jj]].z, 2));
        for (tt=0; tt<200; tt++) {
          temp_radius2 = sqrt(pow(ionisedpart[ii].x - ionisedpart[ionisedpart[current_paricle_nn[jj]].neighstor[tt]].x, 2) +
                              pow(ionisedpart[ii].y - ionisedpart[ionisedpart[current_paricle_nn[jj]].neighstor[tt]].y, 2) +
                              pow(ionisedpart[ii].z - ionisedpart[ionisedpart[current_paricle_nn[jj]].neighstor[tt]].z, 2));
          if (temp_radius > temp_radius2) {
            ionisedpart[current_paricle_nn[jj]].neighstor[tt] = ii;
            break;
          }
        }
      }

      for (pp=0; pp<nos; pp++) {

        // Work out the distances for both test and candidate particle
        distanceii = sqrt(pow(ionisedpart[ii].x - ionisedpart[sinkid[pp]].x, 2) +
                          pow(ionisedpart[ii].y - ionisedpart[sinkid[pp]].y, 2) +
                          pow(ionisedpart[ii].z - ionisedpart[sinkid[pp]].z, 2));
        distancejj = sqrt(pow(ionisedpart[current_paricle_nn[jj]].x - ionisedpart[sinkid[pp]].x, 2) +
                          pow(ionisedpart[current_paricle_nn[jj]].y - ionisedpart[sinkid[pp]].y, 2) +
                          pow(ionisedpart[current_paricle_nn[jj]].z - ionisedpart[sinkid[pp]].z, 2));

        // If the candidate particle is closer than the test particle it is a candidate (Also has controle so a particle cant be its own neighbour)
        if (distancejj < distanceii && ii != current_paricle_nn[jj]) {

          // Use the dot product to work out the angle between the conneting line and the neighbour particle
          dot = ((ionisedpart[ii].x - ionisedpart[sinkid[pp]].x)*(ionisedpart[current_paricle_nn[jj]].x - ionisedpart[sinkid[pp]].x)) +
                ((ionisedpart[ii].y - ionisedpart[sinkid[pp]].y)*(ionisedpart[current_paricle_nn[jj]].y - ionisedpart[sinkid[pp]].y)) +
                ((ionisedpart[ii].z - ionisedpart[sinkid[pp]].z)*(ionisedpart[current_paricle_nn[jj]].z - ionisedpart[sinkid[pp]].z));
          mag = (sqrt(pow(ionisedpart[ii].x - ionisedpart[sinkid[pp]].x, 2) +
                      pow(ionisedpart[ii].y - ionisedpart[sinkid[pp]].y, 2) +
                      pow(ionisedpart[ii].z - ionisedpart[sinkid[pp]].z, 2)))*
                (sqrt(pow(ionisedpart[current_paricle_nn[jj]].x - ionisedpart[sinkid[pp]].x, 2) +
                      pow(ionisedpart[current_paricle_nn[jj]].y - ionisedpart[sinkid[pp]].y, 2) +
                      pow(ionisedpart[current_paricle_nn[jj]].z - ionisedpart[sinkid[pp]].z, 2)));
          angletest = acos(dot/mag);

          // If the partcle is a sink set angle to be max (Stops non-relavant sinks becoming a neighbour)
          if (ionisedpart[current_paricle_nn[jj]].sink == 1) angletest = 2.*pi;

          // If the neighbour is the relavant source set angletest to be neg hence particle will be neighbour
          if (current_paricle_nn[jj] == sinkid[pp]) angletest = -1e50;

          // If Neighbour is closest so far to the connecting line then set current neighbour to be the neighbour
          if (angletest<ionisedpart[ii].angle[pp]) {
            ionisedpart[ii].angle[pp] = angletest;  //Set new comparison angle to be that of the neighbour
            ionisedpart[ii].neigh[pp] = current_paricle_nn[jj];  //Write particle id to neigh array
          }
        }
      }
    }
  }

  // Debug message
  if (debug == 1) {
    gettimeofday(&end, NULL);
    delta = (((end.tv_sec  - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec) / 1.e6)-delta;
    cout << delta << "s to neigbour step one compleate" << endl;
  }


  // Loop over all particles
#pragma omp parallel for private(current_paricle_nn,dot,mag,angletest,jj,pp,distanceii,distancejj,ii,tt,Nneigb,temp_radius,temp_radius2)
  for (ii=0;ii<N;ii++) {

    // Loop over all sources
    for (pp=0; pp<nos; pp++) {

      // For each particle that considers me a neighbour
      for (jj=0; jj<ionisedpart[ii].neighstorcont; jj++) {

        //Work out the distances for both the test and candidate particle
        distanceii = sqrt(pow(ionisedpart[ii].x - ionisedpart[sinkid[pp]].x, 2) +
                          pow(ionisedpart[ii].y - ionisedpart[sinkid[pp]].y, 2) +
                          pow(ionisedpart[ii].z - ionisedpart[sinkid[pp]].z, 2));
        distancejj = sqrt(pow(ionisedpart[ionisedpart[ii].neighstor[jj]].x - ionisedpart[sinkid[pp]].x, 2) +
                          pow(ionisedpart[ionisedpart[ii].neighstor[jj]].y - ionisedpart[sinkid[pp]].y, 2) +
                          pow(ionisedpart[ionisedpart[ii].neighstor[jj]].z - ionisedpart[sinkid[pp]].z, 2));

        //If the candidate particle is closer than the test particle it is a candidate.
        // (Also has control so a particle cant be its own neighbour)
        if (distancejj < distanceii and ii != ionisedpart[ii].neighstor[jj]) {

          // Use the dot product to work out the angle between the conneting line and the candidate particle
          dot = ((ionisedpart[ii].x - ionisedpart[sinkid[pp]].x)*(ionisedpart[ionisedpart[ii].neighstor[jj]].x - ionisedpart[sinkid[pp]].x)) +
                ((ionisedpart[ii].y - ionisedpart[sinkid[pp]].y)*(ionisedpart[ionisedpart[ii].neighstor[jj]].y - ionisedpart[sinkid[pp]].y)) +
                ((ionisedpart[ii].z - ionisedpart[sinkid[pp]].z)*(ionisedpart[ionisedpart[ii].neighstor[jj]].z - ionisedpart[sinkid[pp]].z));
          mag = (sqrt(pow(ionisedpart[ii].x - ionisedpart[sinkid[pp]].x, 2) +
                      pow(ionisedpart[ii].y - ionisedpart[sinkid[pp]].y, 2) +
                      pow(ionisedpart[ii].z - ionisedpart[sinkid[pp]].z, 2)))*
                (sqrt(pow(ionisedpart[ionisedpart[ii].neighstor[jj]].x - ionisedpart[sinkid[pp]].x, 2) +
                      pow(ionisedpart[ionisedpart[ii].neighstor[jj]].y - ionisedpart[sinkid[pp]].y, 2) +
                      pow(ionisedpart[ionisedpart[ii].neighstor[jj]].z - ionisedpart[sinkid[pp]].z, 2)));
          angletest = acos(dot/mag);

          // If the partcle is a sink set angle to be max (Stops non-relavant sinks becoming a neighbour)
          if (ionisedpart[ionisedpart[ii].neighstor[jj]].sink == 1) angletest = 2.*pi;

          // If the neighbour is the relavant source set angletest to be neg hence
          // particle will be neighbour (Ensures link to sources)
          if (ionisedpart[ii].neighstor[jj] == sinkid[pp]) angletest = -1e50;

          // If candidate is closest so far to the connecting line then set current candidate to be the neighbour
          if (angletest<ionisedpart[ii].angle[pp]) {
            ionisedpart[ii].angle[pp] = angletest;                        //Set new comparison angle to be that of the neighbour
            ionisedpart[ii].neigh[pp] = ionisedpart[ii].neighstor[jj];    //Write particle id to neigh array
          }
        }
      }
    }
  }


  if (debug == 1) {
    gettimeofday(&end, NULL);
    delta = (((end.tv_sec  - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec) / 1.e6)-delta;
    cout << delta << "s to All neighbours found" << endl;
  }


  // Begin working out if particle are ionised
  int change = N;                           //Variable for checking convergence
  int finalcheck = 0;                       //Finalcheck is a controle that initials the final check of the solution.

  //loop until no changes are made (We have converged)
  //-----------------------------------------------------------------------------------------------
  while (change!=0 or finalcheck==0) {

    if (debug == 1) cout<<"The number of adjustments made is " << change << endl;

    // Set state of finalcheck (If no changes made in last iteration then carry out final check)
    if (change == 0) finalcheck = 1;

    // Set change to zero (This will increase as changes are made in the calculation)
    change = 0;

    // Open mp initalised.  Loop over each particle
#pragma omp parallel for schedule(dynamic)  //Begin open mp in dynamic mode i.e one argument is passed at a time rather than a block
    for (ii=1; ii<N; ii++) {

      // Call fucntion to determine if test particle is ionised
      photoncount(ionisedpart, sinkid, ndot, N, nos, ii, change);
    }

#pragma omp parallel for private(ii,pp)
    for (ii=1; ii<N; ii++) {
      for(pp=0; pp<nos; pp++) {
        ionisedpart[ii].checked[pp]=0;
      }
    }
  }
  //-----------------------------------------------------------------------------------------------


  if (debug==1){
    gettimeofday(&end, NULL);
    delta = (((end.tv_sec  - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec) / 1.e6)-delta;
    cout << delta << "s to Iterations compleate" << endl;
  }


  // Smooth the temperature for all particles
  //-----------------------------------------------------------------------------------------------
  DOUBLE rad,s,w,invmu;
  if (smoothing) {

#pragma omp parallel for private(jj,ii,rad,s,w,invmu,current_paricle_nn,Nneigb) //Initalise openmp
    for (ii=0; ii<N; ii++) {

      if (ionisedpart[ii].fionised == 1) {

        //Find NNnumber nearest neighbours
        Nneigb = neib->GetGatherNeighbourList(partdata[ii].r, partdata[ii].h*3., partgen,
                                              N, maxneigh, current_paricle_nn);

#ifdef RAD_OPTIMISE
        while (Nneigb == -1) {
          maxneigh=maxneigh*2;
          delete [] current_paricle_nn;
          current_paricle_nn = new int[maxneigh];
          Nneigb = neib->GetGatherNeighbourList(partdata[ii].r,partdata[ii].h*3.,partgen,N,maxneigh,current_paricle_nn); //Find NNnumber nearest neighbours
        }
#endif

        // For each of the neighbours
        for (jj=0; jj<Nneigb; jj++) {

          if (ionisedpart[current_paricle_nn[jj]].fionised == 0) {

            rad = (sqrt(pow(ionisedpart[current_paricle_nn[jj]].x - ionisedpart[ii].x, 2) +
                        pow(ionisedpart[current_paricle_nn[jj]].y - ionisedpart[ii].y, 2) +
                        pow(ionisedpart[current_paricle_nn[jj]].z - ionisedpart[ii].z, 2)));
            s   = rad/(ionisedpart[ii].h*1.5);          //Work out s for smoothing kernal

            // Work out w for the kernal
            if (s < 1) w = 1 - (3./2.)*pow(s,2.) + (3./4.)*pow(s,3.);
            else if (s < 2) w = (1./4.)*pow(2-s,3.);
            else w = 0;

            if (ionisedpart[current_paricle_nn[jj]].t < ti*w) {
              ionisedpart[current_paricle_nn[jj]].t = ti*w;
            }
          }
          else {
            ionisedpart[current_paricle_nn[jj]].t = ti;
          }
        }
      }
    }
  }
  //-----------------------------------------------------------------------------------------------


  if (debug == 1) {
    gettimeofday(&end, NULL);
    delta = (((end.tv_sec  - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec) / 1.e6)-delta;
    cout << delta << "s to Temperatures smoothed" << endl;
  }


  DOUBLE theta,thi,photon_acceleration;
  //-----------------------------------------------------------------------------------------------
#pragma omp parallel for private(jj,ii,rad,s,w,invmu,current_paricle_nn,Nneigb,theta,thi,photon_acceleration) //Initalise openmp
  for (ii=0; ii<N; ii++) {

    // If the particle is ionised then its temperature must be the ionised temperature
    if (ionisedpart[ii].fionised == 1) ionisedpart[ii].t = ti;

    // If the particles temp is less than the neutral temp because of smoothing set it back to the neutral temp
    if (ionisedpart[ii].t < tn) ionisedpart[ii].t = tn;

    // Work out corrected inverted mean gas particle mass to calculate the internal energy
    invmu=(((ionisedpart[ii].t - tn)/mu_ion) + ((ti - ionisedpart[ii].t)/mu_bar))/(ti - tn);
    partdata[ii].mu_bar = 1/invmu;
    ionisedpart[ii].u=ionisedpart[ii].t/tempscale/gammam1*invmu;

    // Set particle ionisation state
    if (ionisedpart[ii].t == tn) {
      partdata[ii].ionstate = 0;
    }
    else if (ionisedpart[ii].fionised == 1) {
      partdata[ii].ionstate = 2;
    }
    else  {
      partdata[ii].ionstate = 1;
    }
    partdata[ii].u = ionisedpart[ii].u;

  //Working out radiation pressure (NOT COMPLEATE)
  //photon_acceleration=3.4455561764e-34*rad_cont*ionisedpart[ii].rho/(mu_bar*mu_bar);

  //for(jj=0;jj<nos;jj++)
    //{
    //Copying ionised state over to holding array
    //ionisation_fraction[ii][newnosid[jj]]=ionisedpart[ii].ionised[jj];
    //if(ionisedpart[ii].fionised==1)
      //{
      //theta=atan((ionisedpart[ii].y-ionisedpart[sinkid[jj]].y)/(ionisedpart[ii].z-ionisedpart[sinkid[jj]].z));
      //thi=atan((ionisedpart[ii].y-ionisedpart[sinkid[jj]].y)/(ionisedpart[ii].x-ionisedpart[sinkid[jj]].x));
      //ionisedpart[ii].rad_pre_acc[0]=ionisedpart[ii].rad_pre_acc[0]+ionisedpart[ii].prob[jj]*photon_acceleration*sin(theta)*cos(thi);
      //ionisedpart[ii].rad_pre_acc[1]=ionisedpart[ii].rad_pre_acc[1]+ionisedpart[ii].prob[jj]*photon_acceleration*sin(theta)*sin(thi);
      //ionisedpart[ii].rad_pre_acc[2]=ionisedpart[ii].rad_pre_acc[2]+ionisedpart[ii].prob[jj]*photon_acceleration*cos(theta);
     // }
   // }

  //partdata[ii].rad_pres[0]=0;//ionisedpart[ii].rad_pre_acc[0];
  //partdata[ii].rad_pres[1]=0;//ionisedpart[ii].rad_pre_acc[1];
  //partdata[ii].rad_pres[2]=0;//ionisedpart[ii].rad_pre_acc[2];

  }
  //-----------------------------------------------------------------------------------------------


#ifndef RAD_OPTIMISE
//Memory De-allocation
  delete [] sinkid;
  delete [] ndot;
  delete [] newnosid;

#pragma omp parallel for private(ii)
  for(ii=0; ii<N; ii++) {
  delete [] ionisedpart[ii].angle;
  delete [] ionisedpart[ii].neigh;
  delete [] ionisedpart[ii].photons;
  delete [] ionisedpart[ii].neighstor;
  delete [] ionisedpart[ii].checked;
  delete [] ionisedpart[ii].prob;
  delete [] ionisedpart[ii].ionised;
  delete [] ionisedpart[ii].rad_pre_acc;
  }
#else
  delete [] current_paricle_nn;
#endif
  delete [] ionisedpart;


  gettimeofday(&end, NULL);
  delta = ((end.tv_sec  - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec) / 1.e6;
  cout << "The time taken to calculate ionisation temperatures = " << delta << " s" << endl;

  if (debug==1 or debug==2) {
    ofstream myfile;
    myfile.open ("timing.dat",ios::app);
    myfile <<delta<<"\n";
  }

  return;
}


template class MultipleSourceIonisation<3,MeshlessFVParticle>;
template class MultipleSourceIonisation<3,GradhSphParticle>;
template class MultipleSourceIonisation<3,SM2012SphParticle>;
