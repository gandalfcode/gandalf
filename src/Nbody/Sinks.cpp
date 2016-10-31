//=================================================================================================
//  Sinks.cpp
//  All routines for creating new sinks and accreting gas and updating all
//  sink particle propterties.
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "Parameters.h"
#include "Sph.h"
#include "Nbody.h"
#include "Sinks.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
using namespace std;



//=================================================================================================
//  Sinks::Sinks()
/// Sinks class constructor
//=================================================================================================
template <int ndim>
Sinks<ndim>::Sinks(NeighbourSearch<ndim>* _neibsearch)
{
  neibsearch       = _neibsearch;
  allocated_memory = false;
  Nsink            = 0;
  Nsinkmax         = 0;
}



//=================================================================================================
//  Sinks::~Sinks()
/// Sinks class destructor
//=================================================================================================
template <int ndim>
Sinks<ndim>::~Sinks()
{
  DeallocateMemory();
}



//=================================================================================================
//  Sinks::AllocateMemory
/// Allocate all memory required for storing sink particle data.
//=================================================================================================
template <int ndim>
void Sinks<ndim>::AllocateMemory(int N)
{
  if (N > Nsinkmax) {
    if (allocated_memory) DeallocateMemory();
    Nsinkmax = N;
    sink = new SinkParticle<ndim>[Nsinkmax];
    allocated_memory = true;
  }

  return;
}



//=================================================================================================
//  Sinks::DeallocateMemory
/// Deallocate all sink particle arrays
//=================================================================================================
template <int ndim>
void Sinks<ndim>::DeallocateMemory(void)
{
  if (allocated_memory) delete[] sink;
  allocated_memory = false;

  return;
}



//=================================================================================================
//  Sinks::SearchForNewSinkParticles
/// Searches through all SPH particles for new sink particle candidates, and
/// if a particle satisfies all tests, then a sink is created.
//=================================================================================================
template <int ndim>
void Sinks<ndim>::SearchForNewSinkParticles
 (const int n,                         ///< [in] Current integer time
  const FLOAT t,                       ///< [in] Current time
  Hydrodynamics<ndim> *hydro,          ///< [inout] Object containing SPH ptcls
  Nbody<ndim> *nbody)                  ///< [inout] Object containing star ptcls
{
  bool sink_flag;                      // Flag if particle is to become a sink
  int i;                               // Particle counter
  int isink;                           // i.d. of hydro particle to form sink from
  int k;                               // Dimension counter
  int s;                               // Sink counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT rho_max;                       // Maximum density of sink candidates

  debug2("[Sinks::SearchForNewSinkParticles]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SEARCH_NEW_SINKS");


  // Continuous loop to search for new sinks.  If a new sink is found, then repeat
  // entire process to search for other sinks on current timestep.
  // If no sinks are found, then exit and return to main program.
  //===============================================================================================
  do {
    isink = -1;
    rho_max = (FLOAT) 0.0;

    // Loop over all SPH particles finding the particle with the highest
    // density that obeys all of the formation criteria, if any do.
    //---------------------------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      sink_flag = true;
      Particle<ndim>& part = hydro->GetParticlePointer(i);

      // Make sure we don't include dead particles
      if (part.flags.is_dead()) continue;

      // Only consider hydro particles located at a local potential minimum
      if (!part.flags.check_flag(potmin)) continue;

      // If density of a hydro particle is too low, skip to next particle
      if (part.rho < rho_sink) continue;

      // Make sure candidate particle is at the end of its current timestep
      if (n%part.nstep != 0) continue;

      // If hydro particle neighbours a nearby sink, skip to next particle
      for (s=0; s<Nsink; s++) {
        for (k=0; k<ndim; k++) dr[k] = part.r[k] - sink[s].star->r[k];
        drsqd = DotProduct(dr, dr, ndim);
        if (drsqd < pow(sink_radius*part.h + sink[s].radius, 2)) sink_flag = false;
      }
      if (!sink_flag) continue;

      // If candidate particle has passed all the tests, then check if it is
      // the most dense candidate.  If yes, record the particle id and density
      if (sink_flag && part.rho > rho_max) {
        isink = i;
        rho_max = part.rho;
      }

    }
    //---------------------------------------------------------------------------------------------

#if defined MPI_PARALLEL
    // We need to know what the other processors have found
    vector<FLOAT> rho_maxs(mpicontrol->Nmpi);
    MPI_Allgather(&rho_max, 1, GANDALF_MPI_FLOAT,
                  &rho_maxs[0], 1, GANDALF_MPI_FLOAT, MPI_COMM_WORLD);
    const int proc_max = std::max_element(rho_maxs.begin(), rho_maxs.end()) - rho_maxs.begin();
    const FLOAT global_rho_max = rho_maxs[proc_max];
    if (global_rho_max > 0.0) {
      // A sink is being created, but not on this processor - mark it
      if (rho_max < global_rho_max) isink = -2;
    }
#endif


    // If all conditions have been met, then create a new sink particle.
    // Also, set minimum sink smoothing lengths
    if (isink >= 0) {
      Particle<ndim>& part_sink = hydro->GetParticlePointer(isink);
      hydro->hmin_sink = min(hydro->hmin_sink, part_sink.h);
      CreateNewSinkParticle(isink, t, part_sink, hydro, nbody);
    }
#if defined MPI_PARALLEL
    if (isink != -1) {
      // The owner of the new sink broadcasts it to everyone
      MPI_Bcast(&sink[Nsink], sizeof(SinkParticle<ndim>), MPI_BYTE,proc_max, MPI_COMM_WORLD);
      MPI_Bcast(&nbody->stardata[nbody->Nstar], sizeof(StarParticle<ndim>),
                MPI_BYTE, proc_max, MPI_COMM_WORLD);
      if (isink == -2) {
        sink[Nsink].star = &(nbody->stardata[nbody->Nstar]);
        nbody->nbodydata[nbody->Nnbody] = &(nbody->stardata[nbody->Nstar]);
      }
    }
#endif

    // Calculate total mass inside sink (direct sum for now since this is not computed that often).
    if (isink != -1) {
      sink[Nsink].mmax = (FLOAT) 0.0;
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        if (part.flags.is_dead()) continue;
        for (k=0; k<ndim; k++) dr[k] = sink[Nsink].star->r[k] - part.r[k];
        drsqd = DotProduct(dr,dr,ndim);
        if (drsqd < pow(sink[Nsink].radius,2)) sink[Nsink].mmax += part.m;
      }
#if defined MPI_PARALLEL
      MPI_Allreduce(MPI_IN_PLACE, &(sink[Nsink].mmax), 1,
                    GANDALF_MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif
    }

    if (isink >= 0) {
      Particle<ndim>& part_sink = hydro->GetParticlePointer(isink);
      cout << "--------------------------------------------------------------------------" << endl;
      cout << "Created new sink particle : " << isink << "     Nsink : " << Nsink+1  << endl;
      cout << "radius : " << sink[Nsink].radius << "    " << part_sink.h << endl;
      cout << "m : " << sink[Nsink].star->m << "    mmax : " << sink[Nsink].mmax << endl;
      cout << "r : " << sink[Nsink].star->r[0] << "   " << sink[Nsink].star->r[0] << endl;
      cout << "--------------------------------------------------------------------------" << endl;
    }

    if (isink != -1) {
      nbody->Nstar++;
      nbody->Nnbody++;
      Nsink++;
    }

  } while (isink != -1);
  //===============================================================================================


  return;
}



//=================================================================================================
//  Sinks::CreateNewSinkParticle
/// Create a new sink particle from specified SPH particle 'isink' and then
/// removes particle from main arrays.
//=================================================================================================
template <int ndim>
void Sinks<ndim>::CreateNewSinkParticle
 (const int isink,                     ///< [in]    i.d. of the above SPH particle
  const FLOAT t,                       ///< [in]    Current time
  Particle<ndim> &part,                ///< [inout] SPH particle to be turned in a sink particle
  Hydrodynamics<ndim> *hydro,          ///< [inout] Object containing SPH ptcls
  Nbody<ndim> *nbody)                  ///< [inout] Object containing star ptcls
{
  int k;                               // Dimension counter

  debug2("[Sinks::CreateNewSinkParticle]");

  // If we've reached the maximum number of sinks, then throw exception
  if (Nsink == Nsinkmax || nbody->Nstar == nbody->Nstarmax) {
    cout << "Run out of memory : " << Nsink << "    " << Nsinkmax << endl;
    ExceptionHandler::getIstance().raise("Error : run out of memory for new sinks");
  }

  // First create new star and set N-body pointer to star
  sink[Nsink].star = &(nbody->stardata[nbody->Nstar]);
  nbody->nbodydata[nbody->Nnbody] = &(nbody->stardata[nbody->Nstar]);

  // Calculate new sink radius depending on chosen sink parameter
  if (sink_radius_mode == "fixed") {
    sink[Nsink].radius = sink_radius;
  }
  else if (sink_radius_mode == "hmult") {
    sink[Nsink].radius = sink_radius*part.h;
  }
  else {
    sink[Nsink].radius = hydro->kernp->kernrange*part.h;
  }


  // Calculate all other sink properties based on radius and SPH particle properties
  sink[Nsink].star->h            = hydro->kernp->invkernrange*sink[Nsink].radius;
  sink[Nsink].star->invh         = (FLOAT) 1.0/part.h;
  sink[Nsink].star->radius       = sink[Nsink].radius;
  //sink[Nsink].star->hfactor      = pow(sink[Nsink].star->invh,ndim);
  sink[Nsink].star->m            = part.m;
  sink[Nsink].star->gpot         = part.gpot;
  sink[Nsink].star->gpe_internal = (FLOAT) 0.0;
  sink[Nsink].star->dt           = part.dt;
  sink[Nsink].star->tlast        = t;
  sink[Nsink].star->nstep        = part.nstep;
  sink[Nsink].star->nlast        = part.nlast;
  sink[Nsink].star->level        = part.level;
  sink[Nsink].star->active       = part.flags.check_flag(active);
  sink[Nsink].star->Ncomp        = 1;
  for (k=0; k<ndim; k++) sink[Nsink].star->r[k]     = part.r[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->v[k]     = part.v[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->a[k]     = part.a[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->adot[k]  = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) sink[Nsink].star->a2dot[k] = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) sink[Nsink].star->a3dot[k] = (FLOAT) 0.0;
  for (k=0; k<ndim; k++) sink[Nsink].star->r0[k]    = part.r0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->v0[k]    = part.v0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->a0[k]    = part.a0[k];
  for (k=0; k<ndim; k++) sink[Nsink].star->adot0[k] = (FLOAT) 0.0; //part.adot0[k];
  for (k=0; k<3; k++) sink[Nsink].angmom[k]         = (FLOAT) 0.0;


  // Remove SPH particle from main arrays
  part.m      = (FLOAT) 0.0;
  part.flags.unset_flag(active);
  part.flags.set_flag(dead);

  return;
}



//=================================================================================================
//  Sinks::AcceteMassToSinks
/// Identify all SPH particles inside sinks and accrete some fraction (or all)
/// of the gas mass to the sinks if selected accretion criteria are satisfied.
//=================================================================================================
template <int ndim>
void Sinks<ndim>::AccreteMassToSinks
 (const int n,                         ///< [in] Integer timestep
  const FLOAT timestep,                ///< [in] Minimum timestep size
  Hydrodynamics<ndim> *hydro,          ///< [inout] Object containing SPH ptcls
  Nbody<ndim> *nbody)                  ///< [inout] Object containing star ptcls
{
  debug2("[Sinks::AccreteMassToSinks]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SINK_ACCRETE_MASS");

  Particle<ndim> *partdata = hydro->GetParticleArray();

  // Allocate local memory and initialise values
  for (int i=0; i<hydro->Ntot; i++) hydro->GetParticlePointer(i).sinkid = -1;
  for (int s=0; s<Nsinkmax; s++) sink[s].Ngas = 0;

#ifdef MPI_PARALLEL
  Box<ndim> mydomain = mpicontrol->MyDomain();
  list<int> ghosts_accreted;
#endif


  // Set-up all parallel threads for computing sink accretion
  //===============================================================================================
#if defined MPI_PARALLEL
#pragma omp parallel default(none) shared(ghosts_accreted,hydro,mydomain,nbody,partdata)
#else
#pragma omp parallel default(none) shared(hydro,nbody,partdata)
#endif
  {
    int i,j,k;                               // Particle and dimension counters
    int Nlist;                               // Max. no of gas particles inside sink
    int Nneib;                               // No. of particles inside sink
    int Nneibmax = 128;                      // Max. no. of particles inside sink
    int s;                                   // Sink counter
    //  int saux;                                // Aux. sink i.d.
    FLOAT asqd;                              // Acceleration squared
    FLOAT dr[ndim];                          // Relative position vector
    FLOAT drmag;                             // Distance
    FLOAT drsqd;                             // Distance squared
    FLOAT dt;                                // Sink/star timestep
    FLOAT dv[ndim];                          // Relative velocity vector
    FLOAT dvtang[ndim];                      // Relative tangential velocity vector
    FLOAT efrac;                             // Energy fraction
    FLOAT macc;                              // Accreted mass
    FLOAT macc_temp;                         // Temp. accreted mass variable
    FLOAT mold;                              // Old mass
    FLOAT mtemp;                             // Aux. mass variable
    //FLOAT rsqdmin;                           // Distance (sqd) to closest sink
    FLOAT rold[ndim];                        // Old sink position
    FLOAT vold[ndim];                        // Old sink velocity
    FLOAT wnorm;                             // Kernel normalisation factor
    int *ilist = new int[Nneibmax];          // ..
    int *neiblist = new int[Nneibmax];       // List of particle ids
    FLOAT *rsqdlist= new FLOAT[Nneibmax];    // Array of particle-sink distances


    // Determine which sink each SPH particle accretes to.  If none, flag -1
    // (note we should really use the tree to compute this)
    //---------------------------------------------------------------------------------------------
#pragma omp for schedule(dynamic,1)
    for (s=0; s<Nsink; s++) {

      // For MPI simulations, only consider sink particles owned by this domain
#if defined MPI_PARALLEL
      if (!ParticleInBox(*(sink[s].star), mydomain)) continue;
#endif

      // Find the list of particles inside
      do {
        Nlist = neibsearch->GetGatherNeighbourList
         (sink[s].star->r, sink[s].radius, partdata, hydro->Nhydro, Nneibmax, neiblist);

        // If there are too many neighbours so the buffers are filled,
        // reallocate the arrays and recompute the neighbour lists.
        if (Nlist == -1) {
          delete[] rsqdlist;
          delete[] neiblist;
          delete[] ilist;
          Nneibmax *= 2;
          ilist = new int[Nneibmax];
          neiblist = new int[Nneibmax];
          rsqdlist = new FLOAT[Nneibmax];
        };

      } while (Nlist == -1);

      // Loop over all potential sink neighbours
      for (j=0; j<Nlist; j++) {
        i = neiblist[j];
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        if (part.flags.is_dead()) continue;

        for (k=0; k<ndim; k++) dr[k] = part.r[k] - sink[s].star->r[k];
        drsqd = DotProduct(dr, dr, ndim);

        if (drsqd <= sink[s].radius*sink[s].radius) {
#pragma omp critical
          part.sinkid = s;
          sink[s].Ngas++;
        }

      }

    }
    //---------------------------------------------------------------------------------------------


    // Calculate the accretion timescale and the total mass accreted from all ptcls for each sink.
    //---------------------------------------------------------------------------------------------
#pragma omp for schedule(dynamic,1)
    for (s=0; s<Nsink; s++) {

      // Only accrete from local sinks
#if defined MPI_PARALLEL
      if (!ParticleInBox(*(sink[s].star), mydomain)) continue;
#endif

      /*cout << "Accreting?? : " << s << "   " << Nsink << "   " << sink[s].Ngas << "    " << n
           << "    " << sink[s].star->nlast << endl;
      cout << "r0 : " << sink[s].star->r[0] << "    " << sink[s].star->r[1] << "    " << sink[s].star->r[2] << endl;
      cout << "v0 : " << sink[s].star->v[0] << "    " << sink[s].star->v[1] << "    " << sink[s].star->v[2] << endl;
      cout << "a0 : " << sink[s].star->a[0] << "    " << sink[s].star->a[1] << "    " << sink[s].star->a[2] << endl;
  */
      // Skip sink if it contains no gas, or unless it's at the beginning of its current step.
      //if (sink[s].Ngas == 0 || !sink[s].star->active) continue;
      //if (sink[s].Ngas == 0 || n%sink[s].star->nstep != 0) continue;
      //if (sink[s].Ngas == 0 || n%sink[s].star->nstep != sink[s].star->nstep/2) continue;

      //if (sink[s].Ngas == 0 || sink[s].star->nlast != n) continue;
      if (sink[s].Ngas == 0 || n%sink[s].star->nstep != 0) continue;


      // Initialise all variables for current sink
      Nneib = 0;
      wnorm = (FLOAT) 0.0;
      sink[s].menc     = (FLOAT) 0.0;
      sink[s].trad     = (FLOAT) 0.0;
      sink[s].tvisc    = (FLOAT) 1.0;
      sink[s].ketot    = (FLOAT) 0.0;
      sink[s].rotketot = (FLOAT) 0.0;
      sink[s].gpetot   = (FLOAT) 0.0;

      Nlist = neibsearch->GetGatherNeighbourList
       (sink[s].star->r, sink[s].radius, partdata, hydro->Nhydro, Nneibmax, neiblist);

      // Loop over all potential sink neighbours
      for (j=0; j<Nlist; j++) {
        i = neiblist[j];
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        if (part.flags.is_dead()) continue;

        if (part.sinkid == s) {
          for (k=0; k<ndim; k++) dr[k] = part.r[k] - sink[s].star->r[k];
          drsqd = DotProduct(dr, dr, ndim);
          if (drsqd > sink[s].radius*sink[s].radius) continue;
          ilist[Nneib]    = i;
          rsqdlist[Nneib] = drsqd;  //*part.m;
          Nneib++;
          part.levelneib = max(part.levelneib, sink[s].star->level);
        }

      }

      // Sort particle ids by increasing distance from the sink
      InsertionSortIds(Nneib, ilist, rsqdlist);


      // Calculate all important quantities (e.g. energy contributions) due to
      // all particles inside the sink
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nneib; j++) {
        i = ilist[j];
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        if (part.flags.is_dead()) continue;

        for (k=0; k<ndim; k++) dr[k] = part.r[k] - sink[s].star->r[k];
        drsqd = DotProduct(dr,dr,ndim);
        drmag = sqrt(drsqd) + small_number;
        for (k=0; k<ndim; k++) dr[k] /= drmag;

        sink[s].menc += part.m;
        wnorm += part.m*hydro->kernp->w0(drmag*sink[s].star->invh)*
          pow(sink[s].star->invh,ndim)/part.rho;

        // Sum total grav. potential energy of all particles inside sink
        sink[s].gpetot += (FLOAT) 0.5*part.m*(sink[s].star->m + sink[s].menc)*
          sink[s].star->invh*hydro->kernp->wpot(drmag*sink[s].star->invh);

        // Compute rotational component of kinetic energy
        for (k=0; k<ndim; k++) dv[k] = part.v[k] - sink[s].star->v[k];
        for (k=0; k<ndim; k++) dvtang[k] = dv[k] - DotProduct(dv,dr,ndim)*dr[k];

        // Compute total and rotational kinetic energies
        sink[s].ketot += part.m*DotProduct(dv,dv,ndim)*
          hydro->kernp->w0(drmag*sink[s].star->invh)*pow(sink[s].star->invh,ndim)/part.rho;
        sink[s].rotketot += part.m*DotProduct(dvtang,dvtang,ndim)*
          hydro->kernp->w0(drmag*sink[s].star->invh)*pow(sink[s].star->invh,ndim)/part.rho;

        // Add contributions to average timescales from particles
        sink[s].tvisc *= pow(sqrt(drmag)/part.sound/part.sound,part.m);
        sink[s].trad += fabs((FLOAT) 4.0*pi*drsqd*part.m*DotProduct(dv,dr,ndim)*
                             hydro->kernp->w0(drmag*sink[s].star->invh)*pow(sink[s].star->invh,ndim));

      }
      //-------------------------------------------------------------------------------------------


      // Normalise SPH sums correctly
      sink[s].ketot *= (FLOAT) 0.5*sink[s].menc/wnorm;
      sink[s].rotketot *= (FLOAT) 0.5*sink[s].menc/wnorm;


      // Calculate the sink accretion timescale and the total amount of mass accreted by sink s
      // this timestep.  If the contained mass is greater than the maximum allowed, accrete the
      // excess mass.  Otherwise, accrete a small amount based on freefall/viscous timescale.
      //-------------------------------------------------------------------------------------------
      if (smooth_accretion == 1) {
        efrac = min((FLOAT) 2.0*sink[s].rotketot/sink[s].gpetot,(FLOAT) 1.0);
        sink[s].tvisc = (sqrt(sink[s].star->m + sink[s].menc)*
                              pow(sink[s].tvisc,(FLOAT) 1.0/sink[s].menc))/alpha_ss;
        sink[s].trad  = sink[s].menc / sink[s].trad;
        sink[s].trot  = twopi*sqrt(pow(sink[s].radius,3)/(sink[s].menc + sink[s].star->m));

        // Finally calculate accretion timescale and mass accreted
        // If there's too much mass inside the sink, artificially increase accretion rate to
        // restore equilibrium (between mass entering sink and that being accreted) quicker.
        sink[s].taccrete = pow(sink[s].trad, (FLOAT) 1.0 - efrac)*pow(sink[s].tvisc,efrac);
        if (sink[s].mmax > small_number && sink[s].menc > sink[s].mmax) {
          sink[s].taccrete *= pow(sink[s].mmax/sink[s].menc,2);
        }
        dt = (FLOAT) sink[s].star->nstep*timestep;
        macc = sink[s].menc*max((FLOAT) 1.0 - (FLOAT) exp(-dt/sink[s].taccrete), (FLOAT) 0.0);

        /*cout << "efrac : " << efrac << "    taccrete : " << sink[s].taccrete << "   "
             << sink[s].tvisc << "    " << sink[s].trad << "    " << sink[s].trot << endl;
        cout << "energy : " << sink[s].ketot << "   " << sink[s].rotketot << "    " << sink[s].gpetot << endl;
        cout << "macc : " << macc << "     macc/mmean : " << macc/hydro->mmean
             << "    " << sink[s].menc << "    mmax : " << sink[s].mmax
             << "    mmax/mmean : " << sink[s].mmax/hydro->mmean << "     dmdt : " << macc/dt << endl;
        */
      }
      else {
        macc = sink[s].menc;
      }


      // Now accrete SPH particles to sink
      //-------------------------------------------------------------------------------------------
      macc_temp = macc;
      for (k=0; k<ndim; k++) rold[k] = sink[s].star->r[k];
      for (k=0; k<ndim; k++) vold[k] = sink[s].star->v[k];
      mold = sink[s].star->m;

      for (k=0; k<ndim; k++) sink[s].star->r[k] *= sink[s].star->m;
      for (k=0; k<ndim; k++) sink[s].star->v[k] *= sink[s].star->m;
      for (k=0; k<ndim; k++) sink[s].star->a[k] *= sink[s].star->m;
      //for (k=0; k<ndim; k++) sink[s].star->adot[k] *= sink[s].star->m;


      // Loop over all neighbouring particles
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nneib; j++) {
        i = ilist[j];

        Particle<ndim>& part = hydro->GetParticlePointer(i);
        if (part.flags.is_dead()) continue;

        mtemp = min(part.m, macc_temp);
        dt = part.dt;

        // Special conditions for total particle accretion
        if (smooth_accretion == 0 || part.m - mtemp < smooth_accrete_frac*hydro->mmean ||
            dt < smooth_accrete_dt*sink[s].trot) {
          mtemp = part.m;
        }
        macc_temp -= mtemp;

        // Now accrete COM quantities to sink particle
        sink[s].star->m += mtemp;
        for (k=0; k<ndim; k++) sink[s].star->r[k] += mtemp*part.r[k];
        for (k=0; k<ndim; k++) sink[s].star->v[k] += mtemp*part.v[k];
        for (k=0; k<ndim; k++) sink[s].star->a[k] += mtemp*part.a[k];
        //for (k=0; k<ndim; k++) sink[s].star->adot[k] += mtemp*part.adot[k];
        sink[s].utot += mtemp*part.u;

        // If we've reached/exceeded the mass limit, do not include more ptcls
        if (macc_temp < small_number) break;
      }
      //-------------------------------------------------------------------------------------------


      // Normalise COM quantities
      for (k=0; k<ndim; k++) sink[s].star->r[k] /= sink[s].star->m;
      for (k=0; k<ndim; k++) sink[s].star->v[k] /= sink[s].star->m;
      for (k=0; k<ndim; k++) sink[s].star->a[k] /= sink[s].star->m;
      //for (k=0; k<ndim; k++) sink[s].star->adot[k] /= sink[s].star->m;

      //if (n%sink[s].star->nstep == 0) {
      for (k=0; k<ndim; k++) sink[s].star->r0[k] = sink[s].star->r[k];
      for (k=0; k<ndim; k++) sink[s].star->v0[k] = sink[s].star->v[k];
      for (k=0; k<ndim; k++) sink[s].star->a0[k] = sink[s].star->a[k];
      //for (k=0; k<ndim; k++) sink[s].star->adot0[k] = sink[s].star->adot[k];
      //}

      // Calculate angular momentum of old COM around new COM
      for (k=0; k<ndim; k++) dr[k] = rold[k] - sink[s].star->r[k];
      for (k=0; k<ndim; k++) dv[k] = vold[k] - sink[s].star->v[k];
      if (ndim == 3) {
        sink[s].angmom[0] += mold*(dr[1]*dv[2] - dr[2]*dv[1]);
        sink[s].angmom[1] += mold*(dr[2]*dv[0] - dr[0]*dv[2]);
        sink[s].angmom[2] += mold*(dr[0]*dv[1] - dr[1]*dv[0]);
      }
      else if (ndim == 2) {
        sink[s].angmom[2] += mold*(dr[0]*dv[1] - dr[1]*dv[0]);
      }


      // Now add angular momentum contribution of individual SPH particles
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nneib; j++) {
        i = ilist[j];

        Particle<ndim>& part = hydro->GetParticlePointer(i);
        if (part.flags.is_dead()) continue;

        mtemp = min(part.m, macc);
        dt = part.dt;

        // Special conditions for total particle accretion
        if (smooth_accretion == 0 || part.m - mtemp < smooth_accrete_frac*hydro->mmean ||
            dt < smooth_accrete_dt*sink[s].trot) {
          mtemp       = part.m;
          part.m      = (FLOAT) 0.0;
          part.flags.set_flag(dead);
          part.flags.unset_flag(active);
        }
        else {
          //part.m -= mtemp;
          hydro->AccreteMassFromParticle(mtemp, part);
        }
#if defined MPI_PARALLEL
        if (i > hydro->Nhydro) {
          // We are accreting a MPI ghost, so we need to record that to transmit it to the owner
#pragma omp critical (ghost_accreted)
          ghosts_accreted.push_back(i);
        }
#endif
        macc -= mtemp;


        // Calculate angular momentum of old COM around new COM
        for (k=0; k<ndim; k++) dr[k] = part.r[k] - sink[s].star->r[k];
        for (k=0; k<ndim; k++) dv[k] = part.v[k] - sink[s].star->v[k];
        if (ndim == 3) {
          sink[s].angmom[0] += mtemp*(dr[1]*dv[2] - dr[2]*dv[1]);
          sink[s].angmom[1] += mtemp*(dr[2]*dv[0] - dr[0]*dv[2]);
          sink[s].angmom[2] += mtemp*(dr[0]*dv[1] - dr[1]*dv[0]);
        }
        else if (ndim == 2) {
          sink[s].angmom[2] += mtemp*(dr[0]*dv[1] - dr[1]*dv[0]);
        }

        // If we've reached/exceeded the mass limit, do not include more ptcls
        if (macc < small_number) break;
      }
      //-------------------------------------------------------------------------------------------


      // Calculate internal sink timestep here
      asqd = DotProduct(sink[s].star->a,sink[s].star->a,ndim);
      sink[s].star->dt_internal = (FLOAT) 0.4*sqrt(sink[s].radius/(sqrt(asqd) + small_number));


    }
    //---------------------------------------------------------------------------------------------


    // Free local thread memory
    delete[] rsqdlist;
    delete[] ilist;
    delete[] neiblist;

  }
  //===============================================================================================


#if defined MPI_PARALLEL
  mpicontrol->UpdateMpiGhostParents(ghosts_accreted, hydro);
  mpicontrol->UpdateSinksAfterAccretion(this);
#endif

  // Quick sanity-check for accreted particles
  for (int i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    if (!(part.flags.is_dead() || part.m > 0.0)) {
      cout << "Accretion problem? : " << i << "   " << part.flags.get() << "   " << part.m
           << "   " << part.h << "   " << part.sinkid << "   " << hydro->mmean << endl;
      ExceptionHandler::getIstance().raise("Error : sink accreting dead or zero-mass particles");
    }
    assert(part.flags.is_dead() || part.m > 0.0);
  }

  return;
}



// Create template class instances of the main SphSimulation object for
// each dimension used (1, 2 and 3)
template class Sinks<1>;
template class Sinks<2>;
template class Sinks<3>;
