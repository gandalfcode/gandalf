//=================================================================================================
//  MeshlessFVTree.cpp
//  Contains all functions for building, stocking and walking for the
//  binary KD tree for Meshless-FV particles.
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


#include <cstdlib>
#include <cassert>
#include <iostream>
#include <string>
#include <math.h>
#include <utility>
#include "Precision.h"
#include "Exception.h"
#include "MfvNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Debug.h"
#include "NeighbourManager.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=================================================================================================
//  MeshlessFVTree::MeshlessFVTree
/// MeshlessFVTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
MeshlessFVTree<ndim,ParticleType>::MeshlessFVTree
 (string tree_type,
  int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, multipole_method _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing,
  ParticleTypeRegister& types):
 HydroTree<ndim,ParticleType>
  (tree_type, _Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing, types)
{

}



//=================================================================================================
//  MeshlessFVTree::~MeshlessFVTree
/// MeshlessFVTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
MeshlessFVTree<ndim,ParticleType>::~MeshlessFVTree()
{
}



//=================================================================================================
//  MeshlessFVTree::UpdateAllSphProperties
/// Update all gather SPH properties (e.g. rho, div_v) for all active particles in domain.
/// Loops over all cells containing active particles, performs a tree walk for all particles in
/// the cell, and then calls SPH class routine to compute properties from neighbours.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MeshlessFVTree<ndim,ParticleType>::UpdateAllProperties
 (MeshlessFV<ndim> *mfv,                   ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                      ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox)                 ///< [in] Simuation box
{
  int cactive;                             // No. of active tree cells
  vector<TreeCellBase<ndim> > celllist;    // List of active cells
  MeshlessFVParticle<ndim> *mfvdata = mfv->GetMeshlessFVParticleArray();
#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[MeshlessFVTree::UpdateAllProperties]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFV_PROPERTIES");

  // Make sure we have enough neibmanagers
  for (int t = neibmanagerbufdens.size(); t < Nthreads; ++t) {
    neibmanagerbufdens.push_back(NeighbourManagerDensity(mfv, simbox));
  }

  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) return;


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,nbody,mfv,mfvdata)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int celldone;                              // Flag if cell is done
    int cc;                                    // Aux. cell counter
    int j;                                     // Aux. particle counter
    int Nactive;                               // No. of active particles in cell
    int okflag;                                // Flag if particle is done
    FLOAT hrangesqd;                           // Kernel extent
    FLOAT hmax;                                // Maximum smoothing length
    int* activelist = activelistbuf[ithread];   // Local array of active particle-ids
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // Local array of active particles
    NeighbourManager<ndim,DensityParticle>& neibmanager = neibmanagerbufdens[ithread];



    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim> cell = celllist[cc];
      celldone = 1;
      hmax = cell.hmax;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,mfvdata,activelist);

      // Skip particles that have an up-to-date density estimate
      for (j=0; j<Nactive; j++) {
        if (mfvdata[activelist[j]].flags.check(update_density)) {
          activepart[j] = mfvdata[activelist[j]];
        }
        else {
          activelist[j] = activelist[--Nactive];
          j-- ;
        }
      }

      if (Nactive == 0) continue;

      // If hmax is too small so the neighbour lists are invalid, make hmax
      // larger and then recompute for the current active cell.
      //-------------------------------------------------------------------------------------------
      do {
        hmax = (FLOAT) 1.05*hmax;
        cell.hmax = hmax;
        celldone = 1;


        // Compute neighbour list for cell from particles on all trees
        neibmanager.set_target_cell(cell) ;
        tree->ComputeGatherNeighbourList(cell,mfvdata,hmax,neibmanager);
        ghosttree->ComputeGatherNeighbourList(cell,mfvdata,hmax,neibmanager);
#ifdef MPI_PARALLEL
        mpighosttree->ComputeGatherNeighbourList(cell,mfvdata,hmax,neibmanager);
#endif
        neibmanager.EndSearchGather(cell, mfvdata);



        // Loop over all active particles in the cell
        //-----------------------------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernrangesqd*hmax*hmax;

          Typemask densmask = mfv->types[activepart[j].ptype].hmask;

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernrangesqd*hmax*hmax;

          NeighbourList<DensityParticle> neiblist =
              neibmanager.GetParticleNeibGather(activepart[j],densmask,hrangesqd);

          // Compute smoothing length and other gather properties for ptcl i
          okflag = mfv->ComputeH(activepart[j], hmax, neiblist, nbody);

          // If h-computation is invalid, then break from loop and recompute larger neighbour lists
          if (okflag == 0) {
            celldone = 0;
            break;
          }

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          neibmanager.VerifyNeighbourList(activelist[j], mfv->Ntot, mfvdata, "gather");
          neibmanager.VerifyReducedNeighbourList(activelist[j], neiblist, mfv->Ntot,
                                                 mfvdata, densmask, "gather");
#endif
        }
        //-----------------------------------------------------------------------------------------

      } while (celldone == 0);
      //-------------------------------------------------------------------------------------------

      // Once cell is finished, copy all active particles back to main memory and record that we
      // have done a density update
      for (j=0; j<Nactive; j++) {
        activepart[j].flags.unset(update_density) ;
        mfvdata[activelist[j]] = activepart[j];
      }

      tree->UpdateHmaxLeaf(cell, mfvdata) ;

    }
    //=============================================================================================


  }
  //===============================================================================================

  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->RunningTime() - twork;
  int Nactivetot=0;
  tree->AddWorkCost(celllist, twork, Nactivetot) ;
#ifdef OUTPUT_ALL
  cout << "Time computing smoothing lengths : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif

  // Update tree smoothing length values here
  timer.EndTiming();
  CodeTiming::BlockTimer timer2 = timing->StartNewTimer("UPDATE_HMAX");

  tree->UpdateAllHmaxValues(mfvdata, false);

  return;
}



//=================================================================================================
//  MeshlessFVTree::UpdateGradientMatrices
/// Compute hydro forces for all active SPH particles.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MeshlessFVTree<ndim,ParticleType>::UpdateGradientMatrices
 (MeshlessFV<ndim> *mfv,                   ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                      ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox)                 ///< [in] Simulation domain box
{
  int cactive;                             // No. of active cells
  vector<TreeCellBase<ndim> > celllist;            // List of active cells
#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[MeshlessFVTree::UpdateGradientMatrices]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFV_UPDATE_GRADIENTS");

  // Make sure we have enough neibmanagers
  for (int t = neibmanagerbufgradient.size(); t < Nthreads; ++t)
    neibmanagerbufgradient.push_back(NeighbourManagerGradient(mfv, simbox));

  int Ntot = mfv->Ntot;
  MeshlessFVParticle<ndim> *mfvdata = mfv->GetMeshlessFVParticleArray();


  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    return;
  }

  // Update ghost tree smoothing length values here
#ifdef MPI_PARALLEL
  if (mfv->Nmpighost > 0) mpighosttree->UpdateAllHmaxValues(mfvdata);
#endif

  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,nbody,mfv,mfvdata,Ntot,simbox)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                        // Aux. cell counter
    int i;                                         // Particle id
    int j;                                         // Aux. particle counter
    int k;                                         // Dimension counter
    int Nactive;                                   // ..
    int* activelist = activelistbuf[ithread];      // ..
    int* levelneib  = levelneibbuf[ithread];       // ..
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..
    NeighbourManagerGradient& neibmanager = neibmanagerbufgradient[ithread];
    for (i=0; i<Ntot; i++) levelneib[i] = 0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim>& cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, mfvdata, activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = mfvdata[activelist[j]];

      // Compute neighbour list for cell from real and periodic ghost particles
      neibmanager.set_target_cell(cell);
      tree->ComputeNeighbourAndGhostList(cell, neibmanager);
#ifdef MPI_PARALLEL
      mpighosttree->ComputeNeighbourList(cell,neibmanager);
#endif
      neibmanager.EndSearch(cell,mfvdata);

      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        // If particle is NOT a hydro particle (and therefore doesn't need gradients), skip to next
        if (!mfv->types[activepart[j].ptype].hydro_forces) continue;

        // Make local copy of hmask for active particle
        Typemask hmask = mfv->types[activepart[j].ptype].hmask;

        activepart[j].levelneib = 0;

        const bool do_pair_once=false;
        NeighbourList<GradientParticle> neiblist =
            neibmanager.GetParticleNeib(activepart[j],hmask,do_pair_once);

#if defined(VERIFY_ALL)
        neibmanager.VerifyNeighbourList(i, mfv->Nhydro, mfvdata, "all");
        neibmanager.VerifyReducedNeighbourList(i, neiblist, mfv->Nhydro, mfvdata, hmask, "all");
#endif
        // Compute all neighbour contributions to gradients
        mfv->ComputeGradients(activepart[j], neiblist);

      }
      //-------------------------------------------------------------------------------------------

      // Update levelneib for neighbours
      const int Nneib_cell = neibmanager.GetNumAllNeib();
      for (int jj=0; jj<Nneib_cell; jj++) {
    	std::pair<int,GradientParticle*> neighbour=neibmanager.GetNeibI(jj);
    	const int i=neighbour.first;
    	GradientParticle& neibpart=*(neighbour.second);
    	levelneib[i]=max(levelneib[i],neibpart.levelneib);
       }


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) {
          for (int kk=0; kk<ndim; kk++) mfvdata[i].B[k][kk] = activepart[j].B[k][kk];
        }
        for (int var=0; var<ndim+2; var++) {
          for (k=0; k<ndim; k++) mfvdata[i].grad[var][k] = activepart[j].grad[var][k];
          mfvdata[i].alpha_slope[var] = activepart[j].alpha_slope[var];
        }
        mfvdata[i].vsig_max = activepart[j].vsig_max;
        mfvdata[i].levelneib = activepart[j].levelneib;
        mfvdata[i].flags = activepart[j].flags;
      }


    }
    //=============================================================================================

    // Copy back the neighbour levels
#ifdef _OPENMP
    int Nthreads = omp_get_num_threads() ;
#else
    int Nthreads = 1 ;
#endif
#pragma omp for schedule(static)
      for(i=0; i<Ntot; ++i) {
        for (k=0; k<Nthreads; k++)
          mfvdata[i].levelneib = max(mfvdata[i].levelneib, levelneibbuf[k][i]);
      }

  }
  //===============================================================================================


  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->RunningTime() - twork;
  int Nactivetot=0;
  tree->AddWorkCost(celllist, twork, Nactivetot) ;
#ifdef OUTPUT_ALL
  cout << "Time computing gradients : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif

  return;
}



//=================================================================================================
//  MeshlessFVTree::UpdateGodunovFluxes
/// Compute hydro forces for all active SPH particles.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MeshlessFVTree<ndim,ParticleType>::UpdateGodunovFluxes
 (FLOAT timestep,                          ///< [in] Lowest timestep value
  MeshlessFV<ndim> *mfv,                   ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                      ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox)                 ///< [in] Simulation domain box
{
  int cactive;                             // No. of active cells
  vector<TreeCellBase<ndim> > celllist;            // List of active cells
#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[MeshlessFVTree::UpdateGodunovFluxes]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFV_UPDATE_FLUXES");

  // Make sure we have enough neibmanagers
  for (int t = neibmanagerbufflux.size(); t < Nthreads; ++t)
    neibmanagerbufflux.push_back(NeighbourManagerFlux(mfv, simbox));

  int Ntot = mfv->Ntot;
  MeshlessFVParticle<ndim> *mfvdata = mfv->GetMeshlessFVParticleArray();

  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    return;
  }

   typedef FLOAT (*fluxArray)[ndim+2];
   typedef FLOAT (*rdmdtArray)[ndim];

    fluxArray* dQBufferGlob = new fluxArray[Nthreads];
    fluxArray* fluxBufferGlob = new fluxArray[Nthreads];
    rdmdtArray* rdmdtBufferGlob = new rdmdtArray[Nthreads];

  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,mfv,mfvdata,Ntot,cout,dQBufferGlob,fluxBufferGlob,rdmdtBufferGlob,timestep)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int i;                                         // Particle id
    int k;                                         // Dimension counter
    int Nactive;                                   // ..
    int* activelist   = activelistbuf[ithread];    // ..
    fluxArray dQBuffer      = new FLOAT[Ntot][ndim+2];  // ..
    dQBufferGlob[ithread] = dQBuffer;
    fluxArray fluxBuffer   = new FLOAT[Ntot][ndim+2];  // ..
    fluxBufferGlob[ithread] = fluxBuffer;
    rdmdtArray rdmdtBuffer    = new FLOAT[Ntot][ndim];    // ..
    rdmdtBufferGlob[ithread] = rdmdtBuffer;
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..
    NeighbourManager<ndim,FluxParticle>& neibmanager = neibmanagerbufflux[ithread];

    for (int i=0; i<Ntot; i++) {
      for (int k=0; k<ndim+2; k++) fluxBuffer[i][k] = (FLOAT) 0.0;
      for (int k=0; k<ndim+2; k++) dQBuffer[i][k] = (FLOAT) 0.0;
      for (int k=0; k<ndim; k++) rdmdtBuffer[i][k] = (FLOAT) 0.0;
    }


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (int cc=0; cc<cactive; cc++) {
     TreeCellBase<ndim>& cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,mfvdata,activelist);

      // Make local copies of active particles
      for (int j=0; j<Nactive; j++) {
        activepart[j] = mfvdata[activelist[j]];
        for (int k=0; k<ndim+2; k++) activepart[j].dQ[k]   = (FLOAT) 0.0;
        for (int k=0; k<ndim+2; k++) activepart[j].dQdt[k] = (FLOAT) 0.0;
        for (int k=0; k<ndim; k++) activepart[j].rdmdt[k]  = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell from real and periodic ghost particles
      neibmanager.set_target_cell(cell);
      tree->ComputeNeighbourAndGhostList(cell, neibmanager);
      neibmanager.EndSearch(cell,mfvdata);


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (int j=0; j<Nactive; j++) {
        i = activelist[j];

        // If particle is not a hydro particle (e.g. cdm), then skip to next active particle
        if (mfv->types[activepart[j].ptype].hydro_forces == false) continue;

        // Make a local copy of the hydro neighbour mask
        Typemask hydromask = mfv->types[activepart[j].ptype].hydromask;

        bool do_pair_once=true;
        NeighbourList<FluxParticle> neiblist =
            neibmanager.GetParticleNeib(activepart[j],hydromask,do_pair_once);

#if defined(VERIFY_ALL)
        neibmanager.VerifyNeighbourList(i, mfv->Nhydro, mfvdata, "all");
#endif

        // Compute all neighbour contributions to hydro fluxes
        mfv->ComputeGodunovFlux(activepart[j], neiblist, timestep);

      }
      //-------------------------------------------------------------------------------------------

      const int Nneib_cell = neibmanager.GetNumAllNeib();
      // Accumulate fluxes for neighbours
      for (int jj=0; jj<Nneib_cell; jj++) {
    	std::pair<int,FluxParticle*> neighbour=neibmanager.GetNeibI(jj);
    	const int i=neighbour.first;
    	FluxParticle& neibpart=*(neighbour.second);
        if (!neibpart.flags.is_mirror()) {
	        if (neibpart.flags.check(active))
	          for (k=0; k<ndim+2; k++) fluxBuffer[i][k] += neibpart.dQdt[k];
          for (k=0; k<ndim+2; k++) dQBuffer[i][k] += neibpart.dQ[k];
          for (k=0; k<ndim; k++) rdmdtBuffer[i][k] += neibpart.rdmdt[k];
        }
      }
      // Add all active particles contributions to main array
      for (int j=0; j<Nactive; j++) {
        i = activelist[j];
        for (int k=0; k<ndim; k++) rdmdtBuffer[i][k] += activepart[j].rdmdt[k];
        for (int k=0; k<ndim+2; k++) fluxBuffer[i][k] += activepart[j].dQdt[k];
        for (int k=0; k<ndim+2; k++) dQBuffer[i][k] += activepart[j].dQ[k];
      }

    }
    //=============================================================================================


    // Add all buffers back to main arrays
#pragma omp for schedule(static)
      for (i=0; i<Ntot; i++) {
	for (int ithread=0; ithread<Nthreads; ithread++) {
        if (mfvdata[i].flags.check(active)) {
		  for (int k=0; k<ndim; k++) mfvdata[i].rdmdt[k] += rdmdtBufferGlob[ithread][i][k];
		  for (int k=0; k<ndim+2; k++) mfvdata[i].dQdt[k] += fluxBufferGlob[ithread][i][k];
		}
		for (int k=0; k<ndim+2; k++) mfvdata[i].dQ[k] += dQBufferGlob[ithread][i][k];
	}
      }

    delete[] rdmdtBuffer;
    delete[] fluxBuffer;
    delete[] dQBuffer;


  }
  //===============================================================================================


  delete[] rdmdtBufferGlob;
  delete[] fluxBufferGlob;
  delete[] dQBufferGlob;

  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->RunningTime() - twork;
  int Nactivetot=0;
  tree->AddWorkCost(celllist, twork, Nactivetot) ;
#ifdef OUTPUT_ALL
  cout << "Time computing fluxes : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif


  return;
}



//=================================================================================================
//  MeshlessFVTree::UpdateAllGravForces
/// Compute all local 'gather' properties of currently active particles, and
/// then compute each particle's contribution to its (active) neighbour
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to
/// construct local neighbour lists for all particles  inside the cell.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MeshlessFVTree<ndim,ParticleType>::UpdateAllGravForces
 (MeshlessFV<ndim> *mfv,               ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int cactive;                         // No. of active cells
  vector<TreeCellBase<ndim> > celllist;            // List of active cells
  //ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[MeshlessFVTree::UpdateAllGravForces]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MFV_GRAV_FORCES");

  // Make sure we have enough neibmanagers
  for (int t = neibmanagerbufgrav.size(); t < Nthreads; ++t)
    neibmanagerbufgrav.push_back(NeighbourManagerGrav(mfv, simbox));

#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  MeshlessFVParticle<ndim> *partdata = mfv->GetMeshlessFVParticleArray();

  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    return;
  }

  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(celllist,cactive,ewald,mfv,nbody,partdata,simbox,cout)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                      // Aux. cell counter
    int Nactive;                                 // ..
    FLOAT aperiodic[ndim];                       // ..
    FLOAT draux[ndim];                           // Aux. relative position vector
    FLOAT potperiodic;                           // ..
    int *activelist  = activelistbuf[ithread];   // ..
    ParticleType<ndim>* activepart  = activepartbuf[ithread];   // ..
    Typemask gravmask = mfv->types.gravmask;
    NeighbourManagerGrav neibmanager = neibmanagerbufgrav[ithread];

    neibmanager.set_multipole_type(multipole) ;


    Typemask hydromask ;
    // This creates a mask which is always false. The purpose is that in this way no particle will be added to the
    // list of hydro neighbours
    for (int k=0; k< Ntypes; ++k){
    	hydromask[k] = false;
    }

    bool self_gravity =  mfv->self_gravity == 1 ;

    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim>& cell = celllist[cc];
      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, partdata, activelist);

      // Make local copies of active particles
      for (int j=0; j<Nactive; j++) activepart[j] = partdata[activelist[j]];

      // Zero/initialise all summation variables for active particles
      for (int j=0; j<Nactive; j++)
        for (int k=0; k<ndim; k++) activepart[j].atree[k] = (FLOAT) 0.0;

      // Do the self-gravity contribution, or just the stars
      if (self_gravity) {

        // Include self-term for potential
        // TODO:
        //   Check: Is this now accounted for in the neighbour list?
        for (int j=0; j<Nactive; j++)
          activepart[j].gpot = (activepart[j].m/activepart[j].h)*mfv->kernp->wpot((FLOAT) 0.0);


        // Compute neighbour list for cell depending on physics options
        neibmanager.set_target_cell(cell);
        tree->ComputeGravityInteractionAndGhostList(cell, neibmanager);
        neibmanager.EndSearchGravity(cell,partdata);

        MultipoleMoment<ndim>* gravcell;
        int Ngravcell = neibmanager.GetGravCell(&gravcell);

        // Loop over all active particles in the cell
        //-------------------------------------------------------------------------------------------
        for (int j=0; j<Nactive; j++) {

          // Only calculate gravity for active particle types that have self-gravity activated
          if (mfv->types[activepart[j].ptype].self_gravity) {

            GravityNeighbourLists<GravParticle> neiblists =
              neibmanager.GetParticleNeibGravity(activepart[j],hydromask);

            // Compute forces with hydro neighbours
            mfv->ComputeSmoothedGravForces(activepart[j], neiblists.neiblist);

            // Compute forces with non-hydro neighbours
            mfv->ComputeSmoothedGravForces(activepart[j], neiblists.smooth_gravlist);

            // Compute direct gravity forces between distant particles
            mfv->ComputeDirectGravForces(activepart[j], neiblists.directlist);

            // Compute gravitational force due to distant cells
            if (multipole == monopole) {
              ComputeCellMonopoleForces(activepart[j].gpot, activepart[j].atree,
                                        activepart[j].r, Ngravcell, gravcell);
            }
            else if (multipole == quadrupole) {
              ComputeCellQuadrupoleForces(activepart[j].gpot, activepart[j].atree,
                                          activepart[j].r, Ngravcell, gravcell);
            }

            // Add the periodic correction force for SPH and direct-sum neighbours
            if (simbox.PeriodicGravity) {
              int Ntotneib = neibmanager.GetNumAllNeib();

              for (int jj=0; jj<Ntotneib; jj++) {
                if (!gravmask[neibmanager[jj].ptype]) continue;
                for (int k=0; k<ndim; k++) draux[k] = neibmanager[jj].r[k] - activepart[j].r[k];
                ewald->CalculatePeriodicCorrection(neibmanager[jj].m, draux, aperiodic, potperiodic);
                for (int k=0; k<ndim; k++) activepart[j].atree[k] += aperiodic[k];
                activepart[j].gpot += potperiodic;
              }

              // Now add the periodic correction force for all cell COMs
              for (int jj=0; jj<Ngravcell; jj++) {
                for (int k=0; k<ndim; k++) draux[k] = gravcell[jj].r[k] - activepart[j].r[k];
                ewald->CalculatePeriodicCorrection(gravcell[jj].m, draux, aperiodic, potperiodic);
                for (int k=0; k<ndim; k++) activepart[j].atree[k] += aperiodic[k];
                activepart[j].gpot += potperiodic;
              }
            }
          }
        }
        //-------------------------------------------------------------------------------------------


        // Compute 'fast' multipole terms here
        if (multipole == fast_monopole || multipole == fast_quadrupole) {
          neibmanager.ComputeFastMultipoleForces(Nactive, activepart, mfv->types) ;
        }
      } // End of self-gravity for this cell

      // Compute all star forces for active particles
      for (int j=0; j<Nactive; j++) {
        if (activelist[j] < mfv->Nhydro) {
          mfv->ComputeStarGravForces(nbody->Nnbody, nbody->nbodydata, activepart[j]);
        }
      }

      // Add all active particles contributions to main array
      for (int j=0; j<Nactive; j++) {
        const int i = activelist[j];
        for (int k=0; k<ndim; k++) partdata[i].a[k]     += activepart[j].atree[k];
        for (int k=0; k<ndim; k++) partdata[i].atree[k] += activepart[j].atree[k];
        partdata[i].gpot  += activepart[j].gpot;
      }

    }
    //=============================================================================================



  }
  //===============================================================================================

#ifdef MPI_PARALLEL
  twork = timing->RunningTime() - twork;
  int Nactivetot=0;
  tree->AddWorkCost(celllist, twork, Nactivetot) ;
#ifdef OUTPUT_ALL
  cout << "Time computing fluxes : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif

  return;
}


template class MeshlessFVTree<1,MeshlessFVParticle>;
template class MeshlessFVTree<2,MeshlessFVParticle>;
template class MeshlessFVTree<3,MeshlessFVParticle>;
