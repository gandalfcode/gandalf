//=================================================================================================
//  GradhSphTree.cpp
//  Contains all functions for building, stocking and walking for the
//  binary KD tree for SPH particles.
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
#include <vector>
#include "Precision.h"
#include "Exception.h"
#include "SphNeighbourSearch.h"
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
//  GradhSphTree::GradhSphTree
/// GradhSphTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
GradhSphTree<ndim,ParticleType>::GradhSphTree
 (string tree_type,
  int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, string _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing,
  ParticleTypeRegister& types):
 SphTree<ndim,ParticleType>
  (tree_type, _Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing, types)
{
}



//=================================================================================================
//  GradhSphTree::~GradhSphTreeNeighbourManagerHydro
/// GradhSphTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
GradhSphTree<ndim,ParticleType>::~GradhSphTree()
{
}



//=================================================================================================
//  GradhSphTree::UpdateAllSphProperties
/// Update all gather SPH properties (e.g. rho, div_v) for all active particles in domain.
/// Loops over all cells containing active particles, performs a tree walk for all particles in
/// the cell, and then calls SPH class routine to compute properties from neighbours.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphTree<ndim,ParticleType>::UpdateAllSphProperties
 (Sph<ndim> *sph,                          ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                      ///< [in] Pointer to N-body object
{
  int cactive;                             // No. of active tree cells
  vector<TreeCellBase<ndim> > celllist;		   // List of active tree cells
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*> (sph->GetSphParticleArray());
#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[GradhSphTree::UpdateAllSphProperties]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_PROPERTIES");

  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);
  assert(cactive <= tree->gtot);

  // If there are no active cells, return to main loop
  if (cactive == 0) return;


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,nbody,sph,sphdata)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int celldone;                              // Flag if cell is done
    int cc;                                    // Aux. cell counter
    int i;                                     // Particle id
    int j;                                     // Aux. particle counter
    int jj;                                    // Aux. particle counter
    int k;                                     // Dimension counter
    int Nactive;                               // No. of active particles in cell
    int Ngather;                               // No. of gather neighbours
    int Nneib;                                 // No. of neighbours from tree-walk
    int okflag;                                // Flag if particle is done
    FLOAT draux[ndim];                         // Aux. relative position vector var
    FLOAT drsqdaux;                            // Distance squared
    FLOAT hrangesqd;                           // Kernel extent
    FLOAT hmax;                                // Maximum smoothing length
    FLOAT rp[ndim];                            // Local copy of particle position
    FLOAT *mu2 = 0;                            // Trimmed array (dummy for grad-h)
    int Nneibmax = Nneibmaxbuf[ithread];       // Local copy of neighbour buffer size
    int* activelist = activelistbuf[ithread];  // Local array of active particle ids
    int* neiblist = new int[Nneibmax];         // Local array of neighbour particle ids
    int* ptype    = new int[Nneibmax];         // Local array of particle types
    vector<typename Sph<ndim>::DensityParticle> ngb;           // Local array of neighbour data
    vector<typename Sph<ndim>::DensityParticle> ngb2;          // Local array of reduced neighbour data
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // Local array of active particles

    ngb.resize(Nneibmax) ; ngb2.reserve(Nneibmax);


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim>& cell = celllist[cc];

      celldone = 1;
      hmax = cell.hmax;


      // If hmax is too small so the neighbour lists are invalid, make hmax
      // larger and then recompute for the current active cell.
      //-------------------------------------------------------------------------------------------
      do {
        hmax = (FLOAT) 1.05*hmax;
        celldone = 1;

        // Find list of active particles in current cell
        Nactive = tree->ComputeActiveParticleList(cell, sphdata, activelist);
        for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

        // Compute neighbour list for cell from particles on all trees
        Nneib = 0;
        Nneib = tree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
        Nneib = ghosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
#ifdef MPI_PARALLEL
        Nneib = mpighosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
#endif

        // If there are too many neighbours so the buffers are filled,
        // reallocate the arrays and recompute the neighbour lists.
        while (Nneib == -1) {
          delete[] neiblist;
          Nneibmax = 2*Nneibmax;
          neiblist = new int[Nneibmax];
          ngb.resize(Nneibmax);

          Nneib = 0;
          Nneib = tree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
          Nneib = ghosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,Nneibmax,Nneib,neiblist);
#ifdef MPI_PARALLEL
          Nneib = mpighosttree->ComputeGatherNeighbourList(cell,sphdata,hmax,
                                                           Nneibmax,Nneib,neiblist);
#endif
        };

        // Make local copies of important neib information (mass and position)
        for (jj=0; jj<Nneib; jj++) {
          j         = neiblist[jj];
          ngb[jj] = typename Sph<ndim>::DensityParticle(sphdata[j]) ;
        }

        // Loop over all active particles in the cell
        //-----------------------------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
          i = activelist[j];
          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernrangesqd*hmax*hmax;
          ngb2.resize(0);

          // Compute distance (squared) to all
          //---------------------------------------------------------------------------------------
          for (jj=0; jj<Nneib; jj++) {

            // Only include particles of appropriate types in density calculation
            if (!sph->types[activepart[j].ptype].hmask[ngb[jj].ptype]) continue;


            for (k=0; k<ndim; k++) draux[k] = ngb[jj].r[k] - rp[k];
            drsqdaux = DotProduct(draux,draux,ndim) + small_number;

            // Record distance squared and masses for all potential gather neighbours
            if (drsqdaux <= hrangesqd) {
              ngb2.push_back(ngb[jj]) ;
            }

          }
          //---------------------------------------------------------------------------------------

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) this->CheckValidNeighbourList(i, Ntot, Nneib, neiblist, sphdata, "gather");
#endif

          // Compute smoothing length and other gather properties for ptcl i
          okflag = sph->ComputeH(activepart[j], hmax, ngb2, nbody);

          // If h-computation is invalid, then break from loop and recompute
          // larger neighbour lists
          if (okflag == 0) {
            celldone = 0;
            break;
          }

        }
        //-----------------------------------------------------------------------------------------

      } while (celldone == 0);
      //-------------------------------------------------------------------------------------------

      // Once cell is finished, copy all active particles back to main memory
      for (j=0; j<Nactive; j++) sphdata[activelist[j]] = activepart[j];


    }
    //=============================================================================================

    // Free-up all memory
    delete[] neiblist;

  }
  //===============================================================================================

  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->RunningTime() - twork;
  int Nactivetot=0;
  tree->AddWorkCost(celllist, twork, Nactivetot);
#ifdef OUTPUT_ALL
  cout << "Time computing smoothing lengths : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif


  // Update tree smoothing length values here
  tree->UpdateAllHmaxValues(sphdata);

  return;
}



//=================================================================================================
//  GradhSphTree::UpdateAllSphHydroForces
/// Compute hydro forces for all active SPH particles.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphTree<ndim,ParticleType>::UpdateAllSphHydroForces
 (Sph<ndim> *sph,                          ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                      ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox)                 ///< [in] Simulation domain box
{
  int cactive;                             // No. of active cells
  vector<TreeCellBase<ndim> > celllist;    // List of active tree cells
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*> (sph->GetSphParticleArray());
#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[GradhSphTree::UpdateAllSphHydroForces]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_HYDRO_FORCES");

  // Make sure we have enough neibmanagers
  for (int t = neibmanagerbufhydro.size(); t < Nthreads; ++t) {
    neibmanagerbufhydro.push_back(NeighbourManagerHydro(sph, simbox));
  }

  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) return;


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,nbody,simbox,sph,sphdata)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int* activelist   = activelistbuf[ithread];    // ..
    int* levelneib    = levelneibbuf[ithread];     // ..
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..
    NeighbourManager<ndim,HydroParticle>& neibmanager = neibmanagerbufhydro[ithread];

    for (int i=0; i<sph->Ntot; i++) levelneib[i] = 0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (int cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim>& cell = celllist[cc];

      // Find list of active particles in current cell
      const int Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

      // Make local copies of active particles
      for (int j=0; j<Nactive; j++) {
        activepart[j] = sphdata[activelist[j]];
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].dalphadt  = (FLOAT) 0.0;
        activepart[j].gpot      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        for (int k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell from real and periodic ghost particles


      neibmanager.clear();
      tree->ComputeNeighbourAndGhostList(cell, neibmanager);
      neibmanager.EndSearch(cell,sphdata);

      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (int j=0; j<Nactive; j++) {
        bool do_hydro = sph->types[activepart[j].ptype].hydro_forces;

        if (do_hydro) {
          Typemask hydromask  = sph->types[activepart[j].ptype].hydromask;
          const bool do_pair_once = false;
          NeighbourList<HydroParticle> neiblist =
              neibmanager.GetParticleNeib(activepart[j],hydromask,do_pair_once);

#if defined(VERIFY_ALL)
          neibmanager.VerifyNeighbourList(activelist[j], sph->Nhydro, sphdata, "all");
          neibmanager.VerifyReducedNeighbourList(activelist[j], neiblist, sph->Nhydro,
                                                 sphdata, hydromask, "all");
#endif

          // Compute all neighbour contributions to hydro forces
          typename ParticleType<ndim>::HydroMethod* method = (typename ParticleType<ndim>::HydroMethod*) sph;
          method->ComputeSphHydroForces(activepart[j],neiblist);
        }
      }
      //-------------------------------------------------------------------------------------------

      // Update levelneib for neighbours
      const int Nneib_cell = neibmanager.GetNumAllNeib();
      for (int jj=0; jj<Nneib_cell; jj++) {
        std::pair<int,HydroParticle*> neighbour=neibmanager.GetNeibI(jj);
        const int i=neighbour.first;
        HydroParticle& neibpart=*(neighbour.second);
        levelneib[i]=max(levelneib[i],neibpart.levelneib);
      }


      // Compute all star forces for active particles
      if (nbody->Nnbody > 0) {
        for (int j=0; j<Nactive; j++) {
          if (activelist[j] < sph->Nhydro) {
            sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,activepart[j]);
          }
        }
      }


      // Add all active particles contributions to main array
      for (int j=0; j<Nactive; j++) {
        const int i = activelist[j];
        for (int k=0; k<ndim; k++) sphdata[i].a[k]     += activepart[j].a[k];
        for (int k=0; k<ndim; k++) sphdata[i].a[k]     += activepart[j].atree[k];
        for (int k=0; k<ndim; k++) sphdata[i].atree[k] += activepart[j].atree[k];
        sphdata[i].gpot     += activepart[j].gpot;
        sphdata[i].dudt     += activepart[j].dudt;
        sphdata[i].dalphadt += activepart[j].dalphadt;
        sphdata[i].div_v    += activepart[j].div_v;
        levelneib[i]        = max(levelneib[i], activepart[j].levelneib);
      }

    }
    //=============================================================================================


    // Propagate the changes in levelneib to the main array
#pragma omp for
    for (int i=0; i<sph->Ntot; i++) {
      for (int ithread=0; ithread<Nthreads; ithread++)
        sphdata[i].levelneib = max(sphdata[i].levelneib, levelneibbuf[ithread][i]);
    }


  }
  //===============================================================================================


  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->RunningTime() - twork;
  int Nactivetot=0;
  tree->AddWorkCost(celllist, twork, Nactivetot) ;
#ifdef OUTPUT_ALL
  cout << "Time computing forces : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif

  return;
}



//=================================================================================================
//  GradhSphTree::UpdateAllSphForces
/// Compute all forces on active SPH particles (hydro + gravity) for periodic boundary conditions.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphTree<ndim,ParticleType>::UpdateAllSphForces
 (Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int cactive;                         // No. of active cells
  vector<TreeCellBase<ndim> > celllist;            // List of active cells
  ParticleType<ndim>* sphdata = reinterpret_cast<ParticleType<ndim>*> (sph->GetSphParticleArray());
#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[GradhSphTree::UpdateAllSphForces]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_ALL_FORCES");

  // Make sure we have enough neibmanagers
  for (int t = neibmanagerbufhydro.size(); t < Nthreads; ++t) {
    neibmanagerbufhydro.push_back(NeighbourManagerHydro(sph, simbox));
  }

  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) return;


  // Set-up all OMP threads
  //===============================================================================================
  #pragma omp parallel default(none) shared(celllist,cactive,ewald,nbody,simbox,sph,sphdata,cout)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                      // Aux. cell counter
    int Nactive;                                 // ..
    FLOAT aperiodic[ndim];                       // ..
    FLOAT potperiodic;                           // ..
    int *activelist  = activelistbuf[ithread];   // ..
    int *levelneib   = levelneibbuf[ithread];    // ..
    ParticleType<ndim>* activepart  = activepartbuf[ithread];   // ..
    Typemask gravmask = sph->types.gravmask;
    NeighbourManager<ndim,HydroParticle>& neibmanager = neibmanagerbufhydro[ithread];

    // Zero timestep level array
    for (int i=0; i<sph->Ntot; i++) levelneib[i] = 0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim> &cell = celllist[cc];

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, sphdata, activelist);

      // Make local copies of active particles
      for (int j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Zero/initialise all summation variables for active particles
      for (int j=0; j<Nactive; j++) {
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        activepart[j].gpot      = (activepart[j].m/activepart[j].h)*sph->kernp->wpot(0.0);
        for (int k=0; k<ndim; k++) activepart[j].a[k]     = (FLOAT) 0.0;
        for (int k=0; k<ndim; k++) activepart[j].atree[k] = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      neibmanager.clear();
      tree->ComputeGravityInteractionAndGhostList(cell, neibmanager);
      neibmanager.EndSearchGravity(cell,sphdata);

      MultipoleMoment<ndim>* gravcell;
      const int Ngravcell = neibmanager.GetGravCell(&gravcell);

      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (int j=0; j<Nactive; j++) {
        const bool do_grav  = sph->types[activepart[j].ptype].self_gravity ;
        Typemask hydromask = sph->types[activepart[j].ptype].hydromask ;
        GravityNeighbourLists<HydroParticle> neiblists =
            neibmanager.GetParticleNeibGravity(activepart[j],hydromask);

        // Compute forces between SPH neighbours (hydro and gravity)
        typename ParticleType<ndim>::HydroMethod* method = (typename ParticleType<ndim>::HydroMethod*) sph;

        if (neiblists.neiblist.size() > 0) {
          method->ComputeSphHydroGravForces(activepart[j], neiblists.neiblist);
        }

        if (do_grav) {

          // Compute soften grav forces between non-SPH neighbours (hydro and gravity)
          method->ComputeSphGravForces(activepart[j], neiblists.smooth_gravlist);

          // Compute direct gravity forces between distant particles
          method->ComputeDirectGravForces(activepart[j], neiblists.directlist);

          // Compute gravitational force due to distant cells
          if (multipole == "monopole") {
            ComputeCellMonopoleForces(activepart[j].gpot, activepart[j].atree,
                                      activepart[j].r, Ngravcell, gravcell);
          }
          else if (multipole == "quadrupole") {
            ComputeCellQuadrupoleForces(activepart[j].gpot, activepart[j].atree,
                                        activepart[j].r, Ngravcell, gravcell);
         }

          // Add the periodic correction force for SPH and direct-sum neighbours
          if (simbox.PeriodicGravity) {
            int Ntotneib = neibmanager.GetNumAllNeib();

            for (int jj=0; jj< Ntotneib; jj++) {
              if (!gravmask[neibmanager[jj].ptype]) continue;
              FLOAT draux[ndim];
              for (int k=0; k<ndim; k++) draux[k] = neibmanager[jj].r[k] - activepart[j].r[k];
              ewald->CalculatePeriodicCorrection(neibmanager[jj].m, draux, aperiodic, potperiodic);
              for (int k=0; k<ndim; k++) activepart[j].atree[k] += aperiodic[k];
              activepart[j].gpot += potperiodic;
            }


            // Now add the periodic correction force for all cell COMs
            for (int jj=0; jj<Ngravcell; jj++) {
              FLOAT draux[ndim];
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
      if (multipole == "fast_monopole") {
        ComputeFastMonopoleForces(Nactive, Ngravcell, gravcell, cell, activepart, sph->types);
      }
      else if (multipole == "fast_quadrupole") {
        ComputeFastQuadrupoleForces(Nactive, Ngravcell, gravcell, cell, activepart, sph->types);
      }

      // Set gpot_hydro (RadWS) before sink contribution
      for (int j=0; j<Nactive; j++) {
        activepart[j].gpot_hydro = activepart[j].gpot;
      }

      // Compute all star forces for active particles
      if (nbody->Nnbody > 0) {
        for (int j=0; j<Nactive; j++) {
          if (activelist[j] < sph->Nhydro) {
            sph->ComputeStarGravForces(nbody->Nnbody, nbody->nbodydata, activepart[j]);
          }
        }
      }

      // Add all active particles contributions to main array
      for (int j=0; j<Nactive; j++) {
        int i = activelist[j];
        for (int k=0; k<ndim; k++) sphdata[i].a[k]     += activepart[j].a[k];
        for (int k=0; k<ndim; k++) sphdata[i].a[k]     += activepart[j].atree[k];
        for (int k=0; k<ndim; k++) sphdata[i].atree[k] += activepart[j].atree[k];
        sphdata[i].gpot  += activepart[j].gpot;
        sphdata[i].gpot_hydro += activepart[j].gpot_hydro;
        sphdata[i].dudt  += activepart[j].dudt;
        sphdata[i].div_v += activepart[j].div_v;
        levelneib[i]      = max(levelneib[i],activepart[j].levelneib);
      }

      // Update levelneib for neighbours
      const int Nneib_cell = neibmanager.GetNumAllNeib();
      for (int jj=0; jj<Nneib_cell; jj++) {
        std::pair<int,HydroParticle*> neighbour=neibmanager.GetNeibI(jj);
        const int i=neighbour.first;
        HydroParticle& neibpart=*(neighbour.second);
        levelneib[i]=max(levelneib[i],neibpart.levelneib);
      }

    }
    //=============================================================================================


    // Propagate the changes in levelneib to the main array
#pragma omp for
    for (int i=0; i<sph->Ntot; i++) {
      for (int ithread=0; ithread<Nthreads; ithread++)
        sphdata[i].levelneib = max(sphdata[i].levelneib, levelneibbuf[ithread][i]);
    }

  }
  //===============================================================================================

  // Compute time spent in routine and in each cell for load balancing
#ifdef MPI_PARALLEL
  twork = timing->RunningTime() - twork;
  int Nactivetot=0;
  tree->AddWorkCost(celllist, twork, Nactivetot) ;
#ifdef OUTPUT_ALL
  cout << "Time computing forces : " << twork << "     Nactivetot : " << Nactivetot << endl;
#endif
#endif


  return;
}



template class GradhSphTree<1,GradhSphParticle>;
template class GradhSphTree<2,GradhSphParticle>;
template class GradhSphTree<3,GradhSphParticle>;
