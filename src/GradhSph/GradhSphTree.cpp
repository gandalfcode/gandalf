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
#include "Precision.h"
#include "Exception.h"
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Debug.h"
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
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing, ParticleTypeRegister& types):
 SphTree<ndim,ParticleType>
  (tree_type, _Nleafmax, _Nmpi, _pruning_level_min, _pruning_level_max, _thetamaxsqd,
   _kernrange, _macerror, _gravity_mac, _multipole, _box, _kern, _timing, types)
{ }



//=================================================================================================
//  GradhSphTree::~GradhSphTree
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
 (int Nhydro,                              ///< [in] No. of SPH particles
  int Ntot,                                ///< [in] No. of SPH + ghost particles
  Sph<ndim> *sph,                          ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody)                      ///< [in] Pointer to N-body object
{
  int cactive;                             // No. of active tree cells
  vector<TreeCellBase<ndim> > celllist;		   // List of active tree cells
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[GradhSphTree::UpdateAllSphProperties]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_PROPERTIES");

  ParticleType<ndim>* sphdata =
      reinterpret_cast<ParticleType<ndim>*> (sph->GetSphParticleArray());

  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);
  assert(cactive <= tree->gtot);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    return;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,cout,nbody,sph,sphdata,Ntot)
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
    FLOAT* gpot   = new FLOAT[Nneibmax];       // Local array of particle potentials
    FLOAT* gpot2  = new FLOAT[Nneibmax];       // Local reduced array of neighbour potentials
    FLOAT* drsqd  = new FLOAT[Nneibmax];       // Local array of distances (squared)
    FLOAT* m      = new FLOAT[Nneibmax];       // Local array of particle masses
    FLOAT* m2     = new FLOAT[Nneibmax];       // Local reduced array of neighbour masses
    FLOAT* r      = new FLOAT[Nneibmax*ndim];  // Local array of particle positions
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // Local array of active particles


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
          delete[] r;
          delete[] m2;
          delete[] m;
          delete[] drsqd;
          delete[] gpot2;
          delete[] gpot;
          delete[] ptype;
          delete[] neiblist;
          Nneibmax = 2*Nneibmax;
          neiblist = new int[Nneibmax];
          ptype    = new int[Nneibmax];
          gpot     = new FLOAT[Nneibmax];
          gpot2    = new FLOAT[Nneibmax];
          drsqd    = new FLOAT[Nneibmax];
          m        = new FLOAT[Nneibmax];
          m2       = new FLOAT[Nneibmax];
          r        = new FLOAT[Nneibmax*ndim];
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
          gpot[jj]  = sphdata[j].gpot;
          m[jj]     = sphdata[j].m;
          ptype[jj] = sphdata[j].ptype;
          for (k=0; k<ndim; k++) r[ndim*jj + k] = sphdata[j].r[k];
        }

        // Loop over all active particles in the cell
        //-----------------------------------------------------------------------------------------
        for (j=0; j<Nactive; j++) {
          i = activelist[j];
          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];

          // Set gather range as current h multiplied by some tolerance factor
          hrangesqd = kernrangesqd*hmax*hmax;
          Ngather = 0;

          // Compute distance (squared) to all
          //---------------------------------------------------------------------------------------
          for (jj=0; jj<Nneib; jj++) {

            // Only include particles of appropriate types in density calculation
            if (!sph->types[activepart[j].ptype].hmask[ptype[jj]]) continue ;

            for (k=0; k<ndim; k++) draux[k] = r[ndim*jj + k] - rp[k];
            drsqdaux = DotProduct(draux,draux,ndim) + small_number;

            // Record distance squared and masses for all potential gather neighbours
            if (drsqdaux <= hrangesqd) {
              gpot[Ngather]  = gpot[jj];
              drsqd[Ngather] = drsqdaux;
              m2[Ngather]    = m[jj];
              Ngather++;
            }

          }
          //---------------------------------------------------------------------------------------

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) this->CheckValidNeighbourList(i, Ntot, Nneib, neiblist, sphdata, "gather");
#endif

          // Compute smoothing length and other gather properties for ptcl i
          okflag = sph->ComputeH(i, Ngather, hmax, m2, mu2, drsqd, gpot, activepart[j], nbody);

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
    delete[] r;
    delete[] m2;
    delete[] m;
    delete[] drsqd;
    delete[] gpot2;
    delete[] gpot;
    delete[] ptype;
    delete[] neiblist;

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
  tree->UpdateAllHmaxValues(sphdata);

  return;
}



//=================================================================================================
//  GradhSphTree::UpdateAllSphHydroForces
/// Compute hydro forces for all active SPH particles.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphTree<ndim,ParticleType>::UpdateAllSphHydroForces
 (int Nhydro,                              ///< [in] No. of SPH particles
  int Ntot,                                ///< [in] No. of SPH + ghost particles
  Sph<ndim> *sph,                          ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                      ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox)                 ///< [in] Simulation domain box
{
  int cactive;                             // No. of active cells
  vector<TreeCellBase<ndim> > celllist;                // List of active tree cells
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[GradhSphTree::UpdateAllSphHydroForces]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_HYDRO_FORCES");

  ParticleType<ndim>* sphdata =
      reinterpret_cast<ParticleType<ndim>*> (sph->GetSphParticleArray());


  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    return;
  }

  // Update ghost tree smoothing length values here
  tree->UpdateAllHmaxValues(sphdata);
  //if (ghosttree->Ntot > 0) ghosttree->UpdateHmaxValues(ghosttree->celldata[0],sphdata);


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(cactive,celllist,nbody,simbox,sph,sphdata)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                        // Aux. cell counter
    int i;                                         // Particle id
    int j;                                         // Aux. particle counter
    int jj;                                        // Aux. particle counter
    int k;                                         // Dimension counter
    int Nactive;                                   // ..
    int Nneib;                                     // ..
    int Nhydroaux;                                 // ..
    FLOAT draux[ndim];                             // Aux. relative position vector
    FLOAT drsqd;                                   // Distance squared
    FLOAT hrangesqdi;                              // Kernel gather extent
    FLOAT rp[ndim];                                // Local copy of particle position
    Typemask hydromask;                            // Mask for computing hydro forces
    int Nneibmax      = Nneibmaxbuf[ithread];      // ..
    int* activelist   = activelistbuf[ithread];    // ..
    int* levelneib    = levelneibbuf[ithread];     // ..
    int* neiblist     = new int[Nneibmax];         // ..
    int* sphlist      = new int[Nneibmax];         // ..
    FLOAT* dr         = new FLOAT[Nneibmax*ndim];  // ..
    FLOAT* drmag      = new FLOAT[Nneibmax];       // ..
    FLOAT* invdrmag   = new FLOAT[Nneibmax];       //..
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];     // ..

    for (i=0; i<sph->Ntot; i++) levelneib[i] = 0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim>& cell = celllist[cc];
      //assert (cell.id>=0);


      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) {
        activepart[j] = sphdata[activelist[j]];
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].dalphadt  = (FLOAT) 0.0;
        activepart[j].gpot      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        for (k=0; k<ndim; k++) activepart[j].a[k] = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell from real and periodic ghost particles
      Nneib = 0;
      Nneib = tree->ComputeNeighbourAndGhostList
        (cell, sphdata, Nneibmax, Nneib, neiblist, neibpart);
      //Nneib = ghosttree->ComputeNeighbourList(cell,sphdata,Nneibmax,Nneib,neiblist,neibpart);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour list.
      while (Nneib == -1) {
        delete[] neibpartbuf[ithread];
        delete[] invdrmag;
        delete[] drmag;
        delete[] dr;
        delete[] sphlist;
        delete[] neiblist;
        Nneibmax                  = 2*Nneibmax;
        Nneibmaxbuf[ithread]      = Nneibmax;
        Ngravcellmaxbuf[ithread] *= 2;
        neiblist                  = new int[Nneibmax];
        sphlist                   = new int[Nneibmax];
        dr                        = new FLOAT[Nneibmax*ndim];
        drmag                     = new FLOAT[Nneibmax];
        invdrmag                  = new FLOAT[Nneibmax];
        neibpartbuf[ithread]      = new ParticleType<ndim>[Nneibmax];
        neibpart                  = neibpartbuf[ithread];
        Nneib = 0;
        Nneib = tree->ComputeNeighbourAndGhostList
          (cell, sphdata, Nneibmax, Nneib, neiblist, neibpart);
      };


      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        bool do_hydro = sph->types[activepart[j].ptype].hydro_forces ;
        if (do_hydro){
          Typemask hydromask  = sph->types[activepart[j].ptype].hydromask;

          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
          hrangesqdi = activepart[j].hrangesqd;
          Nhydroaux = 0;

          // Validate that gather neighbour list is correct
#if defined(VERIFY_ALL)
          if (neibcheck) this->CheckValidNeighbourList(i, sph->Nhydro + sph->NPeriodicGhost,
                                                       Nneib, neiblist, sphdata, "all");
#endif

          // Compute distances and the inverse between the current particle and all neighbours here,
          // for both gather and inactive scatter neibs.  Only consider particles with j > i to
          // compute pair forces once unless particle j is inactive.
          //-----------------------------------------------------------------------------------------
          for (jj=0; jj<Nneib; jj++) {

            // Skip non-hydro particles and the current active particle.
            if (hydromask[neibpart[jj].ptype] == false) continue;

            for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - rp[k];
            drsqd = DotProduct(draux, draux, ndim) + small_number;

            if (drsqd <= small_number) continue ;

            // Compute relative position and distance quantities for pair
            if (drsqd <= hrangesqdi || drsqd <= neibpart[jj].hrangesqd) {
              drmag[Nhydroaux] = sqrt(drsqd);
              invdrmag[Nhydroaux] = (FLOAT) 1.0/drmag[Nhydroaux];
              for (k=0; k<ndim; k++) dr[Nhydroaux*ndim + k] = draux[k]*invdrmag[Nhydroaux];
              levelneib[neiblist[jj]] = max(levelneib[neiblist[jj]],activepart[j].level);
              sphlist[Nhydroaux] = jj;
              Nhydroaux++;
            }
          }
          //-----------------------------------------------------------------------------------------

          // Compute all neighbour contributions to hydro forces
          sph->ComputeSphHydroForces(i,Nhydroaux,sphlist,drmag,invdrmag,dr,activepart[j],neibpart);
        }
      }
      //-------------------------------------------------------------------------------------------


      // Compute all star forces for active particles
      if (nbody->Nnbody > 0) {
        for (j=0; j<Nactive; j++) {
          if (activelist[j] < sph->Nhydro) {
            sph->ComputeStarGravForces(nbody->Nnbody,nbody->nbodydata,activepart[j]);
          }
        }
      }


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k] += activepart[j].a[k];
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
    for (i=0; i<sph->Ntot; i++) {
      for (int ithread=0; ithread<Nthreads; ithread++)
        sphdata[i].levelneib = max(sphdata[i].levelneib, levelneibbuf[ithread][i]);
    }


    // Free-up local memory for OpenMP thread
    delete[] invdrmag;
    delete[] drmag;
    delete[] dr;
    delete[] sphlist;
    delete[] neiblist;

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
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int cactive;                         // No. of active cells
  vector<TreeCellBase<ndim> > celllist;            // List of active cells
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[GradhSphTree::UpdateAllSphForces]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_ALL_FORCES");

  ParticleType<ndim>* sphdata =
      reinterpret_cast<ParticleType<ndim>*> (sph->GetSphParticleArray());

  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    return;
  }


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
    int i;                                       // Particle id
    int j;                                       // Aux. particle counter
    int jj;                                      // Aux. particle counter
    int k;                                       // Dimension counter
    int okflag;                                  // Flag if h-rho iteration is valid
    int Nactive;                                 // ..
    int Ndirect;                                 // ..
    int Ndirectaux;                              // ..
    int Ngrav;
    int Ngravcell;                               // No. of gravity cells
    int Nneib;                                   // No. of neighbours
    int Nhydroaux;                               // ..
    int Nhydroneib;                              // ..
    FLOAT aperiodic[ndim];                       // ..
    FLOAT macfactor;                             // Gravity MAC factor
    FLOAT draux[ndim];                           // Aux. relative position vector
    FLOAT drsqd;                                 // Distance squared
    FLOAT hrangesqdi;                            // Kernel gather extent
    FLOAT potperiodic;                           // ..
    FLOAT rp[ndim];                              // ..
    int Nneibmax     = Nneibmaxbuf[ithread];     // ..
    int Ngravcellmax = Ngravcellmaxbuf[ithread]; // ..
    int *activelist  = activelistbuf[ithread];   // ..
    int *neiblist    = new int[Nneibmax];        // ..
    int *sphlist     = new int[Nneibmax];        // ..
    int *sphauxlist  = new int[Nneibmax];        // ..
    int *directlist  = new int[Nneibmax];        // ..
    int *gravlist    = new int[Nneibmax];
    int *levelneib   = levelneibbuf[ithread];    // ..
    ParticleType<ndim>* activepart  = activepartbuf[ithread];   // ..
    ParticleType<ndim>* neibpart    = neibpartbuf[ithread];     // ..
    MultipoleMoment<ndim>* gravcell = cellbuf[ithread];         // ..

    // Zero timestep level array
    for (i=0; i<sph->Ntot; i++) levelneib[i] = 0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim> &cell = celllist[cc];
      macfactor = (FLOAT) 0.0;


      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, sphdata, activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
        for (j=0; j<Nactive; j++) macfactor =
          max(macfactor,pow((FLOAT) 1.0/activepart[j].gpot,twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        activepart[j].gpot      = (activepart[j].m/activepart[j].h)*sph->kernp->wpot(0.0);
        for (k=0; k<ndim; k++) activepart[j].a[k]     = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeGravityInteractionAndGhostList
        (cell, sphdata, macfactor, Nneibmax, Ngravcellmax, Nneib, Nhydroneib,
         Ndirect, Ngravcell, neiblist, sphlist, directlist, gravcell, neibpart);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour lists.
      while (okflag < 0 || Nneib > Nneibmax) {
        delete[] neibpartbuf[ithread];
        delete[] cellbuf[ithread];
        delete[] directlist;
        delete[] sphauxlist;
        delete[] sphlist;
        delete[] neiblist;
        delete[] gravlist;
        Nneibmax                 = 2*Nneibmax;
        Ngravcellmax             = 2*Ngravcellmax;
        Nneibmaxbuf[ithread]     = Nneibmax;
        Ngravcellmaxbuf[ithread] = Ngravcellmax;
        neiblist                 = new int[Nneibmax];
        sphlist                  = new int[Nneibmax];
        sphauxlist               = new int[Nneibmax];
        directlist               = new int[Nneibmax];
        gravlist                 = new int[Nneibmax];
        neibpartbuf[ithread]     = new ParticleType<ndim>[Nneibmax];
        cellbuf[ithread]         = new MultipoleMoment<ndim>[Ngravcellmax];
        neibpart                 = neibpartbuf[ithread];
        gravcell                 = cellbuf[ithread];
        okflag = tree->ComputeGravityInteractionAndGhostList
          (cell, sphdata, macfactor, Nneibmax, Ngravcellmax, Nneib, Nhydroneib,
           Ndirect, Ngravcell, neiblist, sphlist, directlist, gravcell, neibpart);
      };

      // Prune the directlist of non-gravitating particles
      Typemask gravmask;
      gravmask = sph->types.gravmask;

      for (j=0, i=0; j<Ndirect; j++)
      	if (gravmask[neibpart[directlist[j]].ptype]) {
      	  if (i != j) directlist[i] = directlist[j] ;
            i++ ;
      	}
      Ndirect = i ;

      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        //bool do_hydro = sph->types[activepart[j].ptype].hydro_forces ;
        bool do_grav  = sph->types[activepart[j].ptype].self_gravity ;

        Typemask hydromask ;
        for (k=0; k< Ntypes; ++k){
        	hydromask[k] = sph->types[activepart[j].ptype].hydromask[k] ;
        }

        Nhydroaux = 0;
        Ndirectaux = Ndirect;
        Ngrav = 0;
        for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
        hrangesqdi = activepart[j].hrangesqd;

        //-----------------------------------------------------------------------------------------
        for (jj=0; jj < Nhydroneib; jj++) {
          int ii = sphlist[jj];

          // Compute relative position and distance quantities for pair
          for (k=0; k<ndim; k++) draux[k] = neibpart[ii].r[k] - rp[k];
          drsqd = DotProduct(draux,draux,ndim) + small_number;

          if (drsqd <= small_number) continue ;

          // Record if neighbour is direct-sum or and SPH neighbour.
          // If SPH neighbour, also record max. timestep level for neighbour
          if (drsqd > hrangesqdi && drsqd >= neibpart[ii].hrangesqd && do_grav) {
            if (gravmask[neibpart[ii].ptype]) directlist[Ndirectaux++] = ii;
          }
          else {
            if (hydromask[neibpart[ii].ptype]){
              sphauxlist[Nhydroaux++] = ii;
              levelneib[neiblist[ii]] = max(levelneib[neiblist[ii]], activepart[j].level);
            }
            else if (gravmask[neibpart[ii].ptype] && do_grav){
              gravlist[Ngrav++] = ii;
            }
          }
        }
        //-----------------------------------------------------------------------------------------


        // Compute forces between SPH neighbours (hydro and gravity)
        if (Nhydroaux > 0)
          sph->ComputeSphHydroGravForces(i, Nhydroaux, sphauxlist, activepart[j], neibpart);

        if (do_grav){
          // Compute soften grav forces between non-SPH neighbours (hydro and gravity)
          sph->ComputeSphGravForces(i, Ngrav, gravlist, activepart[j], neibpart);


          // Compute direct gravity forces between distant particles
          sph->ComputeDirectGravForces(i, Ndirectaux, directlist, activepart[j], neibpart);

          // Compute gravitational force due to distant cells
          if (multipole == "monopole") {
            ComputeCellMonopoleForces(activepart[j].gpot, activepart[j].a,
                                      activepart[j].r, Ngravcell, gravcell);
          }
          else if (multipole == "quadrupole") {
            ComputeCellQuadrupoleForces(activepart[j].gpot, activepart[j].a,
                                        activepart[j].r, Ngravcell, gravcell);
         }

          // Add the periodic correction force for SPH and direct-sum neighbours
          if (simbox.PeriodicGravity){
            for (jj=0; jj<Nneib; jj++) {

        	  if (!gravmask[neibpart[jj].ptype]) continue ;

              for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - activepart[j].r[k];
              ewald->CalculatePeriodicCorrection(neibpart[jj].m, draux, aperiodic, potperiodic);
              for (k=0; k<ndim; k++) activepart[j].a[k] += aperiodic[k];
              activepart[j].gpot += potperiodic;
            }

            // Now add the periodic correction force for all cell COMs
            for (jj=0; jj<Ngravcell; jj++) {
              for (k=0; k<ndim; k++) draux[k] = gravcell[jj].r[k] - activepart[j].r[k];
              ewald->CalculatePeriodicCorrection(gravcell[jj].m, draux, aperiodic, potperiodic);
              for (k=0; k<ndim; k++) activepart[j].a[k] += aperiodic[k];
              activepart[j].gpot += potperiodic;
            }
          }
        }

      }
      //-------------------------------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole") {
        ComputeFastMonopoleForces(Nactive, Ngravcell, gravcell, cell, activepart);
      }
      else if (multipole == "fast_quadrupole") {
        ComputeFastQuadrupoleForces(Nactive, Ngravcell, gravcell, cell, activepart);
      }



      // Compute all star forces for active particles
      for (j=0; j<Nactive; j++) {
        sph->ComputeStarGravForces(nbody->Nnbody, nbody->nbodydata, activepart[j]);
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k] += activepart[j].a[k];
        sphdata[i].gpot  += activepart[j].gpot;
        sphdata[i].dudt  += activepart[j].dudt;
        sphdata[i].div_v += activepart[j].div_v;
        levelneib[i]      = max(levelneib[i],activepart[j].levelneib);
      }

    }
    //=============================================================================================


    // Propagate the changes in levelneib to the main array
#pragma omp for
    for (i=0; i<sph->Ntot; i++) {
      for (int ithread=0; ithread<Nthreads; ithread++)
        sphdata[i].levelneib = max(sphdata[i].levelneib, levelneibbuf[ithread][i]);
    }

    // Free-up local memory for OpenMP thread
    delete[] gravlist;
    delete[] directlist;
    delete[] sphauxlist;
    delete[] sphlist;
    delete[] neiblist;

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
//  GradhSphTree::UpdateAllSphGravForces
/// Compute all gravitational forces on active SPH particles for periodic boundary conditions.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void GradhSphTree<ndim,ParticleType>::UpdateAllSphGravForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] No. of SPH + ghost particles
  Sph<ndim> *sph,                      ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int cactive;                         // No. of active cells
  vector<TreeCellBase<ndim> > celllist;            // List of active cells
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

#ifdef MPI_PARALLEL
  double twork = timing->RunningTime();  // Start time (for load balancing)
#endif

  debug2("[GradhSphTree::UpdateAllSphGravForces]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SPH_ALL_GRAV_FORCES");

  ParticleType<ndim>* sphdata =
      reinterpret_cast<ParticleType<ndim>*> (sph->GetSphParticleArray());

  // Find list of all cells that contain active particles
  cactive = tree->ComputeActiveCellList(celllist);

  // If there are no active cells, return to main loop
  if (cactive == 0) {
    return;
  }


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
    int i;                                       // Particle id
    int j;                                       // Aux. particle counter
    int jj;                                      // Aux. particle counter
    int k;                                       // Dimension counter
    int okflag;                                  // Flag if h-rho iteration is valid
    int Nactive;                                 // ..
    int Ndirect;                                 // ..
    int Ndirectaux;                              // ..
    int Ngravcell;                               // No. of gravity cells
    int Nneib;                                   // No. of neighbours
    int Nhydroaux;                               // ..
    int Nhydroneib;                              // ..
    FLOAT aperiodic[ndim];                       // ..
    FLOAT macfactor;                             // Gravity MAC factor
    FLOAT draux[ndim];                           // Aux. relative position vector
    FLOAT drsqd;                                 // Distance squared
    FLOAT hrangesqdi;                            // Kernel gather extent
    FLOAT potperiodic;                           // ..
    FLOAT rp[ndim];                              // ..
    int Nneibmax     = Nneibmaxbuf[ithread];     // ..
    int Ngravcellmax = Ngravcellmaxbuf[ithread]; // ..
    int *activelist  = activelistbuf[ithread];   // ..
    int *neiblist    = new int[Nneibmax];        // ..
    int *sphlist     = new int[Nneibmax];        // ..
    int *sphauxlist  = new int[Nneibmax];        // ..
    int *directlist  = new int[Nneibmax];        // ..
    int	*gravlist    = new int[Nneibmax];        // ..
    int *levelneib   = levelneibbuf[ithread];    // ..
    ParticleType<ndim>* activepart = activepartbuf[ithread];   // ..
    ParticleType<ndim>* neibpart   = neibpartbuf[ithread];     // ..
    MultipoleMoment<ndim>* gravcell       = cellbuf[ithread];         // ..

    // Zero timestep level array
    for (i=0; i<sph->Ntot; i++) levelneib[i] = 0;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim> &cell = celllist[cc];
      macfactor = (FLOAT) 0.0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
        for (j=0; j<Nactive; j++) macfactor =
          max(macfactor,pow((FLOAT) 1.0/activepart[j].gpot,twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].div_v     = (FLOAT) 0.0;
        activepart[j].dudt      = (FLOAT) 0.0;
        activepart[j].levelneib = 0;
        activepart[j].gpot      = (activepart[j].m/activepart[j].h)*sph->kernp->wpot(0.0);
        for (k=0; k<ndim; k++) activepart[j].a[k]     = (FLOAT) 0.0;
      }

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeGravityInteractionAndGhostList
        (cell, sphdata, macfactor, Nneibmax, Ngravcellmax, Nneib, Nhydroneib,
         Ndirect, Ngravcell, neiblist, sphlist, directlist, gravcell, neibpart);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour lists.
      while (okflag < 0 || Nneib > Nneibmax) {
        delete[] neibpartbuf[ithread];
        delete[] cellbuf[ithread];
        delete[] gravlist;
        delete[] directlist;
        delete[] sphauxlist;
        delete[] sphlist;
        delete[] neiblist;
        Nneibmax                 = 2*Nneibmax;
        Ngravcellmax             = 2*Ngravcellmax;
        Nneibmaxbuf[ithread]     = Nneibmax;
        Ngravcellmaxbuf[ithread] = Ngravcellmax;
        neiblist                 = new int[Nneibmax];
        sphlist                  = new int[Nneibmax];
        sphauxlist               = new int[Nneibmax];
        directlist               = new int[Nneibmax];
        gravlist                 = new int[Nneibmax];
        neibpartbuf[ithread]     = new ParticleType<ndim>[Nneibmax];
        cellbuf[ithread]         = new MultipoleMoment<ndim>[Ngravcellmax];
        neibpart                 = neibpartbuf[ithread];
        gravcell                 = cellbuf[ithread];
        okflag = tree->ComputeGravityInteractionAndGhostList
          (cell, sphdata, macfactor, Nneibmax, Ngravcellmax, Nneib, Nhydroneib,
           Ndirect, Ngravcell, neiblist, sphlist, directlist, gravcell, neibpart);
      };

      // Prune the directlist of non-gravitating particles
      Typemask gravmask;
      gravmask = sph->types.gravmask;

      for (j=0, i=0; j<Ndirect; j++)
      	if (gravmask[neibpart[directlist[j]].ptype]) {
      	  if (i != j) directlist[i] = directlist[j] ;
            i++ ;
      	}
      Ndirect = i ;

      // Loop over all active particles in the cell
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        if (sph->types[activepart[j].ptype].self_gravity){

          Nhydroaux = 0;
          Ndirectaux = Ndirect;
          for (k=0; k<ndim; k++) rp[k] = activepart[j].r[k];
          hrangesqdi = activepart[j].hrangesqd;

          //-----------------------------------------------------------------------------------------
          for (jj=0; jj<Nhydroneib; jj++) {
            int ii = sphlist[jj];

            // Skip non-gravitating particles and the current active particle.
            if (gravmask[neibpart[jj].ptype] == false) continue ;

            // Compute relative position and distance quantities for pair
            for (k=0; k<ndim; k++) draux[k] = neibpart[ii].r[k] - rp[k];
            drsqd = DotProduct(draux, draux, ndim) + small_number;

            if (drsqd <= small_number) continue ;

            // Record if neighbour is direct-sum or and SPH neighbour.
            // If SPH neighbour, also record max. timestep level for neighbour
            if (drsqd > hrangesqdi && drsqd >= neibpart[ii].hrangesqd) {
              directlist[Ndirectaux++] = ii;
           }
            else {
              sphauxlist[Nhydroaux++] = ii;
              levelneib[neiblist[ii]] = max(levelneib[neiblist[ii]], activepart[j].level);
            }
          }

          //-----------------------------------------------------------------------------------------


          // Compute forces between SPH neighbours (hydro and gravity)
          sph->ComputeSphGravForces(i, Nhydroaux, sphauxlist, activepart[j], neibpart);

          // Compute direct gravity forces between distant particles
          sph->ComputeDirectGravForces(i, Ndirectaux, directlist, activepart[j], neibpart);

          // Compute gravitational force due to distant cells
          if (multipole == "monopole") {
            ComputeCellMonopoleForces(activepart[j].gpot, activepart[j].a,
                                      activepart[j].r, Ngravcell, gravcell);
          }
          else if (multipole == "quadrupole") {
            ComputeCellQuadrupoleForces(activepart[j].gpot, activepart[j].a,
                                        activepart[j].r, Ngravcell, gravcell);
          }

          // Add the periodic correction force for SPH and direct-sum neighbours
          if (simbox.PeriodicGravity){
            for (jj=0; jj<Nneib; jj++) {
              for (k=0; k<ndim; k++) draux[k] = neibpart[jj].r[k] - activepart[j].r[k];
              ewald->CalculatePeriodicCorrection(neibpart[jj].m, draux, aperiodic, potperiodic);
              for (k=0; k<ndim; k++) activepart[j].a[k] += aperiodic[k];
              activepart[j].gpot += potperiodic;
            }

            // Now add the periodic correction force for all cell COMs
            for (jj=0; jj<Ngravcell; jj++) {
              for (k=0; k<ndim; k++) draux[k] = gravcell[jj].r[k] - activepart[j].r[k];
              ewald->CalculatePeriodicCorrection(gravcell[jj].m, draux, aperiodic, potperiodic);
              for (k=0; k<ndim; k++) activepart[j].a[k] += aperiodic[k];
              activepart[j].gpot += potperiodic;
            }
          }
        }
      }
      //-------------------------------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole") {
        ComputeFastMonopoleForces(Nactive, Ngravcell, gravcell, cell, activepart);
      }
      else if (multipole == "fast_quadrupole") {
        ComputeFastQuadrupoleForces(Nactive, Ngravcell, gravcell, cell, activepart);
      }

      // Compute all star forces for active particles
      for (j=0; j<Nactive; j++) {
        sph->ComputeStarGravForces(nbody->Nnbody, nbody->nbodydata, activepart[j]);
      }

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k] += activepart[j].a[k];
        sphdata[i].gpot  += activepart[j].gpot;
        sphdata[i].dudt  += activepart[j].dudt;
        sphdata[i].div_v += activepart[j].div_v;
        levelneib[i]      = max(levelneib[i],activepart[j].levelneib);
      }

    }
    //=============================================================================================


    // Propagate the changes in levelneib to the main array
#pragma omp for
    for (i=0; i<sph->Ntot; i++) {
      for (int ithread=0; ithread<Nthreads; ithread++)
        sphdata[i].levelneib = max(sphdata[i].levelneib, levelneibbuf[ithread][i]);
    }

    // Free-up local memory for OpenMP thread
    delete[] gravlist;
    delete[] directlist;
    delete[] sphauxlist;
    delete[] sphlist;
    delete[] neiblist;

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

