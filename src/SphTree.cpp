//=============================================================================
//  SphTree.cpp
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
//=============================================================================


#include <cstdlib>
#include <cassert>
#include <iostream>
#include <numeric>
#include <string>
#include <math.h>
#include <vector>
#include "Precision.h"
#include "Exception.h"
#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "SphParticle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
using namespace std;



//=============================================================================
//  SphTree::SphTree
/// SphTree constructor.  Initialises various variables.
//=============================================================================
template <int ndim, template<int> class ParticleType>
SphTree<ndim,ParticleType>::SphTree
(int Nleafmaxaux,
 int Nmpiaux,
 FLOAT thetamaxsqdaux,
 FLOAT kernrangeaux,
 FLOAT macerroraux,
 string gravity_mac_aux,
 string multipole_aux,
 DomainBox<ndim> *boxaux,
 SphKernel<ndim> *kernaux,
 CodeTiming *timingaux):
  SphNeighbourSearch<ndim>(kernrangeaux,boxaux,kernaux,timingaux),
  Nleafmax(Nleafmaxaux),
  Nmpi(Nmpiaux),
  thetamaxsqd(thetamaxsqdaux),
  invthetamaxsqd(1.0/thetamaxsqdaux),
  gravity_mac(gravity_mac_aux),
  macerror(macerroraux),
  multipole(multipole_aux)
{
  allocated_buffer = false;
  neibcheck        = true;
  Ntot             = 0;
  Ntotmax          = 0;
  Ntotmaxold       = 0;
#if defined _OPENMP
  Nthreads         = omp_get_max_threads();
#else
  Nthreads         = 1;
#endif

  // Set-up main tree object
  tree = new KDTree<ndim,ParticleType>(Nleafmaxaux, thetamaxsqdaux,
                                       kernrangeaux, macerroraux,
                                       gravity_mac_aux, multipole_aux);

  // Set-up ghost-particle tree object
  ghosttree = new KDTree<ndim,ParticleType>(Nleafmaxaux, thetamaxsqdaux,
                                            kernrangeaux, macerroraux,
                                            gravity_mac_aux, multipole_aux);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = new KDTree<ndim,ParticleType>(Nleafmaxaux, thetamaxsqdaux,
                                               kernrangeaux, macerroraux,
                                               gravity_mac_aux, multipole_aux);

  // Set-up multiple pruned trees, one for each MPI process
  prunedtree = new KDTree<ndim,ParticleType>*[Nmpi];
  for (int j=0; j<Nmpi; j++) {
    prunedtree[j] = new KDTree<ndim,ParticleType>
      (Nleafmaxaux, thetamaxsqdaux, kernrangeaux, macerroraux,
       gravity_mac_aux, multipole_aux);
  }

  Ncellexport = new int[Nmpi];
  Npartexport = new int[Nmpi];
  cellexportlist = new int*[Nmpi];
  for (int j=0; j<Nmpi; j++) cellexportlist[j] = new int[1];

  ids_sent_particles.resize(Nmpi);
#endif

}



//=============================================================================
//  SphTree::~SphTree
/// SphTree destructor.  Deallocates tree memory upon object destruction.
//=============================================================================
template <int ndim, template<int> class ParticleType>
SphTree<ndim,ParticleType>::~SphTree()
{
  if (tree->allocated_tree) {
    DeallocateMemory();
    tree->DeallocateTreeMemory();
  }
}



//=============================================================================
//  SphTree::AllocateMemory
/// Allocate memory for tree as requested.  If more memory is required
/// than currently allocated, tree is deallocated and reallocated here.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::AllocateMemory(Sph<ndim> *sph)
{
  int ithread;                      // Thread id number

  debug2("[SphTree::AllocateMemory]");

  if (!allocated_buffer) {

    Nneibmaxbuf = new int[Nthreads];
    Ndirectmaxbuf = new int[Nthreads];
    Ngravcellmaxbuf = new int[Nthreads];
    levelneibbuf = new int*[Nthreads];
    activelistbuf = new int*[Nthreads];
    activepartbuf = new ParticleType<ndim>*[Nthreads];
    neibpartbuf = new ParticleType<ndim>*[Nthreads];

    for (ithread=0; ithread<Nthreads; ithread++) {
      Nneibmaxbuf[ithread] = max(1,4*sph->Ngather);
      Ndirectmaxbuf[ithread] = max(1,4*sph->Ngather);
      Ngravcellmaxbuf[ithread] = max(1,4*sph->Ngather);
      levelneibbuf[ithread] = new int[Ntotmax];
      activelistbuf[ithread] = new int[Nleafmax];
      activepartbuf[ithread] = new ParticleType<ndim>[Nleafmax];
      neibpartbuf[ithread] = new ParticleType<ndim>[Nneibmaxbuf[ithread]];
    }
    allocated_buffer = true;

  }

  return;
}



//=============================================================================
//  SphTree::DeallocateTreeMemory
/// Deallocates all binary tree memory
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::DeallocateMemory(void)
{
  int ithread;                      // Thread id number

  debug2("[SphTree::DeallocateTreeMemory]");

  if (allocated_buffer) {

    for (ithread=0; ithread<Nthreads; ithread++) {
      delete[] neibpartbuf[ithread];
      delete[] activepartbuf[ithread];
      delete[] levelneibbuf[ithread];
    }
    delete[] neibpartbuf;
    delete[] activepartbuf;
    delete[] activelistbuf;
    delete[] levelneibbuf;
    delete[] Ngravcellmaxbuf;
    delete[] Ndirectmaxbuf;
    delete[] Nneibmaxbuf;

  }

  return;
}



//=============================================================================
//  SphTree::BuildTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::BuildTree
(bool rebuild_tree,                 ///< Flag to rebuild tree
 int n,                             ///< Integer time
 int ntreebuildstep,                ///< Tree build frequency
 int ntreestockstep,                ///< Tree stocking frequency
 int Npart,                         ///< No. of particles
 int Npartmax,                      ///< Max. no. of particles
 SphParticle<ndim> *sph_gen,        ///< Particle data array
 Sph<ndim> *sph,                    ///< Pointer to SPH object
 FLOAT timestep)                    ///< Smallest physical timestep
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphTree::BuildTree]");
  timing->StartTimingSection("BUILD_TREE",2);

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif


  // For tree rebuild steps
  //---------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    // Delete any dead particles from main SPH arrays before we re-build tree
    sph->DeleteDeadParticles();

    Ntotold = Ntot;
    Ntot = sph->Ntot;
    Ntotmaxold = Ntotmax;
    Ntotmax = max(Ntotmax,Ntot);
    Ntotmax = max(Ntotmax,sph->Nsphmax);

    tree->ifirst = 0;
    tree->ilast = sph->Nsph - 1;
    tree->Ntot = sph->Nsph;
    tree->Ntotmaxold = tree->Ntotmax;
    tree->Ntotmax = max(tree->Ntotmax,tree->Ntot);
    tree->Ntotmax = max(tree->Ntotmax,sph->Nsphmax);
    tree->BuildTree(Npart,Npartmax,sphdata,timestep);

    AllocateMemory(sph);
#ifdef MPI_PARALLEL
    if (tree->Ntotmax > tree->Ntotmaxold) {
      for (i=Nmpi-1; i>=0; i--) delete[] cellexportlist[i];
      for (i=0; i<Nmpi; i++) cellexportlist[i] = new int[tree->gmax];
    }
#endif

  }

  // Else stock the tree
  //---------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    tree->StockTree(tree->kdcell[0],sphdata);

  }

  // Otherwise simply extrapolate tree cell properties
  //---------------------------------------------------------------------------
  else {

    tree->ExtrapolateCellProperties(timestep);

  }
  //---------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  timing->EndTimingSection("BUILD_TREE");

  return;
}



//=============================================================================
//  SphTree::BuildGhostTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::BuildGhostTree
(bool rebuild_tree,                 ///< Flag to rebuild tree
 int n,                             ///< Integer time
 int ntreebuildstep,                ///< Tree build frequency
 int ntreestockstep,                ///< Tree stocking frequency
 int Npart,                         ///< No. of particles
 int Npartmax,                      ///< Max. no. of particles
 SphParticle<ndim> *sph_gen,        ///< Particle data array
 Sph<ndim> *sph,                    ///< Pointer to SPH object
 FLOAT timestep)                    ///< Smallest physical timestep
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  // If no periodic ghosts exist, do not build tree
  //if (sph->NPeriodicGhost == 0) return;

  debug2("[SphTree::BuildGhostTree]");
  timing->StartTimingSection("BUILD_GHOST_TREE",2);

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif


  // For tree rebuild steps
  //---------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    ghosttree->ifirst = sph->Nsph;
    ghosttree->ilast = sph->Nsph + sph->NPeriodicGhost - 1;
    ghosttree->Ntot = sph->NPeriodicGhost;
    ghosttree->Ntotmaxold = ghosttree->Ntotmax;
    ghosttree->Ntotmax = max(ghosttree->Ntotmax,ghosttree->Ntot);
    ghosttree->Ntotmax = max(ghosttree->Ntotmax,sph->Nsphmax);
    ghosttree->BuildTree(ghosttree->Ntot,ghosttree->Ntotmax,sphdata,timestep);

  }

  // Else stock the tree
  //---------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    ghosttree->StockTree(ghosttree->kdcell[0],sphdata);

  }

  // Otherwise simply extrapolate tree cell properties
  //---------------------------------------------------------------------------
  else {

    //ExtrapolateCellProperties(kdcell[0],timestep);
    ghosttree->ExtrapolateCellProperties(timestep);

  }
  //---------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  timing->EndTimingSection("BUILD_GHOST_TREE");


  return;
}



//=============================================================================
//  SphTree::GetGatherNeighbourList
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
int SphTree<ndim,ParticleType>::GetGatherNeighbourList
(FLOAT rp[ndim],                    ///< Position vector
 FLOAT rsearch,                     ///< Gather search radius
 SphParticle<ndim> *sph_gen,        ///< Pointer to SPH particle array
 int Nsph,                          ///< No. of SPH particles
 int Nneibmax,                      ///< Max. no. of neighbours
 int *neiblist)                     ///< List of neighbouring particles
{
  int Nneib = 0;                    // No. of (non-dead) neighbours
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[SphTree::GetGatherNeighbourList]");

  Nneib = tree->ComputeGatherNeighbourList(sphdata,rp,rsearch,
                                           Nneibmax,neiblist);

  return Nneib;
}



//=============================================================================
//  SphTree::UpdateAllSphDerivatives
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateAllSphDerivatives
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph)                    ///< [in] Pointer to SPH object
{
  return;
}



//=============================================================================
//  SphTree::UpdateAllSphDudt
/// Compute all local 'gather' properties of currently active particles, and
/// then compute each particle's contribution to its (active) neighbour
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateAllSphDudt
(int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph)                    ///< [in] Pointer to SPH object
{
  return;
}


//=============================================================================
//  SphTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in the tree and update all active
/// particle counters.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateActiveParticleCounters(SphParticle<ndim> * sphdata_gen,
                                            Sph<ndim> * sph) {

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sphdata_gen);
  tree->UpdateActiveParticleCounters(sphdata);
}



//=============================================================================
//  SphTree::SearchBoundaryGhostParticles
/// Search domain to create any required ghost particles near any boundaries.
/// Currently only searches to create periodic or mirror ghost particles.
//=============================================================================
template <int ndim, template <int> class ParticleType>
void SphTree<ndim, ParticleType >::SearchBoundaryGhostParticles
(FLOAT tghost,                      ///< Ghost particle 'lifetime'
 DomainBox<ndim> simbox,            ///< Simulation box structure
 Sph<ndim> *sph)                    ///< Sph object pointer
{
  int c;
  int i;                            // Particle counter
  const FLOAT grange = ghost_range*kernrange;
  KDTreeCell<ndim> *cell;           // ..
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray());

  // Set all relevant particle counters
  sph->Nghost         = 0;
  sph->NPeriodicGhost = 0;
  sph->Nmpighost      = 0;
  sph->Nghostmax      = sph->Nsphmax - sph->Nsph;
  sph->Ntot           = sph->Nsph;


  // If all boundaries are open, immediately return to main loop
  if (simbox.x_boundary_lhs == "open" && simbox.x_boundary_rhs == "open" &&
      simbox.y_boundary_lhs == "open" && simbox.y_boundary_rhs == "open" &&
      simbox.z_boundary_lhs == "open" && simbox.z_boundary_rhs == "open")
    return;


  debug2("[SphTree::SearchBoundaryGhostParticles]");


  // Create ghost particles in x-dimension
  //===========================================================================
  if ((simbox.x_boundary_lhs == "open" &&
       simbox.x_boundary_rhs == "open") == 0) {

    // Start from root-cell
    c = 0;

    //-------------------------------------------------------------------------
    while (c < tree->Ncell) {
      cell = &(tree->kdcell[c]);

      // If x-bounding box overlaps edge of x-domain, open cell
      //-----------------------------------------------------------------------
      if (cell->bbmin[0] + min(0.0,cell->v[0]*tghost) <
          simbox.boxmin[0] + grange*cell->hmax ||
          cell->bbmax[0] + max(0.0,cell->v[0]*tghost) >
          simbox.boxmax[0] - grange*cell->hmax) {

        // If not a leaf-cell, then open cell to first child cell
        if (cell->level != tree->ltot)
          c++;

        else if (cell->N == 0)
          c = cell->cnext;

        // If leaf-cell, check through particles in turn to find ghosts
        else if (cell->level == tree->ltot) {
          i = cell->ifirst;
    	    while (i != -1) {
            sph->CheckXBoundaryGhostParticle(i,tghost,simbox);
            if (i == cell->ilast) break;
    	      i = tree->inext[i];
          };
          c = cell->cnext;
        }
      }

      // If not in range, then open next cell
      //-----------------------------------------------------------------------
      else
        c = cell->cnext;

    }
    //-------------------------------------------------------------------------

    sph->Ntot = sph->Nsph + sph->Nghost;
  }


  // Create ghost particles in y-dimension
  //===========================================================================
  if (ndim >= 2 && (simbox.y_boundary_lhs == "open" &&
		    simbox.y_boundary_rhs == "open") == 0) {

    // Start from root-cell
    c = 0;

    //-------------------------------------------------------------------------
    while (c < tree->Ncell) {
      cell = &(tree->kdcell[c]);

      // If x-bounding box overlaps edge of x-domain, open cell
      //-----------------------------------------------------------------------
      if (cell->bbmin[1] + min(0.0,cell->v[1]*tghost) <
          simbox.boxmin[1] + grange*cell->hmax ||
          cell->bbmax[1] + max(0.0,cell->v[1]*tghost) >
          simbox.boxmax[1] - grange*cell->hmax) {

        // If not a leaf-cell, then open cell to first child cell
        if (cell->level != tree->ltot)
          c++;

        else if (cell->N == 0)
          c = cell->cnext;

        // If leaf-cell, check through particles in turn to find ghosts
        else if (cell->level == tree->ltot) {
          i = cell->ifirst;
    	    while (i != -1) {
            sph->CheckYBoundaryGhostParticle(i,tghost,simbox);
            if (i == cell->ilast) break;
    	      i = tree->inext[i];
          };
          c = cell->cnext;
        }
      }

      // If not in range, then open next cell
      //-----------------------------------------------------------------------
      else
        c = cell->cnext;

    }
    //-------------------------------------------------------------------------


    // Check x-ghosts (which are not part of tree) by direct-sum
    for (i=sph->Nsph; i<sph->Ntot; i++) {
      sph->CheckYBoundaryGhostParticle(i,tghost,simbox);
    }

    sph->Ntot = sph->Nsph + sph->Nghost;
  }


  // Create ghost particles in z-dimension
  //===========================================================================
  if (ndim == 3 && (simbox.z_boundary_lhs == "open" &&
		    simbox.z_boundary_rhs == "open") == 0) {

    // Start from root-cell
    c = 0;

    //-------------------------------------------------------------------------
    while (c < tree->Ncell) {
      cell = &(tree->kdcell[c]);

      // If x-bounding box overlaps edge of x-domain, open cell
      //-----------------------------------------------------------------------
      if (cell->bbmin[2] + min(0.0,cell->v[2]*tghost) <
          simbox.boxmin[2] + grange*cell->hmax ||
          cell->bbmax[2] + max(0.0,cell->v[2]*tghost) >
          simbox.boxmax[2] - grange*cell->hmax) {

        // If not a leaf-cell, then open cell to first child cell
        if (cell->level != tree->ltot)
          c++;

        else if (cell->N == 0)
          c = cell->cnext;

        // If leaf-cell, check through particles in turn to find ghosts
        else if (cell->level == tree->ltot) {
          i = cell->ifirst;
    	    while (i != -1) {
            sph->CheckZBoundaryGhostParticle(i,tghost,simbox);
            if (i == cell->ilast) break;
    	      i = tree->inext[i];
          };
          c = cell->cnext;
        }
      }

      // If not in range, then open next cell
      //-----------------------------------------------------------------------
      else
        c = cell->cnext;

    }
    //-------------------------------------------------------------------------


    // Check x- and y-ghosts (which are not part of tree) by direct-sum
    for (i=sph->Nsph; i<sph->Ntot; i++) {
      sph->CheckZBoundaryGhostParticle(i,tghost,simbox);
    }


    sph->Ntot = sph->Nsph + sph->Nghost;
  }


  // Quit here if we've run out of memory for ghosts
  if (sph->Ntot > sph->Nsphmax) {
    string message="Not enough memory for ghost particles";
    ExceptionHandler::getIstance().raise(message);
  }

  sph->NPeriodicGhost = sph->Nghost;

  return;
}



#ifdef MPI_PARALLEL
//=============================================================================
//  GradhSphTree::UpdateGravityExportList
/// Compute all local 'gather' properties of currently active particles, and
/// then compute each particle's contribution to its (active) neighbour
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateGravityExportList
(int rank,                          ///< [in] ..
 int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int ithread;                      // ..
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  int k;                            // Dimension counter
  int okflag;                       // Flag if h-rho iteration is valid
  int Nactive;                      // ..
  int Ngravcell;                    // No. of gravity cells
  int Ngravcellmax;                 // ..
  int Ngravcelltemp;                // ..
  FLOAT macfactor;                  // ..

  int *activelist;                  // ..
  KDTreeCell<ndim> *cell;           // Pointer to binary tree cell
  KDTreeCell<ndim> **celllist;      // List of pointers to binary tree cells
  KDTreeCell<ndim> **gravcelllist;  // List of pointers to grav. cells
  ParticleType<ndim> *activepart;   // ..
  ParticleType<ndim> *neibpart;     // ..
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  debug2("[GradhSphTree::UpdateDistantSphForces]");
  timing->StartTimingSection("SPH_DISTANT_FORCES",2);


  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[2*tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);

  // Reset all export lists
  for (j=0; j<Nmpi; j++) Ncellexport[j] = 0;


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) shared(celllist,cactive,sph,sphdata,cout) \
  private(activepart,activelist,cc,cell,directlist,draux,drsqd,gravcelllist) \
  private(hrangesqdi,i,interactlist,ithread,j,jj,k,levelneib,macfactor)	\
  private(neiblist,neibpart,Nactive,Ndirect,Ndirectaux,Ndirectmax) \
  private(Ngravcell,Ngravcellmax,Ninteract,Nneib,Nneibmax,okflag,rp)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif
    Ngravcellmax = Nprunedcellmax;
    activelist = activelistbuf[ithread];
    activepart = activepartbuf[ithread];
    gravcelllist = new KDTreeCell<ndim>*[Ngravcellmax];


    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];
      macfactor = 0.0;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell,sphdata,activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = sphdata[activelist[j]];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") {
        for (j=0; j<Nactive; j++)
          macfactor = max(macfactor,pow(1.0/activepart[j].gpot,twothirds));
      }

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].gpot = activepart[j].m*activepart[j].invh*sph->kernp->wpot(0.0);
      }


      // Loop over all distant pruned trees and compute list of cells.
      // If pruned tree is too close, record cell id for exporting
      //-----------------------------------------------------------------------
      for (j=0; j<Nmpi; j++) {
        if (j == rank) continue;

        Ngravcelltemp = prunedtree[j]->ComputeDistantGravityInteractionList
          (cell,macfactor,Ngravcellmax,Ngravcell,gravcelllist);

        // If pruned tree is too close (flagged by -1), then record cell id
	      // for exporting to other MPI processes
        if (Ngravcelltemp == -1)
          cellexportlist[j][Ncellexport[j]++] = cell->id;
        else
	        Ngravcell = Ngravcelltemp;

      }
      //-----------------------------------------------------------------------


      // Loop over all active particles in the cell
      //-----------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {
        i = activelist[j];

        // Compute gravitational force due to distant cells
        if (multipole == "monopole")
          tree->ComputeCellMonopoleForces(activepart[j].gpot,
                                          activepart[j].agrav,
                                          activepart[j].r,
                                          Ngravcell,gravcelllist);
        else if (multipole == "quadrupole")
          tree->ComputeCellQuadrupoleForces(activepart[j].gpot,
                                            activepart[j].agrav,
                                            activepart[j].r,
                                            Ngravcell,gravcelllist);

      }
      //-----------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == "fast_monopole")
        tree->ComputeFastMonopoleForces(Nactive,Ngravcell,gravcelllist,
				                                cell,activepart);

      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
    	i = activelist[j];
        for (k=0; k<ndim; k++) sphdata[i].a[k] = activepart[j].a[k];
        for (k=0; k<ndim; k++) sphdata[i].agrav[k] = activepart[j].agrav[k];
        for (k=0; k<ndim; k++) sphdata[i].a[k] += sphdata[i].agrav[k];
        sphdata[i].gpot = activepart[j].gpot;
      }

    }
    //=========================================================================


    // Free-up local memory for OpenMP thread
    delete[] gravcelllist;

  }
  //===========================================================================

  delete[] celllist;

  timing->EndTimingSection("SPH_DISTANT_FORCES");

  return;
}



//=============================================================================
//  SphTree::UpdateHydroExportList
/// Compute all local 'gather' properties of currently active particles, and
/// then compute each particle's contribution to its (active) neighbour
/// neighbour hydro forces.  Optimises the algorithm by using grid-cells to
/// construct local neighbour lists for all particles  inside the cell.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UpdateHydroExportList
(int rank,                          ///< [in] ..
 int Nsph,                          ///< [in] No. of SPH particles
 int Ntot,                          ///< [in] No. of SPH + ghost particles
 SphParticle<ndim> *sph_gen,        ///< [inout] Pointer to SPH ptcl array
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 Nbody<ndim> *nbody)                ///< [in] Pointer to N-body object
{
  bool overlapflag;                 // ..
  int cactive;                      // No. of active cells
  int cc;                           // Aux. cell counter
  int i;                            // Particle id
  int ithread;                      // ..
  int j;                            // Aux. particle counter
  int jj;                           // Aux. particle counter
  KDTreeCell<ndim> *cell;           // Pointer to binary tree cell
  KDTreeCell<ndim> **celllist;      // List of pointers to binary tree cells


  debug2("[SphTree::UpdateDistantSphForces]");
  timing->StartTimingSection("MPI_HYDRO_EXPORT",2);


  // Find list of all cells that contain active particles
  celllist = new KDTreeCell<ndim>*[2*tree->gtot];
  cactive = tree->ComputeActiveCellList(celllist);

  // Reset all export lists
  for (j=0; j<Nmpi; j++) Ncellexport[j] = 0;


  // Set-up all OMP threads
  //===========================================================================
#pragma omp parallel default(none) shared(celllist,cactive,cout) \
  private(cc,cell,i,ithread,j,jj,overlapflag)
  {
#if defined _OPENMP
    ithread = omp_get_thread_num();
#else
    ithread = 0;
#endif

    // Loop over all active cells
    //=========================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      cell = celllist[cc];

      // Loop over all distant pruned trees and compute list of cells.
      // If pruned tree is too close, record cell id for exporting
      //-----------------------------------------------------------------------
      for (j=0; j<Nmpi; j++) {
        if (j == rank) continue;

        overlapflag = prunedtree[j]->ComputeHydroTreeCellOverlap(cell);

        // If pruned tree is too close (flagged by -1), then record cell id
        // for exporting to other MPI processes
        if (overlapflag) cellexportlist[j][Ncellexport[j]++] = cell->id;

      }
      //-----------------------------------------------------------------------

    }
    //=========================================================================

  }
  //===========================================================================

  delete[] celllist;

  timing->EndTimingSection("MPI_HYDRO_EXPORT");

  return;
}



//=============================================================================
//  SphTree::BuildPrunedTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::BuildPrunedTree
(int pruning_level,
 int rank)
{
  int c;                            // Cell counter
  int cnew;                         // ..
  int cnext;                        // ..
  int c1;                           // ..
  int c2;                           // ..
  int i;                            // Particle counter
  int k;                            // Dimension counter

  debug2("[SphTree::BuildPrunedTree]");
  timing->StartTimingSection("BUILD_PRUNED_TREE",2);

  cnew = 0;
  Nprunedcellmax = 0;


  // Set level at which tree will be pruned (for all trees)
  //---------------------------------------------------------------------------
  for (i=0; i<Nmpi; i++) {
    prunedtree[i]->ltot_old = prunedtree[i]->ltot;
    prunedtree[i]->ltot     = pruning_level;
    prunedtree[i]->gmax     = pow(2,pruning_level);
    prunedtree[i]->Ncellmax = 2*prunedtree[i]->gmax - 1;
    prunedtree[i]->Ncell    = 2*prunedtree[i]->gmax - 1;
    Nprunedcellmax += prunedtree[i]->Ncellmax;

    // Allocate (or reallocate if needed) all tree memory
    prunedtree[i]->AllocateTreeMemory();

    // If the number of levels in the tree has changed (due to destruction or
    // creation of new particles) then re-create tree data structure
    // including linked lists and cell pointers
    if (prunedtree[i]->ltot != prunedtree[i]->ltot_old)
      prunedtree[i]->CreateTreeStructure();

  }
  //---------------------------------------------------------------------------


  // Now walk through main tree cell-by-cell and copy all important data
  // to pruned tree cells
  //---------------------------------------------------------------------------
  for (c=0; c<tree->Ncell; c++) {

    // If cell is on a lower level, skip over
    if (tree->kdcell[c].level > pruning_level) continue;

    // Otherwise, record all data from cell, except for cell links which
    // are maintained to ensure a valid tree
    c1 = prunedtree[rank]->kdcell[cnew].c1;
    c2 = prunedtree[rank]->kdcell[cnew].c2;
    cnext = prunedtree[rank]->kdcell[cnew].cnext;
    prunedtree[rank]->kdcell[cnew] = tree->kdcell[c];
    prunedtree[rank]->kdcell[cnew].c1 = c1;
    prunedtree[rank]->kdcell[cnew].c2 = c2;
    prunedtree[rank]->kdcell[cnew].cnext = cnext;

    cnew++;

  }
  //---------------------------------------------------------------------------

  cout << "Pruned tree size : " << prunedtree[rank]->Ncell << "    "
       << prunedtree[rank]->ltot << "    " << pruning_level << endl;
  for (c=0; c<prunedtree[rank]->Ncell; c++) {
    cout << "bb[" << c << "] : " << prunedtree[rank]->kdcell[c].bbmin[0]
         << "    " << prunedtree[rank]->kdcell[c].bbmax[0]
         << "    level : " << prunedtree[rank]->kdcell[c].level << endl;
  }

  return;
}



//=============================================================================
//  SphTree::BuildMpiGhostTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::BuildMpiGhostTree
(bool rebuild_tree,                 ///< Flag to rebuild tree
 int n,                             ///< Integer time
 int ntreebuildstep,                ///< Tree build frequency
 int ntreestockstep,                ///< Tree stocking frequency
 int Npart,                         ///< No. of particles
 int Npartmax,                      ///< Max. no. of particles
 SphParticle<ndim> *sph_gen,        ///< Particle data array
 Sph<ndim> *sph,                    ///< Pointer to SPH object
 FLOAT timestep)                    ///< Smallest physical timestep
{
  int i;                            // Particle counter
  int k;                            // Dimension counter
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph_gen);

  // If no MPI ghosts exist, do not build tree
  //if (sph->Nmpighost == 0) return;

  debug2("[SphTree::BuildGhostTree]");
  timing->StartTimingSection("BUILD_MPIGHOST_TREE",2);

  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif

  cout << "BUILDING TREE WITH " << sph->Nmpighost << " MPI GHOSTS!!" << endl;

  // For tree rebuild steps
  //---------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    mpighosttree->ifirst = sph->Nsph + sph->NPeriodicGhost;
    mpighosttree->ilast = sph->Nsph + sph->NPeriodicGhost +sph->Nmpighost - 1;
    mpighosttree->Ntot = sph->Nmpighost;
    mpighosttree->Ntotmaxold = mpighosttree->Ntotmax;
    mpighosttree->Ntotmax = max(mpighosttree->Ntotmax,mpighosttree->Ntot);
    mpighosttree->Ntotmax = max(mpighosttree->Ntotmax,sph->Nsphmax);
    mpighosttree->BuildTree(mpighosttree->Ntot,mpighosttree->Ntotmax,
                            sphdata,timestep);

  }

  // Else stock the tree
  //---------------------------------------------------------------------------
  else if (n%ntreestockstep == 0) {

    mpighosttree->StockTree(mpighosttree->kdcell[0],sphdata);

  }

  // Otherwise simply extrapolate tree cell properties
  //---------------------------------------------------------------------------
  else {

    mpighosttree->ExtrapolateCellProperties(timestep);

  }
  //---------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  timing->EndTimingSection("BUILD_MPIGHOST_TREE");


  return;
}



//=============================================================================
//  SphTree::SearchMpiGhostParticles
/// ...
//=============================================================================
template <int ndim, template<int> class ParticleType>
int SphTree<ndim,ParticleType>::SearchMpiGhostParticles
(const FLOAT tghost,                ///< [in] Expected ghost life-time
 const Box<ndim> &mpibox,           ///< [in] Bounding box of MPI domain
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 vector<int> &export_list)          ///< [out] List of particle ids
{
  int c;                            // Cell counter
  int i;
  int k;
  int Nexport = 0;                  // No. of MPI ghosts to export
  FLOAT scattermin[ndim];
  FLOAT scattermax[ndim];
  const FLOAT grange = ghost_range*kernrange;
  KDTreeCell<ndim> *cell;


  // Start from root-cell
  c = 0;

  //---------------------------------------------------------------------------
  while (c < tree->Ncell) {
    cell = &(tree->kdcell[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] =
        cell->bbmin[k] + min(0.0,cell->v[k]*tghost) - grange*cell->hmax;
      scattermax[k] =
        cell->bbmax[k] + max(0.0,cell->v[k]*tghost) + grange*cell->hmax;
    }


    // If maximum cell scatter box overlaps MPI domain, open cell
    //-------------------------------------------------------------------------
    if (tree->BoxOverlap(scattermin,scattermax,mpibox.boxmin,mpibox.boxmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (cell->level != tree->ltot)
        c++;

      else if (cell->N == 0)
        c = cell->cnext;

      // If leaf-cell, check through particles in turn to find ghosts and
      // add to list to be exported
      else if (cell->level == tree->ltot) {
        i = cell->ifirst;
        while (i != -1) {
          export_list.push_back(i);
          Nexport++;
          if (i == cell->ilast) break;
          i = tree->inext[i];
        };
        c = cell->cnext;
      }
    }

    // If not in range, then open next cell
    //-------------------------------------------------------------------------
    else
      c = cell->cnext;

  }
  //---------------------------------------------------------------------------



  // Start from root-cell of ghost-tree
  c = 0;

  //---------------------------------------------------------------------------
  while (c < ghosttree->Ncell) {
    cell = &(ghosttree->kdcell[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] =
        cell->bbmin[k] + min(0.0,cell->v[k]*tghost) - grange*cell->hmax;
      scattermax[k] =
        cell->bbmax[k] + max(0.0,cell->v[k]*tghost) + grange*cell->hmax;
    }

    // If maximum cell scatter box overlaps MPI domain, open cell
    //-------------------------------------------------------------------------
    if (ghosttree->BoxOverlap(scattermin,scattermax,
        mpibox.boxmin,mpibox.boxmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (cell->level != ghosttree->ltot)
        c++;

      else if (cell->N == 0)
        c = cell->cnext;

      // If leaf-cell, check through particles in turn to find ghosts and
      // add to list to be exported
      else if (cell->level == ghosttree->ltot) {
        i = cell->ifirst;
        while (i != -1) {
          export_list.push_back(i);
          Nexport++;
          if (i == cell->ilast) break;
          i = ghosttree->inext[i];
        };
        c = cell->cnext;
      }
    }

    // If not in range, then open next cell
    //-------------------------------------------------------------------------
    else
      c = cell->cnext;

  }
  //---------------------------------------------------------------------------


  return Nexport;
}



//=============================================================================
//  SphTree::SearchHydroExportParticles
/// ...
//=============================================================================
template <int ndim, template<int> class ParticleType>
int SphTree<ndim,ParticleType>::SearchHydroExportParticles
(const Box<ndim> &mpibox,           ///< [in] Bounding box of MPI domain
 Sph<ndim> *sph,                    ///< [in] Pointer to SPH object
 vector<KDTreeCell<ndim> *> &cell_list)          ///< [out] List of particle ids
{
  int c;                            // Cell counter
  int i;
  int k;
  int Nexport = 0;                  // No. of MPI ghosts to export
  FLOAT scattermin[ndim];
  FLOAT scattermax[ndim];
  const FLOAT grange = ghost_range*kernrange;
  KDTreeCell<ndim> *cell;
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray());


  // Start from root-cell
  c = 0;

  //---------------------------------------------------------------------------
  while (c < tree->Ncell) {
    cell = &(tree->kdcell[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = cell->bbmin[k] - grange*cell->hmax;
      scattermax[k] = cell->bbmax[k] + grange*cell->hmax;
    }


    // If maximum cell scatter box overlaps MPI domain, open cell
    //-------------------------------------------------------------------------
    if (tree->BoxOverlap(scattermin,scattermax,mpibox.boxmin,mpibox.boxmax)) {

      // If not a leaf-cell, then open cell to first child cell
      if (cell->level != tree->ltot)
        c++;

      else if (cell->N == 0)
        c = cell->cnext;

      // If leaf-cell and active, add the cell to the list of cells being exported
      else if (cell->level == tree->ltot) {
        if (cell->Nactive>0) {
          Nexport += cell->Nactive;
          cell_list.push_back(cell);
        }
        c = cell->cnext;
      }
    }

    // If not in range, then open next cell
    //-----------------------------------------------------------------------
    else
      c = cell->cnext;

  }
  //-------------------------------------------------------------------------

  return Nexport;

}



//=============================================================================
//  SphTree::FindMpiTransferParticles
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::FindMpiTransferParticles
(Sph<ndim>* sph,                            ///< [in] Pointer to sph class
 vector<vector<int> >& particles_to_export, ///< [inout] Vector that for each
                                            ///< node gives the list of particles to export
 vector<int>& all_particles_to_export,      ///< [inout] Vector containing all the particles that will be exported by this processor
 const vector<int>& potential_nodes,        ///< [in] Vector containing the potential nodes we might be sending particles to
 MpiNode<ndim>* mpinodes)                   ///< [in] Array of other mpi nodes
{
  int c;
  int i;
  int inode;
  int node_number;
  KDTreeCell<ndim> *cell;
  ParticleType<ndim> *sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray());


  // Loop over potential domains and walk the tree for each bounding box
  //-------------------------------------------------------------------------
  for (inode=0; inode<potential_nodes.size(); inode++) {
    node_number = potential_nodes[inode];

    Box<ndim>& nodebox = mpinodes[node_number].domain;


    // Start from root-cell
    c = 0;

    //-------------------------------------------------------------------------
    while (c < Ncell) {
      cell = &(tree->kdcell[c]);

      // If maximum cell scatter box overlaps MPI domain, open cell
      //-----------------------------------------------------------------------
      if (tree->BoxOverlap(cell->bbmin,cell->bbmax,
                           nodebox.boxmin,nodebox.boxmax)) {

	// If not a leaf-cell, then open cell to first child cell
	if (cell->level != tree->ltot)
	  c++;

	else if (cell->N == 0)
	  c = cell->cnext;

	// If leaf-cell, check through particles in turn to find ghosts and
	// add to list to be exported
	else if (cell->level == tree->ltot) {
	  i = cell->ifirst;
	  while (i != -1) {
	    if (ParticleInBox(sphdata[i],mpinodes[node_number].domain)) {
	      particles_to_export[node_number].push_back(i);
	      all_particles_to_export.push_back(i);
	    }
	    if (i == cell->ilast) break;
	    i = tree->inext[i];
	  };
	  c = cell->cnext;
	}
      }

      // If not in range, then open next cell
      //-----------------------------------------------------------------------
      else
	c = cell->cnext;

    }
    //-------------------------------------------------------------------------

  }
  //---------------------------------------------------------------------------

  return;
}



//=============================================================================
//  SphTree::GetExportInfo
/// Get the array with the information that needs to be exported to the given
/// processor (note: Nproc is ignored at the moment, as we always need to export
/// all particles to the other processors)
//=============================================================================
template <int ndim, template<int> class ParticleType>
int SphTree<ndim,ParticleType>::GetExportInfo (
    int Nproc,        ///< [in] Number of processor we want to send the information to
    Sph<ndim>*  sph,  ///< [in] Pointer to sph object
    vector<char >& send_buffer,   ///< [inout] Vector where the particles to export will be stored
    MpiNode<ndim>& mpinode,
    int rank,
    int Nmpi)       ///< [in] Array with information for the other mpi nodes
    {

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray() );

  assert(tree->Nimportedcell==0);

  const bool first_proc = (Nproc==0) || (rank==0 && Nproc==1);
  const bool hydro_only = !sph->self_gravity && sph->hydro_forces;

  int Nactive=0, cactive;

  // Get active cells and their number (so that we know how much memory to allocate)
  vector<KDTreeCell<ndim>*> celllist;
  celllist.reserve(tree->gtot);
  if (hydro_only) {
    Nactive += SearchHydroExportParticles(mpinode.domain,sph,celllist);
    cactive = celllist.size();
  }
  else {

    for (int i=0; i<sph->Nsph; i++) {
      if (sphdata[i].active)
        Nactive++;
    }

    celllist.resize(tree->gtot);
    cactive = tree->ComputeActiveCellList(&celllist[0]);
  }

  // Work out size of the information we are sending
  //Header consists of number of particles and number of cells
  const int size_header = 2*sizeof(int);
  const int size_particles = Nactive*sizeof(ParticleType<ndim>);
  const int size_cells = cactive*sizeof(KDTreeCell<ndim>);
  const int old_size = send_buffer.size();
  send_buffer.resize(size_particles+size_cells+size_header+old_size);

  //Write the header
  copy(&send_buffer[old_size],&Nactive);
  copy(&send_buffer[old_size+sizeof(int)],&cactive);

  //Clear the array needed for bookkeeping (which active particles we sent to which processor)
  vector<int>& ids_active_particles = ids_sent_particles[Nproc];
  ids_active_particles.clear(); ids_active_particles.reserve(Nactive);

  //Copy cells to export inside array
  int offset = size_header+old_size;
  int activelist[Nleafmax];
  int exported_particles = 0;
  for (int i=0; i<cactive; i++) {
    copy(&send_buffer[offset],celllist[i]);
    const int Nactive_cell = tree->ComputeActiveParticleList(celllist[i],sphdata,activelist);
    // Update the ifirst and ilast pointers in the cell
    KDTreeCell<ndim>* exported_cell = reinterpret_cast<KDTreeCell<ndim>*> (&send_buffer[offset]);
    exported_cell->ifirst = exported_particles;
    exported_cell->ilast = exported_particles+Nactive_cell-1;
    offset += sizeof(KDTreeCell<ndim>);
    // Copy active particles
    for (int iparticle=0; iparticle<Nactive_cell; iparticle++) {
      ids_active_particles.push_back(activelist[iparticle]);
      copy(&send_buffer[offset],&sphdata[activelist[iparticle]]);
      offset += sizeof(ParticleType<ndim>);
    }
    exported_particles += Nactive_cell;
  }
  assert(exported_particles == Nactive);
  assert(offset == send_buffer.size());

  return size_particles+size_cells+size_header;

}


//=============================================================================
//  SphTree::UnpackExported
/// Unpack the information exported from the other processors, contaning the particles
/// that were exported and
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UnpackExported (
    vector<char >& received_array,
    vector<int>& Nbytes_exported_from_proc,
    Sph<ndim>* sph) {


  int offset = 0;


  assert(sph->NImportedParticles==0);
  tree->Nimportedcell=0;
  tree->Ncelltot=tree->Ncell;

  N_imported_part_per_proc.resize(Nbytes_exported_from_proc.size());

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray() );

  for (int Nproc = 0; Nproc<Nbytes_exported_from_proc.size(); Nproc++) {

    int N_received_bytes = Nbytes_exported_from_proc[Nproc];
    int N_received_particles; int N_received_cells;

    if (N_received_bytes == 0) {
      N_imported_part_per_proc[Nproc]=0;
      continue;
    }

    copy(&N_received_particles,&received_array[offset]);
    N_imported_part_per_proc[Nproc]=N_received_particles;
    copy(&N_received_cells,&received_array[offset+sizeof(int)]);

    //Ensure there is enough memory
    if (sph->Ntot + N_received_particles > sph->Nsphmax) {
      ExceptionHandler::getIstance().raise("Error while receiving imported particles: not enough memory!");
    }
    if (tree->Ncelltot + N_received_cells > tree->Ncellmax) {
      ExceptionHandler::getIstance().raise("Error while receiving imported cells: not enough memory!");
    }


    //Copy received particles inside SPH main arrays and received cells inside the tree array
    //Also update the linked list
    int particle_index = sph->Ntot;
    offset += 2*sizeof(int);
    for (int icell=0; icell<N_received_cells; icell++) {
      KDTreeCell<ndim>& dest_cell = tree->kdcell[icell+tree->Ncelltot];
      copy(&dest_cell,&received_array[offset]);
      offset += sizeof(KDTreeCell<ndim>);
      // Offset the ifirst and ilast pointers
      dest_cell.ifirst += sph->Ntot;
      dest_cell.ilast += sph->Ntot;
      // Now copy the received particles inside the SPH main arrays
      for (int iparticle=0;iparticle<dest_cell.Nactive;iparticle++) {
        copy( &sphdata[particle_index] , &received_array[offset]);
        tree->inext[particle_index] = particle_index + 1;
        particle_index++; offset+=sizeof(ParticleType<ndim>);
      }

    }

    //Update the SPH counters
    sph->Ntot += N_received_particles;
    sph->NImportedParticles += N_received_particles;

    //Update the tree counters
    tree->Nimportedcell += N_received_cells;
    tree->Ncelltot      += N_received_cells;
    tree->Ntot          = sph->Ntot;

  }

  assert (offset == std::accumulate(Nbytes_exported_from_proc.begin(), Nbytes_exported_from_proc.end(),0));


}


//=============================================================================
//  SphTree::GetBackExportInfo
/// Return the data to transmit back to the other processors (particle acceleration etc.)
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::GetBackExportInfo (
    vector<char >& send_buffer, ///< [inout] These arrays will be overwritten with the information to send
    vector<int>& Nbytes_exported_from_proc,
    vector<int>& Nbytes_to_each_proc,
    Sph<ndim>* sph,   ///< [in] Pointer to the SPH object
    int rank
    ) {

  int InitialNImportedParticles = sph->NImportedParticles;

  //loop over the processors, removing particles as we go
  int removed_particles=0;
  send_buffer.resize(sph->NImportedParticles * sizeof(ParticleType<ndim>));
  for (int Nproc=0 ; Nproc < N_imported_part_per_proc.size(); Nproc++ ) {

    const int N_received_particles = N_imported_part_per_proc[Nproc];

//    //Copy the accelerations and gravitational potential of the particles
//    ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray() );
//    int j=0;
//    for (int i=sph->Nsph - N_received_particles; i<sph->Nsph; i++) {
//      for (int k=0; k<ndim; k++)
//        send_buffer[removed_particles+j].a[k] = sphdata[i].a[k];
//      send_buffer[removed_particles+j].gpot = sphdata[i].gpot;
//      j++;
//    }

    //Copy the particles inside the send buffer
    ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray() );
    int j=0;
    const int start_index = sph->Nsph + sph->Nghost + removed_particles;
    for (int i=start_index; i<start_index + N_received_particles; i++) {
      copy (&send_buffer[(removed_particles+j)*sizeof(ParticleType<ndim>)],&sphdata[i]);
      j++;
    }

    assert(j==N_received_particles);

    removed_particles += j;

    //Decrease the particle counter
    sph->Ntot -= N_received_particles;
    sph->NImportedParticles -= N_received_particles;

    //Update the information about how much data we are sending
    Nbytes_exported_from_proc[Nproc] = N_received_particles*sizeof(ParticleType<ndim>);

    //Update the information with how much data we are receiving
    Nbytes_to_each_proc[Nproc] = ids_sent_particles[Nproc].size()*sizeof(ParticleType<ndim>);

  }

  tree->Ncelltot = tree->Ncell;
  tree->Nimportedcell=0;
  tree->Ntot          -= InitialNImportedParticles;

  assert(sph->NImportedParticles == 0);
  assert(sph->Ntot == sph->Nsph + sph->Nghost);
  assert(send_buffer.size() == removed_particles*sizeof(ParticleType<ndim>));

}



//=============================================================================
//  SphTree::UnpackReturnedExportInfo
/// Unpack the data that was returned by the other processors, summing the accelerations to the particles
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::UnpackReturnedExportInfo (
    vector<char >& received_information,
    vector<int>& recv_displs,
    Sph<ndim>* sph,
    int rank
    ) {

  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray() );

  //For each processor, sum up the received quantities for each particle
  for(int Nproc=0; Nproc<recv_displs.size(); Nproc++ ) {

    if (rank==Nproc)
      continue;

    const vector<int>& ids_active_particles = ids_sent_particles[Nproc];

    for (int i=0; i< ids_active_particles.size(); i++ ) {
    const int j = ids_active_particles[i];

      ParticleType<ndim>* received_particle = reinterpret_cast<ParticleType<ndim>*>
      (&received_information[i * sizeof(ParticleType<ndim>) + recv_displs[Nproc] ]);

      assert(sphdata[j].iorig == received_particle->iorig);

      for (int k=0; k<ndim; k++) {
        sphdata[j].a[k] += received_particle->a[k];
        sphdata[j].agrav[k] += received_particle->agrav[k];
      }
      sphdata[j].gpot += received_particle->gpot;
      sphdata[j].gpe += received_particle->gpe;
      sphdata[j].dudt += received_particle->dudt;
      sphdata[j].div_v += received_particle->div_v;
      sphdata[j].levelneib = max(sphdata[j].levelneib, received_particle->levelneib);

    }

  }

}



//=============================================================================
//  SphTree::CommunicatePrunedTrees
/// ..
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::CommunicatePrunedTrees (
    vector<int>& my_matches,
    int rank ) {

  //---------------------------------------------------------------------------
  for (int iturn=0; iturn<my_matches.size(); iturn++) {

    //See who we need to communicate with
    int inode = my_matches[iturn];
    cout << "inode : " << rank << "   " << inode << "    " << iturn << endl;

    //Decide if sending or receiving first
    bool send_turn;
    if (rank<inode) {
      send_turn = true;
    }
    else {
      send_turn = false;
    }

    //Do the actual communication
    for (int i=0; i<2; i++) {
      if (send_turn) {
        KDTree<ndim, ParticleType>* tree = prunedtree[rank];
        MPI_Send(tree->kdcell,tree->Ncell*sizeof(KDTreeCell<ndim>),
                 MPI_CHAR,inode,3,MPI_COMM_WORLD);
        send_turn=false;
      }
      else {
        KDTree<ndim, ParticleType>* tree = prunedtree[inode];
        MPI_Status status;
        MPI_Recv(tree->kdcell,tree->Ncell*sizeof(KDTreeCell<ndim>),
                 MPI_CHAR,inode,3,MPI_COMM_WORLD,&status);
        send_turn=true;
      }
    }

  }
  //---------------------------------------------------------------------------

  MPI_Barrier(MPI_COMM_WORLD);
  for (int i=0; i<Nmpi; i++) {
    cout << "Writing pruned tree " << i << " for process " << rank << endl;
    cout << "Ncell : " << prunedtree[i]->Ncell << endl;
    cout << "r : " << prunedtree[i]->kdcell[0].r[0]
         << "   box : " << prunedtree[i]->kdcell[0].bbmin[0]
         << "   " << prunedtree[i]->kdcell[0].bbmax[0] << endl;
  }

}


#endif



#if defined(VERIFY_ALL)
//=============================================================================
//  SphNeighbourSearch::CheckValidNeighbourList
/// Checks that the neighbour list generated by the grid is valid in that it
/// (i) does include all true neighbours, and
/// (ii) all true neigbours are only included once and once only.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void SphTree<ndim,ParticleType>::CheckValidNeighbourList
(int i,                             ///< [in] Particle i.d.
 int Ntot,                          ///< [in] Total no. of particles
 int Nneib,                         ///< [in] No. of potential neighbours
 int *neiblist,                     ///< [in] List of potential neighbour i.d.s
 ParticleType<ndim> *partdata,      ///< [in] Array of particle data
 string neibtype)                   ///< [in] Neighbour search type
{
  int count;                        // Valid neighbour counter
  int j;                            // Neighbour particle counter
  int k;                            // Dimension counter
  int Ntrueneib = 0;                // No. of 'true' neighbours
  int *trueneiblist;                // List of true neighbour ids
  FLOAT drsqd;                      // Distance squared
  FLOAT dr[ndim];                   // Relative position vector

  // Allocate array to store local copy of potential neighbour ids
  trueneiblist = new int[Ntot];

  // First, create list of 'true' neighbours by looping over all particles
  if (neibtype == "gather") {
    for (j=0; j<Ntot; j++) {
      for (k=0; k<ndim; k++)
	dr[k] = partdata[j].r[k] - partdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd <= kernrangesqd*partdata[i].h*partdata[i].h)
	trueneiblist[Ntrueneib++] = j;
    }
  }
  else if (neibtype == "all") {
    for (j=0; j<Ntot; j++) {
      for (k=0; k<ndim; k++)
        dr[k] = partdata[j].r[k] - partdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (drsqd < kernrangesqd*partdata[i].h*partdata[i].h ||
    	  drsqd < kernrangesqd*partdata[j].h*partdata[j].h)
 	     trueneiblist[Ntrueneib++] = j;
     }
   }

  // Now compare each given neighbour with true neighbour list for validation
  for (j=0; j<Ntrueneib; j++) {
    count = 0;
    for (k=0; k<Nneib; k++)
      if (neiblist[k] == trueneiblist[j]) count++;

    // If the true neighbour is not in the list, or included multiple times,
    // then output to screen and terminate program
    if (count != 1) {
      InsertionSortIds(Nneib,neiblist);
      cout << "Problem with neighbour lists : " << i << "  " << j << "   "
	   << count << "   "
	   << partdata[i].r[0] << "   " << partdata[i].h << endl;
      cout << "Nneib : " << Nneib << "   Ntrueneib : " << Ntrueneib << endl;
      PrintArray("neiblist     : ",Nneib,neiblist);
      PrintArray("trueneiblist : ",Ntrueneib,trueneiblist);
      string message = "Problem with neighbour lists in grid search";
      ExceptionHandler::getIstance().raise(message);
    }
  }

  delete[] trueneiblist;

  return;
}
#endif



template class SphTree<1,GradhSphParticle>;
template class SphTree<2,GradhSphParticle>;
template class SphTree<3,GradhSphParticle>;
template class SphTree<1,SM2012SphParticle>;
template class SphTree<2,SM2012SphParticle>;
template class SphTree<3,SM2012SphParticle>;
template class SphTree<1,GodunovSphParticle>;
template class SphTree<2,GodunovSphParticle>;
template class SphTree<3,GodunovSphParticle>;
