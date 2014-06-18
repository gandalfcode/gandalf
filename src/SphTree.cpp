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
/// Allocate memory for binary tree as requested.  If more memory is required 
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



#ifdef MPI_PARALLEL
//=============================================================================
//  SphTree::BuildGhostTree
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
#endif



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



#if defined MPI_PARALLEL
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

  //cout << "SEARCHING FOR MPI GHOSTS : " << mpibox.boxmin[0] << "   " << mpibox.boxmax[0] << "    " << Ncell << endl;

  // Start from root-cell
  c = 0;
  
  //-------------------------------------------------------------------------
  while (c < tree->Ncell) {
    cell = &(tree->kdcell[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = 
	cell->bbmin[k] + min(0.0,cell->v[k]*tghost) - grange*cell->hmax;
      scattermax[k] = 
	cell->bbmax[k] + max(0.0,cell->v[k]*tghost) + grange*cell->hmax;
    }

    //cout << "Checking MPI boxes : " << Nexport << "   " << c << "   " 
    // << scattermin[0] << "   " << scattermax[0] << "   "
    // << tree->BoxOverlap(scattermin,scattermax,mpibox.boxmin,mpibox.boxmax)
    // << endl;

    
    // If maximum cell scatter box overlaps MPI domain, open cell
    //-----------------------------------------------------------------------
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
	  //cout << "Found new MPI ghost : " << i << "   " << Nexport << endl;
	  export_list.push_back(i);
	  Nexport++;
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



  // Start from root-cell of ghost-tree
  c = 0;
  
  //-------------------------------------------------------------------------
  while (c < ghosttree->Ncell) {
    cell = &(ghosttree->kdcell[c]);

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = 
	cell->bbmin[k] + min(0.0,cell->v[k]*tghost) - grange*cell->hmax;
      scattermax[k] = 
	cell->bbmax[k] + max(0.0,cell->v[k]*tghost) + grange*cell->hmax;
    }

    //cout << "Checking MPI boxes : " << Nexport << "   " << c << "   " 
    // << scattermin[0] << "   " << scattermax[0] << "   "
    // << ghosttree->BoxOverlap(scattermin,scattermax,mpibox.boxmin,mpibox.boxmax)
    // << endl;

    
    // If maximum cell scatter box overlaps MPI domain, open cell
    //-----------------------------------------------------------------------
    if (ghosttree->BoxOverlap(scattermin,scattermax,mpibox.boxmin,mpibox.boxmax)) {

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
	  //cout << "Found new MPI ghost : " << i << "   " << Nexport << endl;
	  export_list.push_back(i);
	  Nexport++;
	  if (i == cell->ilast) break;
	  i = ghosttree->inext[i];
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
            
      // If leaf-cell, check through particles in turn to find ghosts and
      // add to list to be exported
      else if (cell->level == tree->ltot) {
        i = cell->ifirst;
        while (i != -1) {
          if (sphdata[i].active) {
            export_list.push_back(i);
            Nexport++;
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
