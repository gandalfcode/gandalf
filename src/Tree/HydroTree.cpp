//=================================================================================================
//  HydroTree.cpp
//  Contains all functions for managing the tree for hydrodynamical particles.
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
#include <numeric>
#include <string>
#include <math.h>
#include <vector>
#include "Precision.h"
#include "Exception.h"
#include "Hydrodynamics.h"
#include "NeighbourSearch.h"
#include "Parameters.h"
#include "InlineFuncs.h"
#include "Particle.h"
#include "Debug.h"
#if defined _OPENMP
#include <omp.h>
#endif
#if defined MPI_PARALLEL
#include "CommunicationHandler.h"
#endif
using namespace std;




//=================================================================================================
//  HydroTree::HydroTree
/// HydroTree constructor.  Initialises various variables.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
HydroTree<ndim,ParticleType>::HydroTree
 (string tree_type,
  int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max, FLOAT _thetamaxsqd,
  FLOAT _kernrange, FLOAT _macerror, string _gravity_mac, multipole_method _multipole,
  DomainBox<ndim>* _box, SmoothingKernel<ndim>* _kern, CodeTiming* _timing,
  ParticleTypeRegister& types):
  neibcheck(true),
  kernrange(_kernrange),
  kernrangesqd(_kernrange*_kernrange),
  Nleafmax(_Nleafmax),
  Nmpi(_Nmpi),
  pruning_level_min(_pruning_level_min),
  pruning_level_max(_pruning_level_max),
  thetamaxsqd(_thetamaxsqd),
  invthetamaxsqd((FLOAT) 1.0/_thetamaxsqd),
  macerror(_macerror),
  gravity_mac(_gravity_mac),
  multipole(_multipole),
  timing(_timing),
  box(_box),
  kernp(_kern)
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
#ifdef MPI_PARALLEL
  Npartexport = new int[Nmpi];
  cellexportlist.resize(Nmpi);
  ids_sent_particles.resize(Nmpi);
  ids_sent_cells.resize(Nmpi);
  N_imported_part_per_proc.resize(Nmpi);
  N_imported_cells_per_proc.resize(Nmpi);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif


  // Set-up main tree object
  tree = CreateTree(tree_type, _Nleafmax, _thetamaxsqd, _kernrange,
                    _macerror, _gravity_mac, _multipole, *_box, types,false);

  // Set-up ghost-particle tree object
  ghosttree = CreateTree(tree_type, _Nleafmax, _thetamaxsqd, _kernrange,
      _macerror, _gravity_mac, _multipole, *_box, types,false);

#ifdef MPI_PARALLEL
  // Set-up ghost-particle tree object
  mpighosttree = CreateTree(tree_type, _Nleafmax, _thetamaxsqd, _kernrange,
                            _macerror, _gravity_mac, _multipole, *_box, types,false);

  // Set-up multiple pruned trees, one for each MPI process
  prunedtree = new TreeBase<ndim>*[Nmpi] ;
  sendprunedtree =  new TreeBase<ndim>*[Nmpi] ;

  for (int i=0; i<Nmpi; i++) {
    prunedtree[i] = CreateTree(tree_type, _Nleafmax, _thetamaxsqd, _kernrange,
                               _macerror, _gravity_mac, _multipole, *_box, types,true);
    sendprunedtree[i] = CreateTree(tree_type, _Nleafmax, _thetamaxsqd, _kernrange,
                                   _macerror, _gravity_mac, _multipole, *_box, types,true);
  }
#endif
}



//=================================================================================================
//  HydroTree::~HydroTree
/// HydroTree destructor.  Deallocates tree memory upon object destruction.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
HydroTree<ndim,ParticleType>::~HydroTree()
{
  // Free up the trees that we own.
  delete tree ;
  delete ghosttree ;
#ifdef MPI_PARALLEL
  for (int i=0; i<Nmpi; i++) {
    delete prunedtree[i];
    delete sendprunedtree[i];
  }
  delete[] prunedtree ;
  delete[] sendprunedtree;
#endif
}

//=================================================================================================
//  HydroTree::CreateTree
/// Construct the required tree type based on parameters
//=================================================================================================
template <int ndim, template <int> class ParticleType>
TreeBase<ndim>* HydroTree<ndim,ParticleType>::CreateTree
(string tree_type,
 int Nleafmax, FLOAT thetamaxsqd,
 FLOAT kernrange, FLOAT macerror,
 string gravity_mac, multipole_method multipole,
 const DomainBox<ndim>& domain,
 const ParticleTypeRegister& reg,
 const bool IAmPruned
)
{
  TreeBase<ndim> * t = NULL;
  if (tree_type == "bruteforce") {
    typedef BruteForceTree<ndim,ParticleType, BruteForceTreeCell> __tree ;

    t = new __tree(Nleafmax, thetamaxsqd, kernrange, macerror,
                   gravity_mac,  multipole,domain, reg, IAmPruned);
  }
  else if (tree_type == "kdtree") {
    typedef KDTree<ndim,ParticleType, KDTreeCell>  __tree ;

    t = new __tree(Nleafmax, thetamaxsqd, kernrange, macerror,
                   gravity_mac,  multipole,domain, reg, IAmPruned);
  }
  else if (tree_type == "octtree") {
    typedef OctTree<ndim,ParticleType, OctTreeCell>  __tree ;

    t = new __tree(Nleafmax, thetamaxsqd, kernrange, macerror,
                   gravity_mac,  multipole,domain, reg, IAmPruned);
  }
  else if (tree_type == "treeray") {
    typedef OctTree<ndim,ParticleType, TreeRayCell> __tree ;

    t = new __tree(Nleafmax, thetamaxsqd, kernrange, macerror,
                   gravity_mac,  multipole,domain, reg, IAmPruned);
  }
  else {
    string message = "Tree Type not recognised: " + tree_type ;
    ExceptionHandler::getIstance().raise(message);
  }

  return t ;
}


template <int ndim, template <int> class ParticleType>
MAC_Type HydroTree<ndim,ParticleType>::GetOpeningCriterion() const {
  return tree->GetMacType() ;
}

template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::SetOpeningCriterion(MAC_Type value) {
  tree->SetMacType(value) ;
  ghosttree->SetMacType(value) ;
#ifdef MPI_PARALLEL
  for (int i=0; i<Nmpi; i++) {
    prunedtree[i]->SetMacType(value) ;
    sendprunedtree[i]->SetMacType(value) ;
  }
#endif
}



//=================================================================================================
//  HydroTree::AllocateMemory
/// Allocate memory for tree as requested.  If more memory is required
/// than currently allocated, tree is deallocated and reallocated here.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::AllocateMemory
 (const int Ngather)                   ///< [in] Average no. of gather neighbours
{
  debug2("[HydroTree::AllocateMemory]");

  if (!allocated_buffer) {

    Nneibmaxbuf     = new int[Nthreads];
    Ngravcellmaxbuf = new int[Nthreads];
    activelistbuf   = new int*[Nthreads];
    levelneibbuf    = new int*[Nthreads];
    neiblistbuf     = new int*[Nthreads];
    ptypebuf        = new int*[Nthreads];
    activepartbuf   = new ParticleType<ndim>*[Nthreads];

    for (int ithread=0; ithread<Nthreads; ithread++) {
      Nneibmaxbuf[ithread]     = max(1, 8*Ngather);
      Ngravcellmaxbuf[ithread] = max(1, 16*Ngather);
      activelistbuf[ithread]   = new int[Nleafmax];
      levelneibbuf[ithread]    = new int[Ntotmax];
      neiblistbuf[ithread]     = new int[Nneibmaxbuf[ithread]];
      ptypebuf[ithread]        = new int[Nneibmaxbuf[ithread]];
      activepartbuf[ithread]   = new ParticleType<ndim>[Nleafmax];
    }
    allocated_buffer = true;

  }

  return;
}


//=================================================================================================
//  HydroTree::ReallocateMemory
/// Reallocate memory for tree when the number of particles has changed
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::ReallocateMemory
 ()
{
    for (int ithread=0; ithread<Nthreads; ithread++) {

    	delete[] levelneibbuf[ithread];
    	levelneibbuf[ithread] = new int[Ntotmax];
    	assert (Ntot <= Ntotmax);
    }
}



//=================================================================================================
//  HydroTree::DeallocateTreeMemory
/// Deallocates all binary tree memory
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::DeallocateMemory(void)
{
  int ithread;                         // Thread id number

  debug2("[HydroTree::DeallocateTreeMemory]");

  if (allocated_buffer) {

    for (ithread=0; ithread<Nthreads; ithread++) {
      delete[] levelneibbuf[ithread];
      delete[] activepartbuf[ithread];
      delete[] ptypebuf[ithread];
      delete[] neiblistbuf[ithread];
      delete[] levelneibbuf[ithread];
      delete[] activelistbuf[ithread];
    }
    delete[] levelneibbuf;
    delete[] activepartbuf;
    delete[] activelistbuf;
    delete[] Ngravcellmaxbuf;
    delete[] Nneibmaxbuf;
  }

  return;
}



//=================================================================================================
//  HydroTree::BuildTree
/// Main routine to control how the tree is built, re-stocked and interpolated during each step.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::BuildTree
 (const bool rebuild_tree,             ///< [in] Flag to rebuild tree
  const int n,                         ///< [in] Integer time
  const int ntreebuildstep,            ///< [in] Tree build frequency
  const int ntreestockstep,            ///< [in] Tree stocking frequency
  const FLOAT timestep,                ///< [in] Smallest physical timestep
  Hydrodynamics<ndim> *hydro)          ///< [inout] Pointer to Hydrodynamics object
{
  ParticleType<ndim> *partdata = hydro->template GetParticleArray<ParticleType>();
  CodeTiming::BlockTimer timer = timing->StartNewTimer("BUILD_TREE");

  debug2("[HydroTree::BuildTree]");


  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif


  // For tree rebuild steps
  //-----------------------------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    // Delete any dead particles from main Hydrodynamics arrays before we re-build tree
    hydro->DeleteDeadParticles();

    Ntotold    = Ntot;
    Ntot       = hydro->Ntot;
    Ntotmaxold = Ntotmax;
    Ntotmax    = max(Ntotmax, Ntot);
    Ntotmax    = max(Ntotmax, hydro->Nhydromax);
    tree->Ntot = hydro->Nhydro;
    tree->BuildTree(0, hydro->Nhydro-1, hydro->Nhydro, hydro->Nhydromax, timestep, partdata);

    AllocateMemory(hydro->Ngather);
    if (Ntotmaxold < Ntotmax) ReallocateMemory();

  }

  // Else stock the tree
  //-----------------------------------------------------------------------------------------------
  else {
    
    tree->StockTree(partdata,true);

  }
  //-----------------------------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  return;
}



//=================================================================================================
//  HydroTree::BuildGhostTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::BuildGhostTree
 (const bool rebuild_tree,             ///< [in] Flag to rebuild tree
  const int n,                         ///< [in] Integer time
  const int ntreebuildstep,            ///< [in] Tree build frequency
  const int ntreestockstep,            ///< [in] Tree stocking frequency
  const FLOAT timestep,                ///< [in] Smallest physical timestep
  Hydrodynamics<ndim> *hydro)          ///< [inout] Pointer to Hydrodynamics object
{
  ParticleType<ndim> *partdata = hydro->template GetParticleArray<ParticleType>();
  CodeTiming::BlockTimer timer = timing->StartNewTimer("BUILD_GHOST_TREE");

  debug2("[HydroTree::BuildGhostTree]");


  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif

  if (hydro->Ntot > Ntotmax) {
	  Ntotmax = hydro->Ntot;
	  Ntot = hydro->Ntot;
	  ReallocateMemory();
  }


  // For tree rebuild steps
  //-----------------------------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    ghosttree->Ntot       = hydro->NPeriodicGhost;
    int max_particles    = max(ghosttree->Ntot, hydro->Nhydromax);
    ghosttree->BuildTree(hydro->Nhydro, hydro->Nhydro + hydro->NPeriodicGhost - 1,
                         ghosttree->Ntot, max_particles, timestep, partdata);

  }

  // Else stock the tree
  //-----------------------------------------------------------------------------------------------
  else {

    ghosttree->StockTree(partdata,true);

  }
  //-----------------------------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif

  return;
}



//=================================================================================================
//  HydroTree::GetGatherNeighbourList
/// Compute the gather neighbour list at the point 'rp' of all particles within a search radius
/// of 'rsearch'.  Searches through real and ghost neighbour trees.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int HydroTree<ndim,ParticleType>::GetGatherNeighbourList
 (FLOAT rp[ndim],                      ///< [in] Position vector
  FLOAT rsearch,                       ///< [in] Gather search radius
  Particle<ndim> *part_gen,            ///< [in] Pointer to Hydrodynamics particle array
  int Nhydro,                          ///< [in] No. of hydro particles
  int Nneibmax,                        ///< [in] Max. no. of neighbours
  int *neiblist)                       ///< [out] List of neighbouring particles
{
  int Nneib = 0;                       // No. of (non-dead) neighbours
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[HydroTree::GetGatherNeighbourList]");

  Nneib = tree->ComputeGatherNeighbourList(partdata, rp, rsearch, Nneibmax, Nneib, neiblist);
  Nneib = ghosttree->ComputeGatherNeighbourList(partdata, rp, rsearch, Nneibmax, Nneib, neiblist);
#ifdef MPI_PARALLEL
  Nneib = mpighosttree->ComputeGatherNeighbourList(partdata, rp, rsearch, Nneibmax, Nneib, neiblist);
#endif

  return Nneib;
}



//=================================================================================================
//  HydroTree::UpdateActiveParticleCounters
/// Loop through all leaf cells in the tree and update all active particle counters.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::UpdateActiveParticleCounters
 (Hydrodynamics<ndim> *hydro)          ///< [in] Pointer to hydrodynamics object
{
  debug2("[HydroTree::UpdateActiveParticleCounters]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("HYDRO_TREE_UPDATE_ACTIVE_COUNTERS");

  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();
  tree->UpdateActiveParticleCounters(partdata);
}



//=================================================================================================
//  HydroTree::SearchBoundaryGhostParticles
/// Search domain to create any required ghost particles near any boundaries.
/// Currently only searches to create periodic or mirror ghost particles.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::SearchBoundaryGhostParticles
 (FLOAT tghost,                                ///< [in] Ghost particle 'lifetime'
  const DomainBox<ndim> &simbox,               ///< [in] Simulation box structure
  Hydrodynamics<ndim> *hydro)                  ///< [inout] Hydrodynamics object pointer
{
  int i;                                       // Particle counter
  const FLOAT grange = ghost_range*kernrange;  // Range of ghost particles (in terms of h)

  // Set all relevant particle counters
  hydro->Nghost         = 0;
  hydro->NPeriodicGhost = 0;
  hydro->Nmpighost      = 0;
  hydro->Ntot           = hydro->Nhydro;


  // If all boundaries are open, immediately return to main loop
  if (simbox.boundary_lhs[0] == openBoundary && simbox.boundary_rhs[0] == openBoundary &&
      simbox.boundary_lhs[1] == openBoundary && simbox.boundary_rhs[1] == openBoundary &&
      simbox.boundary_lhs[2] == openBoundary && simbox.boundary_rhs[2] == openBoundary) return;


  debug2("[HydroTree::SearchBoundaryGhostParticles]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("SEARCH_BOUNDARY_GHOSTS");


  // Loop over all dimensions and create ghost particles
  for (int k=0; k<ndim; k++) {
    if ((simbox.boundary_lhs[k] != openBoundary || simbox.boundary_rhs[k] != openBoundary)) {

      // Do the real particles using the tree
      tree->GenerateBoundaryGhostParticles(tghost, grange, k, simbox, hydro) ;

      // Include ghosts-of-ghosts by doing ghosts explicitly.
      if (k > 0) {
        for (i=hydro->Nhydro; i<hydro->Ntot; i++) {
          hydro->CheckBoundaryGhostParticle(i, k, tghost,simbox);
        }
      }

      hydro->Ntot = hydro->Nhydro + hydro->Nghost;
    }
  }
  hydro->NPeriodicGhost = hydro->Nghost;


  if (hydro->Ntot > Ntotmax) {
	  Ntotmax = hydro->Ntot;
	  ReallocateMemory();
  }

  return;
}


//=================================================================================================
//  HydroTree::UpdateAllStarGasForces
/// Calculate the gravitational acceleration on all star particles due to
/// all gas particles via the tree.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::UpdateAllStarGasForces
 (Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to SPH object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  DomainBox<ndim> &simbox,             ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald gravity object pointer
{
  int Nactive;                         // No. of active particles in cell
  int *activelist;                     // List of active particle ids
  NbodyParticle<ndim> *star;           // Pointer to star particle

  int Ntot = hydro->Ntot;
  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();


  debug2("[GradhSphTree::UpdateAllStarGasForces]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("STAR_GAS_GRAV_FORCES");

  // Make list of all active stars
  Nactive = 0;
  activelist = new int[nbody->Nstar];
  for (int i=0; i<nbody->Nstar; i++) {
    if (nbody->nbodydata[i]->flags.check(active)) activelist[Nactive++] = i;
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) private(star)\
  shared(activelist,ewald,hydro,Nactive,Ntot,nbody,partdata,simbox,cout)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int i;                                       // Particle id
    int j;                                       // Aux. particle counter
    int okflag;                                  // Flag if h-rho iteration is valid
    int Ndirect;                                 // No. of direct-sum gravity particles
    int Ngravcell;                               // No. of gravity cells
    int Nneib;                                   // No. of neighbours
    int Nneibmax = Ntot; //Nneibmaxbuf[ithread];
    int Ngravcellmax = Ngravcellmaxbuf[ithread]; // ..
    FLOAT macfactor;                             // Gravity MAC factor
    int* neiblist = new int[Nneibmax];           // ..
    int* directlist = new int[Nneibmax];         // ..
    MultipoleMoment<ndim>* gravcell = new MultipoleMoment<ndim>[Ngravcellmax];   // ..


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(dynamic)
    for (j=0; j<Nactive; j++) {
      i = activelist[j];
      star = nbody->nbodydata[i];

      // Compute average/maximum term for computing gravity MAC
      if (gravity_mac == "eigenmac") macfactor = pow((FLOAT) 1.0/star->gpot,twothirds);
      else macfactor = (FLOAT) 0.0;

      // Compute neighbour list for cell depending on physics options
      okflag = tree->ComputeStarGravityInteractionList
       (star, macfactor, Nneibmax, Nneibmax, Ngravcellmax, Nneib,
        Ndirect, Ngravcell, neiblist, directlist, gravcell, partdata);

      // If there are too many neighbours, reallocate the arrays and recompute the neighbour lists.
      while (okflag == -1) {
        delete[] gravcell;
        Ngravcellmax = 2*Ngravcellmax;
        gravcell = new MultipoleMoment<ndim>[Ngravcellmax];
        okflag = tree->ComputeStarGravityInteractionList
         (star, macfactor, Nneibmax, Nneibmax, Ngravcellmax, Nneib,
          Ndirect, Ngravcell, neiblist, directlist, gravcell, partdata);
      };

      // Compute contributions to star force from nearby hydro particles
      nbody->CalculateDirectHydroForces(star, Nneib, Ndirect, neiblist, directlist, hydro, simbox, ewald);

      // Compute gravitational force due to distant cells
      if (multipole == monopole || multipole == fast_monopole) {
        ComputeCellMonopoleForces(star->gpot, star->a, star->r, Ngravcell, gravcell);
      }
      else if (multipole == quadrupole || multipole == fast_quadrupole) {
        ComputeCellQuadrupoleForces(star->gpot, star->a, star->r, Ngravcell, gravcell);
      }


    }
    //=============================================================================================


    // Free-up local memory for OpenMP thread
    delete[] gravcell;
    delete[] directlist;
    delete[] neiblist;

  }
  //===============================================================================================

  delete[] activelist;

  return;
}
//=================================================================================================
// HydroTree::GetMaximumSmoothingLength
/// Returns the maximum smoothing length over all particles in the simulation
//=================================================================================================

template <int ndim, template <int> class ParticleType>
double HydroTree<ndim,ParticleType>::GetMaximumSmoothingLength() const
{
  assert(tree != NULL) ;
  double hmax = tree->GetMaximumSmoothingLength() ;

#if defined MPI_PARALLEL
  for (int n = 0; n < Nmpi; n++)
	hmax = std::max(hmax, prunedtree[n]->GetMaximumSmoothingLength()) ;
#endif
  return hmax ;
}

template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::UpdateTimestepsLimitsFromDistantParticles
(Hydrodynamics<ndim>* hydro,                     ///<[inout] Pointer to Hydrodynamics object
 const bool only_imported_particles)								///<[in] Wheter we need to loop only over imported particles (relevant only for MPI)
 {
  int cactive = 0;                      // No. of active cells
  vector<TreeCellBase<ndim> > celllist; // List of active tree cells
  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();

  debug2("[HydroTree::UpdateTimestepsLimitsFromDistantParticles]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("HYDRO_DISTANT_TIMESTEPS");

  // Find list of all cells that contain active particles
  if (only_imported_particles) {
#ifdef MPI_PARALLEL
	  cactive = tree->ComputeImportedCellList(celllist);
#else
	  string message = "This should not happen - bug in UpdateTimestepsLimitsFromDistantParticles!";
	  ExceptionHandler::getIstance().raise(message);
#endif
  }
  else {
	  cactive = tree->ComputeActiveCellList(celllist);
  }

#ifdef MPI_PARALLEL
    // Reset all export lists
    if (!only_imported_particles) {
		for (int j=0; j<Nmpi; j++) {
		  Npartexport[j] = 0;
		  cellexportlist[j].clear();
		  cellexportlist[j].reserve(tree->gmax);
		}
    }
#endif

#pragma omp parallel default(none) shared(celllist,cactive,cout,partdata)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                    // Aux. cell counter
    int Nactive;                               // No. of active particles in current cell
    int *activelist                = activelistbuf[ithread];
    ParticleType<ndim> *activepart = activepartbuf[ithread];



    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim>& cell = celllist[cc] ;

      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, partdata, activelist);
      if (Nactive == 0) continue ;

      // Make local copies of active particles & update vsig_max
      for (int j=0; j<Nactive; j++) activepart[j] = partdata[activelist[j]];

      tree->ComputeSignalVelocityFromDistantInteractions(cell, Nactive, activepart, partdata);

#ifdef MPI_PARALLEL
      if (!only_imported_particles) {
		  // Loop over pruned trees as well
		  // Doing it after the local tree walk is more efficient as vsig_max might have been already increased
		  for (int i=0; i<Nmpi; i++) {
			  if (i==rank) continue;

			  bool to_export=prunedtree[i]->ComputeSignalVelocityFromDistantInteractions(cell, Nactive, activepart, partdata);

			  // If the previous function returned true, we need to export this particle
			  if (to_export) {
	#pragma omp critical
			  {
				  cellexportlist[i].push_back(cell.id);
			  }

	#pragma omp atomic
			  Npartexport[i] += Nactive;
			  }

		  }
      }
#endif

      for (int j=0; j<Nactive; j++) {
        partdata[activelist[j]].vsig_max = activepart[j].vsig_max;
      }
    }
  }
 }

#ifdef MPI_PARALLEL
//=================================================================================================
//  HydroTree::UpdateGravityExportList
/// Compute gravitational forces due to other MPI node trees (using pruned trees).
/// If the other domains are too close (so the pruned trees are not adequate), then flag
/// cell to be exported to that MPI node for a full local tree-walk.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::UpdateGravityExportList
 (int rank,                            ///< [in] MPI rank
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to Hydrodynamics object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  const DomainBox<ndim> &simbox,       ///< [in] Simulation domain box
  Ewald<ndim> *ewald)                  ///< [in] Ewald object
{
  int cactive;                          // No. of active cells
  vector<TreeCellBase<ndim> > celllist; // List of active tree cells
  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();

  debug2("[GradhHydroTree::UpdateGravityExportForces]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("HYDRO_DISTANT_FORCES");



  // Find list of all cells that contain active particles
  assert(tree->Nimportedcell==0);
  cactive = tree->ComputeActiveCellList(celllist);


  // Reset all export lists
  for (int j=0; j<Nmpi; j++) {
    Npartexport[j] = 0;
    cellexportlist[j].clear();
    cellexportlist[j].reserve(tree->gmax);
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(rank,simbox,celllist,cactive,cout,hydro,partdata,ewald)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    int cc;                                    // Aux. cell counter
    int i;                                     // Particle id
    int j;                                     // Aux. particle counter
    int k;                                     // Dimension counter
    int Nactive;                               // No. of active particles in current cell
    int Ngravcelltemp;                         // Aux. gravity cell counter
    int *activelist                = activelistbuf[ithread];
    ParticleType<ndim> *activepart = activepartbuf[ithread];
    vector<MultipoleMoment<ndim> > gravcelllist;


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim>& cell = celllist[cc] ;
      gravcelllist.clear();


      // Find list of active particles in current cell
      Nactive = tree->ComputeActiveParticleList(cell, partdata, activelist);

      // Make local copies of active particles
      for (j=0; j<Nactive; j++) activepart[j] = partdata[activelist[j]];

      // Zero/initialise all summation variables for active particles
      for (j=0; j<Nactive; j++) {
        activepart[j].gpot = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].a[k]     = (FLOAT) 0.0;
        for (k=0; k<ndim; k++) activepart[j].atree[k] = (FLOAT) 0.0;
      }


      // Loop over all distant pruned trees and compute list of cells.
      // If pruned tree is too close, record cell id for exporting
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nmpi; j++) {
        if (j == rank) continue;
        const int Ngravcellold = gravcelllist.size();
        Ngravcelltemp = prunedtree[j]->ComputeDistantGravityInteractionList
          (cell, simbox ,gravcelllist);

        // If pruned tree is too close to be used (flagged by -1), then record cell id
        // for exporting those particles to other MPI processes
        if (Ngravcelltemp == -1) {
#pragma omp critical
          {
        	  cellexportlist[j].push_back(cell.id);
          }
#pragma omp atomic
          Npartexport[j] += Nactive;
          gravcelllist.resize(Ngravcellold);
        }

      }
      //-------------------------------------------------------------------------------------------


      // Compute gravitational force due to distant cells
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nactive; j++) {

        if (multipole == monopole) {
          ComputeCellMonopoleForces(activepart[j].gpot, activepart[j].atree, activepart[j].r,
                                    gravcelllist.size(), &gravcelllist[0]);
        }
        else if (multipole == quadrupole) {
          ComputeCellQuadrupoleForces(activepart[j].gpot, activepart[j].atree, activepart[j].r,
                                      gravcelllist.size(), &gravcelllist[0]);
        }

        // Add the periodic correction force the cells
        if (simbox.PeriodicGravity){
          int Ngravcell = gravcelllist.size() ;

          for (int jj=0; jj<Ngravcell; jj++) {
            MultipoleMoment<ndim> & gravcell = gravcelllist[jj] ;

            FLOAT draux[ndim], aperiodic[ndim], potperiodic ;
            for (int k=0; k<ndim; k++) draux[k] = gravcell.r[k] - activepart[j].r[k];

            ewald->CalculatePeriodicCorrection(gravcell.m, draux, aperiodic, potperiodic);

            for (int k=0; k<ndim; k++) activepart[j].atree[k] += aperiodic[k];
            activepart[j].gpot += potperiodic;
          }
        }

      }
      //-------------------------------------------------------------------------------------------


      // Compute 'fast' multipole terms here
      if (multipole == fast_monopole) {
        ComputeFastMonopoleForces(Nactive, gravcelllist.size(), &gravcelllist[0], cell,
                                 activepart, hydro->types);
      }
      if (multipole == fast_quadrupole) {
        ComputeFastQuadrupoleForces(Nactive, gravcelllist.size(), &gravcelllist[0], cell,
                                    activepart, hydro->types);
      }


      // Add all active particles contributions to main array
      for (j=0; j<Nactive; j++) {
        i = activelist[j];
        for (k=0; k<ndim; k++) partdata[i].a[k]     += activepart[j].atree[k];
        for (k=0; k<ndim; k++) partdata[i].atree[k] += activepart[j].atree[k];
        partdata[i].gpot += activepart[j].gpot;
      }

    }
    //=============================================================================================

  }
  //===============================================================================================

  return;
}



//=================================================================================================
//  HydroTree::UpdateHydroExportList
/// Check if any local particles need to be exported to other MPI nodes by comparing the
/// smoothing length box overlaps.  If so, then flag for exporting.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::UpdateHydroExportList
 (int rank,                            ///< [in] MPI rank
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to Hydrodynamics object
  Nbody<ndim> *nbody,                  ///< [in] Pointer to N-body object
  const DomainBox<ndim> &simbox)       ///< [in] Simulation domain box
{
  int cactive;                         // No. of active cells
  TreeCellBase<ndim> **celllist;           // List of pointers to binary tree cells
  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();

  debug2("[HydroTree::UpdateHydroExportList]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("MPI_HYDRO_EXPORT");


  // Find list of all cells that contain active particles
  celllist = new TreeCellBase<ndim>*[2*tree->gtot];
  cactive = tree->ComputeActiveCellPointers(celllist);

  // Reset all export lists
  for (int j=0; j<Nmpi; j++) {
    Npartexport[j] = 0;
    cellexportlist[j].clear();
    cellexportlist[j].reserve(tree->gmax);
  }


  // Set-up all OMP threads
  //===============================================================================================
#pragma omp parallel default(none) shared(simbox,celllist,cactive,cout,rank,partdata)
  {
#if defined _OPENMP
    const int ithread = omp_get_thread_num();
#else
    const int ithread = 0;
#endif
    bool overlapflag;                            // Flag if cells overlap
    int cc;                                      // Aux. cell counter
    int j;                                       // Aux. particle counter
    int *activelist = activelistbuf[ithread];    // List of active particle ids


    // Loop over all active cells
    //=============================================================================================
#pragma omp for schedule(guided)
    for (cc=0; cc<cactive; cc++) {
      TreeCellBase<ndim> *cellptr = celllist[cc];
      TreeCellBase<ndim>& cell = *cellptr;

      // Loop over all distant pruned trees and compute list of cells.
      // If pruned tree is too close, record cell id for exporting
      //-------------------------------------------------------------------------------------------
      for (j=0; j<Nmpi; j++) {
        if (j == rank) continue;

        overlapflag = prunedtree[j]->ComputeHydroTreeCellOverlap(cellptr, simbox);

        // If pruned tree is too close, then record cell id
        // for exporting to other MPI processes
        if (overlapflag) {

#pragma omp critical
          {
        	  cellexportlist[j].push_back(cell.id);
          }
          const int Nactive = tree->ComputeActiveParticleList(cell, partdata, activelist);

#pragma omp atomic
          Npartexport[j] += Nactive;
          // assert(Ncellexport[j] <= tree->gmax);
        }

      }
      //-------------------------------------------------------------------------------------------

    }
    //=============================================================================================

#ifdef VERIFY_ALL
    //PrintArray("Ncellexport : ",Nmpi,Ncellexport);
    PrintArray("Npartexport : ",Nmpi,Npartexport);
#endif


  }
  //===============================================================================================

  delete[] celllist;


  return;
}



//=================================================================================================
//  HydroTree::BuildPrunedTree
/// Constructs a pruned version of the local tree ready to be exported to other MPI processes.
/// Copies all levels up to and including 'pruning_level'.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::BuildPrunedTree
 (const int rank,                      ///< [in] Rank of local MPI node
  const DomainBox<ndim> &simbox,       ///< [in] Simulation domain box object
  const MpiNode<ndim> *mpinode,        ///< [in] Pointer to MPI node array
  Hydrodynamics<ndim> *hydro)          ///< [in] Pointer to Hydrodynamics object
  {
  bool localNode;                              // Is this pruned tree for the local node?
  int i;                                       // Particle counter
  TreeBase<ndim> *treeptr;                     // Pointer to tree object in question

  debug2("[HydroTree::BuildPrunedTree]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("BUILD_PRUNED_TREE");

  Nprunedcellmax = 0;
#if defined(OUTPUT_ALL)
  cout << "Levels : " << pruning_level_min << "    " << tree->ltot << endl;
#endif

  // Update all work counters in the tree for load-balancing purposes
  tree->UpdateWorkCounters();

  MPI_Request req[Nmpi-1];
  MPI_Status status[Nmpi-1];
  int j=0;
  int size_not_known=0;
  // Post all the receives
  for (int i=0; i<Nmpi; i++) {

	  if (i==rank)
		  continue;

	  TreeBase<ndim>* treeptr = prunedtree[i];

	  // Guess the maximum number of cells and allocate memory
	  const int max_cells = treeptr->GetMaxCellNumber(pruning_level_max);

	  if (max_cells != -1) {
	    // In this case we have an upper limit on how much data we are receiving
	    // Post the receive!

        treeptr->AllocateTreeMemory(0,max_cells,false);

        int cell_size = treeptr->GetTreeCellSize() ;
        MPI_Irecv(treeptr->GetCellDataPointer(),max_cells*cell_size,MPI_CHAR,i,3,MPI_COMM_WORLD,&req[j] );
	  }
	  else {
	    // In this case, no idea on how much data we are receiving. Start to send, and we will come back to this problem later
	    size_not_known=1;
	  }


	  j++;
  }


  MPI_Request send_req[Nmpi-1];
  j=0;
  // Set level at which tree will be pruned (for all trees)
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nmpi; i++) {

    // If creating pruned tree for local node, then set pointer to correct array.
    // Otherwise, set pointer so send buffer is filled with pruned tree
    if (i == rank) {
      localNode = true;
      treeptr = prunedtree[i];
    }
    else {
      localNode = false;
      treeptr = sendprunedtree[i];
    }

    treeptr->ltot     = pruning_level_max;
    treeptr->Ntot     = 0;
    treeptr->gmax     = 0;
    int max_cells = max(1,treeptr->GetMaxCellNumber(pruning_level_max));
    treeptr->AllocateTreeMemory(0,max_cells,false);


    treeptr->Ncell = tree->CreatePrunedTreeForMpiNode
      (mpinode[i], simbox, (FLOAT) 0.0, localNode, pruning_level_min, pruning_level_max,
       max_cells, treeptr);


    // If insufficient memory was allocated, then re-allocate larger array and repeat.
    //---------------------------------------------------------------------------------------------
    while (treeptr->Ncell == -1) {

     max_cells *= 2;

      treeptr->AllocateTreeMemory(0,max_cells,false);
      treeptr->Ncell = tree->CreatePrunedTreeForMpiNode
        (mpinode[i], simbox, (FLOAT) 0.0, localNode, pruning_level_min, pruning_level_max,
         max_cells, treeptr);

    }
    //---------------------------------------------------------------------------------------------


    // Allocate (or reallocate if needed) all tree memory
    Nprunedcellmax += treeptr->Ncell;

    // Now that the tree is ready, start the sending
    if (i != rank) {
      int cell_size = treeptr->GetTreeCellSize() ;
      MPI_Isend(treeptr->GetCellDataPointer(),treeptr->Ncell*cell_size,MPI_CHAR,i,3,MPI_COMM_WORLD,&send_req[j]);
      j++;
    }

  }
  //-----------------------------------------------------------------------------------------------

  // If we didn't know the size of the receives, we still need to post the receives!


  if (size_not_known) {

      vector<int> flags(Nmpi-1);
      int Ncompleted=0;
      while (Ncompleted<Nmpi-1) {
        int j=0;
        for (int i=0; i<Nmpi; i++) {

            if (i==rank)
                continue;

            if (flags[j]) {
              j++;
              continue;
            }

            TreeBase<ndim>* treeptr = prunedtree[i];

            // See how much stuff we have received
            MPI_Status status;
            MPI_Iprobe(i,3,MPI_COMM_WORLD,&flags[j],&status);

            if (flags[j]) {
              // We know how much stuff we are receiving
              // Allocate memory and post the receive

              int Nbytes_received;
              int cell_size = treeptr->GetTreeCellSize() ;
              MPI_Get_count(&status,MPI_CHAR,&Nbytes_received);
              const int max_cells = Nbytes_received/cell_size;

              treeptr->AllocateTreeMemory(0,max_cells,false);

              MPI_Irecv(treeptr->GetCellDataPointer(),max_cells*cell_size,MPI_CHAR,i,3,MPI_COMM_WORLD,&req[j] );
              Ncompleted++;
            }


            j++;
        }
    }
  }

  // Now wait for all sends to be completed
  MPI_Waitall(Nmpi-1,send_req,MPI_STATUSES_IGNORE);

  // Wait for all receives to be completed
  MPI_Waitall(Nmpi-1,req,status);

  // Set the number of cells
  j=0;
  for (int i=0; i<Nmpi; i++) {
	  if (i==rank)
		  continue;

	  int count;
	  MPI_Get_count(&status[j], MPI_CHAR, &count);

	  treeptr = prunedtree[i];

	  int cell_size = treeptr->GetTreeCellSize() ;
	  treeptr->Ncell= count/cell_size;
	  treeptr->first_stock = true;

#ifdef OUTPUT_ALL
	  cout << "on rank " << rank << " we received a pruned tree with " << treeptr->Ncell << " from " << i << endl;
#endif

	  j++;

  }

  return;
}


//=================================================================================================
//  HydroTree::StockPrunedTree
/// Stock the pruned trees without doing a full rebuild
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::StockPrunedTree
(const int rank,				     ///< [in] Rank of local MPI node
 Hydrodynamics<ndim>* hydro		)    ///< [inout] Pointer to hydrodynamics object
{
	  debug2("[HydroTree::StockPrunedTree]");
	  CodeTiming::BlockTimer timer = timing->StartNewTimer("STOCK_PRUNED_TREE");

	  ParticleType<ndim> *partdata = hydro->template GetParticleArray<ParticleType>();
	  // Update all work counters in the tree for load-balancing purposes
	  tree->UpdateWorkCounters();

	  // First thing to do is posting the receives
	  MPI_Request req[Nmpi-1];
	  int which_completed[Nmpi-1];
	  int req_completed=0;
	  MPI_Status status[Nmpi-1];
	  int j=0;

	  vector< vector <char> > receive_buffer(Nmpi);

	  for (int i=0; i<Nmpi; i++) {

		  // No need to send anything to ourselves
		  if (i==rank)
			  continue;

		  TreeBase<ndim>* treeptr = prunedtree[i];

	      const int cell_size = treeptr->GetTreeCellSize() ;

	      // Note that this is an overestimate; we are only receiving Nleaf cells
	      const int max_size = cell_size*treeptr->Ncell;

	      // Allocate buffer
		  receive_buffer[i].resize(max_size);

		  // Post the receive
		  MPI_Irecv(&receive_buffer[i][0],max_size,MPI_CHAR,i,4,MPI_COMM_WORLD,&req[j]);

		  j++;


	  }

	  MPI_Request send_req[Nmpi-1];
	  vector<vector<char> > send_buffer(Nmpi);
	  j=0;
	  // Loop over all the local pruned trees. Stock them and send them one by one as they get ready
	  for (int i=0; i<Nmpi; i++) {

		  TreeBase<ndim>* treeptr;
		  if (i==rank) {
			  treeptr=prunedtree[i];
		  }
		  else {
			  treeptr = sendprunedtree[i];
		  }

		  // Copy the stocked cells - easier than recomputing
		  treeptr->UpdateLeafCells(tree);

		  if (i != rank) {
			  // Post the send
			  const int Nleaf = treeptr->GetNLeafCells();
			  const int cell_size = treeptr->GetTreeCellSize();
			  const int size_send = cell_size*Nleaf;
			  send_buffer[i].reserve(size_send);

			  treeptr->CopyLeafCells(send_buffer[i],TreeBase<ndim>::to_buffer);

			  assert(send_buffer[i].size() == size_send);

			  MPI_Isend(&send_buffer[i][0],size_send,MPI_CHAR,i,4,MPI_COMM_WORLD,&send_req[j]);
			  j++;

		  }
		  else {
			  // Stock the tree (this is used for load balancing, even if we have a full version of this tree)
			  treeptr->StockTree(partdata,false);
		  }
	  }

	  // Now loop to see which receives are finished
	  while (req_completed<Nmpi-1) {
		  int completed_now;
		  MPI_Waitsome(Nmpi-1,req,&completed_now,which_completed,status);

		  // Unpack the information
		  for (int i=0; i<completed_now; i++) {
			  const int j=which_completed[i];
			  int iproc=j;
			  if (j >= rank)
				  iproc += 1;
			  TreeBase<ndim>* treeptr = prunedtree[iproc];
			  // See how much information we have actually received
			  int count;
			  MPI_Get_count(&status[i], MPI_CHAR, &count);
			  receive_buffer[iproc].resize(count);
			  // Copy the information into the leaf cells
			  treeptr->CopyLeafCells(receive_buffer[iproc],TreeBase<ndim>::from_buffer);
			  // Stock the pruned tree (but not the leaf cells!)
			  treeptr->StockTree(partdata,false);
		  }

		  req_completed += completed_now;

	  }

	  // Give a chance to the send requests to complete (even if they should already have) to free their resources
	  MPI_Waitall(Nmpi-1,send_req,MPI_STATUSES_IGNORE);


}



//=================================================================================================
//  HydroTree::BuildMpiGhostTree
/// Main routine to control how the tree is built, re-stocked and interpolated
/// during each timestep.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::BuildMpiGhostTree
 (const bool rebuild_tree,             ///< Flag to rebuild tree
  const int n,                         ///< Integer time
  const int ntreebuildstep,            ///< Tree build frequency
  const int ntreestockstep,            ///< Tree stocking frequency
  const FLOAT timestep,                ///< Smallest physical timestep
  Hydrodynamics<ndim> *hydro)          ///< Pointer to Hydrodynamics object
{
  ParticleType<ndim> *partdata = hydro->template GetParticleArray<ParticleType>();
  CodeTiming::BlockTimer timer = timing->StartNewTimer("BUILD_MPIGHOST_TREE");

  debug2("[HydroTree::BuildMpiGhostTree]");


  // Activate nested parallelism for tree building routines
#ifdef _OPENMP
  omp_set_nested(1);
#endif

#ifdef OUTPUT_ALL
  cout << "BUILDING TREE WITH " << hydro->Nmpighost << " MPI GHOSTS!!" << endl;
#endif

  if (hydro->Ntot > Ntotmax) {
	  Ntotmax = hydro->Ntot;
	  Ntot = hydro->Ntot;
	  ReallocateMemory();
  }

  // For tree rebuild steps
  //-----------------------------------------------------------------------------------------------
  if (n%ntreebuildstep == 0 || rebuild_tree) {

    mpighosttree->Ntot       = hydro->Nmpighost;
    const int max_particles    = max(mpighosttree->Ntot,hydro->Nhydromax);
    mpighosttree->BuildTree(hydro->Nhydro + hydro->NPeriodicGhost,
                            hydro->Nhydro + hydro->NPeriodicGhost + hydro->Nmpighost - 1,
                            mpighosttree->Ntot, max_particles, timestep, partdata);

  }

  // Else stock the tree
  //-----------------------------------------------------------------------------------------------
  else {

    mpighosttree->StockTree(partdata,true);

  }
  //-----------------------------------------------------------------------------------------------

#ifdef _OPENMP
  omp_set_nested(0);
#endif


  return;
}



//=================================================================================================
//  HydroTree::SearchMpiGhostParticles
/// Search through local domain for any MPI ghost particles that should be sent to the given
/// domain by checking the predicted position over a time tghost.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int HydroTree<ndim,ParticleType>::SearchMpiGhostParticles
 (const FLOAT tghost,                  ///< [in] Expected ghost life-time
  const Box<ndim> &mpibox,             ///< [in] Bounding box of MPI domain
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to Hydrodynamics object
  vector<int> &export_list)            ///< [out] List of particle ids
{

  const FLOAT grange = 2.0*ghost_range*kernrange;

  // Walk both trees.
  //-----------------------------------------------------------------------------------------------
  int Nexport = tree->FindBoxGhostParticles(tghost, grange, mpibox, export_list) ;
  Nexport += ghosttree->FindBoxGhostParticles(tghost, grange, mpibox, export_list) ;

  return Nexport;
}



//=================================================================================================
//  HydroTree::FindMpiTransferParticles
/// This function finds particles that this thread needs to export to other MPI threads. It
/// fills two vectors: 1) the list of particles to go to each processor and 2) the combined list
/// of all particles to send.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::FindMpiTransferParticles
 (Hydrodynamics<ndim> *hydro,                ///< [in] Pointer to Hydrodynamics class
  vector<vector<int> >& particles_to_export, ///< [inout] Vector that for each
                                             ///< node gives the list of particles to export
  vector<int>& all_particles_to_export,      ///< [inout] Vector containing all the particles
                                             ///<         that will be exported by this processor
  const vector<int>& potential_nodes,        ///< [in] Vector containing the potential nodes we
                                             ///<      might be sending particles to
  MpiNode<ndim>* mpinodes)                   ///< [in] Array of other mpi nodes
{
  int jnode;                                 // Aux. node counter
  int inode;                                 // MPI node id



  // Loop over potential domains and walk the tree for each bounding box
  //-----------------------------------------------------------------------------------------------
  int num_pot_nodes = potential_nodes.size();
  for (jnode=0; jnode<num_pot_nodes; jnode++) {

    inode = potential_nodes[jnode];
    Box<ndim>& nodebox = mpinodes[inode].domain;

    // Start from root-cell
    tree->FindBoxOverlapParticles(nodebox, particles_to_export[inode],
                                                     hydro->template GetParticleArray<ParticleType>()) ;

    // Copy particles to per processor array
    all_particles_to_export.insert(all_particles_to_export.end(),
                                        particles_to_export[inode].begin(),
                                        particles_to_export[inode].end());
  }
  //-----------------------------------------------------------------------------------------------

#ifndef NDEBUG
  VerifyUniqueIds(all_particles_to_export.size(),hydro->Nhydro,&all_particles_to_export[0]);
  vector<int> temp(all_particles_to_export); std::sort(temp.begin(),temp.end());
  for (int i=0; i<hydro->Nhydro;i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    Box<ndim>& domainbox = mpinodes[rank].domain;
    if (ParticleInBox(part,domainbox)) {
      const bool WillExport = std::binary_search(temp.begin(),temp.end(),i);
      if (WillExport) {
        // Deal with edge case when the particle is exactly on the boundary
        // (in which case exporting it is not wrong)
        bool edge=false;
        for (int k=0; k<ndim; k++) {
          if (part.r[k] == domainbox.min[k] || part.r[k] == domainbox.max[k]) {
            edge =true;
          }
        }
        if (!edge) assert(!std::binary_search(temp.begin(),temp.end(),i));
      }
    }
    else {
      //Edge case when the particle is at the right edge of the domain
      bool edge = false;
      for (int k=0; k<ndim; k++) if (part.r[k]==domainbox.max[k]) edge=true;
      if (!edge) assert(std::binary_search(temp.begin(),temp.end(),i));
      //if (edge) assert(!std::binary_search(temp.begin(),temp.end(),i));
      for (int jnode=0; jnode<potential_nodes.size(); jnode++) {
          inode=potential_nodes[jnode];
          vector<int>::iterator it = std::find(particles_to_export[inode].begin(),particles_to_export[inode].end(),i);
          const bool InDomain = ParticleInBox(part,mpinodes[inode].domain);
          if (InDomain) {
              assert(it != particles_to_export[inode].end());
          }
          else {
              assert( it == particles_to_export[inode].end());
          }
	  if (it != particles_to_export[inode].end()) {
              assert(InDomain);
          }
          else {
              assert(!InDomain);
          }
      }
    }
  }
#endif

  return;
}



//=================================================================================================
//  HydroTree::FindLoadBalancingDivision
/// Find the predicted cell-cell dividing point that approximately balances the CPU work load
/// in order to achieve load balancing amongst all MPI nodes.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
FLOAT HydroTree<ndim,ParticleType>::FindLoadBalancingDivision
 (int k_divide,                        ///< Dimension of cell division
  FLOAT r_old,                         ///< Old position of cell division
  FLOAT boxmin[ndim],                  ///< Minimum extent of parent MPI tree cell
  FLOAT boxmax[ndim])                  ///< Maximum extent of parent MPI tree cell
{
  int i;                               // MPI node counter
  int k;                               // Dimension counter
  FLOAT r_divide = r_old;              // Cell division location
  FLOAT r_max = boxmax[k_divide];      // Max. for bisection iteration of division
  FLOAT r_min = boxmin[k_divide];      // Min. for bisection iteration of division
  FLOAT workleft;                      // Work computed on LHS of division
  FLOAT workright;                     // Work computed on RHS of division
  FLOAT workfrac;                      // Fraction of work on LHS
  FLOAT worktol = (FLOAT) 0.001;        // Work balance tolerance for iteration
  FLOAT boxleftmin[ndim];              // Min. extent of left-box division
  FLOAT boxleftmax[ndim];              // Max. extent of left-box division
  FLOAT boxrightmin[ndim];             // Min. extent of right-box division
  FLOAT boxrightmax[ndim];             // Max. extent of right-box division

  // Set box extents from MPI tree node extent
  for (k=0; k<ndim; k++) {
    boxleftmin[k] = boxmin[k];
    boxleftmax[k] = boxmax[k];
    boxrightmin[k] = boxmin[k];
    boxrightmax[k] = boxmax[k];
  }


  // Find the work-balance position through bisection iteration
  //-----------------------------------------------------------------------------------------------
  do {
    boxleftmax[k_divide] = r_divide;
    boxrightmin[k_divide] = r_divide;
    workleft = (FLOAT) 0.0;
    workright = (FLOAT) 0.0;

    // Compute work included in left-hand side from pruned trees of all MPI domains
    for (i=0; i<Nmpi; i++) {
      workleft += prunedtree[i]->ComputeWorkInBox(boxleftmin, boxleftmax);
    }

    // Compute work included in right-hand side from pruned trees of all MPI domains
    for (i=0; i<Nmpi; i++) {
      workright += prunedtree[i]->ComputeWorkInBox(boxrightmin, boxrightmax);
    }

#ifdef OUTPUT_ALL
    cout << "workleft : " << workleft << "     workright : " << workright
         << "       workfrac : " << workleft / (workleft + workright)
         << "    rold : " << r_old << "     r_new : " << r_divide << endl;
#endif

    // If fraction of work on either side of division is too inbalanced, calculate new position
    // of division and perform another iteration.  Otherwise exit iteration loop.
    workfrac = workleft / (workleft + workright);
    if (workfrac < 0.5 - worktol) r_min = r_divide;
    else if (workfrac > 0.5 + worktol) r_max = r_divide;
    else break;

    r_divide = 0.5*(r_min + r_max);

  } while (fabs(workfrac - 0.5) > worktol);
  //-----------------------------------------------------------------------------------------------

  return r_divide;
}



//=================================================================================================
//  HydroTree::FindParticlesToTransfer
/// Compute on behalf of the MpiControl class the particles that are outside
/// the domain after a load balancing step and need to be transferred to other nodes
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::FindParticlesToTransfer
 (Hydrodynamics<ndim> *hydro,                ///< [in] Pointer to sph class
  vector<vector<int> >& id_export_buffers,   ///< [inout] List of ids to export for each node
  vector<int>& all_ids_export_buffer,        ///< [inout] List of all ids to export from proc
  const vector<int>& potential_nodes,        ///< [in] Potential nodes we might send particles to
  MpiNode<ndim>* mpinodes)                   ///< [in] Array of other mpi nodes
{
  // Loop over particles and prepare the ones to export
  int nodes_size = potential_nodes.size() ;
  for (int i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i) ;

    // Loop over potential domains and see if we need to transfer this particle to them
    for (int inode=0; inode<nodes_size; inode++) {
      int node_number = potential_nodes[inode];

      if (ParticleInBox(part, mpinodes[node_number].domain)) {
        id_export_buffers[node_number].push_back(i);
        all_ids_export_buffer.push_back(i);

        // The particle can belong only to one domain, so we can break from this loop
        break;
      }
    }
  }

  return;
}



//=================================================================================================
//  HydroTree::GetExportInfo
/// Get the array with the information that needs to be exported to the given processor
//=================================================================================================
template <int ndim, template <int> class ParticleType>
int HydroTree<ndim,ParticleType>::GetExportInfo
 (int iproc,                           ///< [in] No. of processor we want to send data to
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to GetParticleArray object
  vector<char >& send_buffer,          ///< [inout] Vector where the ptcls to export will be stored
  MpiNode<ndim>& mpinode,              ///< ..
  int rank,                            ///< ..
  int Nmpi)                            ///< [in] Array with information for the other mpi nodes
{
  int cactive = cellexportlist[iproc].size();
  int Nactive = Npartexport[iproc];

  const int size_particles  = tree->GetSizeOfExportedParticleData(Nactive) ;
  const int size_cells      = tree->GetSizeOfExportedCellData(cactive) ;
  const int old_size        = send_buffer.size();

  vector<int>& celllist = cellexportlist[iproc];


  assert(tree->Nimportedcell == 0);

  // Clear the array needed for bookkeeping (which active particles we sent to which processor)
  vector<int>& ids_active_particles = ids_sent_particles[iproc];
  ids_active_particles.clear();
  ids_active_particles.reserve(Nactive);

  // Correspondingly clear also the same kind of information for the active cells
  vector<int>& ids_active_cells = ids_sent_cells[iproc];
  ids_active_cells.clear();
  ids_active_cells.reserve(cactive);

  // Work out size of the information we are sending
  // Header consists of number of particles and number of cells
  send_buffer.reserve(size_particles + size_cells + old_size);


  // Write the packed data
  int exported_particles =
      tree->PackParticlesAndCellsForMPITransfer(cellexportlist[iproc],
                                                ids_active_cells,
                                                ids_active_particles,
                                                send_buffer,
                                                hydro->template GetParticleArray<ParticleType>()) ;

  assert(exported_particles == Nactive);
  assert(ids_active_particles.size() == static_cast<unsigned int>(Nactive)) ;
//  assert(ids_active_cells.size() == static_cast<unsigned int>(cactive)) ;

  return size_particles + size_cells;
}



//=================================================================================================
//  HydroTree::UnpackExported
/// Unpack the information exported from the other processors, contaning the particles
/// that were exported and
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::UnpackExported
 (vector<char >& received_array,
  Hydrodynamics<ndim> *hydro,
  const int iproc,
  vector< vector<char> >& receive_header,
  const int rank,
  const bool first_unpack)
{

  if (first_unpack) {
    tree->Nimportedcell = 0;
    assert(hydro->NImportedParticles == 0);
  }

  // Gather information about how many cells and particles we have received from each processor
  vector<int> imported_part_from_j(Nmpi);
  vector<int> imported_cell_from_j(Nmpi);
  for (int j=0; j<Nmpi-1; j++) {
    int i=j;
    if (i>= rank)
      i += 1;

    std::vector<char>::const_iterator iter = receive_header[j].begin() ;

    unsigned int read_size ;
    unpack_bytes(&read_size, iter) ;
    unpack_bytes(&imported_part_from_j[i], iter) ;
    unpack_bytes(&imported_cell_from_j[i], iter) ;

    assert(i != iproc || read_size == received_array.size()) ;
    assert(iter == receive_header[j].end()) ;
  }


  // Ensure there is enough memory
  if (first_unpack) {
    const int N_received_cells_total =
        std::accumulate(imported_cell_from_j.begin(),imported_cell_from_j.end(),0);
    const int N_received_part_total  =
        std::accumulate(imported_part_from_j.begin(),imported_part_from_j.end(),0);

	  hydro->AllocateMemory(hydro->Ntot + N_received_part_total);

	  if (hydro->Ntot + N_received_part_total > Ntotmax) {
		  Ntotmax = hydro->Ntot + N_received_part_total;
		  ReallocateMemory();
	  }

	  tree->ReallocateMemory(hydro->Ntot + N_received_part_total,tree->Ncell+N_received_cells_total);
  }

  // Get the info about how many particles/cells we received from this thread.
  int N_received_bytes     = received_array.size();
  int N_received_cells     = imported_cell_from_j[iproc];
  int N_received_particles = imported_part_from_j[iproc];

  if (N_received_bytes == 0) {
    N_imported_part_per_proc[iproc] = 0;
    N_imported_cells_per_proc[iproc] = 0;
    return;
  }

  N_imported_part_per_proc[iproc] = N_received_particles;
  N_imported_cells_per_proc[iproc] = N_received_cells;

  // Copy received particles inside main arrays and received cells inside the tree array
  // Also update the linked list
  const vector<int>::iterator nth_part = imported_part_from_j.begin() + iproc;
  const vector<int>::iterator nth_cell = imported_cell_from_j.begin() + iproc;
  const int offset_parts = hydro->Nhydro + hydro->Nghost +
      std::accumulate(imported_part_from_j.begin(),nth_part,0);
  const int offset_cells = tree->Ncell +
      std::accumulate(imported_cell_from_j.begin(),nth_cell,0);

  tree->UnpackParticlesAndCellsFromMPITransfer(offset_parts, N_received_particles,
                                               offset_cells, N_received_cells,
                                               received_array,
                                               hydro);

  //-----------------------------------------------------------------------------------------------


  // Update the hydro counters
  hydro->Ntot += N_received_particles;
  hydro->NImportedParticles += N_received_particles;

  // Update the tree counters
  tree->Nimportedcell += N_received_cells;
  tree->Ntot = hydro->Ntot;
  Ntot = hydro->Ntot;

#ifdef OUTPUT_ALL
  cout << "Importing " << N_received_particles << " from rank " << iproc << " on rank " << rank << endl;
#endif

  return;
}


//=================================================================================================
//  HydroTree::GetBackExportInfo
/// Return the data to transmit back to the other processors (particle acceleration etc.)
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::GetBackExportInfo
 (vector<char >& send_buffer,              ///< [inout] These arrays will be overwritten with the information to send
  Hydrodynamics<ndim> *hydro,              ///< [in] Pointer to the GetParticleArray object
  const int rank,						   ///< [in] Our rank
  const int iproc)                         ///< [in] Rank that we are sending to
{
  const int N_received_particles = N_imported_part_per_proc[iproc];
  const int N_received_cells     = N_imported_cells_per_proc[iproc];
  const int size_imp_part  = tree->GetSizeOfReturnedParticleData(N_received_particles);
  const int size_imp_cells = tree->GetSizeOfReturnedCellData(N_received_cells);


  // Work out which particles we need
  const vector<int>::iterator nth = N_imported_part_per_proc.begin()+iproc;
  const int part_start_index = hydro->Nhydro + hydro->Nghost +
      std::accumulate(N_imported_part_per_proc.begin(), nth, 0);

  const vector<int>::iterator nth_cell = N_imported_cells_per_proc.begin()+iproc;
  const int cell_start_index = tree->Ncell +
      std::accumulate(N_imported_cells_per_proc.begin(), nth_cell, 0);


  // Pack the particles for sending
  send_buffer.clear();
  send_buffer.reserve(size_imp_part+size_imp_cells);
  tree->PackParticlesAndCellsForMPIReturn(part_start_index, N_received_particles,
                                          cell_start_index, N_received_cells,
                                          send_buffer,
                                          hydro->template GetParticleArray<ParticleType>());


  assert(send_buffer.size() == static_cast<unsigned int>(size_imp_part+size_imp_cells)) ;

  return;
}



//=================================================================================================
//  HydroTree::UnpackReturnedExportInfo
/// Unpack the data that was returned by the other processors, summing the accelerations to the particles
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void HydroTree<ndim,ParticleType>::UnpackReturnedExportInfo
 (vector<char >& received_information,   ///< ..
  Hydrodynamics<ndim> *hydro,            ///< ..
  const int rank,						 ///< Our local rank
  const int iproc)                       ///< Remote processor
{
  // Unpack the info for this processor.

	const vector<int>& ids_active_particles = ids_sent_particles[iproc];
	const vector<int>& ids_active_cells = ids_sent_cells[iproc];

	tree->UnpackParticlesAndCellsForMPIReturn(ids_active_particles, ids_active_cells,
	                                          received_information, hydro) ;

  //-----------------------------------------------------------------------------------------------

  return;
}
#endif



#if defined(VERIFY_ALL)
//=================================================================================================
//  NeighbourSearch::CheckValidNeighbourList
/// Checks that the neighbour list generated by the grid is valid in that it
/// (i) does include all true neighbours, and
/// (ii) all true neigbours are only included once and once only.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void HydroTree<ndim,ParticleType>::CheckValidNeighbourList
 (int i,                               ///< [in] Particle i.d.
  int Ntot,                            ///< [in] Total no. of particles
  int Nneib,                           ///< [in] No. of potential neighbours
  int *neiblist,                       ///< [in] List of potential neighbour i.d.s
  ParticleType<ndim> *partdata,        ///< [in] Array of particle data
  string neibtype)                     ///< [in] Neighbour search type
{
  bool invalid_flag = false;           // Flag if neighbour list is invalid
  int count;                           // Valid neighbour counter
  int j;                               // Neighbour particle counter
  int k;                               // Dimension counter
  int Ntrueneib = 0;                   // No. of 'true' neighbours
  int *trueneiblist;                   // List of true neighbour ids
  FLOAT drsqd;                         // Distance squared
  FLOAT dr[ndim];                      // Relative position vector

  // Allocate array to store local copy of potential neighbour ids
  trueneiblist = new int[Ntot];


  // First, create list of 'true' neighbours by looping over all particles
  if (neibtype == "gather") {
    for (j=0; j<Ntot; j++) {
      for (k=0; k<ndim; k++) dr[k] = partdata[j].r[k] - partdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (partdata[j].flags.is_dead()) continue;
      if (drsqd <= kernrangesqd*partdata[i].h*partdata[i].h) trueneiblist[Ntrueneib++] = j;
    }
  }
  else if (neibtype == "all") {
    for (j=0; j<Ntot; j++) {
      for (k=0; k<ndim; k++) dr[k] = partdata[j].r[k] - partdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      if (partdata[j].flags.is_dead()) continue;
      if (drsqd < kernrangesqd*partdata[i].h*partdata[i].h ||
          drsqd < kernrangesqd*partdata[j].h*partdata[j].h) trueneiblist[Ntrueneib++] = j;
    }
  }


  // Now compare each given neighbour with true neighbour list for validation
  for (j=0; j<Ntrueneib; j++) {
    count = 0;
    for (k=0; k<Nneib; k++) {
      if (neiblist[k] == trueneiblist[j]) {
        count++;
      }
    }

    // If the true neighbour is not in the list, or included multiple times,
    // then output to screen and terminate program
    if (count != 1) {
      for (k=0; k<ndim; k++) dr[k] = partdata[trueneiblist[j]].r[k] - partdata[i].r[k];
      drsqd = DotProduct(dr,dr,ndim);
      cout << "Could not find neighbour " << j << "   " << trueneiblist[j] << "     " << i
           << "      " << sqrt(drsqd)/kernrange/partdata[i].h << "     "
           << sqrt(drsqd)/kernrange/partdata[trueneiblist[j]].h << "    "
           << partdata[trueneiblist[j]].r[0] << "   type : "
           << partdata[trueneiblist[j]].ptype << endl;
      invalid_flag = true;
    }

  }


  // If the true neighbour is not in the list, or included multiple times,
  // then output to screen and terminate program
  if (invalid_flag) {
    cout << "Problem with neighbour lists : " << i << "  " << j << "   "
         << count << "   " << partdata[i].r[0] << "   " << partdata[i].h << endl;
    cout << "Nneib : " << Nneib << "   Ntrueneib : " << Ntrueneib
         << "    neibtype : " << neibtype << endl;
    InsertionSort(Nneib,neiblist);
    PrintArray("neiblist     : ",Nneib,neiblist);
    PrintArray("trueneiblist : ",Ntrueneib,trueneiblist);
    string message = "Problem with neighbour lists in tree search";
    ExceptionHandler::getIstance().raise(message);
  }


  delete[] trueneiblist;

  return;
}
#endif



template class HydroTree<1,GradhSphParticle>;
template class HydroTree<2,GradhSphParticle>;
template class HydroTree<3,GradhSphParticle>;

template class HydroTree<1,SM2012SphParticle>;
template class HydroTree<2,SM2012SphParticle>;
template class HydroTree<3,SM2012SphParticle>;

template class HydroTree<1,MeshlessFVParticle>;
template class HydroTree<2,MeshlessFVParticle>;
template class HydroTree<3,MeshlessFVParticle>;
