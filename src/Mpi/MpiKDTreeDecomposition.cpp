//=================================================================================================
//  MpiKDTreeDecomposition.cpp
//  Contains functions for Main MPI class which controls the distribution of
//  work amongst all MPI tasks for the current simulation, including load
//  balancing and moving and copying particles between nodes.
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
#include <numeric>
#include "Constants.h"
#include "Precision.h"
#include "SmoothingKernel.h"
#include "DomainBox.h"
#include "Diagnostics.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
#include "MpiControl.h"
#include "MpiTree.h"
#include "Sinks.h"
using namespace std;



//=================================================================================================
//  MpiKDTreeDecomposition::MpiKDTreeDecomposition
/// MpiKDTreeDecomposition constructor (currently only calls base MpiControlType constructor)
//=================================================================================================
template <int ndim, template<int> class ParticleType>
MpiKDTreeDecomposition<ndim, ParticleType>::MpiKDTreeDecomposition() :
  MpiControlType<ndim,ParticleType>()
{
}



//=================================================================================================
//  MpiKDTreeDecomposition::CreateInitialDomainDecomposition
/// Creates a binary tree containing all particles in order to determine how to distribute the
/// particles across all MPI nodes with an equal amount of CPU work per MPI node.  If creating the
/// initial partition (i.e. before we have calculated the timestep), we give the particles equal
/// weighting and therefore each node will have equal numbers of particles.  For later steps
/// (i.e. when we know the timesteps and work information), we split the domains to give each MPI
/// node equal amounts of work.  This routine should only be called for the root process.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MpiKDTreeDecomposition<ndim, ParticleType>::CreateInitialDomainDecomposition
 (Hydrodynamics<ndim> *hydro,          ///< Pointer to main Hydrodynamics object
  Nbody<ndim> *nbody,                  ///< Pointer to main N-body object
  Parameters *simparams,               ///< Simulation parameters
  DomainBox<ndim> simbox,              ///< Simulation domain box
  bool &initial_h_provided)            ///< Receives from root whether or not initial h was provided
{
  int i;                               // Particle counter
  int initial_h = initial_h_provided;  // ..
  int inode;                           // Node counter
  int k;                               // Dimension counter
  ParticleType<ndim> *partbuffer;      // ..

  debug2("[MpiKDTreeDecomposition::CreateInitialDomainDecomposition]");

  // Create MPI tree object (for all nodes)
  mpitree = new MpiTree<ndim,ParticleType>(Nmpi);

  // Broadcast whether or not the initial h was provided
  MPI_Bcast(&initial_h, 1, MPI_INT, 0, MPI_COMM_WORLD);
  initial_h_provided = initial_h;


  // For periodic simulations, set bounding box of root node to be the periodic box size.
  // Otherwise, set to extend to infinity.
  if (simbox.boundary_lhs[0] == openBoundary) mpibox.boxmin[0] = -big_number;
  else mpibox.boxmin[0] = simbox.boxmin[0];
  if (simbox.boundary_rhs[0] == openBoundary) mpibox.boxmax[0] = big_number;
  else mpibox.boxmax[0] = simbox.boxmax[0];
  if (ndim > 1) {
    if (simbox.boundary_lhs[1] == openBoundary) mpibox.boxmin[1] = -big_number;
    else mpibox.boxmin[1] = simbox.boxmin[1];
    if (simbox.boundary_rhs[1] == openBoundary) mpibox.boxmax[1] = big_number;
    else mpibox.boxmax[1] = simbox.boxmax[1];
  }
  if (ndim == 3) {
    if (simbox.boundary_lhs[2] == openBoundary) mpibox.boxmin[2] = -big_number;
    else mpibox.boxmin[2] = simbox.boxmin[2];
    if (simbox.boundary_rhs[2] == openBoundary) mpibox.boxmax[2] = big_number;
    else mpibox.boxmax[2] = simbox.boxmax[2];
  }
  //mpitree->box = &mpibox;

#ifdef OUTPUT_ALL
  cout << "Simulation bounding box" << endl;
  for (k=0; k<ndim; k++) {
    cout << "r[" << k << "]  :  " << mpibox.boxmin[k] << "   " << mpibox.boxmax[k] << endl;
  }
#endif


  // For main process, create load balancing tree, transmit information to all
  // other nodes including particle data
  //===============================================================================================
  if (rank == 0) {

    // Set number of tree members to total no. of hydro particles (inc. ghosts)
    mpitree->Nhydro  = hydro->Nhydro;
    mpitree->Ntot    = hydro->Nhydro;
    mpitree->Ntotmax = hydro->Nhydromax;

    // Create all other MPI node objects
    this->AllocateMemory(mpitree->Ntotmax);

    // Get pointer to hydro particles and cast it to the right type
    ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

    // Compute the size of all tree-related arrays now we know number of points
    mpitree->ComputeTreeSize();

    // Allocate (or reallocate if needed) all tree memory
    mpitree->AllocateMemory();

    // Create tree data structure including linked lists and cell pointers
    mpitree->CreateTreeStructure(mpinode);

    // Set properties for root cell before constructing tree
    mpitree->tree[0].N      = hydro->Nhydro;
    mpitree->tree[0].ifirst = 0;
    mpitree->tree[0].ilast  = hydro->Nhydro - 1;
    for (k=0; k<ndim; k++) mpitree->tree[0].boxmin[k] = mpibox.boxmin[k];
    for (k=0; k<ndim; k++) mpitree->tree[0].boxmax[k] = mpibox.boxmax[k];
    for (i=0; i<hydro->Nhydro; i++) mpitree->inext[i] = -1;
    for (i=0; i<hydro->Nhydro-1; i++) mpitree->inext[i] = i + 1;
    for (i=0; i<hydro->Nhydro; i++) mpitree->ids[i] = i;

    // Recursively divide tree up until we've reached bottom level
    mpitree->DivideTreeCell(0, mpitree->Ntot-1, partdata, mpitree->tree[0]);

#ifdef OUTPUT_ALL
    cout << "Tree[" << rank << "] : " << mpitree->ltot << "   " << mpitree->Ncell << endl;
#endif

    // Broadcast MPI tree to all other nodes
    MPI_Bcast(&mpitree->Nhydro, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mpitree->Ntot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mpitree->Ntotmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(mpitree->tree, mpitree->Ncell*sizeof(MpiTreeCell<ndim>),
              MPI_CHAR, 0, MPI_COMM_WORLD);


    // Update all MPI node bounding boxes
    //---------------------------------------------------------------------------------------------
    for (inode=0; inode<Nmpi; inode++) {
      int icell = mpitree->g2c[inode];

      // Create bounding boxes containing particles in each sub-tree
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmin[k] = mpitree->tree[icell].boxmin[k];
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmax[k] = mpitree->tree[icell].boxmax[k];

#ifdef OUTPUT_ALL
      cout << "CHECKING MPITREE : " << inode << "   " << icell << "   "
           << &mpitree->tree[icell] << endl;
#endif

      // Copy particle ids to node lists
      mpinode[inode].Nhydro = 0;
      i = mpitree->tree[icell].ifirst;

      while (i != -1) {
        mpinode[inode].ids[mpinode[inode].Nhydro++] = i;
        if (i == mpitree->tree[icell].ilast) break;
        i = mpitree->inext[i];
      };
      mpinode[inode].Ntot = mpinode[inode].Nhydro;

#ifdef OUTPUT_ALL
      cout << "MPIDOMAIN : " << inode << "   Nhydro : " << mpinode[inode].Nhydro << "    box : "
           << mpinode[inode].domain.boxmin[0] << "    " << mpinode[inode].domain.boxmax[0] << endl;
#endif

    }
    //---------------------------------------------------------------------------------------------

#ifdef OUTPUT_ALL
    cout << "CHECKING SENDPARTICLES : " << mpinode[inode].ids << "   " << partdata << endl;
#endif

    // Send particles to all other domains
    for (inode=1; inode<Nmpi; inode++) {
      this->SendParticles(inode, mpinode[inode].Nhydro, mpinode[inode].ids, partdata);
#ifdef OUTPUT_ALL
      cout << "Sent " << mpinode[inode].Nhydro << " particles to node " << inode << endl;
#endif
    }
#ifdef OUTPUT_ALL
    cout << "Sent all particles to other processes" << endl;
#endif

    // Delete all other particles from local domain
    hydro->Nhydro = mpinode[0].Nhydro;
    partbuffer = new ParticleType<ndim>[hydro->Nhydro];
    for (i=0; i<hydro->Nhydro; i++) partbuffer[i] = partdata[mpinode[0].ids[i]];
    for (i=0; i<hydro->Nhydro; i++) partdata[i] = partbuffer[i];
    delete[] partbuffer;
#ifdef OUTPUT_ALL
    cout << "Deleted all other particles from root node" << endl;
#endif

  }

  // For other nodes, receive all bounding box and particle data once transmitted by main process.
  //===============================================================================================
  else {

    // Receive all broadcasts of MPI tree
    MPI_Bcast(&mpitree->Nhydro, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mpitree->Ntot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mpitree->Ntotmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Allocate all memory and prepare important variables for tree
    this->AllocateMemory(mpitree->Ntotmax);
    mpitree->ComputeTreeSize();
    mpitree->AllocateMemory();
    mpitree->CreateTreeStructure(mpinode);

#ifdef OUTPUT_ALL
    cout << "Tree[" << rank << "] : " << mpitree->ltot << "    " << mpitree->Ncell << endl;
#endif

    // Now receive all data from tree nodes
    MPI_Bcast(mpitree->tree, mpitree->Ncell*sizeof(MpiTreeCell<ndim>), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Update all MPI node bounding boxes
    for (inode=0; inode<Nmpi; inode++) {
      int icell = mpitree->g2c[inode];

      // Create bounding boxes containing particles in each sub-tree
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmin[k] = mpitree->tree[icell].boxmin[k]; //bbmin[k];
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmax[k] = mpitree->tree[icell].boxmax[k]; //bbmax[k];
    }

#ifdef OUTPUT_ALL
    //cout << "CHECKING SENDPARTICLES : " << mpinode[inode].ids << "   " << partdata << endl;
    cout << "Memory allocated?    Nhydromax : " << hydro->Nhydromax << endl;
#endif

    // Now, receive particles form main process and copy to local main array
    this->ReceiveParticles(0, (hydro->Nhydro), &partbuffer);

    hydro->AllocateMemory(hydro->Nhydro);
    mpinode[rank].Nhydro = hydro->Nhydro;

    // Get pointer to hydro particles and cast it to the right type
    ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());


    // Update all MPI node bounding boxes
#ifdef OUTPUT_ALL
    //---------------------------------------------------------------------------------------------
    for (inode=0; inode<Nmpi; inode++) {
      int icell = mpitree->g2c[inode];
      cout << "CHECKING MPITREE : " << inode << "   " << icell << "   "
           << &mpitree->tree[icell] << "    N : " << mpitree->tree[icell].N << endl;
      cout << "MPIDOMAIN : " << inode << "   Nhydro : " << mpinode[inode].Nhydro << "    box : "
           << mpinode[inode].domain.boxmin[0] << "    " << mpinode[inode].domain.boxmax[0] << endl;
    }
    //---------------------------------------------------------------------------------------------

    cout << "Received particles on node " << rank << "   Nhydro : " << hydro->Nhydro << endl;
#endif

    for (i=0; i<hydro->Nhydro; i++) partdata[i] = partbuffer[i];
    delete[] partbuffer;
#ifdef OUTPUT_ALL
    cout << "Deallocated partbuffer" << endl;
#endif

  }
  //===============================================================================================


  // Update all bounding boxes
  this->UpdateAllBoundingBoxes(hydro->Nhydro, hydro, hydro->kernp);

  // Share the stars with all other domains
  MPI_Bcast(&(nbody->Nstar), 1, MPI_INT, 0, MPI_COMM_WORLD);
  nbody->AllocateMemory(nbody->Nstar);

  if (nbody->Nstar > 0) {
    MPI_Bcast(nbody->stardata, sizeof(StarParticle<ndim>)*nbody->Nstar,
              MPI_BYTE, 0, MPI_COMM_WORLD);
  }

  return;
}



//=================================================================================================
//  MpiKDTreeDecomposition::LoadBalancing
/// If we are on a load balancing step, then determine which level of the binary partition we are
/// adjusting for load balancing.  Next, adjust the domain boundaries at that level (and for all
/// child domains).  Then broadcast the new domain boundaries to all other nodes to determine
/// which particles should be transfered to new nodes.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MpiKDTreeDecomposition<ndim, ParticleType >::LoadBalancing
 (Hydrodynamics<ndim> *hydro,          ///< [inout] Pointer to main Hydrodynamics object
  Nbody<ndim> *nbody)                  ///< [inout] Pointer to main N-body object
{
  int c;                               // MPI tree cell counter
  int c2;                              // i.d. of second child cell
  int inode;                           // MPI node counter
  int k;                               // Dimension counter
  int l;                               // MPI tree level counter
  int lbalance = 0;                    // Load balance level (always top level for now)
  FLOAT rold;                          // Position of previous load balancing division
  DOUBLE worktot = 0.0;                // Total work on all nodes

  // If running on only one MPI node, return immediately
  if (Nmpi == 1) return;
  debug2("[MpiKDTreeDecomposition::LoadBalancing]");
  timing->StartTimingSection("MPI_LOAD_BALANCING");

  //Get pointer to sph particles and cast it to the right type
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

  // Sum-up total work on all MPI nodes
  for (inode=0; inode<Nmpi; inode++) worktot += 0.0;


  // Starting with the highest MpiTree division, start adjusting divisional positions to achieve
  // equal amounts of work on each side of the divide.  Use the extrapolated cells from each
  // pruned tree to compute work done.
  //-----------------------------------------------------------------------------------------------
  for (l=lbalance; l<mpitree->ltot; l++) {

    // Loop over all MPI tree cells on current balancing level
    for (c=0; c<mpitree->Ncell; c++) {
      if (mpitree->tree[c].level != l) continue;
      c2   = mpitree->tree[c].c2;
      k    = mpitree->tree[c].k_divide;
      rold = mpitree->tree[c].r_divide;

#ifdef OUTPUT_ALL
      cout << "Previous load balancing division for " << c << "    rold : " << rold
           << "     bb : " << mpitree->tree[c].boxmin[k]<< "   "
           << mpitree->tree[c].boxmax[k] << "   k_divide : " << k << endl;
#endif

      // In case of extreme movement of the load balancing positions (perhaps due to latency or
      // large movement of active particles between nodes), then set some arbitary position of
      // the first division guess in order to search for correct division
      if (mpitree->tree[c].boxmin[k] > mpitree->tree[c].r_divide ||
          mpitree->tree[c].boxmax[k] < mpitree->tree[c].r_divide) {
        mpitree->tree[c].r_divide = 0.5*(mpitree->tree[c].boxmin[k] + mpitree->tree[c].boxmax[k]);
      }
      assert(mpitree->tree[c].boxmin[k] < mpitree->tree[c].r_divide);
      assert(mpitree->tree[c].boxmax[k] > mpitree->tree[c].r_divide);

      // Now find new division between child cells that is load-balanced
      mpitree->tree[c].r_divide = neibsearch->FindLoadBalancingDivision
        (mpitree->tree[c].k_divide, mpitree->tree[c].r_divide,
         mpitree->tree[c].boxmin, mpitree->tree[c].boxmax);
#ifdef OUTPUT_ALL
      cout << "Moved load balancing division for " << c << "    rold : " << rold
           << "     rnew : " << mpitree->tree[c].r_divide << "    k_divide : "
           << mpitree->tree[c].k_divide << endl;
#endif
    }

    // Update all cell bounding boxes now new divisions have been computed
    mpitree->UpdateBoundingBoxes();

    // Finally update all MPI node bounding boxes
    for (inode=0; inode<Nmpi; inode++) {
      int icell = mpitree->g2c[inode];

      // Create bounding boxes containing particles in each sub-tree
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmin[k] = mpitree->tree[icell].boxmin[k];
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmax[k] = mpitree->tree[icell].boxmax[k];
    }

  }
  //-----------------------------------------------------------------------------------------------

  // Write out 'final' mpinode bounding boxes
#ifdef OUTPUT_ALL
  for (int inode=0; inode<Nmpi; inode++) {
    cout << "Node " << inode << "     box : " << mpinode[inode].domain.boxmin[0] << "     "
         << mpinode[inode].domain.boxmax[0] << endl;
  }
#endif


  // Update all node bounding boxes now domains have been reset
  this->UpdateAllBoundingBoxes(hydro->Nhydro, hydro, hydro->kernp);


  // Prepare lists of particles that now occupy other processor domains that need to be transfered
  // First construct the list of nodes that we might be sending particles to
  vector<int> potential_nodes;
  potential_nodes.reserve(Nmpi);
  for (int inode=0; inode<Nmpi; inode++) {
    if (inode == rank) continue;
    if (BoxesOverlap(mpinode[inode].domain, mpinode[rank].rbox)) potential_nodes.push_back(inode);
  }


  // Find the ptcls that need to be transferred - delegate to NeighbourSearch
  vector<vector<int> > particles_to_transfer(Nmpi);
  vector<int> all_particles_to_export;
  // Send and receive particles from/to all other nodes
  vector<ParticleType<ndim> > sendbuffer, recvbuffer;
  neibsearch->FindParticlesToTransfer(hydro, particles_to_transfer, all_particles_to_export,
                                      potential_nodes, mpinode);


  //-----------------------------------------------------------------------------------------------
  for (int iturn = 0; iturn<my_matches.size(); iturn++) {
    int inode = my_matches[iturn];

    int N_to_transfer = particles_to_transfer[inode].size();
    sendbuffer.clear();
    sendbuffer.resize(N_to_transfer);
    recvbuffer.clear();

    // Copy particles into send buffer
    for (int ipart=0; ipart < N_to_transfer; ipart++) {
      int index = particles_to_transfer[inode][ipart];
      sendbuffer[ipart] = partdata[index];
    }

    // Decide if we have to send or receive first
    bool send_turn;
    if (rank < inode) send_turn = true;
    else send_turn = false;

    // Do the actual communication, sending and receiving in the right order
    for (int i=0; i<2; i++) {
      if (send_turn) {
        //cout << "Sending " << N_to_transfer << " from " << rank << " to " << inode << endl;
        MPI_Send(&sendbuffer[0], N_to_transfer, particle_type, inode, tag_bal, MPI_COMM_WORLD);
        send_turn = false;
#ifdef OUTPUT_ALL
        cout << "TRANSFERING " << N_to_transfer << " particles from " << rank << " to node " << inode << endl;
#endif
      }
      else {
        int N_to_receive;
        MPI_Status status;
        MPI_Probe(inode, tag_bal, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, particle_type, &N_to_receive);
        recvbuffer.resize(N_to_receive);
#ifdef OUTPUT_ALL
        cout << "Rank " << rank << " receiving " << N_to_receive << " from " << inode << endl;
#endif
        if (hydro->Nhydro + N_to_receive > hydro->Nhydromax) {
          cout << "Memory problem : " << rank << " " << hydro->Nhydro
               << " " << N_to_receive << " " << hydro->Nhydromax <<endl;
          string message = "Not enough memory for transfering particles";
          ExceptionHandler::getIstance().raise(message);
        }
        MPI_Recv(&recvbuffer[0], N_to_receive, particle_type, inode,
                 tag_bal, MPI_COMM_WORLD, &status);
        send_turn = true;
      }
    }

    // Copy particles from receive buffer to main arrays
    int running_counter = hydro->Nhydro;
    // TODO: check we have enough memory
    for (int i=0; i< recvbuffer.size(); i++) {
      partdata[running_counter] = recvbuffer[i];
      running_counter++;
    }
    hydro->Nhydro = running_counter;

  }
  //-----------------------------------------------------------------------------------------------


  // Remove transferred particles
  for (int i=0; i<all_particles_to_export.size(); i++) {
    partdata[all_particles_to_export[i]].itype = dead;
  }
  hydro->DeleteDeadParticles();

  timing->EndTimingSection("MPI_LOAD_BALANCING");

  return;
}



// Template class for each particle type
template class MpiKDTreeDecomposition<1, GradhSphParticle>;
template class MpiKDTreeDecomposition<2, GradhSphParticle>;
template class MpiKDTreeDecomposition<3, GradhSphParticle>;
template class MpiKDTreeDecomposition<1, SM2012SphParticle>;
template class MpiKDTreeDecomposition<2, SM2012SphParticle>;
template class MpiKDTreeDecomposition<3, SM2012SphParticle>;
