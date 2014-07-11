//=============================================================================
//  MpiControl.cpp
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
//=============================================================================


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <numeric>
#include "Constants.h"
#include "Precision.h"
#include "SphKernel.h"
#include "DomainBox.h"
#include "Diagnostics.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
#include "MpiControl.h"
#include "MpiTree.h"
using namespace std;



//=============================================================================
//  MpiControl::MpiControl()
/// MPI control class constructor.  Initialises all MPI control variables,
/// plus calls MPI routines to find out rank and number of processors.
//=============================================================================
template <int ndim>
MpiControl<ndim>::MpiControl()
{
  int len;                          // Length of host processor name string
  Box<ndim> dummy;                  // Dummy box variable

  allocated_mpi = false;
  balance_level = 0;

  // Find local processor rank, total no. of processors and host processor name
  MPI_Comm_size(MPI_COMM_WORLD,&Nmpi);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Get_processor_name(hostname,&len);

  // Create diagnostics data structure in database
  diagnostics_type = Diagnostics<ndim>::CreateMpiDataType();
  MPI_Type_commit(&diagnostics_type);

  // Create and commit the box datatype
  box_type = CreateBoxType(dummy);
  MPI_Type_commit(&box_type);

  // Create and commit the exported particle datatype
  ExportParticleType = ExportParticleInfo<ndim>::CreateMpiDataType();
  MPI_Type_commit(&ExportParticleType);

  // Create and commit the exported back particle datatype
  ExportBackParticleType = ExportBackParticleInfo<ndim>::CreateMpiDataType();
  MPI_Type_commit(&ExportBackParticleType);

  // Allocate buffer to send and receive boxes
  boxes_buffer.resize(Nmpi);

  // Allocate the buffers needed to send and receive particles
  num_particles_export_per_node.resize(Nmpi);
  displacements_send.resize(Nmpi);
  num_particles_to_be_received.resize(Nmpi);
  receive_displs.resize(Nmpi);
  Nbytes_exported_from_proc.resize(Nmpi);
  Nbytes_to_each_proc.resize(Nmpi);

  CreateLeagueCalendar();

#ifdef VERIFY_ALL
  if (this->rank == 0)
    printf("MPI working.  Nmpi : %d   rank : %d   hostname : %s\n",
           Nmpi,rank,hostname);
  else
    printf("%d is running too!!\n",this->rank);

  if (Nmpi > 1) {
    if (rank ==0) {
      SphParticle<ndim> particle;
      particle.gradrho[ndim-1]=-1;
      MPI_Send(&particle,1,particle_type,1,0,MPI_COMM_WORLD);
    }
    else if (rank ==1) {
      SphParticle<ndim> particle;
      MPI_Status status;
      MPI_Recv(&particle,1,particle_type,0,0,MPI_COMM_WORLD,&status);
      if (particle.gradrho[ndim-1]!=-1)
        cerr << "Error in transmitting particles: the last field has not been received correctly!" << endl;
    }
  }
#endif

}



//=============================================================================
//  MpiControl::~MpiControl()
/// MPI node class destructor.
//=============================================================================
template <int ndim, template<int> class ParticleType>
MpiControlType<ndim, ParticleType>::MpiControlType() :
MpiControl<ndim>() {
  //Allocate buffers
  particles_to_export_per_node.resize(Nmpi);

  // Create and commit the particle datatype
  particle_type = ParticleType<ndim>::CreateMpiDataType();
  MPI_Type_commit(&particle_type);

  bruteforce = NULL;

}



//=============================================================================
//  MpiControl::~MpiControl()
/// MPI node class destructor.
//=============================================================================
template <int ndim>
MpiControl<ndim>::~MpiControl()
{
  //MPI_Type_free(&particle_type);
  MPI_Type_free(&box_type);
  MPI_Type_free(&diagnostics_type);
}



//=============================================================================
//  MpiControl::AllocateMemory
/// Allocate all memory for MPI control class.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::AllocateMemory(int _Ntot)
{
  mpinode = new MpiNode<ndim>[Nmpi];
  for (int inode=0; inode<Nmpi; inode++) {
    mpinode[inode].Ntotmax = (2*_Ntot)/Nmpi;
    mpinode[inode].ids = new int[mpinode[inode].Ntotmax];
    mpinode[inode].worksent = new FLOAT[Nmpi];
    mpinode[inode].workreceived = new FLOAT[Nmpi];
  }

  return;
}



//=============================================================================
//  MpiControl::DeallocateMemory
/// Deallocate all MPI control class memory.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::DeallocateMemory(void)
{
  for (int inode=0; inode<Nmpi; inode++) {
    delete[] mpinode[inode].worksent;
    delete[] mpinode[inode].workreceived;
  }
  delete[] mpinode;

  return;
}



//=============================================================================
//  MpiControl::CreateLeagueCalendar
/// Create the calendar for the League. In this analogy with soccer,
/// each communication between two nodes is like a football match.
/// Since everybody neads to speak with everybody, this is effectively
/// like a league. We use the scheduling/Berger algorithm
/// (http://en.wikipedia.org/wiki/Round-robin_tournament#Scheduling_algorithm,
/// http://it.wikipedia.org/wiki/Algoritmo_di_Berger) to organize the event.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::CreateLeagueCalendar(void)
{

  //First check that number of processes is even
  if (! Nmpi%2) {
    std::string error = "The number of MPI processes must be even!";
    ExceptionHandler::getIstance().raise(error);
  }

  int Nturns = Nmpi-1;

  //Allocate memory for the calendar
  my_matches.resize(Nturns);


  if (rank==0) {

    //Create a vector containing the full calendar
    //First index is process
    std::vector<std::vector<int> > calendar(Nmpi);
    //And second is turn
    for (int iteam=0; iteam< calendar.size(); iteam++) {
      std::vector<int>& calendar_team = calendar[iteam];
      calendar_team.resize(Nturns);
    }

    // Create pairs table
    std::vector<std::vector<int> > pairs (Nturns);
    for (int iturn=0; iturn< pairs.size(); iturn++) {
      std::vector<int>& pairs_turn = pairs[iturn];
      pairs_turn.resize(Nmpi);
      // Fill in the pairs table
      pairs_turn[0]=Nturns;
      for (int i=1; i<pairs_turn.size(); i++) {
        pairs_turn[i] = (i+iturn) % (Nmpi-1);
      }
      // Can now fill in the calendar
      for (int istep =0; istep<Nmpi/2; istep++) {
        int first_team = pairs_turn[istep];
        int size = pairs_turn.size()-1;
        int second_team = pairs_turn[size-istep];
        calendar[first_team][iturn]=second_team;
        calendar[second_team][iturn]=first_team;
      }
    }


#if defined VERIFY_ALL
    //Validate the calendar
    std::vector<bool> other_teams(Nmpi-1);

    for (int iteam=0; iteam<calendar.size();iteam++) {

      std::fill(other_teams.begin(), other_teams.end(), false);

      for (int iturn=0; iturn<calendar[iteam].size(); iturn++) {
        int opponent = calendar[iteam][iturn];

        //1st check: verify the matches are correctly reported in both locations
        if (calendar[opponent][iturn] != iteam) {
          string msg = "Error 1 in validating the calendar!";
          ExceptionHandler::getIstance().raise(msg);
        }

        //2nd check: verify each team is played against only once
        int index_other_teams = opponent>=iteam ? opponent-1 : opponent;

        if (other_teams[index_other_teams]) {
          string msg = "Error 2 in validating the calendar!";
          ExceptionHandler::getIstance().raise(msg);
        }
        other_teams[index_other_teams] = true;
      }

      for (int jteam; jteam<other_teams.size(); jteam++) {
        if (!other_teams[jteam]) {
          string msg = "Error 3 in validating the calendar!";
          ExceptionHandler::getIstance().raise(msg);
        }
      }
    }

    cout << "Calendar validated!" << endl;

#endif


    //Copy our calendar to the vector
    for (int iturn=0; iturn<Nturns; iturn++) {
      my_matches[iturn] = calendar[0][iturn];
    }


    //Now transmit the calendar to the other nodes
    for (int inode=1; inode < Nmpi; inode++) {
      MPI_Send(&calendar[inode][0], Nturns, MPI_INT, inode, tag_league, MPI_COMM_WORLD);
    }

  }

  else {
    MPI_Status status;
    MPI_Recv(&my_matches[0], Nturns, MPI_INT, 0, tag_league, MPI_COMM_WORLD, &status);

  }


  return;
}



//=============================================================================
//  MpiControl::CreateInitialDomainDecomposition
/// Creates a binary tree containing all particles in order to determine how
/// to distribute the particles across all MPI nodes with an equal amount of
/// CPU work per MPI node.  If creating the initial partition (i.e. before
/// we have calculated the timestep), we give the particles equal weighting
/// and therefore each node will have equal numbers of particles.  For later
/// steps (i.e. when we know the timesteps and work information), we split
/// the domains to give each MPI node equal amounts of work.  This routine
/// should only be called for the root process.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim, ParticleType>::CreateInitialDomainDecomposition
(Sph<ndim> *sph,                    ///< Pointer to main SPH object
 Nbody<ndim> *nbody,                ///< Pointer to main N-body object
 Parameters *simparams,             ///< Simulation parameters
 DomainBox<ndim> simbox,            ///< Simulation domain box
 bool& initial_h_provided)          ///< Receives from root whether or not initial h was provided
{
  int i;                            // Particle counter
  int inode;                        // Node counter
  int k;                            // Dimension counter
  int okflag;                       // ..
  FLOAT boxbuffer[2*ndim*Nmpi];     // Bounding box buffer
  ParticleType<ndim> *partbuffer;   // ..

  //Broadcast whether or not the initial h was provided
  int initial_h = initial_h_provided;
  MPI_Bcast(&initial_h,1,MPI_INT,0,MPI_COMM_WORLD);
  initial_h_provided=initial_h;


  // For main process, create load balancing tree, transmit information to all
  // other nodes including particle data
  //===========================================================================
  if (rank == 0) {

    debug2("[MpiControl::CreateInitialDomainDecomposition]");

    // Create MPI binary tree for organising domain decomposition
    mpitree = new MpiTree<ndim,ParticleType>(Nmpi);

    // Set number of tree members to total no. of SPH particles (inc. ghosts)
    mpitree->Nsph = sph->Nsph;
    mpitree->Ntot = sph->Nsph;
    mpitree->Ntotmax = sph->Nsphmax;


    // Create all other MPI node objects
    this->AllocateMemory(mpitree->Ntotmax);

    // Get pointer to sph particles and cast it to the right type
    ParticleType<ndim>* sphdata =
      static_cast<ParticleType<ndim>* > (sph->GetParticlesArray());


    // For periodic simulations, set bounding box of root node to be the
    // periodic box size.  Otherwise, set to be the particle bounding box.
    if (simbox.x_boundary_lhs == "open") mpibox.boxmin[0] = -big_number;
    else mpibox.boxmin[0] = simbox.boxmin[0];
    if (simbox.x_boundary_rhs == "open") mpibox.boxmax[0] = big_number;
    else mpibox.boxmax[0] = simbox.boxmax[0];
    if (ndim > 1) {
      if (simbox.y_boundary_lhs == "open") mpibox.boxmin[1] = -big_number;
      else mpibox.boxmin[1] = simbox.boxmin[1];
      if (simbox.y_boundary_rhs == "open") mpibox.boxmax[1] = big_number;
      else mpibox.boxmax[1] = simbox.boxmax[1];
    }
    if (ndim == 3) {
      if (simbox.z_boundary_lhs == "open") mpibox.boxmin[2] = -big_number;
      else mpibox.boxmin[2] = simbox.boxmin[2];
      if (simbox.z_boundary_rhs == "open") mpibox.boxmax[2] = big_number;
      else mpibox.boxmax[2] = simbox.boxmax[2];
    }
    //mpitree->box = &mpibox;

    cout << "Simulation bounding box" << endl;
    for (k=0; k<ndim; k++)
      cout << "r[" << k << "]  :  " << mpibox.boxmin[k] << "   "
           << mpibox.boxmax[k] << endl;


    // Compute the size of all tree-related arrays now we know number of points
    mpitree->ComputeTreeSize();

    // Allocate (or reallocate if needed) all tree memory
    mpitree->AllocateMemory();

    // Create tree data structure including linked lists and cell pointers
    mpitree->CreateTreeStructure(mpinode);

    // Set properties for root cell before constructing tree
    mpitree->tree[0].N = sph->Nsph;
    mpitree->tree[0].ifirst = 0;
    mpitree->tree[0].ilast = sph->Nsph - 1;
    for (k=0; k<ndim; k++) mpitree->tree[0].bbmin[k] = mpibox.boxmin[k];
    for (k=0; k<ndim; k++) mpitree->tree[0].bbmax[k] = mpibox.boxmax[k];
    for (i=0; i<sph->Nsph; i++) mpitree->inext[i] = -1;
    for (i=0; i<sph->Nsph-1; i++) mpitree->inext[i] = i + 1;
    for (i=0; i<sph->Nsph; i++) mpitree->ids[i] = i;

    // Recursively divide tree up until we've reached bottom level
    mpitree->DivideTreeCell(0,mpitree->Ntot-1,sphdata,mpitree->tree[0]);


    // Copy details from each tree leaf cell into
    //-------------------------------------------------------------------------
    for (inode=0; inode<Nmpi; inode++) {
      int icell = mpitree->g2c[inode];

      // Create bounding boxes containing particles in each sub-tree
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmin[k] =
        mpitree->tree[icell].bbmin[k];
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmax[k] =
        mpitree->tree[icell].bbmax[k];

      cout << "CHECKING MPITREE : " << inode << "   " << icell << "   "
           << &mpitree->tree[icell] << endl;

      // Copy particle ids to node lists
      mpinode[inode].Nsph = 0;
      i = mpitree->tree[icell].ifirst;

      while (i != -1) {
        mpinode[inode].ids[mpinode[inode].Nsph++] = i;
        if (i == mpitree->tree[icell].ilast) break;
        i = mpitree->inext[i];
      };
      mpinode[inode].Ntot = mpinode[inode].Nsph;

      cout << "MPIDOMAIN : " << inode << "   Nsph : " << mpinode[inode].Nsph
	         << "    box : " << mpinode[inode].domain.boxmin[0]
	         << "     " << mpinode[inode].domain.boxmax[0] << endl;
    }
    //-------------------------------------------------------------------------


    // Pack all bounding box data into single array
    for (inode=0; inode<Nmpi; inode++) {
      for (k=0; k<ndim; k++)
        boxbuffer[2*ndim*inode + k] = mpinode[inode].domain.boxmin[k];
      for (k=0; k<ndim; k++)
        boxbuffer[2*ndim*inode + ndim + k] = mpinode[inode].domain.boxmax[k];
    }

    // Now broadcast all bounding boxes to other processes
    MPI_Bcast(boxbuffer,2*ndim*Nmpi,MPI_DOUBLE,0,MPI_COMM_WORLD);


    cout << "CHECKING SENDPARTICLES : " << mpinode[inode].ids
         << "   " << sphdata << endl;


    // Send particles to all other domains
    for (inode=1; inode<Nmpi; inode++) {
      SendParticles(inode,mpinode[inode].Nsph,mpinode[inode].ids,sphdata);
      cout << "Sent " << mpinode[inode].Nsph
           << " particles to node " << inode << endl;
    }

    cout << "Sent all particles to other processes" << endl;

    // Delete all other particles from local domain
    sph->Nsph = mpinode[0].Nsph;
    partbuffer = new ParticleType<ndim>[sph->Nsph];
    for (i=0; i<sph->Nsph; i++) partbuffer[i] = sphdata[mpinode[0].ids[i]];
    for (i=0; i<sph->Nsph; i++) sphdata[i] = partbuffer[i];
    delete[] partbuffer;
    cout << "Deleted all other particles from root node" << endl;

  }

  // For other nodes, receive all bounding box and particle data once
  // transmitted by main process.
  //===========================================================================
  else {

    // Create MPI node objects
    this->AllocateMemory(sph->Nsph);

    // Receive bounding box data for domain and unpack data
    MPI_Bcast(boxbuffer,2*ndim*Nmpi,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // Unpack all bounding box data
    for (inode=0; inode<Nmpi; inode++) {
      for (k=0; k<ndim; k++)
        mpinode[inode].domain.boxmin[k] = boxbuffer[2*ndim*inode + k];
      for (k=0; k<ndim; k++)
        mpinode[inode].domain.boxmax[k] = boxbuffer[2*ndim*inode + ndim + k];
      if (rank == 1) {
        cout << "Node " << inode << "    rank : " << rank << endl;
        cout << "xbox : " << mpinode[inode].domain.boxmin[0]
             << "    " << mpinode[inode].domain.boxmax[0] << endl;
        if (ndim > 1)
          cout << "ybox : " << mpinode[inode].domain.boxmin[1] << "    "
               << mpinode[inode].domain.boxmax[1] << endl;
        if (ndim == 3)
          cout << "zbox : " << mpinode[inode].domain.boxmin[2] << "    "
               << mpinode[inode].domain.boxmax[2] << endl;
      }
    }

    // Now, receive particles form main process and copy to local main array
    ReceiveParticles(0, (sph->Nsph), &partbuffer);

    sph->AllocateMemory(sph->Nsph);
    mpinode[rank].Nsph = sph->Nsph;

    // Get pointer to sph particles and cast it to the right type
    ParticleType<ndim>* sphdata =
      static_cast<ParticleType<ndim>* > (sph->GetParticlesArray());


    cout << "Received particles on node " << rank
         << "   Nsph : " << sph->Nsph << endl;

    for (i=0; i<sph->Nsph; i++) sphdata[i] = partbuffer[i];
    delete[] partbuffer;
    cout << "Deallocated partbuffer" << endl;

  }
  //===========================================================================


  return;
}



//=============================================================================
//  MpiControl::UpdateAllBoundingBoxes
/// Update local copy of all bounding boxes from all other MPI domains.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::UpdateAllBoundingBoxes
(int Npart,                         ///< No. of SPH particles
 Sph<ndim> *sph,        ///< Pointer to SPH data
 SphKernel<ndim> *kernptr)          ///< Pointer to kernel object
{
  int inode;                        // MPI node counter

  if (rank == 0) debug2("[MpiControl::UpdateAllBoundingBoxes]");

  // Update local bounding boxes
  mpinode[rank].UpdateBoundingBoxData(Npart,sph,kernptr);

  // Do an all_gather to receive the new array
  MPI_Allgather(&mpinode[rank].hbox,1,box_type,&boxes_buffer[0],
                1,box_type,MPI_COMM_WORLD);

  // Save the information inside the nodes
  for (inode=0; inode<Nmpi; inode++) {
    mpinode[inode].hbox = boxes_buffer[inode];
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Do an all_gather to receive the new array
  MPI_Allgather(&mpinode[rank].rbox,1,box_type,&boxes_buffer[0],
                1,box_type,MPI_COMM_WORLD);

  // Save the information inside the nodes
  for (inode=0; inode<Nmpi; inode++) {
    mpinode[inode].rbox = boxes_buffer[inode];
  }

  return;
}



//=============================================================================
//  MpiControl::LoadBalancing
/// If we are on a load balancing step, then determine which level of
/// the binary partition we are adjusting for load balancing.  Next, adjust
/// the domain boundaries at that level (and for all child domains).
/// Then broadcast the new domain boundaries to all other nodes to determine
/// which particles should be transfered to new nodes.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim, ParticleType >::LoadBalancing
(Sph<ndim> *sph,                    ///< Pointer to main SPH object
 Nbody<ndim> *nbody)                ///< Pointer to main N-body object
{
  int c;                            // MPI tree cell counter
  int c1;                           // ..
  int c2;                           // ..
  int i;                            // Particle counter
  int inode;                        // MPI node counter
  int k;                            // Dimension counter
  int kk;                           // ..
  int l;                            // ..
  int okflag;                       // Successful communication flag
  FLOAT rnew;                       // New boundary position for load balancing
  FLOAT boxbuffer[2*ndim*Nmpi];     // Bounding box buffer
  FLOAT workbuffer[1+ndim+Nmpi];    // Node work information buffer
  DOUBLE worktot = 0.0;             // Total work on all nodes
  MPI_Status status;                // MPI status flag

  // If running on only one MPI node, return immediately
  if (Nmpi == 1) return;

  //Get pointer to sph particles and cast it to the right type
  ParticleType<ndim>* sphdata =
    static_cast<ParticleType<ndim>* > (sph->GetParticlesArray());


  // Sum-up total work on all MPI nodes
  worktot = 0.0;
  for (inode=0; inode<Nmpi; inode++) worktot += 0.0;


  // Loop over all non-leaf levels (starting from the root) and compute the
  // work on each MPI-tree division, adjusting the divide to balance the
  // work on both sides.
  //---------------------------------------------------------------------------
  for (l=0; l<mpitree->ltot-1; l++) {

    // Loop over all MPI tree cells on current balancing level
    for (c=0; c<mpitree->Ncell; c++) {
      if (mpitree->tree[c].level != l) continue;
      c2 = mpitree->tree[c].c2;

      // Now find new division between child cells that is load-balanced
      mpitree->tree[c].r_divide = neibsearch->FindLoadBalancingDivision
        (mpitree->tree[c].k_divide,mpitree->tree[c].r_divide,
         mpitree->tree[c].bbmin,mpitree->tree[c].bbmax,
         mpitree->tree[c+1].nodes,mpitree->tree[c2].nodes,mpinode);

    }

    // Update all cell bounding boxes now new divisions have been computed
    mpitree->UpdateBoundingBoxes();

    // Finally update all MPI node bounding boxes
    for (inode=0; inode<Nmpi; inode++) {
      int icell = mpitree->g2c[inode];

      // Create bounding boxes containing particles in each sub-tree
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmin[k] =
        mpitree->tree[icell].bbmin[k];
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmax[k] =
        mpitree->tree[icell].bbmax[k];
    }

  }
  //---------------------------------------------------------------------------

  // Prepare lists of particles that now occupy other processor domains that
  // need to be transfered

  // First construct the list of nodes that we might be sending particles to
  std::vector<int> potential_nodes;
  potential_nodes.reserve(Nmpi);
  for (int inode=0; inode<Nmpi; inode++) {
    if (inode == rank) continue;
    if (BoxesOverlap(mpinode[inode].domain,mpinode[rank].rbox)) {
      potential_nodes.push_back(inode);
    }
  }



  // Find the ptcls that need to be transferred - delegate to NeighbourSearch
  std::vector<std::vector<int> > particles_to_transfer (Nmpi);
  std::vector<int> all_particles_to_export;
  BruteForceSearch<ndim,ParticleType> bruteforce(sph->kernp->kernrange,
                                                 &mpibox,sph->kernp,timing);
  bruteforce.FindParticlesToTransfer(sph, particles_to_transfer,
                                     all_particles_to_export,
                                     potential_nodes, mpinode);

  // Send and receive particles from/to all other nodes
  std::vector<ParticleType<ndim> > sendbuffer, recvbuffer;
  for (int iturn = 0; iturn<my_matches.size(); iturn++) {
    int inode = my_matches[iturn];

    int N_to_transfer = particles_to_transfer[inode].size();
    cout << "Transfer!!  Rank : " << rank << "    N_to_transfer : "
         << N_to_transfer << "    dest : " << inode << endl;
    sendbuffer.clear(); sendbuffer.resize(N_to_transfer);
    recvbuffer.clear();

    // Copy particles into send buffer
    for (int ipart=0; ipart < N_to_transfer; ipart++) {
      int index = particles_to_transfer[inode][ipart];
      sendbuffer[ipart] = sphdata[index];
    }

    // Do the actual send and receive

    //Decide if we have to send or receive first
    bool send_turn;
    if (rank < inode) {
      send_turn=true;
    }
    else {
      send_turn=false;
    }

    //Do the actual communication, sending and receiving in the right order
    for (int i=0; i < 2; i++) {
      if (send_turn) {
        cout << "Sending " << N_to_transfer << " from " << rank << " to " << inode << endl;
        MPI_Send(&sendbuffer[0], N_to_transfer, particle_type, inode,
                 tag_bal, MPI_COMM_WORLD);
        send_turn = false;
        cout << "Sent " << N_to_transfer << " from " << rank << " to " << inode << endl;
      }
      else {
        int N_to_receive;
        MPI_Status status;
        MPI_Probe(inode, tag_bal, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, particle_type, &N_to_receive);
        recvbuffer.resize(N_to_receive);
        cout << "Rank " << rank << " receiving " << N_to_receive << " from " << inode << endl;
        if (sph->Nsph+N_to_receive > sph->Nsphmax) {
          cout << "Memory problem : " << rank << " " << sph->Nsph
               << " " << N_to_receive << " " << sph->Nsphmax <<endl;
          string message = "Not enough memory for transfering particles";
          ExceptionHandler::getIstance().raise(message);
        }
        MPI_Recv(&recvbuffer[0], N_to_receive, particle_type, inode,
                 tag_bal, MPI_COMM_WORLD, &status);
        send_turn = true;
        cout << "Rank " << rank << " received " << N_to_receive << " from " << inode << endl;
      }
    }

    // Copy particles from receive buffer to main arrays
    int running_counter = sph->Nsph;
    // TODO: check we have enough memory
    for (int i=0; i< recvbuffer.size(); i++) {
      sphdata[running_counter] = recvbuffer[i];
      running_counter++;
    }
    sph->Nsph = running_counter;

  }


  // Remove transferred particles

  for (int i=0; i<all_particles_to_export.size(); i++) {
    sphdata[all_particles_to_export[i]].itype=dead;
  }

  sph->DeleteDeadParticles();


  return;

}


//=============================================================================
//  MpiControlType::ExportParticlesBeforeForceLoop
/// Export the particles that need force contribution from other processors
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim,ParticleType>::ExportParticlesBeforeForceLoop (Sph<ndim>* sph) {

  //Get the information to send to the other processors
  vector<char> send_buffer;

  for (int Nproc=0; Nproc<Nmpi; Nproc++) {
    //No need to send anything to ourselves
    if (Nproc == rank)
      continue;
    //Append at the end of send_vector the information we are sending and get how much it is
    Nbytes_to_each_proc[Nproc] = neibsearch->GetExportInfo(Nproc, sph, send_buffer, mpinode[Nproc],rank,Nmpi);
  }
  int Nbytes_to_be_exported = send_buffer.size();
  assert(std::accumulate(Nbytes_to_each_proc.begin(),Nbytes_to_each_proc.end(),0)==Nbytes_to_be_exported);

  //First need to know how many bytes each processor is sending
  MPI_Alltoall(&Nbytes_to_each_proc[0],1,MPI_INT,&Nbytes_exported_from_proc[0],1,MPI_INT,MPI_COMM_WORLD);

  //Can now compute the displacements
  vector<int> displs_recv(Nmpi), displs_send(Nmpi);
  int running_counter_recv=0, running_counter_send=0;
  for (int inode=1; inode<Nmpi; inode++) {
    running_counter_recv += Nbytes_exported_from_proc[inode-1];
    running_counter_send += Nbytes_to_each_proc[inode-1];
    displs_recv[inode]=running_counter_recv;
    displs_send[inode]=running_counter_send;
  }

  //Compute total number of particles to be received (including the ones from ourselves)
  int Nbytes_received_exported = std::accumulate(Nbytes_exported_from_proc.begin(),Nbytes_exported_from_proc.end(),0);

  //Allocate memory to receive all the particles
  vector<char > receive_buffer(Nbytes_received_exported);

 //Perform the actual communication
  MPI_Alltoallv(&send_buffer[0],&Nbytes_to_each_proc[0],&displs_send[0],MPI_CHAR,&receive_buffer[0],
      &Nbytes_exported_from_proc[0],&displs_recv[0],MPI_CHAR, MPI_COMM_WORLD);

  //Unpack the received arrays
  neibsearch->UnpackExported(receive_buffer, Nbytes_exported_from_proc, sph);

}


//=============================================================================
//  MpiControlType::GetExportedParticlesAccelerations
/// Get back the information about the exported particles from the other processors
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim,ParticleType>::GetExportedParticlesAccelerations (Sph<ndim>* sph) {

  vector<char> send_buffer;

  //Get the array with the acceleration for every other processor
  //Quite confusingly, note that the role of the two arays with the number of bytes from/to
  //each processor is now reversed (from is to send and to is to receive)
  neibsearch->GetBackExportInfo(send_buffer, Nbytes_exported_from_proc,
                                Nbytes_to_each_proc, sph, rank);

  vector<int> send_displs(Nmpi);
  compute_displs (send_displs, Nbytes_exported_from_proc);

  //Allocate receive buffer
  const int Nbytes_to_receive = std::accumulate(Nbytes_to_each_proc.begin(),
                                                Nbytes_to_each_proc.end(),0);
  vector<char > receive_buffer(Nbytes_to_receive);

  //Prepare receive counts and displacements
  vector<int> recv_displs(Nmpi);
  compute_displs (recv_displs, Nbytes_to_each_proc);

  //Do the actual communication
  MPI_Alltoallv(&send_buffer[0], &Nbytes_exported_from_proc[0],
                &send_displs[0], MPI_CHAR, &receive_buffer[0],
                &Nbytes_to_each_proc[0], &recv_displs[0],
                MPI_CHAR, MPI_COMM_WORLD );

  //Unpack the received information to update accelerations
  neibsearch->UnpackReturnedExportInfo(receive_buffer, recv_displs, sph, rank);

}



//=============================================================================
//  MpiControl::SendReceiveGhosts
/// Compute particles to send to other nodes and receive needed particles
/// from other nodes.  Returns the number of ghosts received. The array with
/// the received particles is stored inside the pointer returned. The caller
/// must NOT delete this array when is finished with it, as the memory is
/// internally managed by this class.
//=============================================================================
template <int ndim, template<int> class ParticleType>
int MpiControlType<ndim, ParticleType>::SendReceiveGhosts
(const FLOAT tghost,                ///< Ghost 'lifetime'
 Sph<ndim>* sph,                    ///< Main SPH object pointer
 ParticleType<ndim>** array)        ///< Main SPH particle array
{
  int i;                            // Particle counter
  int iaux;                         // ..
  int index;                        // ..
  int iparticle;                    // ..
  int j;                            // ..
  int Nexport;                      // ..
  int Npart;                        // ..
  int running_counter;              // ..
  int inode;                        // MPI node index
  vector<int> overlapping_nodes;    // List of nodes that overlap this domain
  vector<int> ghost_export_list;    // List of particles ids to be exported
  ParticleType<ndim>* sphdata = static_cast<ParticleType<ndim>* > (sph->GetParticlesArray());

  if (rank == 0) debug2("[MpiControl::SendReceiveGhosts]");

  // Reserve space for all nodes
  overlapping_nodes.reserve(Nmpi);

  // Loop over domains and find the ones that could overlap us
  for (inode=0; inode<Nmpi; inode++) {
    if (inode == rank) continue;
    if (BoxesOverlap(mpinode[inode].hbox,mpinode[rank].hbox)) {
      overlapping_nodes.push_back(inode);
    }
  }

  // Clear the buffer holding the indexes of the particles we need to export
  for (inode=0; inode<particles_to_export_per_node.size(); inode++) {
    particles_to_export_per_node[inode].clear();
  }


  for (iaux=0; iaux<overlapping_nodes.size(); iaux++) {
    inode = overlapping_nodes[iaux];
    ghost_export_list.clear();
    Nexport = neibsearch->SearchMpiGhostParticles(tghost,mpinode[inode].domain,
  						  sph,ghost_export_list);
    for (j=0; j<ghost_export_list.size(); j++) {
      i = ghost_export_list[j];
      particles_to_export_per_node[inode].push_back(&sphdata[i]);
    }
  }


  //BruteForceSearch<ndim,ParticleType> bruteforce(sph->kernp->kernrange,
  //			                 &mpibox,sph->kernp,timing);
  //bruteforce.FindGhostParticlesToExport(sph,particles_to_export_per_node,
  //                                overlapping_nodes,mpinode);


  // Prepare arrays with no. of particles to export per node and displacements
  fill(num_particles_export_per_node.begin(),
       num_particles_export_per_node.end(),0);

  running_counter = 0;
  for (inode=0; inode<Nmpi; inode++) {
    Npart = particles_to_export_per_node[inode].size();
    num_particles_export_per_node[inode] = Npart;
    displacements_send[inode] = running_counter;
    running_counter += Npart;
  }

  // Compute total number of particles to export
  Nexport = std::accumulate(num_particles_export_per_node.begin(),
                            num_particles_export_per_node.end(),0);

  // Comunicate with all processors the no. of ptcls that everyone is exporting
  vector<int> ones(Nmpi,1);
  vector<int> displs(Nmpi);
  for (inode=0; inode<displs.size(); inode++) {
    displs[inode] = inode;
  }
  MPI_Alltoallv(&num_particles_export_per_node[0], &ones[0], &displs[0],
                MPI_INT, &num_particles_to_be_received[0], &ones[0],
                &displs[0], MPI_INT, MPI_COMM_WORLD);

  // Compute the total number of particles that will be received
  tot_particles_to_receive = accumulate(num_particles_to_be_received.begin(),
                                        num_particles_to_be_received.end(),0);

  // Allocate receive buffer
  particles_receive.clear();
  particles_receive.resize(tot_particles_to_receive);

  cout << "TOT_PARTICLES_TO_EXPORT : " << Nexport << endl;

  // Create vector containing all particles to export
  particles_to_export.resize(Nexport);

  index = 0;
  for (inode=0; inode<Nmpi; inode++) {
    std::vector<ParticleType<ndim>* >& particles_on_this_node = particles_to_export_per_node[inode];
    for (iparticle=0; iparticle<particles_on_this_node.size(); iparticle++) {
      particles_to_export[index] = *particles_on_this_node[iparticle];
      index++;
    }
  }

  assert(index==Nexport);

  // Compute receive displacements
  running_counter = 0;
  for (inode=0; inode<receive_displs.size(); inode++) {
    receive_displs[inode] = running_counter;
    running_counter += num_particles_to_be_received[inode];
  }

  // Send and receive particles
  MPI_Alltoallv(&particles_to_export[0], &num_particles_export_per_node[0],
                &displacements_send[0], particle_type, &particles_receive[0],
                &num_particles_to_be_received[0], &receive_displs[0],
                particle_type, MPI_COMM_WORLD);

  *array = &particles_receive[0];

  return tot_particles_to_receive;

}



//=============================================================================
//  MpiControl::UpdateGhostParticles
/// Update the ghost particles that were previously found. Return the number
/// of Mpi Ghosts; modified the passed pointer with the address of the array
/// of ghosts (note: the caller must NOT call delete on this array, as the
/// memory is internally managed by the MpiControl class).
//=============================================================================
template <int ndim, template<int> class ParticleType>
int MpiControlType<ndim, ParticleType>::UpdateGhostParticles
(ParticleType<ndim>** array)        ///< Main SPH particle array
{
  int index = 0;                    // ..
  int inode;                        // MPI node counter
  int ipart;                        // Particle counter

  // Update the local buffer of particles to send
  for (inode=0; inode<Nmpi; inode++) {
    vector<ParticleType<ndim>* >& particles_on_this_node =
      particles_to_export_per_node[inode];
    for (ipart=0; ipart<particles_on_this_node.size(); ipart++) {
      particles_to_export[index] = *particles_on_this_node[ipart];
      index++;
    }
  }

  // Send and receive particles
  MPI_Alltoallv(&particles_to_export[0], &num_particles_export_per_node[0],
                &displacements_send[0], particle_type, &particles_receive[0],
                &num_particles_to_be_received[0], &receive_displs[0],
                particle_type, MPI_COMM_WORLD);

  *array = &particles_receive[0];

  return tot_particles_to_receive;
}



//=============================================================================
//  MpiControl::SendParticles
/// Given an array of ids and a node, copy particles inside a buffer and
/// send them to the given node.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim, ParticleType>::SendParticles
(int Node,                          ///< ..
 int Nparticles,                    ///< ..
 int* list,                         ///< ..
 ParticleType<ndim>* main_array)    ///< main_array_gen
{
  int i;                            // Particle counter

  //Ensure there is enough memory in the buffer
  sendbuffer.resize(Nparticles);

  //Copy particles from the main arrays to the buffer
  for (i=0; i<Nparticles; i++) {
    sendbuffer[i] = main_array[list[i]];
  }

  MPI_Send(&sendbuffer[0], Nparticles, particle_type, Node,
          tag_srpart, MPI_COMM_WORLD);

  cout << "SENDING : " << Nparticles << "   to " << Node << endl;

  return;
}



//=============================================================================
//  MpiControl::ReceiveParticles
/// Given a node, receive particles from it. Return the number of particles
/// received and a pointer to the array containing the particles. The caller
/// is reponsible to free the array after its usage.
//=============================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim, ParticleType>::ReceiveParticles
(int Node,                          ///< ..
 int& Nparticles,                   ///< ..
 ParticleType<ndim>** array)        ///< pointer_array_gen
{
  const int tag = 1;
  MPI_Status status;

  // "Probe" the message to know how big the message is going to be
  MPI_Probe(Node, tag, MPI_COMM_WORLD, &status);

  // Get the number of elements
  MPI_Get_count(&status, particle_type, &Nparticles);

  cout << "NPARTICLES : " << Nparticles << "   node : " << Node << endl;

  // Allocate enough memory to hold the particles
  *array = new ParticleType<ndim> [Nparticles];

  // Now receive the message
  MPI_Recv(*array, Nparticles, particle_type, Node,
           tag_srpart, MPI_COMM_WORLD, &status);

  return;
}



//=============================================================================
//  MpiControl::CollateDiagnosticsData
/// ..
//=============================================================================
template <int ndim>
void MpiControl<ndim>::CollateDiagnosticsData(Diagnostics<ndim> &diag)
{
  int inode;
  int j;
  int k;
  Diagnostics<ndim> diagaux;
  MPI_Status status;

  //---------------------------------------------------------------------------
  if (rank == 0) {

    // First multiply root node values by mass
    for (k=0; k<ndim; k++) diag.rcom[k] *= diag.mtot;
    for (k=0; k<ndim; k++) diag.vcom[k] *= diag.mtot;

    for (inode=1; inode<Nmpi; inode++) {
      MPI_Recv(&diagaux, 1, diagnostics_type,
               inode, 0, MPI_COMM_WORLD, &status);
      diag.Nsph += diagaux.Nsph;
      diag.Nstar += diagaux.Nstar;
      diag.Etot += diagaux.Etot;
      diag.utot += diagaux.utot;
      diag.ketot += diagaux.ketot;
      diag.gpetot += diagaux.gpetot;
      diag.mtot += diagaux.mtot;
      for (k=0; k<ndim; k++) diag.mom[k] += diagaux.mom[k];
      for (k=0; k<3; k++) diag.angmom[k] += diagaux.angmom[k];
      for (k=0; k<ndim; k++) diag.force[k] += diagaux.force[k];
      for (k=0; k<ndim; k++) diag.force_hydro[k] += diagaux.force_hydro[k];
      for (k=0; k<ndim; k++) diag.force_grav[k] += diagaux.force_grav[k];
      for (k=0; k<ndim; k++) diag.rcom[k] += diagaux.mtot*diagaux.rcom[k];
      for (k=0; k<ndim; k++) diag.vcom[k] += diagaux.mtot*diagaux.vcom[k];
    }

    // Renormalise centre of mass positions and velocities
    for (k=0; k<ndim; k++) diag.rcom[k] /= diag.mtot;
    for (k=0; k<ndim; k++) diag.vcom[k] /= diag.mtot;

  }
  //---------------------------------------------------------------------------
  else {

    MPI_Send(&diag, 1, diagnostics_type, 0, 0, MPI_COMM_WORLD);

  }
  //---------------------------------------------------------------------------

  return;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class MpiControl<1>;
template class MpiControl<2>;
template class MpiControl<3>;

// Template class for each particle type
template class MpiControlType<1, GradhSphParticle>;
template class MpiControlType<2, GradhSphParticle>;
template class MpiControlType<3, GradhSphParticle>;
template class MpiControlType<1, SM2012SphParticle>;
template class MpiControlType<2, SM2012SphParticle>;
template class MpiControlType<3, SM2012SphParticle>;
template class MpiControlType<1, GodunovSphParticle>;
template class MpiControlType<2, GodunovSphParticle>;
template class MpiControlType<3, GodunovSphParticle>;
