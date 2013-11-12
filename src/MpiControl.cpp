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
#include "Constants.h"
#include "Precision.h"
#include "SphKernel.h"
#include "DomainBox.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
#include "MpiControl.h"
using namespace std;



//=============================================================================
//  MpiControl::MpiControl()
/// MPI node class constructor.
//=============================================================================
template <int ndim>
MpiControl<ndim>::MpiControl()
{
  int len;

  allocated_mpi = false;

  MPI_Comm_size(MPI_COMM_WORLD,&Nmpi);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Get_processor_name(hostname, &len);

  if (this->rank == 0)
    printf("MPI working.  Nmpi : %d   rank : %d   hostname : %s\n",rank,Nmpi,hostname);
  else
    printf("%d is running too!!\n",this->rank);

  //Create and commit the particle datatype
  particle_type = SphParticle<ndim>::CreateMpiDataType();
  MPI_Type_commit(&particle_type);

#ifdef VERIFY_ALL
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
template <int ndim>
MpiControl<ndim>::~MpiControl()
{
  MPI_Type_free(&particle_type);
}



//=============================================================================
//  MpiControl::AllocateMemory
/// Allocate all memory for MPI control class.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::AllocateMemory(void)
{

  mpinode = new MpiNode<ndim>[Nmpi];

  return;
}



//=============================================================================
//  MpiControl::DeallocateMemory
/// Deallocate all MPI control class memory.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::DeallocateMemory(void)
{

  delete[] mpinode;

  return;
}



//=============================================================================
//  MpiControl::InitialiseMpiProcess
/// Call all initial MPI routines, to find rank number, no. of processes etc..
//=============================================================================
template <int ndim>
void MpiControl<ndim>::InitialiseMpiProcess(void)
{
  debug2("[MpiControl::InitialiseMpiProcess]");


  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Abort(MPI_COMM_WORLD,0);

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
template <int ndim>
void MpiControl<ndim>::CreateInitialDomainDecomposition
(Sph<ndim> *sph,                   ///< Pointer to main SPH object
 Nbody<ndim> *nbody,               ///< Pointer to main N-body object
 Parameters *simparams,            ///< Simulation parameters
 DomainBox<ndim> simbox)           ///< Simulation domain box
{
  int i;                           // Particle counter
  int k;                           // Dimension counter
  int okflag;                      // ..
  FLOAT boxbuffer[2*ndim*Nmpi];    // Bounding box buffer
  SphParticle<ndim> *partbuffer;   // ..


  // For main process, create load balancing tree, transmit information to all
  // other nodes including particle data
  //---------------------------------------------------------------------------
  if (rank == 0) {

    debug2("[MpiControl::CreateLoadBalancingTree]");

    // Create MPI binary tree for organising domain decomposition
    mpitree = new BinaryTree<ndim>(16,0.1,0.0,"geometric","monopole",1,Nmpi);

    // Create all other MPI node objects
    AllocateMemory();

    // Create binary tree from all SPH particles
    // Set number of tree members to total number of SPH particles (inc. ghosts)
    mpitree->Nsph = sph->Nsph;
    mpitree->Ntot = sph->Nsph;
    mpitree->Ntotmax = max(mpitree->Ntot,mpitree->Ntotmax);
    mpitree->gtot = 0;
    for (i=0; i<sph->Nsph; i++)
      for (k=0; k<ndim; k++) sph->rsph[ndim*i + k] = sph->sphdata[i].r[k];


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
    mpitree->box = &mpibox;


    // Compute the size of all tree-related arrays now we know number of points
    mpitree->ComputeTreeSize();

    // Allocate (or reallocate if needed) all tree memory
    mpitree->AllocateTreeMemory();

    // Create tree data structure including linked lists and cell pointers
    mpitree->CreateTreeStructure();

    // Find ordered list of ptcl positions ready for adding particles to tree
    mpitree->OrderParticlesByCartCoord(sph->sphdata);

    // Now add particles to tree depending on Cartesian coordinates
    mpitree->LoadParticlesToTree(sph->rsph);

    // Create bounding boxes containing particles in each sub-tree
    for (i=0; i<Nmpi; i++) {
      for (k=0; k<ndim; k++) mpinode[i].domain.boxmin[k] = mpitree->subtrees[i]->box.boxmin[k];
      for (k=0; k<ndim; k++) mpinode[i].domain.boxmax[k] = mpitree->subtrees[i]->box.boxmax[k];
    }


    // Pack all bounding box data into single array
    for (i=0; i<Nmpi; i++) {
      for (k=0; k<ndim; k++) boxbuffer[2*ndim*i + k] = mpinode[i].domain.boxmin[k];
      for (k=0; k<ndim; k++) boxbuffer[2*ndim*i + ndim + k] = mpinode[i].domain.boxmax[k];
    }

    // Now broadcast all bounding boxes to other processes
    MPI_Bcast(boxbuffer,2*ndim*Nmpi,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Finally, send particles to all other domains
    for (i=1; i<Nmpi; i++)
      SendParticles(i, mpitree->subtrees[i]->Nsph, mpitree->subtrees[i]->ids, sph->sphdata);

  }

  // For other nodes, receive all bounding box and particle data once
  // transmitted by main process.
  //---------------------------------------------------------------------------
  else {

    // Create MPI node objects
    AllocateMemory();

    // Receive bounding box data for domain and unpack data
    MPI_Bcast(boxbuffer,2*ndim*Nmpi,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Unpack all bounding box data
    for (i=0; i<Nmpi; i++) {
      for (k=0; k<ndim; k++) mpinode[i].domain.boxmin[k] = boxbuffer[2*ndim*i + k];
      for (k=0; k<ndim; k++) mpinode[i].domain.boxmax[k] = boxbuffer[2*ndim*i + ndim + k];
      //if (rank == 1) {
      //cout << "Node " << i << endl;
      //cout << "xbox : " << mpinode[i].domain.boxmin[0] << "    " << mpinode[i].domain.boxmax[0] << endl;
      //cout << "ybox : " << mpinode[i].domain.boxmin[1] << "    " << mpinode[i].domain.boxmax[1] << endl;
      //cout << "zbox : " << mpinode[i].domain.boxmin[2] << "    " << mpinode[i].domain.boxmax[2] << endl;
      //}
    }

    // Now, receive particles form main process
    ReceiveParticles(rank, &mpinode[rank].Nsph, partbuffer);

  }
  //---------------------------------------------------------------------------

  MPI_Barrier(MPI_COMM_WORLD);


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
template <int ndim>
void MpiControl<ndim>::LoadBalancing
(Sph<ndim> *sph,                   ///< Pointer to main SPH object
 Nbody<ndim> *nbody)               ///< Pointer to main N-body object
{
  int i;                           // Particle counter
  int inode;                       // MPI node counter
  int k;                           // Dimension counter
  int okflag;                      // Successful communication flag
  MPI_Status status;               // ..


  // Compute work that will be transmitted to all other domains if using
  // current domain boxes and particle positions
  // (DAVID : Maybe in the future once basic implementation works)

  // Compute total work contained in domain at present
  mpinode[rank].worktot = 0.0;
  for (i=0; i<sph->Nsph; i++) mpinode[rank].worktot += 1.0/sph->sphdata[i].dt;


  // For root process, receive all current node CPU nodes, adjust domain
  // boundaries to balance out work more evenly and then broadcast new domain
  // boxes to all other nodes.
  //---------------------------------------------------------------------------
  if (rank == 0) {


    // Receive all important load balancing information from other nodes
    for (inode=1; inode<Nmpi; inode++) {
      okflag = MPI_Recv(&mpinode[inode].worktot,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
    }


    // Work out tree level at which we are altering the load balancing



    // Adjust bounding box sizes on load balancing level



    // Transmit new bounding box sizes to all other nodes
    MPI_Barrier(MPI_COMM_WORLD);

  }
  //---------------------------------------------------------------------------
  else {


    // Transmit load balancing information to main root node
    okflag = MPI_Send(&mpinode[rank].worktot,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);


    // Receive new bounding box information for all nodes
    MPI_Barrier(MPI_COMM_WORLD);


  }
  //---------------------------------------------------------------------------



  // Prepare lists of particles that now occupy other processor domains that 
  // need to be transfered



  // Send particles to all other nodes



  // Receive particles from all other nodes




  return;
}



//=============================================================================
//  MpiControl::TransferParticlesToNode
/// Once we know the new domain boundaries for all MPI nodes, transfer any 
/// particles that now lie in other domain boxes to those respective MPI 
/// nodes.  Also, receives particles from other domains.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::TransferParticlesToNode(void)
{
  return;
}


//==================================================================================
//  MpiControl::SendParticles
/// Given an array of ids and a node, copy particles inside a buffer and send them to
/// the given node
//==================================================================================
template <int ndim>
void MpiControl<ndim>::SendParticles(int Node, int Nparticles, int* list, SphParticle<ndim>* main_array) {

  const int tag = 1;

  //Ensure there is enough memory in the buffer
  sendbuffer.reserve(Nparticles);

  //Copy particles from the main arrays to the buffer
  for (int i=0; i<Nparticles; i++) {
    sendbuffer[i] = main_array[list[i]];
  }

  MPI_Send (&sendbuffer[0], Nparticles, particle_type, Node, tag, MPI_COMM_WORLD);

}

//==================================================================================
//  MpiControl::ReceiveParticles
/// Given a node, receive particles from it. Return the number of particles received
/// and a pointer to the array containing the particles. The caller is reponsible
/// to free the array after its usage
//==================================================================================
template <int ndim>
void MpiControl<ndim>::ReceiveParticles (int Node, int& Nparticles, SphParticle<ndim>* array) {
  const int tag = 1;
  MPI_Status status;
  //"Probe" the message to know how big the message is going to be
  MPI_Probe(Node, tag, MPI_COMM_WORLD, &status);

  //Get the number of elements
  MPI_Get_count( &status,  particle_type, &Nparticles );

  //Allocate enough memory to hold the particles
  array = new SphParticle<ndim> [Nparticles];

  //Now receive the message
  MPI_Recv(array, Nparticles, particle_type, Node, tag, MPI_COMM_WORLD, &status);

}
// Template class instances for each dimensionality value (1, 2 and 3)
template class MpiControl<1>;
template class MpiControl<2>;
template class MpiControl<3>;
