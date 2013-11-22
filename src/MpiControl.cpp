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
    printf("MPI working.  Nmpi : %d   rank : %d   hostname : %s\n",Nmpi,rank,hostname);
  else
    printf("%d is running too!!\n",this->rank);

  //Create and commit the particle datatype
  particle_type = SphParticle<ndim>::CreateMpiDataType();
  MPI_Type_commit(&particle_type);

  // Create diagnostics data structure in database
  diagnostics_type = Diagnostics<ndim>::CreateMpiDataType();
  MPI_Type_commit(&diagnostics_type);

  //Create and commit the box datatype
  Box<ndim> dummy;
  box_type = CreateBoxType(dummy);
  MPI_Type_commit(&box_type);
  //Allocate buffer to send and receive boxes
  boxes_buffer.resize(Nmpi);

  //Allocate the buffers needed to send and receive particles
  particles_to_export_per_node.resize(Nmpi);
  num_particles_export_per_node.resize(Nmpi);
  particles_to_export.resize(Nmpi);
  displacements_send.resize(Nmpi);
  num_particles_to_be_received.resize(Nmpi);
  receive_displs.resize(Nmpi);

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
  MPI_Type_free(&box_type);
  MPI_Type_free(&diagnostics_type);
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
  int inode;                       // Node counter
  int k;                           // Dimension counter
  int okflag;                      // ..
  FLOAT boxbuffer[2*ndim*Nmpi];    // Bounding box buffer
  SphParticle<ndim> *partbuffer;   // ..


  // For main process, create load balancing tree, transmit information to all
  // other nodes including particle data
  //---------------------------------------------------------------------------
  if (rank == 0) {

    debug2("[MpiControl::CreateInitialDomainDecomposition]");

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
    for (inode=0; inode<Nmpi; inode++) {
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmin[k] =
        mpitree->subtrees[inode]->box.boxmin[k];
      for (k=0; k<ndim; k++) mpinode[inode].domain.boxmax[k] =
        mpitree->subtrees[inode]->box.boxmax[k];
    }


    // Pack all bounding box data into single array
    for (inode=0; inode<Nmpi; inode++) {
      for (k=0; k<ndim; k++) boxbuffer[2*ndim*inode + k] = mpinode[inode].domain.boxmin[k];
      for (k=0; k<ndim; k++) boxbuffer[2*ndim*inode + ndim + k] = mpinode[inode].domain.boxmax[k];
    }

    // Now broadcast all bounding boxes to other processes
    MPI_Bcast(boxbuffer,2*ndim*Nmpi,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Send particles to all other domains
    for (inode=1; inode<Nmpi; inode++) {
      SendParticles(inode, mpitree->subtrees[inode]->Nsph,
                    mpitree->subtrees[inode]->ids, sph->sphdata);
      cout << "Sent " << mpitree->subtrees[inode]->Nsph
           << " particles to node " << inode << endl;
    }

    cout << "Sent all particles to other processes" << endl;

    // Delete all other particles from local domain
    sph->Nsph = mpitree->subtrees[0]->Nsph;
    partbuffer = new SphParticle<ndim>[sph->Nsph];
    for (i=0; i<sph->Nsph; i++) partbuffer[i] = sph->sphdata[mpitree->subtrees[0]->ids[i]];
    for (i=0; i<sph->Nsph; i++) sph->sphdata[i] = partbuffer[i];
    delete[] partbuffer;
    cout << "Deleted all other particles from root node" << endl;

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
    for (inode=0; inode<Nmpi; inode++) {
      for (k=0; k<ndim; k++)
        mpinode[inode].domain.boxmin[k] = boxbuffer[2*ndim*inode + k];
      for (k=0; k<ndim; k++)
        mpinode[inode].domain.boxmax[k] = boxbuffer[2*ndim*inode + ndim + k];
      if (rank == 1) {
      cout << "Node " << i << endl;
      cout << "xbox : " << mpinode[inode].domain.boxmin[0]
           << "    " << mpinode[inode].domain.boxmax[0] << endl;
      cout << "ybox : " << mpinode[inode].domain.boxmin[1] << "    "
           << mpinode[inode].domain.boxmax[1] << endl;
      cout << "zbox : " << mpinode[inode].domain.boxmin[2] << "    "
           << mpinode[inode].domain.boxmax[2] << endl;
      }
    }
    // Now, receive particles form main process and copy to local main array
    ReceiveParticles(0, (sph->Nsph), &partbuffer);

    sph->AllocateMemory(sph->Nsph);
    mpinode[rank].Nsph = sph->Nsph;

    cout << "Received particles on node " << rank << "   Nsph : " << sph->Nsph << endl;

    for (i=0; i<sph->Nsph; i++) sph->sphdata[i] = partbuffer[i];
    delete[] partbuffer;
    cout << "Dellocated partbuffer" << endl;

  }
  //---------------------------------------------------------------------------

  MPI_Barrier(MPI_COMM_WORLD);


  return;
}



//=============================================================================
//  MpiControl::UpdateAllBoundingBoxes
/// Update local copy of all bounding boxes from all other MPI domains.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::UpdateAllBoundingBoxes
(int Npart,                        ///< No. of SPH particles
 SphParticle<ndim> *sphdata,       ///< Pointer to SPH data
 SphKernel<ndim> *kernptr)         ///< Pointer to kernel object
{

  // Update local bounding boxes
  mpinode[rank].UpdateBoundingBoxData(Npart,sphdata,kernptr);

  // Do an all_gather to receive the new array
  MPI_Allgather(&mpinode[rank].hbox,1,box_type,&boxes_buffer[0],
                1,box_type,MPI_COMM_WORLD);

  // Save the information inside the nodes
  for (int i=0; i<Nmpi; i++) {
    mpinode[i].hbox = boxes_buffer[i];
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
template <int ndim>
void MpiControl<ndim>::LoadBalancing
(Sph<ndim> *sph,                   ///< Pointer to main SPH object
 Nbody<ndim> *nbody)               ///< Pointer to main N-body object
{
  int i;                           // Particle counter
  int inode;                       // MPI node counter
  int k;                           // Dimension counter
  int okflag;                      // Successful communication flag
  FLOAT boxbuffer[2*ndim*Nmpi];    // Bounding box buffer
  MPI_Status status;               // MPI status flag


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
      okflag = MPI_Recv(&mpinode[inode].worktot,1,MPI_DOUBLE,
                        0,0,MPI_COMM_WORLD,&status);
    }


    // Work out tree level at which we are altering the load balancing



    // Adjust bounding box sizes on load balancing level



    // Transmit new bounding box sizes to all other nodes
    for (inode=0; inode<Nmpi; inode++) {
      for (k=0; k<ndim; k++)
        boxbuffer[2*ndim*inode + k] = mpinode[inode].domain.boxmin[k];
      for (k=0; k<ndim; k++)
        boxbuffer[2*ndim*inode + ndim + k] = mpinode[inode].domain.boxmax[k];
    }

    // Now broadcast all bounding boxes to other processes
    MPI_Bcast(boxbuffer,2*ndim*Nmpi,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

  }
  //---------------------------------------------------------------------------
  else {


    // Transmit load balancing information to main root node
    okflag = MPI_Send(&mpinode[rank].worktot,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);


    // Receive bounding box data for domain and unpack data
    MPI_Bcast(boxbuffer,2*ndim*Nmpi,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // Unpack all bounding box data
    for (inode=0; inode<Nmpi; inode++) {
      for (k=0; k<ndim; k++) 
        mpinode[inode].domain.boxmin[k] = boxbuffer[2*ndim*inode + k];
      for (k=0; k<ndim; k++) 
        mpinode[inode].domain.boxmax[k] = boxbuffer[2*ndim*inode + ndim + k];
    }
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
//  MpiControl::SendReceiveGhosts
/// Compute particles to send to other nodes and receive needed particles 
/// from other nodes.  Returns the number of ghosts received. The array with 
/// the received particles is stored inside the pointer returned. The caller 
/// must NOT delete this array when is finished with it, as the memory is 
/// internally managed by this class.
//=============================================================================
template <int ndim>
int MpiControl<ndim>::SendReceiveGhosts
(SphParticle<ndim>** array,        ///< Main SPH particle array
 Sph<ndim>* sph)                   ///< Main SPH object pointer
{
  int i;                           // Particle counter
  int index;                       // ..
  int inode;                       // Node counter
  int running_counter;             // ..


  std::vector<int > overlapping_nodes;

  //Reserve space for all nodes
  overlapping_nodes.reserve(Nmpi);

  //Loop over domains and find the ones that could overlap us
  for (inode=0; inode<Nmpi; inode++) {
    if (inode == rank) continue;
    if (BoxesOverlap(mpinode[inode].hbox,mpinode[rank].hbox)) {
      overlapping_nodes.push_back(inode);
    }
  }

  //Clear the buffer holding the indexes of the particles we need to export
  for (inode=0; inode<particles_to_export_per_node.size(); inode++) {
    particles_to_export_per_node[inode].clear();
  }

  //Loop over particles and prepare the ones to export
  for (i=0; i<sph->Ntot; i++) {
    SphParticle<ndim>& part = sph->sphdata[i];

    //Loop over potential domains and see if we need to export this particle to them
    for (inode=0; inode<overlapping_nodes.size(); inode++) {
      int node_number = overlapping_nodes[inode];
      if (ParticleBoxOverlap(part,mpinode[node_number].hbox)) {
        particles_to_export_per_node[node_number].push_back(&part);
      }
    }
  }

  //Prepare arrays with number of particles to export per node and displacements
  std::fill(num_particles_export_per_node.begin(),num_particles_export_per_node.end(),0);

  running_counter = 0;
  for (int inode=0; inode<Nmpi; inode++) {
    int num_particles = particles_to_export_per_node[inode].size();
    num_particles_export_per_node[inode] = num_particles;
    displacements_send[inode]=running_counter;
    running_counter += num_particles;
  }

  //Compute total number of particles to export
  int tot_particles_to_export = std::accumulate(num_particles_export_per_node.begin(),num_particles_export_per_node.end(),0);

  //Comunicate with all the processors the number of particles that everyone is exporting
  std::vector<int> ones (Nmpi,1);
  std::vector<int> displs(Nmpi);
  for (inode=0; inode<displs.size(); inode++) {
    displs[inode] = inode;
  }
  MPI_Alltoallv(&num_particles_export_per_node[0], &ones[0], &displs[0], MPI_INT, &num_particles_to_be_received[0], &ones[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

  //Compute the total number of particles that will be received
  tot_particles_to_receive = std::accumulate(num_particles_to_be_received.begin(),num_particles_to_be_received.end(),0);

  //Allocate receive buffer
  particles_receive.clear();
  particles_receive.resize(tot_particles_to_receive);

  //Create vector containing all particles to export
  particles_to_export.resize(tot_particles_to_export);

  index = 0;
  for (inode=0; inode < Nmpi; inode++) {
    std::vector<SphParticle<ndim>* >& particles_on_this_node = particles_to_export_per_node[inode];
    for (int iparticle=0; iparticle<particles_on_this_node.size(); iparticle++) {
      particles_to_export[index] =  *particles_on_this_node[iparticle];
      index++;
    }
  }

  assert(index==tot_particles_to_export);

  //Compute receive displacements
  running_counter = 0;
  for (inode=0; inode<receive_displs.size(); inode++) {
    receive_displs[inode] = running_counter;
    running_counter += num_particles_to_be_received[inode];
  }

  //Send and receive particles
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
template <int ndim>
int MpiControl<ndim>::UpdateGhostParticles
(SphParticle<ndim>** array)         ///< Main SPH particle array
{
  int index = 0;                    // ..
  int inode;                        // MPI node counter
  int ipart;                        // Particle counter

  //Update the local buffer of particles to send
  for (inode=0; inode<Nmpi; inode++) {
    std::vector<SphParticle<ndim>* >& particles_on_this_node = particles_to_export_per_node[inode];
    for (ipart=0; ipart<particles_on_this_node.size(); ipart++) {
      particles_to_export[index] = *particles_on_this_node[ipart];
      index++;
    }
  }

  //Send and receive particles
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
template <int ndim>
void MpiControl<ndim>::SendParticles
(int Node, 
 int Nparticles, 
 int* list, 
 SphParticle<ndim>* main_array) 
{
  int i;                            // Particle counter
  const int tag = 1;

  //Ensure there is enough memory in the buffer
  sendbuffer.resize(Nparticles);

  //Copy particles from the main arrays to the buffer
  for (i=0; i<Nparticles; i++) {
    sendbuffer[i] = main_array[list[i]];
  }

  MPI_Send(&sendbuffer[0], Nparticles, particle_type, Node, 
           tag, MPI_COMM_WORLD);

  return;
}



//=============================================================================
//  MpiControl::ReceiveParticles
/// Given a node, receive particles from it. Return the number of particles 
/// received and a pointer to the array containing the particles. The caller 
/// is reponsible to free the array after its usage.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::ReceiveParticles
(int Node, 
 int& Nparticles, 
 SphParticle<ndim>** array) 
{

  const int tag = 1;
  MPI_Status status;

  //"Probe" the message to know how big the message is going to be
  MPI_Probe(Node, tag, MPI_COMM_WORLD, &status);

  //Get the number of elements
  MPI_Get_count(&status, particle_type, &Nparticles);

  //Allocate enough memory to hold the particles
  *array = new SphParticle<ndim> [Nparticles];

  //Now receive the message
  MPI_Recv(*array, Nparticles, particle_type, Node, 
           tag, MPI_COMM_WORLD, &status);

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
