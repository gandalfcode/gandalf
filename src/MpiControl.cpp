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
#include "mpi.h"
#include "Constants.h"
#include "Precision.h"
#include "SphKernel.h"
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
  allocated_mpi = false;
}



//=============================================================================
//  MpiControl::~MpiControl()
/// MPI node class destructor.
//=============================================================================
template <int ndim>
MpiControl<ndim>::~MpiControl()
{
}



//=============================================================================
//  MpiControl::AllocateMemory
/// Allocate all memory for MPI control class.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::AllocateMemory(void)
{
  return;
}



//=============================================================================
//  MpiControl::DeallocateMemory
/// Deallocate all MPI control class memory.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::DeallocateMemory(void)
{
  return;
}



//=============================================================================
//  MpiControl::InitialiseMpiProcess
/// Call all initial MPI routines, to find rank number, no. of processes etc..
//=============================================================================
template <int ndim>
void MpiControl<ndim>::InitialiseMpiProcess(void)
{
  int len;
  debug2("[MpiControl::InitialiseMpiProcess]");

  MPI_Comm_size(MPI_COMM_WORLD,&Nmpi);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Get_processor_name(hostname, &len);
  MPI_Barrier(MPI_COMM_WORLD);

  if (this->rank == 0)
    printf("MPI working.  Nmpi : %d   rank : %d   hostname : %s\n",rank,Nmpi,hostname);
  else
    printf("%d is running too!!\n",this->rank);

  //MPI_Finalize();
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Abort(MPI_COMM_WORLD,0);

  return;
}



//=============================================================================
//  MpiControl::CreateLoadBalancingTree
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
void MpiControl<ndim>::CreateLoadBalancingTree(void)
{
  // Create MPI binary tree for organising domain decomposition
  mpitree = new BinaryTree<ndim>(16,0.1,0.0,"geometric","monopole",1,Nmpi);


  // For main process, create load balancing tree, transmit information to all
  // other nodes including particle data
  //---------------------------------------------------------------------------
  if (rank == 0) {

	debug2("[MpiControl::CreateLoadBalancingTree]");


	// Create binary tree from all SPH particles


    // Finally, broadcast all bounding boxes and domain information to all
	// other nodes

  }

  // For other nodes, receive all bounding box and particle data once
  // transmitted by main process.
  //---------------------------------------------------------------------------
  else {


  }
  //---------------------------------------------------------------------------


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
void MpiControl<ndim>::LoadBalancing(void)
{
  return;
}



//=============================================================================
//  MpiControl::TransferParticlesToNode
/// Once we know the new domain boundaries for all MPI nodes, transfer any 
/// particles that now lie in other domain boxes to thos respective MPI 
/// nodes.  Also, receives particles from other domains.
//=============================================================================
template <int ndim>
void MpiControl<ndim>::TransferParticlesToNode(void)
{
  return;
}



// Template class instances for each dimensionality value (1, 2 and 3)
template class MpiControl<1>;
template class MpiControl<2>;
template class MpiControl<3>;
