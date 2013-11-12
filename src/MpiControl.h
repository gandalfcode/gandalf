//=============================================================================
//  MpiControl.h
//  Main MPI class for controlling the distribution of work amongst all MPI 
//  tasks for the current simulation, including load balancing and moving 
//  and copying particles between nodes.
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


#ifndef _MPI_CONTROL_H_
#define _MPI_CONTROL_H_


#include <string>
#include "Precision.h"
#include "MpiNode.h"
#include "Sph.h"
#include "Nbody.h"
#include "SphNeighbourSearch.h"
#include "SphParticle.h"
#include "DomainBox.h"
#if defined MPI_PARALLEL
#include "mpi.h"
#endif
using namespace std;



//=============================================================================
//  Class MpiControl
/// \brief   Main MPI control class for managing MPI simulations.
/// \details Main MPI control class for managing MPI simulations.
/// \author  D. A. Hubber, G. Rosotti
/// \date    09/10/2013
//=============================================================================
template <int ndim>
class MpiControl
{
  void SendParticles(int Node, int Nparticles, int* list, SphParticle<ndim>* );
  void ReceiveParticles (int Node, int& Nparticles, SphParticle<ndim>* array);
 public:

  // Constructor and destructor
  //---------------------------------------------------------------------------
  MpiControl();
  ~MpiControl();


  // Other functions
  //---------------------------------------------------------------------------
  void AllocateMemory(void);
  void DeallocateMemory(void);

  void InitialiseMpiProcess(void);
  void CreateInitialDomainDecomposition(Sph<ndim> *, Nbody<ndim> *, Parameters* , DomainBox<ndim>);
  void LoadBalancing(Sph<ndim> *, Nbody<ndim> *);
  void TransferParticlesToNode(void);


  // MPI control variables
  //---------------------------------------------------------------------------
  bool allocated_mpi;               ///< Flag if memory has been allocated.
  int rank;                         ///< MPI rank of process
  int Nmpi;                         ///< No. of MPI processes
  int Nloadbalance;                 ///< No. of steps between load-balancing

  char hostname[MPI_MAX_PROCESSOR_NAME];
  DomainBox<ndim> mpibox;           ///< ..

  BinaryTree<ndim> *mpitree;        ///< Main MPI load balancing tree
  MpiNode<ndim> *mpinode;           ///< Data for all MPI nodes

  MPI_Datatype particle_type;        ///< Datatype for the particles


};
#endif
