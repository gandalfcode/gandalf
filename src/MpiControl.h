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
#include "MpiTree.h"
#include "SphParticle.h"
#include "DomainBox.h"
#include "Diagnostics.h"
#if defined MPI_PARALLEL
#include "mpi.h"
#endif
using namespace std;


static const int tag_srpart = 1;
static const int tag_league = 2;
static const int tag_bal = 3;


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
 protected:
  //Convenience functions to send and receive particles
  virtual void SendParticles(int Node, int Nparticles, int* list, SphParticle<ndim>* )=0;
  virtual void ReceiveParticles (int Node, int& Nparticles, SphParticle<ndim>** array)=0;


  MPI_Datatype box_type;             ///< Datatype for the box
  MPI_Datatype diagnostics_type;     ///< Datatype for diagnostic info

  //Buffers needed to send and receive particles
  std::vector<int> num_particles_export_per_node;
  std::vector<int> displacements_send;
  std::vector<int> num_particles_to_be_received;
  std::vector<int> receive_displs;
  int tot_particles_to_receive;

  std::vector<Box<ndim> > boxes_buffer;     ///< Buffer needed by the UpdateAllBoundingBoxes routine

  SphNeighbourSearch<ndim>* neibsearch;    ///< Neighbour search class

  void CreateLeagueCalendar();
  std::vector<int> my_matches; ///< List of the matches of this node. For each turn, gives the node we will play with

 public:

  // Constructor and destructor
  //---------------------------------------------------------------------------
  MpiControl();
  ~MpiControl();


  // Other functions
  //---------------------------------------------------------------------------
  void AllocateMemory(int);
  void DeallocateMemory(void);
  void SetNeibSearch(SphNeighbourSearch<ndim>* _neibsearch) {neibsearch=_neibsearch;}

  void CollateDiagnosticsData(Diagnostics<ndim> &);
  virtual void CreateInitialDomainDecomposition(Sph<ndim> *, Nbody<ndim> *, Parameters* , DomainBox<ndim>)=0;
  virtual void LoadBalancing(Sph<ndim> *, Nbody<ndim> *)=0;
  void UpdateAllBoundingBoxes(int, Sph<ndim> *, SphKernel<ndim> *);


  // MPI control variables
  //---------------------------------------------------------------------------
  bool allocated_mpi;               ///< Flag if memory has been allocated.
  int balance_level;                ///< MPI tree level to do load balancing
  int rank;                         ///< MPI rank of process
  int Nmpi;                         ///< No. of MPI processes
  int Nloadbalance;                 ///< No. of steps between load-balancing

  char hostname[MPI_MAX_PROCESSOR_NAME];
  DomainBox<ndim> mpibox;           ///< ..
  MpiNode<ndim> *mpinode;           ///< Data for all MPI nodes
  CodeTiming *timing;               ///< Simulation timing object (pointer)

};


template <int ndim, template<int> class ParticleType>
class MpiControlType : public MpiControl<ndim> {

  using MpiControl<ndim>::balance_level;
  using MpiControl<ndim>::Nmpi;
  using MpiControl<ndim>::mpinode;
  using MpiControl<ndim>::mpibox;
  using MpiControl<ndim>::rank;
  using MpiControl<ndim>::Nloadbalance;
  using MpiControl<ndim>::timing;

  using MpiControl<ndim>::my_matches;
  using MpiControl<ndim>::num_particles_export_per_node;
  using MpiControl<ndim>::displacements_send;
  using MpiControl<ndim>::num_particles_to_be_received;
  using MpiControl<ndim>::receive_displs;
  using MpiControl<ndim>::tot_particles_to_receive;


  virtual void SendParticles(int Node, int Nparticles, int* list, SphParticle<ndim>* );
  virtual void ReceiveParticles (int Node, int& Nparticles, SphParticle<ndim>** array);

  std::vector<ParticleType<ndim> > sendbuffer; ///< Used by the SendParticles routine

  MPI_Datatype particle_type;        ///< Datatype for the particles


  //Buffers needed to send and receive particles
  std::vector<std::vector<ParticleType<ndim>* > > particles_to_export_per_node;
  std::vector<ParticleType<ndim> > particles_to_export;
  std::vector<ParticleType<ndim> > particles_receive;

  MpiTree<ndim,ParticleType> *mpitree;           ///< Main MPI load balancing tree

public:

  MpiControlType ();

  virtual void CreateInitialDomainDecomposition(Sph<ndim> *, Nbody<ndim> *, Parameters* , DomainBox<ndim>);
  virtual void LoadBalancing(Sph<ndim> *, Nbody<ndim> *);
  int SendReceiveGhosts(ParticleType<ndim>** array, Sph<ndim>* sph);
  int UpdateGhostParticles(ParticleType<ndim>** array);
};


#endif
