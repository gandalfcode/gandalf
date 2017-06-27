//=================================================================================================
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
//=================================================================================================


#ifndef _MPI_CONTROL_H_
#define _MPI_CONTROL_H_


#include <string>
#include "Precision.h"
#include "Diagnostics.h"
#include "DomainBox.h"
#include "Hydrodynamics.h"
#include "Nbody.h"
#include "MpiNode.h"
#include "MpiTree.h"
#include "NeighbourSearch.h"
#include "Particle.h"
#if defined MPI_PARALLEL
#include "mpi.h"
#endif
using namespace std;


static const int tag_srpart = 1;
static const int tag_league = 2;
static const int tag_bal = 3;

template <int ndim>
class Sinks;


enum MPI_ghost_update_type  {
    update_sink_parents = 0,
    update_dust_parents = 1
};

//=================================================================================================
//  Class MpiControl
/// \brief   Main MPI control class for managing MPI simulations.
/// \details Main MPI control class for managing MPI simulations.
/// \author  D. A. Hubber, G. Rosotti
/// \date    09/10/2013
//=================================================================================================
template <int ndim>
class MpiControl
{
 protected:

  MPI_Datatype box_type;                   ///< Datatype for the box
  MPI_Datatype diagnostics_type;           ///< Datatype for diagnostic info
  MPI_Datatype ExportParticleType;         ///< Datatype for the information to export
  MPI_Datatype ExportBackParticleType;     ///< Datatype for the information to get back
                                           ///< from exported particles

  // Buffers needed to send and receive particles
  //-----------------------------------------------------------------------------------------------
  int Nreceive_tot;                        ///< Total number of particles to be received
  vector<int> displacements_receive;       ///< ???
  vector<int> displacements_send;          ///< ???
  vector<int> my_matches;                  ///< List of the matches of this node.
                                           ///< For each turn, gives the node we will play with
  vector<int> Nexport_per_node;            ///< No. of particles exported to each node
  vector<int> Nreceive_per_node;           ///< No. of particles received from each node
  vector<Box<ndim> > boxes_buffer;         ///< Buffer needed by UpdateAllBoundingBoxes routine
  NeighbourSearch<ndim>* neibsearch;       ///< Neighbour search class

  void CreateLeagueCalendar();

  DomainBox<ndim>& simbox;                  ///< Copy of the simbox object


 public:

  // MPI control variables
  //-----------------------------------------------------------------------------------------------
  char hostname[MPI_MAX_PROCESSOR_NAME];   ///< ..
  bool allocated_mpi;                      ///< Flag if memory has been allocated.
  int balance_level;                       ///< MPI tree level to do load balancing
  int rank;                                ///< MPI rank of process
  int Nmpi;                                ///< No. of MPI processes
  int Nloadbalance;                        ///< No. of steps between load-balancing
  Box<ndim> partbox;                       ///< Box containing all particles (from all nodes)
  MpiNode<ndim> *mpinode;                  ///< Data for all MPI nodes
  CodeTiming *timing;                      ///< Simulation timing object (pointer)


  // Constructor and destructor
  //-----------------------------------------------------------------------------------------------
  MpiControl(DomainBox<ndim>& simboxaux);
  ~MpiControl();

  void AllocateMemory(int);
  void DeallocateMemory(void);
  void SetNeibSearch(NeighbourSearch<ndim>* _neibsearch) {neibsearch = _neibsearch;}
  void CollateDiagnosticsData(Diagnostics<ndim> &);
  void UpdateAllBoundingBoxes(int, Hydrodynamics<ndim> *, SmoothingKernel<ndim> *,const bool,const bool);
  void ComputeTotalStarGasForces(Nbody<ndim> * nbody);

  //void ExportMpiGhostParticles(FLOAT, DomainBox<ndim>, Hydrodynamics<ndim> *);
  virtual void ExportParticlesBeforeForceLoop(Hydrodynamics<ndim> *) = 0;
  virtual void GetExportedParticlesAccelerations(Hydrodynamics<ndim> *) = 0;
  virtual void CreateInitialDomainDecomposition(Hydrodynamics<ndim> *, Nbody<ndim> *,
                                                Parameters*,  bool&) = 0;
  virtual void LoadBalancing(Hydrodynamics<ndim> *, Nbody<ndim> *) = 0;
  virtual void UpdateMpiGhostParents (list<int>& ids, Hydrodynamics<ndim>* hydro,
                                      MPI_ghost_update_type=update_sink_parents)=0;
  void UpdateSinksAfterAccretion(Sinks<ndim>* sink);

  Box<ndim> MyDomain();

};



//=================================================================================================
//  Class MpiControlType
/// \brief   ..
/// \details ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    09/10/2013
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class MpiControlType : public MpiControl<ndim>
{
public:
  using MpiControl<ndim>::balance_level;
  using MpiControl<ndim>::Nmpi;
  using MpiControl<ndim>::mpinode;
  using MpiControl<ndim>::rank;
  using MpiControl<ndim>::Nloadbalance;
  using MpiControl<ndim>::timing;
  using MpiControl<ndim>::my_matches;
  using MpiControl<ndim>::Nexport_per_node;
  using MpiControl<ndim>::displacements_send;
  using MpiControl<ndim>::Nreceive_per_node;
  using MpiControl<ndim>::displacements_receive;
  using MpiControl<ndim>::Nreceive_tot;
  using MpiControl<ndim>::ExportParticleType;
  using MpiControl<ndim>::ExportBackParticleType;
  using MpiControl<ndim>::neibsearch;
  using MpiControl<ndim>::partbox;


  // Buffers needed to send and receive particles
  //-----------------------------------------------------------------------------------------------
  vector<vector<int> > particles_to_export_per_node; ///< ..
  vector<int> particle_ids_receive;                  ///< Origin location of imported particles
  vector<int> particle_ids_send ;
  vector<ParticleType<ndim> > particles_to_export;   ///< Particles to be exported to other nodes
  vector<ParticleType<ndim> > particles_receive;     ///< Particles to be received from other nodes
  vector<ParticleType<ndim> > sendbuffer;            ///< Used by the SendParticles routine
  MPI_Datatype particle_type;                        ///< Datatype for the particles


  //-----------------------------------------------------------------------------------------------
  MpiControlType(DomainBox<ndim>& simboxaux);
  ~MpiControlType() {};


  virtual void ExportParticlesBeforeForceLoop (Hydrodynamics<ndim> *);
  virtual void GetExportedParticlesAccelerations (Hydrodynamics<ndim> *);
  virtual void UpdateMpiGhostParents (list<int>& ids, Hydrodynamics<ndim>*,MPI_ghost_update_type);

  void SendParticles(int Node, int Nparticles, int* list, ParticleType<ndim>*);
  void ReceiveParticles(int Node, int& Nparticles, ParticleType<ndim>** array);
  int SendReceiveGhosts(const FLOAT, Hydrodynamics<ndim> *,ParticleType<ndim>**);
  int UpdateGhostParticles(ParticleType<ndim>*, ParticleType<ndim>**);

private:
  template<class ReturnDataType> void DoUpdateMpiGhostParents(list<int>& ids, Hydrodynamics<ndim>*);

};


template <int ndim>
inline Box<ndim> MpiControl<ndim>::MyDomain()
{
  if (!allocated_mpi) {
    Box<ndim> result;
    for (int k=0; k<ndim; k++) {
      result.min[k] = -big_number;
      result.max[k] = big_number;
    }
    return result;
  }
  else {
    return mpinode[rank].domain;
  }
}



//=================================================================================================
//  Class MpiKDTreeDecomposition
/// \brief   ..
/// \details ..
/// \author  D. A. Hubber, G. Rosotti
/// \date    23/06/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class MpiKDTreeDecomposition : public MpiControlType<ndim,ParticleType>
{
public:
  using MpiControl<ndim>::balance_level;
  using MpiControl<ndim>::Nmpi;
  using MpiControl<ndim>::mpinode;
  using MpiControl<ndim>::rank;
  using MpiControl<ndim>::Nloadbalance;
  using MpiControl<ndim>::timing;
  using MpiControl<ndim>::my_matches;
  using MpiControl<ndim>::Nexport_per_node;
  using MpiControl<ndim>::displacements_send;
  using MpiControl<ndim>::Nreceive_per_node;
  using MpiControl<ndim>::displacements_receive;
  using MpiControl<ndim>::Nreceive_tot;
  using MpiControl<ndim>::ExportParticleType;
  using MpiControl<ndim>::ExportBackParticleType;
  using MpiControl<ndim>::neibsearch;
  using MpiControl<ndim>::partbox;
  using MpiControl<ndim>::simbox;
  using MpiControlType<ndim,ParticleType>::particles_to_export_per_node;
  using MpiControlType<ndim,ParticleType>::particles_to_export;
  using MpiControlType<ndim,ParticleType>::particles_receive;
  using MpiControlType<ndim,ParticleType>::sendbuffer;
  using MpiControlType<ndim,ParticleType>::particle_type;


  MpiTree<ndim,ParticleType> *mpitree;                 ///< Main MPI load balancing tree


  //-----------------------------------------------------------------------------------------------
  MpiKDTreeDecomposition(DomainBox<ndim>& simboxaux): MpiControlType<ndim,ParticleType>(simboxaux) {};
  //virtual ~MpiKDTreeDecomposition();

  virtual void CreateInitialDomainDecomposition(Hydrodynamics<ndim> *, Nbody<ndim> *,
                                                Parameters*, bool&);
  virtual void LoadBalancing(Hydrodynamics<ndim> *, Nbody<ndim> *);

};
#endif
