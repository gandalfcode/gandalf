//=================================================================================================
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
#include "CommunicationHandler.h"
using namespace std;



//=================================================================================================
//  MpiControl::MpiControl()
/// MPI control class constructor.  Initialises all MPI control variables,
/// plus calls MPI routines to find out rank and number of processors.
//=================================================================================================
template <int ndim>
MpiControl<ndim>::MpiControl(DomainBox<ndim>& simboxaux) : simbox(simboxaux)
{
  int len;                          // Length of host processor name string
  Box<ndim> dummy;                  // Dummy box variable

  allocated_mpi = false;
  balance_level = 0;

  // Find local processor rank, total no. of processors and host processor name
  MPI_Comm_size(MPI_COMM_WORLD, &Nmpi);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(hostname, &len);

  // Verify that number of processes is even
  if (!Nmpi%2) {
    string error = "The number of MPI processes must be even!";
    ExceptionHandler::getIstance().raise(error);
  }

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
  Nexport_per_node.resize(Nmpi);
  displacements_send.resize(Nmpi);
  Nreceive_per_node.resize(Nmpi);
  displacements_receive.resize(Nmpi);

  // Create the list of MPI node pair list for communications
  CreateLeagueCalendar();

#ifdef VERIFY_ALL
  if (this->rank == 0) {
    printf("MPI working.  Nmpi : %d   rank : %d   hostname : %s\n",Nmpi,rank,hostname);
  }
  else {
    printf("%d is running too!!\n",this->rank);
  }

  /*if (Nmpi > 1) {
    if (rank == 0) {
      SphParticle<ndim> particle;
      particle.rho = -1.0;
      MPI_Send(&particle, 1, particle_type, 1, 0, MPI_COMM_WORLD);
    }
    else if (rank ==1) {
      SphParticle<ndim> particle;
      MPI_Status status;
      MPI_Recv(&particle,1,particle_type,0,0,MPI_COMM_WORLD,&status);
      if (particle.rho != -1.0) {
        cerr << "Error in transmitting particles: "
                "the last field has not been received correctly!" << endl;
      }
    }
  }*/
#endif

}



//=================================================================================================
//  MpiControl::~MpiControl()
/// MpiControl destructor.
//=================================================================================================
template <int ndim>
MpiControl<ndim>::~MpiControl()
{
  //MPI_Type_free(&particle_type);
  MPI_Type_free(&box_type);
  MPI_Type_free(&diagnostics_type);
}




//=================================================================================================
//  MpiControlType::MpiControlType()
/// MpiControlType constructor.  Does additional work for specific particle types.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
MpiControlType<ndim, ParticleType>::MpiControlType(DomainBox<ndim>& simboxaux) : MpiControl<ndim>(simboxaux)
{
  // Allocate buffers
  particles_to_export_per_node.resize(Nmpi);

  // Create and commit the particle datatype
  particle_type = ParticleType<ndim>::CreateMpiDataType();
  MPI_Type_commit(&particle_type);
}



//=================================================================================================
//  MpiControl::AllocateMemory
/// Allocate all memory for MPI control class.
//=================================================================================================
template <int ndim>
void MpiControl<ndim>::AllocateMemory(int _Ntot)
{
  mpinode = new MpiNode<ndim>[Nmpi];

  for (int inode=0; inode<Nmpi; inode++) {
    mpinode[inode].Ntotmax = (2*_Ntot)/Nmpi;
    mpinode[inode].ids     = new int[mpinode[inode].Ntotmax];
  }

  allocated_mpi = true;

  return;
}



//=================================================================================================
//  MpiControl::DeallocateMemory
/// Deallocate all MPI control class memory.
//=================================================================================================
template <int ndim>
void MpiControl<ndim>::DeallocateMemory(void)
{
  for (int inode=0; inode<Nmpi; inode++) {
    delete[] mpinode[inode].ids;
  }
  delete[] mpinode;

  return;
}



//=================================================================================================
//  MpiControl::CreateLeagueCalendar
/// Create the calendar for the League. In this analogy with soccer, each communication between
/// two nodes is like a football match.  Since everybody neads to speak with everybody, this is
/// effectively like a league. We use the scheduling/Berger algorithm
/// (http://en.wikipedia.org/wiki/Round-robin_tournament#Scheduling_algorithm,
/// http://it.wikipedia.org/wiki/Algoritmo_di_Berger) to organize the event.
//=================================================================================================
template <int ndim>
void MpiControl<ndim>::CreateLeagueCalendar(void)
{
  int Nturns = Nmpi-1;

  // Allocate memory for the calendar
  my_matches.resize(Nturns);

  //-----------------------------------------------------------------------------------------------
  if (rank == 0) {

    // Create a vector containing the full calendar
    // First index is process
    vector<vector<int> > calendar(Nmpi);
    // And second is turn
    for (unsigned int iteam=0; iteam< calendar.size(); iteam++) {
      vector<int>& calendar_team = calendar[iteam];
      calendar_team.resize(Nturns);
    }

    // Create pairs table
    vector<vector<int> > pairs (Nturns);
    for (unsigned int iturn=0; iturn< pairs.size(); iturn++) {
      vector<int>& pairs_turn = pairs[iturn];
      pairs_turn.resize(Nmpi);
      // Fill in the pairs table
      pairs_turn[0]=Nturns;
      for (unsigned int i=1; i<pairs_turn.size(); i++) {
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
    // Validate the calendar
    vector<bool> other_teams(Nmpi-1);

    for (int iteam=0; iteam<calendar.size();iteam++) {

      std::fill(other_teams.begin(), other_teams.end(), false);

      for (int iturn=0; iturn<calendar[iteam].size(); iturn++) {
        int opponent = calendar[iteam][iturn];

        // 1st check: verify the matches are correctly reported in both locations
        if (calendar[opponent][iturn] != iteam) {
          string msg = "Error 1 in validating the calendar!";
          ExceptionHandler::getIstance().raise(msg);
        }

        // 2nd check: verify each team is played against only once
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

    // Copy our calendar to the vector
    for (int iturn=0; iturn<Nturns; iturn++) {
      my_matches[iturn] = calendar[0][iturn];
    }

    // Now transmit the calendar to the other nodes
    for (int inode=1; inode < Nmpi; inode++) {
      MPI_Send(&calendar[inode][0], Nturns, MPI_INT, inode, tag_league, MPI_COMM_WORLD);
    }

  }
  //-----------------------------------------------------------------------------------------------
  else {
    MPI_Status status;
    MPI_Recv(&my_matches[0], Nturns, MPI_INT, 0, tag_league, MPI_COMM_WORLD, &status);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  MpiControl::UpdateAllBoundingBoxes
/// Update local copy of all bounding boxes from all other MPI domains.
//=================================================================================================
template <int ndim>
void MpiControl<ndim>::UpdateAllBoundingBoxes
 (int Npart,                           ///< [in] No. of hydro particles
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to Hydrodynamics object
  SmoothingKernel<ndim> *kernptr,      ///< [in] Pointer to smoothing kernel object
  const bool UseGhosts,          ///< [in] Wheter or not the bounding box should take periodic ghosts into account
  const bool UseTree)      ///< [in] Wheter or not we should use the tree to compute the bounding box rather than explicitly looping
{
  if (rank == 0) debug2("[MpiControl::UpdateAllBoundingBoxes]");
  CodeTiming::BlockTimer timer = timing->StartNewTimer("UPDATEBOUNDINGBOXES");

  // Update local bounding boxes
  if (UseTree)
    neibsearch->UpdateBoundingBoxData(mpinode[rank],UseGhosts);
  else
    mpinode[rank].UpdateBoundingBoxData(Npart,hydro,kernptr);

  // Do an all_gather to receive the new array
  MPI_Allgather(&(mpinode[rank].hbox), 1, box_type, &boxes_buffer[0], 1, box_type, MPI_COMM_WORLD);

  // Save the information inside the nodes
  for (int inode=0; inode<Nmpi; inode++) {
    mpinode[inode].hbox = boxes_buffer[inode];
  }

  // Do an all_gather to receive the new array
  MPI_Allgather(&(mpinode[rank].rbox), 1, box_type, &boxes_buffer[0], 1, box_type, MPI_COMM_WORLD);

  // Save the information inside the nodes
  for (int inode=0; inode<Nmpi; inode++) {
    mpinode[inode].rbox = boxes_buffer[inode];
  }

  // Now update the global bounding box for the entire domain
  // (needed for adjusting mpitree in load balancing step)
  for (int k=0; k<ndim; k++) partbox.min[k] = simbox.min[k];
  for (int k=0; k<ndim; k++) partbox.max[k] = simbox.max[k];
  for (int inode=0; inode<Nmpi; inode++) {
    for (int k=0; k<ndim; k++) {
      partbox.min[k] = min(partbox.min[k], mpinode[inode].rbox.min[k]);
      partbox.max[k] = max(partbox.max[k], mpinode[inode].rbox.max[k]);
    }
  }

  return;
}



//=================================================================================================
//  MpiControl::ComputeTotalStarGasForces
/// Sums up the forces on the star from each processor to find the total ones
//=================================================================================================
template <int ndim>
void MpiControl<ndim>::ComputeTotalStarGasForces
 (Nbody<ndim> * nbody)
{
  // Need to transmit ndim+1 doubles per star (acceleration+gpot)
  const int number_doubles = (ndim*2+1)*nbody->Nnbody;
  const int size_buffer = sizeof(DOUBLE)*number_doubles;
  NbodyParticle<ndim> **star = nbody->nbodydata;

  DOUBLE* buffer = (DOUBLE*) malloc(size_buffer);

  for (int i=0; i<nbody->Nnbody; i++) {
    buffer[i*(2*ndim+1) + 0] = star[i]->gpot;
    for (int k=0; k<ndim; k++) buffer[i*(2*ndim+1) +1 +k] = star[i]->a[k];
    for (int k=0; k<ndim; k++) buffer[i*(2*ndim+1) +1 +ndim +k] = star[i]->adot[k];
  }

  MPI_Allreduce(MPI_IN_PLACE,buffer,number_doubles,GANDALF_MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  for (int i=0; i<nbody->Nnbody; i++) {
    if (!star[i]->flags.check(active)) continue;
    star[i]->gpot = buffer[i*(2*ndim+1) + 0];
    for (int k=0; k<ndim; k++) star[i]->a[k] = buffer[i*(2*ndim+1) +1 +k];
    for (int k=0; k<ndim; k++) star[i]->adot[k] = buffer[i*(2*ndim+1) +1 +ndim +k];
  }

  free(buffer);

  return;
}



//=================================================================================================
//  MpiControl::UpdateStarsAfterAccretion
/// After the sinks have accreted gas particles, need to communicate their new properties
/// The properties updated are:
/// menc, trad, tvisc, ketot, rotketot, gpetot,
/// taccrete, trot, fhydro, utot, angmom
/// And then for the related star:
/// r, v, a, m, r0, v0, a0, dt_internal
/// TBC: Ngas should not be necessary
//=================================================================================================
template <int ndim>
void MpiControl<ndim>::UpdateSinksAfterAccretion
 (Sinks<ndim>* sink)                             ///< [inout] Pointer to sinks array
{
  const int number_variables = ndim*8 + 11;      // ..
  int local_sinks = 0;                           // ..
  int offset = 0;                                // ..
  Box<ndim> mydomain = this->MyDomain();         // ..
  vector<int> owner(sink->Nsink);                // ..
  vector<int> N_sinks_per_rank(Nmpi);            // ..

  // Find out how many stars we have locally and for each processor; also store owner of each sink
  for (int s=0; s<sink->Nsink; s++) {
    if (ParticleInBox(*(sink->sink[s].star), mydomain)) {
      local_sinks++;
      owner[s] = rank;
    }
  }
  N_sinks_per_rank[rank] = local_sinks;

  // Send around the owner vector
  MPI_Allreduce(MPI_IN_PLACE, &owner[0], sink->Nsink, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // Send around the number of sinks per node
  MPI_Allreduce(MPI_IN_PLACE, &N_sinks_per_rank[0], Nmpi, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // Allocate buffers
  const int size_send_buffer = sizeof(DOUBLE)*local_sinks*number_variables;
  const int size_receive_buffer = sizeof(DOUBLE)*sink->Nsink*number_variables;
  DOUBLE* sendbuffer = (DOUBLE*) malloc(size_send_buffer);
  DOUBLE* receivebuffer = (DOUBLE*) malloc(size_receive_buffer);


  //-----------------------------------------------------------------------------------------------
  for (int s=0; s<sink->Nsink; s++) {
    if (owner[s] == rank) {

      // Fill in the buffer
      sendbuffer[offset++] = sink->sink[s].menc;
      sendbuffer[offset++] = sink->sink[s].trad;
      sendbuffer[offset++] = sink->sink[s].tvisc;
      sendbuffer[offset++] = sink->sink[s].ketot;
      sendbuffer[offset++] = sink->sink[s].rotketot;
      sendbuffer[offset++] = sink->sink[s].gpetot;
      sendbuffer[offset++] = sink->sink[s].taccrete;
      sendbuffer[offset++] = sink->sink[s].trot;
      sendbuffer[offset++] = sink->sink[s].utot;
      for (int k=0; k<3; k++) sendbuffer[offset++] = sink->sink[s].angmom[k];

      StarParticle<ndim>& star = *(sink->sink[s].star);
      for (int k=0; k<ndim; k++) sendbuffer[offset++] = star.r[k];
      for (int k=0; k<ndim; k++) sendbuffer[offset++] = star.v[k];
      for (int k=0; k<ndim; k++) sendbuffer[offset++] = star.a[k];
      sendbuffer[offset++] = star.m;
      for (int k=0; k<ndim; k++) sendbuffer[offset++] = star.r0[k];
      for (int k=0; k<ndim; k++) sendbuffer[offset++] = star.v0[k];
      for (int k=0; k<ndim; k++) sendbuffer[offset++] = star.a0[k];
      sendbuffer[offset++] = star.dt_internal;

    }
  }
  //-----------------------------------------------------------------------------------------------

  vector<int> displ(Nmpi);
  compute_displs(displ, N_sinks_per_rank);
  for (unsigned int i=0; i<N_sinks_per_rank.size(); i++) {
    displ[i] *= number_variables;
    N_sinks_per_rank[i] *= number_variables;
  }

  MPI_Allgatherv(sendbuffer, local_sinks*number_variables, GANDALF_MPI_DOUBLE, receivebuffer,
                 &N_sinks_per_rank[0], &displ[0], GANDALF_MPI_DOUBLE, MPI_COMM_WORLD);

  free(sendbuffer);

  vector<int> consumed(Nmpi);

  // Unpack the received data
  //-----------------------------------------------------------------------------------------------
  for (int s=0; s<sink->Nsink; s++) {
    const int origin = owner[s];

    if (origin == rank) continue;

    // Extract the data from receivebuffer[displ[origin]+consumed[origin];
    int offset = 0;
    int offset_rank        = consumed[origin]*number_variables;
    sink->sink[s].menc     = receivebuffer[displ[origin] + offset_rank + offset++];
    sink->sink[s].trad     = receivebuffer[displ[origin] + offset_rank + offset++];
    sink->sink[s].tvisc    = receivebuffer[displ[origin] + offset_rank + offset++];
    sink->sink[s].ketot    = receivebuffer[displ[origin] + offset_rank + offset++];
    sink->sink[s].rotketot = receivebuffer[displ[origin] + offset_rank + offset++];
    sink->sink[s].gpetot   = receivebuffer[displ[origin] + offset_rank + offset++];
    sink->sink[s].taccrete = receivebuffer[displ[origin] + offset_rank + offset++];
    sink->sink[s].trot     = receivebuffer[displ[origin] + offset_rank + offset++];
    sink->sink[s].utot     = receivebuffer[displ[origin] + offset_rank + offset++];
    for (int k=0; k<3; k++) {
      sink->sink[s].angmom[k] = receivebuffer[displ[origin] + offset_rank + offset++];
    }

    StarParticle<ndim>& star = *(sink->sink[s].star);
    for (int k=0; k<ndim; k++) star.r[k] = receivebuffer[displ[origin] + offset_rank + offset++];
    for (int k=0; k<ndim; k++) star.v[k] = receivebuffer[displ[origin] + offset_rank + offset++];
    for (int k=0; k<ndim; k++) star.a[k] = receivebuffer[displ[origin] + offset_rank + offset++];
    star.m = receivebuffer[displ[origin] + offset_rank + offset++];
    for (int k=0; k<ndim; k++) star.r0[k] = receivebuffer[displ[origin] + offset_rank + offset++];
    for (int k=0; k<ndim; k++) star.v0[k] = receivebuffer[displ[origin] + offset_rank + offset++];
    for (int k=0; k<ndim; k++) star.a0[k] = receivebuffer[displ[origin] + offset_rank + offset++];
    star.dt_internal = receivebuffer[displ[origin] + offset_rank + offset++];

    consumed[origin]++;

  }
  //-----------------------------------------------------------------------------------------------

  free(receivebuffer);

}


//=================================================================================================
// namespace MpiReturnParticle
/// Contains classes for returning data computed for ghosts to their original processes.
//=================================================================================================

namespace MpiReturnParticle {
  //===============================================================================================
  // Sink
  /// Return mass for particles that have been partially accreted and flag fully accreted ones as
  /// dead.
  //===============================================================================================
  template<int ndim>
  class ReturnSink {
  public:
    ReturnSink() : iorig(-1), m(0) {} ;
    ReturnSink(const Particle<ndim>& p) : iorig(p.iorig), m(p.m) {} ;

    void update_received(Particle<ndim>& p) const {
      p.m = m ;
      if (p.m == 0) {
        p.flags.unset(active);
        p.flags.set(dead);
      }
    }

    int iorig ;
  private:
    double m ;
  };

  //===============================================================================================
  // ReturnDust
  /// Return the temperature update from drag forces
  //===============================================================================================
  template<int ndim>
  class ReturnDust {
  public:
    ReturnDust() : iorig(-1), dudt(0) {} ;
    ReturnDust(const Particle<ndim>& p) : iorig(p.iorig), dudt(p.dudt) {
      assert(p.ptype == gas_type);
    };

    void update_received(Particle<ndim>& p) const {
      assert(p.ptype == gas_type) ;
      p.dudt += dudt ;
    }

    int iorig;
  private:
    double dudt;
  };

  //===============================================================================================
  // CreateMPIDataType
  /// Creates the MPI data type for a given return type
  //===============================================================================================
  template<class DataType>
  MPI_Datatype& CreateMPIDataType() {
    static bool first = true ;
    static MPI_Datatype MpiReturnType ;
    if (first) {
      MPI_Datatype types[1] = {MPI_BYTE};
      MPI_Aint offsets[1] = {0};
      int blocklen[1] = {sizeof(DataType)};

      MPI_Type_create_struct(1,blocklen,offsets,types, &MpiReturnType);
      MPI_Type_commit(&MpiReturnType);
      first = false;
    }
    return MpiReturnType ;
  }

} // namespace MpiReturnParticle






//=================================================================================================
//  MpiControl::UpdateMpiGhostParents
/// If the MPI ghosts have been modified locally (e.g. because of accretion), need to propagate
/// the changes to their owner process (where the particles are "real")
//=================================================================================================
template <int ndim, template<int> class ParticleType>
template <class ReturnDataType>
void MpiControlType<ndim, ParticleType>::DoUpdateMpiGhostParents
 (list<int>& ids_ghosts,
  Hydrodynamics<ndim>* hydro)
{
  vector<vector<ReturnDataType> > buffer_proc_data(Nmpi) ;
  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();

  // Work out how many ghosts we need to update per each processor
  vector<int> N_updates_per_proc(Nmpi);

  // This vector holds for each processor where its ghosts start in the main hydro array
  vector<int> i_start_ghost(Nmpi);

  assert(Nreceive_per_node[rank] == 0);
  std::partial_sum(Nreceive_per_node.begin(), Nreceive_per_node.end(), i_start_ghost.begin());

  for (list<int>::iterator it = ids_ghosts.begin(); it != ids_ghosts.end(); it++) {
    const int index = *it;
    assert(index >= hydro->Nhydro + hydro->NPeriodicGhost);
    assert(index <  hydro->Nhydro + hydro->Nghost);

    // Create the return type
    ReturnDataType temp(partdata[index]) ;


    // Find to which processor we should send the particle
    const int proc = std::upper_bound(i_start_ghost.begin(),
                                      i_start_ghost.end(), temp.iorig) - i_start_ghost.begin();

    N_updates_per_proc[proc]++;
    assert(N_updates_per_proc[proc]<=Nreceive_per_node[proc]);

    // Check that iorig/ptype hasn't changed
    assert(temp.iorig == particles_receive[temp.iorig].iorig);
    assert(particles_receive[temp.iorig].ptype == partdata[index].ptype);

    // Tell the target proc where the particle belongs
    temp.iorig = particle_ids_receive[temp.iorig] ;
    buffer_proc_data[proc].push_back(temp) ;
  }

  // Now that we know how many ghosts we need to send to each processor, can fill the buffer
  vector<ReturnDataType> buffer(ids_ghosts.size());

  // Partial sum of the numbers
  vector<int> displs_ghosts(Nmpi);
  compute_displs(displs_ghosts, N_updates_per_proc);

  // Copy from the buffer of each processor to the total buffer
  for (int iproc=0; iproc<Nmpi; iproc++) {
    std::copy(buffer_proc_data[iproc].begin(), buffer_proc_data[iproc].end(),
              buffer.begin() + displs_ghosts[iproc]);
  }

  // Tell everybody how many ghosts we have to update remotely
  vector<int> N_updates_from_proc(Nmpi);
  MPI_Alltoall(&N_updates_per_proc[0], 1, MPI_INT,
               &N_updates_from_proc[0], 1, MPI_INT, MPI_COMM_WORLD);

  vector<int> displs_ghosts_receive(Nmpi);
  compute_displs(displs_ghosts_receive, N_updates_from_proc);

  const int total_ghost_receive = std::accumulate(N_updates_from_proc.begin(),
                                                  N_updates_from_proc.end(), 0);

  vector<ReturnDataType> buffer_receive(total_ghost_receive);

  MPI_Datatype& MpiReturnType =
      MpiReturnParticle::CreateMPIDataType<ReturnDataType>();

  MPI_Alltoallv(&buffer[0], &N_updates_per_proc[0], &displs_ghosts[0], MpiReturnType,
                &buffer_receive[0], &N_updates_from_proc[0], &displs_ghosts_receive[0],
                MpiReturnType, MPI_COMM_WORLD);

  // Update the real particles from the received data
  for (int i=0; i< total_ghost_receive; i++) {
    const ReturnDataType& received = buffer_receive[i] ;

    assert(received.iorig != -1 &&  received.iorig < hydro->Nhydro + hydro->NPeriodicGhost) ;

    // Find the real particle that the ghost corresponds to
    int j = received.iorig ;
    while (partdata[j].flags.check(any_boundary)) {
      assert(j >= hydro->Nhydro);
      j = partdata[j].iorig ;
    }
    assert(j >= 0 && j < hydro->Nhydro);

    received.update_received(partdata[j]) ;
  }

  return;
}

//=================================================================================================
//  MpiControl::UpdateMpiGhostParents
/// If the MPI ghosts have been modified locally (e.g. because of accretion), need to propagate
/// the changes to their owner process (where the particles are "real")
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim, ParticleType>::UpdateMpiGhostParents
 (list<int>& ids_ghosts,
  Hydrodynamics<ndim>* hydro,
  MPI_ghost_update_type update_type)
{
  using namespace MpiReturnParticle ;

  switch(update_type) {
    case update_sink_parents:
      this->DoUpdateMpiGhostParents<ReturnSink<ndim> >(ids_ghosts, hydro) ;
      break;
    case update_dust_parents:
      this->DoUpdateMpiGhostParents<ReturnDust<ndim> >(ids_ghosts, hydro) ;
      break;
    default:
      ExceptionHandler::getIstance().raise("Requested MpiGhost update not recognised");
  }
}

//=================================================================================================
//  MpiControlType::ExportParticlesBeforeForceLoop
/// Export the particles that need force contribution from other processors
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim,ParticleType>::ExportParticlesBeforeForceLoop
 (Hydrodynamics<ndim>* hydro)              ///< Pointer to hydrodynamics object
{
  CodeTiming::BlockTimer timer = timing->StartNewTimer("EXPORT_PARTICLES");
  vector<vector <char> > send_buffer(Nmpi-1);
  vector<vector <char> > receive_buffer(Nmpi-1);
  vector<vector <char> > header_receive(Nmpi-1);

  const int header_size = neibsearch->ExportSize(0, hydro).size();

  // Post the receives for the header
  MPI_Request req_header[Nmpi-1];
  int j=0;
  for (int iproc=0; iproc<Nmpi; iproc++) {

    if (iproc==rank)
      continue;

    header_receive[j].resize(header_size);

    MPI_Irecv(&header_receive[j][0],header_size,MPI_CHAR,iproc,4,MPI_COMM_WORLD,&req_header[j]);
    j++;

  }

  // Send the header to each processor
  j=0;
  MPI_Request sendreq_header[Nmpi-1];
  vector<vector<char> > header(Nmpi-1);
  for (int iproc=0; iproc<Nmpi; iproc++) {

    if (iproc==rank)
      continue;

    header[j] = neibsearch->ExportSize(iproc, hydro);

    MPI_Isend(&header[j][0],header_size,MPI_CHAR,iproc,4,MPI_COMM_WORLD,&sendreq_header[j]);
    j++;

  }

  // Check if all the headers have been received
  int reqheader_completed=0;
  int reqheader_complnow;
  int which_completed[Nmpi-1];
  MPI_Status status_header[Nmpi-1];
  MPI_Request req[Nmpi-1];
  MPI_Status status[Nmpi-1];
  MPI_Testsome(Nmpi-1,req_header,&reqheader_complnow,which_completed,status_header);
  // If some requests have finished, post the receive
  if (reqheader_complnow) {
    for (int i=0; i< reqheader_complnow; i++) {
      const int j=which_completed[i];
      int iproc = j;
      if (j >= rank)
        iproc +=1;

      int size_receive;
      copy(&size_receive,&header_receive[j][0]);
      receive_buffer[j].resize(size_receive);
      MPI_Irecv(&receive_buffer[j][0],size_receive,MPI_CHAR,iproc,5,MPI_COMM_WORLD,&req[j]);

    }
  }
  reqheader_completed = reqheader_complnow;


  // Prepare the information to send
  j=0;
  MPI_Request send_req[Nmpi-1];
  for (int iproc=0; iproc<Nmpi; iproc++) {

    // No need to send anything to ourselves
    if (iproc == rank) continue;

    // Pack the information to send
   neibsearch->GetExportInfo(iproc, hydro, send_buffer[j], mpinode[iproc], rank, Nmpi);

    // Post the send now that we know what to transmit
    MPI_Isend(&send_buffer[j][0],send_buffer[j].size(),MPI_CHAR,iproc,5,MPI_COMM_WORLD,&send_req[j]);

    // If we still don't know how much data we are receiving from some processor, check if the
    // transmission has finished in the meanwhile
    if (reqheader_completed < Nmpi-1) {
      MPI_Testsome(Nmpi-1,req_header,&reqheader_complnow,which_completed,status_header);
      // If some request have finished, post the receive
      if (reqheader_complnow) {
        for (int i=0; i< reqheader_complnow; i++) {
          const int j=which_completed[i];
          int iproc = j;
          if (j >= rank)
            iproc +=1;

          int size_receive;
          copy(&size_receive,&header_receive[j][0]);
          receive_buffer[j].resize(size_receive);
          MPI_Irecv(&receive_buffer[j][0],size_receive,MPI_CHAR,iproc,5,MPI_COMM_WORLD,&req[j]);

        }
        reqheader_completed += reqheader_complnow;
      }
    }


    j++;
  }


  // We have posted all the sends - the only thing to do know is finishing posting all the receives,
  // if that hasn't happened yet
  while (reqheader_completed < Nmpi-1) {
    MPI_Waitsome(Nmpi-1,req_header,&reqheader_complnow,which_completed,status_header);
    // If some request have finished, post the receive
    if (reqheader_complnow) {
      for (int i=0; i< reqheader_complnow; i++) {
        const int j=which_completed[i];
        int iproc = j;
        if (j >= rank)
          iproc +=1;

        int size_receive;
        copy(&size_receive,&header_receive[j][0]);
        receive_buffer[j].resize(size_receive);
        MPI_Irecv(&receive_buffer[j][0],size_receive,MPI_CHAR,iproc,5,MPI_COMM_WORLD,&req[j]);

      }
      reqheader_completed += reqheader_complnow;
    }
  }

  // Now we can only wait to receive all the information
  int req_completed=0;
  bool first_unpack=true;
  while (req_completed<Nmpi-1) {
    int completed_now;
    MPI_Waitsome(Nmpi-1,req,&completed_now,which_completed,status);

    // Unpack the received information into the main arrays
    for (int i=0; i<completed_now; i++) {
        const int j = which_completed[i];
        int iproc = j;
        if (j >= rank)
            iproc +=1;
        neibsearch->UnpackExported(receive_buffer[j], hydro, iproc,
            header_receive, rank, first_unpack);
        first_unpack=false;
    }

    req_completed += completed_now;

  }

  // Sends are probably finished by now, but we do need to wait so that they can be deallocated
  MPI_Waitall(Nmpi-1,sendreq_header,MPI_STATUSES_IGNORE);
  MPI_Waitall(Nmpi-1,send_req,MPI_STATUSES_IGNORE);

  return;
}



//=================================================================================================
//  MpiControlType::GetExportedParticlesAccelerations
/// Get back the information about the exported particles from the other processors
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim,ParticleType>::GetExportedParticlesAccelerations
 (Hydrodynamics<ndim>* hydro)              ///< Pointer to main hydrodynamics object
{
  CodeTiming::BlockTimer timer = timing->StartNewTimer("GET_EXPORTED");
  vector<vector<char> > send_buffer(Nmpi-1);                // ..
  vector<vector<char> > receive_buffer(Nmpi-1);


  // Post all the receives (for performance)
  MPI_Request req[Nmpi-1];
  MPI_Status status[Nmpi-1];
  int j=0;
  for (int i=0; i<Nmpi; i++) {

	  if (i==rank)
		  continue;

	  // We know how much data we are going to receive
	  const int size_receive = neibsearch->ExportInfoSize(i);
	  receive_buffer[j].resize(size_receive);

	  MPI_Irecv(&receive_buffer[j][0],size_receive,MPI_CHAR,i,5,MPI_COMM_WORLD,&req[j]);

	  j++;
  }


  // Get the accelerations to send to the other processors
  // Quite confusingly, note that the role of the two arrays with the number of bytes from/to
  // each processor is now reversed (from is to send and to is to receive)
  MPI_Request send_req[Nmpi-1];
  j=0;
  for (int i=0; i<Nmpi; i++) {

	  if (i==rank)
		  continue;

	  // Prepare the data to send to each processor
	  neibsearch->GetBackExportInfo(send_buffer[j], hydro, rank, i);

	  // Post the send immediately now that we know which data to send
	  MPI_Isend(&send_buffer[j][0],send_buffer[j].size(),MPI_CHAR,i,5,MPI_COMM_WORLD,&send_req[j]);

	  j++;

  }
  //Now that we have extracted all the needed information reset the counters
  neibsearch->ResetCountersExportInfo(hydro);


  // Now wait for all the communication to be finished

  int ncompleted=0;
  int which_finished[Nmpi-1];
  int niter = 0;
  while (ncompleted < Nmpi-1) {
	  int completed_now;
	  MPI_Waitsome(Nmpi-1,req,&completed_now,which_finished,status);

	  // Unpack the received information to update accelerations
	  for (int i=0; i<completed_now; i++) {
		  j = which_finished[i];
		  int iproc = j;
		  if (j >= rank)
			  iproc +=1;
		  neibsearch->UnpackReturnedExportInfo(receive_buffer[j], hydro, rank, iproc);
	  }

	  ncompleted += completed_now;

	  niter++;

  }
  hydro->FinishReturnExport();


  // Vector with sends gets deallocated when the function returns, so need to be sure that the sending request have completed
  // Probably not a big problem - if the receives have finished, very likely also the sends have
  MPI_Waitall(Nmpi-1,send_req,MPI_STATUSES_IGNORE);

  return;
}



//=================================================================================================
//  MpiControl::SendReceiveGhosts
/// Compute particles to send to other nodes and receive needed particles from other nodes.
/// Returns the number of ghosts received. The array with the received particles is stored inside
/// the pointer returned. The caller must NOT delete this array when is finished with it, as the
/// memory is internally managed by this class.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int MpiControlType<ndim, ParticleType>::SendReceiveGhosts
 (const FLOAT tghost,                  ///< Ghost 'lifetime'
  Hydrodynamics<ndim>* hydro,          ///< Main SPH object pointer
  ParticleType<ndim>** array)          ///< Main SPH particle array
{
  int i;                               // Particle counter
  unsigned int iaux;                   // ..
  int index = 0;                       // ..
  int Nexport;                         // ..
  int Npart;                           // ..
  int running_counter;                 // ..
  unsigned int inode;                  // MPI node index
  vector<int> overlapping_nodes;       // List of nodes that overlap this domain
  vector<int> ghost_export_list;       // List of particles ids to be exported
  ParticleType<ndim>* partdata = hydro->template GetParticleArray<ParticleType>();

  if (rank == 0) debug2("[MpiControl::SendReceiveGhosts]");

  // Reserve space for all nodes
  overlapping_nodes.reserve(Nmpi);

  // Loop over domains and find the ones that could overlap us
  for (int inode=0; inode<Nmpi; inode++) {
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
    Nexport = neibsearch->SearchMpiGhostParticles(tghost, mpinode[inode].rbox,
                                                  hydro, ghost_export_list);
    for (unsigned int j=0; j<ghost_export_list.size(); j++) {
      i = ghost_export_list[j];
      particles_to_export_per_node[inode].push_back(i);
      assert(!partdata[i].flags.is_dead());
    }
#ifdef OUTPUT_ALL
    cout << "Nexport : " << Nexport << "     size : " << ghost_export_list.size() << endl;
#endif
    assert((unsigned int) Nexport == ghost_export_list.size());
  }

  // Prepare arrays with no. of particles to export per node and displacements
  fill(Nexport_per_node.begin(), Nexport_per_node.end(), 0);

  running_counter = 0;
  for (int inode=0; inode<Nmpi; inode++) {
    Npart = particles_to_export_per_node[inode].size();
    Nexport_per_node[inode] = Npart;
    displacements_send[inode] = running_counter;
    running_counter += Npart;
  }

  // Compute total number of particles to export
  Nexport = std::accumulate(Nexport_per_node.begin(), Nexport_per_node.end(), 0);

  // Comunicate with all processors the no. of ptcls that everyone is exporting
  vector<int> ones(Nmpi,1);
  vector<int> displs(Nmpi);
  for (inode=0; inode<displs.size(); inode++) {
    displs[inode] = inode;
  }
  MPI_Alltoallv(&Nexport_per_node[0], &ones[0], &displs[0], MPI_INT,
                &Nreceive_per_node[0], &ones[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

  // Compute the total number of particles that will be received
  Nreceive_tot = accumulate(Nreceive_per_node.begin(),
                                        Nreceive_per_node.end(),0);

  // Allocate receive buffer
  particles_receive.clear();
  particles_receive.resize(Nreceive_tot);

#if defined(OUTPUT_ALL)
  cout << "TOT_PARTICLES_TO_EXPORT : " << Nexport << endl;
#endif

  // Create vector containing all particles to export
  particles_to_export.resize(Nexport);
  particle_ids_send.resize(Nexport);

  for (int inode=0; inode<Nmpi; inode++) {
    vector<int >& particles_on_this_node = particles_to_export_per_node[inode];

    for (unsigned int iparticle=0; iparticle<particles_on_this_node.size(); iparticle++) {
      particles_to_export[index] = partdata[particles_on_this_node[iparticle]];
      particle_ids_send[index] = particles_to_export[index].iorig;

      // Record in iorig the location in memory of the particle
      particles_to_export[index].iorig = particles_on_this_node[iparticle];
      assert(!particles_to_export[index].flags.is_dead());
      index++;
    }
  }
  assert(index==Nexport);

  // Compute receive displacements
  running_counter = 0;
  for (inode=0; inode<displacements_receive.size(); inode++) {
    displacements_receive[inode] = running_counter;
    running_counter += Nreceive_per_node[inode];
  }

  // Send and receive particles
  MPI_Alltoallv(&particles_to_export[0], &Nexport_per_node[0], &displacements_send[0],
                particle_type, &particles_receive[0], &Nreceive_per_node[0],
                &displacements_receive[0], particle_type, MPI_COMM_WORLD);

  // Now that we've received the particles, set their iorig to the location in the
  // received array and store the 'original' iorig
  particle_ids_receive.resize(Nreceive_tot);
  for (int i=0; i < Nreceive_tot; i++) {
    particle_ids_receive[i] = particles_receive[i].iorig ;
    particles_receive[i].iorig = i ;
  }

  *array = &particles_receive[0];

  return Nreceive_tot;
}



//=================================================================================================
//  MpiControl::UpdateGhostParticles
/// Update the ghost particles that were previously found. Return the number of Mpi Ghosts;
/// modified the passed pointer with the address of the array of ghosts (note: the caller must
/// NOT call delete on this array, as the memory is internally managed by the MpiControl class).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int MpiControlType<ndim, ParticleType>::UpdateGhostParticles
 (		 ParticleType<ndim>* partdata,			  ///< Main array of particles
		 ParticleType<ndim>** array)          ///< Will contain a pointer to the array of ghosts
{
  int index = 0;                       // ..
  int inode;                           // MPI node counter
  unsigned int ipart;                  // Particle counter

  // Update the local buffer of particles to send
  for (inode=0; inode<Nmpi; inode++) {
    vector<int >& particles_on_this_node = particles_to_export_per_node[inode];
    for (ipart=0; ipart<particles_on_this_node.size(); ipart++) {
      const ParticleType<ndim>& part = partdata[particles_on_this_node[ipart]];
      // Check nobody has swapped the particles
      assert(particle_ids_send[index] == part.iorig);

      particles_to_export[index] = part;
      particles_to_export[index].iorig = particles_on_this_node[ipart];
      index++;
    }
  }

  // Send and receive particles
  MPI_Alltoallv(&particles_to_export[0], &Nexport_per_node[0], &displacements_send[0],
                particle_type, &particles_receive[0], &Nreceive_per_node[0],
                &displacements_receive[0], particle_type, MPI_COMM_WORLD);


  for (int i=0; i < Nreceive_tot; i++) {
    // Check that we received a particle with the expected origin
    assert(particles_receive[i].iorig == particle_ids_receive[i]);
    particles_receive[i].iorig = i ;
  }

  *array = &particles_receive[0];

  return Nreceive_tot;
}



//=================================================================================================
//  MpiControl::SendParticles
/// Given an array of ids and a node, copy particles inside a buffer and
/// send them to the given node.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim, ParticleType>::SendParticles
 (int Node,                            ///< ..
  int Nparticles,                      ///< ..
  int* list,                           ///< ..
  ParticleType<ndim>* main_array)      ///< main_array_gen
{
  // Ensure there is enough memory in the buffer
  sendbuffer.resize(Nparticles);

  //Copy particles from the main arrays to the buffer
  for (int i=0; i<Nparticles; i++) {
    sendbuffer[i] = main_array[list[i]];
  }

  MPI_Send(&sendbuffer[0], Nparticles, particle_type, Node, tag_srpart, MPI_COMM_WORLD);

#ifdef OUTPUT_ALL
  cout << "SENDING : " << Nparticles << "   to " << Node << endl;
#endif

  return;
}



//=================================================================================================
//  MpiControl::ReceiveParticles
/// Given a node, receive particles from it. Return the number of particles received and
/// a pointer to the array containing the particles. The caller is reponsible to free the
/// array after its usage.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void MpiControlType<ndim, ParticleType>::ReceiveParticles
 (int Node,                            ///< ..
  int& Nparticles,                     ///< ..
  ParticleType<ndim>** array)          ///< pointer_array_gen
{
  const int tag = 1;                   // ..
  MPI_Status status;                   // ..

  // "Probe" the message to know how big the message is going to be
  MPI_Probe(Node, tag, MPI_COMM_WORLD, &status);

  // Get the number of elements
  MPI_Get_count(&status, particle_type, &Nparticles);

#ifdef OUTPUT_ALL
  cout << "RECEIVING NPARTICLES : " << Nparticles << "   node : " << Node << endl;
#endif

  // Allocate enough memory to hold the particles
  *array = new ParticleType<ndim> [Nparticles];

  // Now receive the message
  MPI_Recv(*array, Nparticles, particle_type, Node, tag_srpart, MPI_COMM_WORLD, &status);

  return;
}



//=================================================================================================
//  MpiControl::CollateDiagnosticsData
/// ..
//=================================================================================================
template <int ndim>
void MpiControl<ndim>::CollateDiagnosticsData
 (Diagnostics<ndim> &diag)             ///< Main diagnostics object
{
  int inode;                           // MPI node counter
  int k;                               // Dimension counter
  Diagnostics<ndim> diagaux;           // Aux. diagnostics struct
  MPI_Status status;                   // MPI communication status message

  //-----------------------------------------------------------------------------------------------
  if (rank == 0) {

    // First multiply root node values by mass
    for (k=0; k<ndim; k++) diag.rcom[k] *= diag.mtot;
    for (k=0; k<ndim; k++) diag.vcom[k] *= diag.mtot;

    for (inode=1; inode<Nmpi; inode++) {
      MPI_Recv(&diagaux, 1, diagnostics_type, inode, 0, MPI_COMM_WORLD, &status);
      diag.Nhydro += diagaux.Nhydro;
      diag.Etot   += diagaux.Etot;
      diag.utot   += diagaux.utot;
      diag.ketot  += diagaux.ketot;
      diag.gpetot += diagaux.gpetot;
      diag.mtot   += diagaux.mtot;
      for (k=0; k<ndim; k++) diag.mom[k] += diagaux.mom[k];
      for (k=0; k<3; k++) diag.angmom[k] += diagaux.angmom[k];
      for (k=0; k<ndim; k++) diag.force[k] += diagaux.force[k];
      for (k=0; k<ndim; k++) diag.rcom[k] += diagaux.mtot*diagaux.rcom[k];
      for (k=0; k<ndim; k++) diag.vcom[k] += diagaux.mtot*diagaux.vcom[k];
    }

    // Renormalise centre of mass positions and velocities
    for (k=0; k<ndim; k++) diag.rcom[k] /= diag.mtot;
    for (k=0; k<ndim; k++) diag.vcom[k] /= diag.mtot;

  }
  //-----------------------------------------------------------------------------------------------
  else {

    MPI_Send(&diag, 1, diagnostics_type, 0, 0, MPI_COMM_WORLD);

  }
  //-----------------------------------------------------------------------------------------------

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
template class MpiControlType<1, MeshlessFVParticle>;
template class MpiControlType<2, MeshlessFVParticle>;
template class MpiControlType<3, MeshlessFVParticle>;
