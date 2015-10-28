//=================================================================================================
//  BruteForceSearch.cpp
//  Contains all routines for generating SPH neighbour lists using
//  brute-force (i.e. direct summation over all particles).
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



#include <assert.h>
#include <iostream>
#include <math.h>
#include <numeric>
#include "NeighbourSearch.h"
//#include "SphNeighbourSearch.h"
#include "Sph.h"
#include "Parameters.h"
#include "Particle.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "SmoothingKernel.h"
#if defined MPI_PARALLEL
#include "MpiNode.h"
#endif
using namespace std;



//=================================================================================================
//  BruteForceSearch::BuildTree
/// For Brute Force neighbour searching, there is no tree to construct but we delete any
/// dead SPH particles here to be consistent with the tree neighbour search.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::BuildTree
 (const bool rebuild_tree,             ///< [in] Flag to rebuild tree
  const int n,                ///< [in] Integer time
  const int ntreebuildstep,            ///< [in] Tree build frequency
  const int ntreestockstep,            ///< [in] Tree stocking frequency
  const int Npart,                     ///< [in] No. of particles
  const int Npartmax,                  ///< [in] Max. no. of particles
  const FLOAT timestep,                ///< [in] Smallest physical timestep
  Particle<ndim> *part_gen,            ///< [inout] Particle data array
  Hydrodynamics<ndim> *hydro)          ///< [inout] Pointer to Hydrodynamics object
{
  hydro->DeleteDeadParticles();
  return;
}



//=================================================================================================
//  BruteForceSearch::GetGatherNeighbourList
/// Compute gather neighbour list from all hydro particles within a distance rsearch of
/// position rp using brute force.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int BruteForceSearch<ndim,ParticleType>::GetGatherNeighbourList
 (FLOAT rp[ndim],                      ///< [in] Position
  FLOAT rsearch,                       ///< [in] Search radius
  Particle<ndim> *part_gen,            ///< [in] Pointer to particle data array
  int Nhydro,                          ///< [in] No. of SPH particles
  int Nneibmax,                        ///< [in] Max. size of neighbour list
  int *neiblist)                       ///< [out] List of neighbour ids
{
  int i,k;                             // Particle and dimension counters
  int Nneib = 0;                       // No. of (non-dead) neighbours
  FLOAT dr[ndim];                      // Relative distance vector
  FLOAT drsqd;                         // Distance squared
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[BruteForceSearch::GetGatherNeighbourList]");

  // Compute smoothing lengths of all SPH particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<Nhydro; i++) {

    // Skip over dead particles
    if (partdata[i].itype == dead) continue;

    for (k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rp[k];
    drsqd = DotProduct(dr,dr,ndim);

    if (drsqd <= rsearch*rsearch && Nneib < Nneibmax) neiblist[Nneib++] = i;
    else if (drsqd <= rsearch*rsearch && Nneib == Nneibmax) return -1;

  }
  //-----------------------------------------------------------------------------------------------

  return Nneib;
}



//=================================================================================================
//  BruteForceSearch::SearchBoundaryGhostParticles
/// Search domain to create any required ghost particles near any boundaries.
/// Currently only searches to create periodic or mirror ghost particles.
//=================================================================================================
template <int ndim, template <int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::SearchBoundaryGhostParticles
 (FLOAT tghost,                        ///< [in] Ghost particle 'lifetime'
  DomainBox<ndim> &simbox,             ///< [in] Simulation box structure
  Hydrodynamics<ndim> *hydro)          ///< [inout] Sph object pointer
{
  int i;                               // Particle counter

  // Set all relevant particle counters
  hydro->Nghost         = 0;
  hydro->NPeriodicGhost = 0;
  hydro->Nghostmax      = hydro->Nhydromax - hydro->Nhydro;
  hydro->Ntot           = hydro->Nhydro;


  // If all boundaries are open, immediately return to main loop
  if (simbox.boundary_lhs[0] == openBoundary && simbox.boundary_rhs[0] == openBoundary &&
      simbox.boundary_lhs[1] == openBoundary && simbox.boundary_rhs[1] == openBoundary &&
      simbox.boundary_lhs[2] == openBoundary && simbox.boundary_rhs[2] == openBoundary)
    return;


  debug2("[BruteForceSearch::SearchBoundaryGhostParticles]");


  // Create ghost particles in x-dimension
  //-----------------------------------------------------------------------------------------------
  if ((simbox.boundary_lhs[0] == openBoundary && simbox.boundary_rhs[0] == openBoundary) == 0) {

    for (i=0; i<hydro->Ntot; i++) {
      hydro->CheckXBoundaryGhostParticle(i, tghost, simbox);
    }

    hydro->Ntot = hydro->Nhydro + hydro->Nghost;
  }


  // Create ghost particles in y-dimension
  //-----------------------------------------------------------------------------------------------
  if (ndim >= 2 && (simbox.boundary_lhs[1] == openBoundary &&
                    simbox.boundary_rhs[1] == openBoundary) == 0) {

    for (i=0; i<hydro->Ntot; i++) {
      hydro->CheckYBoundaryGhostParticle(i, tghost, simbox);
    }

    hydro->Ntot = hydro->Nhydro + hydro->Nghost;
  }


  // Create ghost particles in z-dimension
  //-----------------------------------------------------------------------------------------------
  if (ndim == 3 && (simbox.boundary_lhs[2] == openBoundary &&
                    simbox.boundary_rhs[2] == openBoundary) == 0) {

    for (i=0; i<hydro->Ntot; i++) {
      hydro->CheckZBoundaryGhostParticle(i, tghost, simbox);
    }

    hydro->Ntot = hydro->Nhydro + hydro->Nghost;
  }


  // Quit here if we've run out of memory for ghosts
  if (hydro->Ntot > hydro->Nhydromax) {
    string message="Not enough memory for ghost particles";
    ExceptionHandler::getIstance().raise(message);
  }

  hydro->NPeriodicGhost = hydro->Nghost;

  return;
}



//=================================================================================================
//  BruteForceSearch::UpdateAllStarGasForces
/// ..
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UpdateAllStarGasForces
 (int Nhydro,                          ///< [in] No. of SPH particles
  int Ntot,                            ///< [in] Total no. of SPH particles (incl. ghosts)
  Particle<ndim> *part_gen,            ///< [in] SPH particle array pointer
  Hydrodynamics<ndim> *hydro,          ///< [inout] Pointer to SPH ptcl array
  Nbody<ndim> *nbody)                  ///< [inout] Pointer to N-body object
{
  int i;                               // Particle and dimension counters
  int Nneib = 0;                       // No. of (non-dead) neighbours
  int *dummy = 0;                      // Dummy var to satisfy function argument
  int *neiblist;                       // Array of (all) particle ids
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (part_gen);

  debug2("[BruteForceSearch::UpdateAllStarGasForces]");

  // Allocate memory for storing neighbour ids and position data
  neiblist = new int[Nhydro];
  for (i=0; i<Nhydro; i++) {
    if (partdata[i].itype != dead) neiblist[Nneib++] = i;
  }

  // Compute smoothing lengths of all SPH particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<nbody->Nnbody; i++) {

    // Skip over inactive particles
    if (!nbody->nbodydata[i]->active) continue;

    // Compute forces between SPH neighbours (hydro and gravity)
    nbody->CalculateDirectHydroForces(nbody->nbodydata[i], Nneib, 0, neiblist, dummy, hydro);

  }
  //-----------------------------------------------------------------------------------------------

  delete[] neiblist;

  return;
}



#if defined MPI_PARALLEL
//=================================================================================================
//  BruteForceSearch::SearchMpiGhostParticles
/// Compute on behalf of the MpiControl class the ghost particles we need to export to other nodes.
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int BruteForceSearch<ndim,ParticleType>::SearchMpiGhostParticles
 (const FLOAT tghost,                  ///< [in] Expected ghost life-time
  const Box<ndim> &mpibox,             ///< [in] Bounding box of MPI domain
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to Hydrodynamics object
  vector<int> &export_list)            ///< [out] List of particle ids
{
  int i;                               // ..
  int k;                               // ..
  int Nexport = 0;                     // No. of MPI ghosts to export
  FLOAT scattermin[ndim];              // ..
  FLOAT scattermax[ndim];              // ..
  const FLOAT grange = ghost_range*kernrange;
  ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

  // Loop over particles and prepare the ones to export
  //-----------------------------------------------------------------------------------------------
  //for (i=0; i<hydro->Nhydro; i++) {
  for (i=0; i<hydro->Nhydro+hydro->NPeriodicGhost; i++) {
    ParticleType<ndim>& part = partdata[i];

    // Construct maximum cell bounding box depending on particle velocities
    for (k=0; k<ndim; k++) {
      scattermin[k] = part.r[k] + min((FLOAT) 0.0, part.v[k]*tghost) - grange*part.h;
      scattermax[k] = part.r[k] + max((FLOAT) 0.0, part.v[k]*tghost) + grange*part.h;
    }

    // If maximum cell scatter box overlaps MPI domain, open cell
    if (BoxOverlap(ndim, scattermin, scattermax, mpibox.boxmin, mpibox.boxmax)) {
      export_list.push_back(i);
      Nexport++;
    }

  }
  //-----------------------------------------------------------------------------------------------

  return Nexport;
}



//=================================================================================================
//  BruteForceSearch::FindMpiTransferParticles
/// Compute on behalf of the MpiControl class the particles that are outside
/// the domain after a load balancing and need to be transferred to other nodes
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::FindMpiTransferParticles
 (Hydrodynamics<ndim> *hydro,              ///< [in] Pointer to sph class
  vector<vector<int> >& id_export_buffers, ///< [inout] Vector that for each node
                                           ///<         gives the list of particles to export
  vector<int>& all_ids_export_buffer,      ///< [inout] Vector containing all the particles that
                                           ///<         will be exported by this processor
  const vector<int>& potential_nodes,      ///< [in] Potential nodes we might send particles to
  MpiNode<ndim>* mpinodes)                 ///< [in] Array of other mpi nodes
{
  int i;
  int inode;
  int node_number;
  ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());


  // Loop over particles and prepare the ones to export
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<hydro->Nhydro; i++) {
    ParticleType<ndim>& part = partdata[i];

    // Loop over potential domains and see if we need to transfer
    // this particle to them
    //---------------------------------------------------------------------------------------------
    for (inode=0; inode<potential_nodes.size(); inode++) {
      node_number = potential_nodes[inode];

      // If particle belongs to this domain, add to vector and break from loop
      if (ParticleInBox(part,mpinodes[node_number].domain)) {
        id_export_buffers[node_number].push_back(i);
        all_ids_export_buffer.push_back(i);
        break;
      }

    }
    //---------------------------------------------------------------------------------------------

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  BruteForceSearch::FindParticlesToTransfer
/// Compute on behalf of the MpiControl class the particles that are outside
/// the domain after a load balancing and need to be transferred to other nodes
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::FindParticlesToTransfer
 (Hydrodynamics<ndim> *hydro,                ///< [in] Pointer to sph class
  vector<vector<int> >& id_export_buffers,   ///< [inout] List of ids to export for each node
  vector<int>& all_ids_export_buffer,        ///< [inout] List of all ids to export from proc
  const vector<int>& potential_nodes,        ///< [in] Potential nodes we might send particles to
  MpiNode<ndim>* mpinodes)                   ///< [in] Array of other mpi nodes
{
  ParticleType<ndim> *partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray());

  // Loop over particles and prepare the ones to export
  for (int i=0; i<hydro->Nhydro; i++) {
    ParticleType<ndim>& part = partdata[i];

    // Loop over potential domains and see if we need to transfer this particle to them
    for (int inode=0; inode<potential_nodes.size(); inode++) {
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
//  BruteForceSearch::GetExportInfo
/// Get the array with the information that needs to be exported to the given processor
/// (NB: Nproc is ignored at the moment, as we must always export all ptcls to other processors).
//=================================================================================================
template <int ndim, template<int> class ParticleType>
int BruteForceSearch<ndim,ParticleType>::GetExportInfo
 (int Nproc,                           ///< [in] No of processor we want to send the information to
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to sph object
  vector<char >& ptcl_export_buffer,   ///< [inout] Buffer containing particles to be exported
  MpiNode<ndim>& mpinode,              ///< ..
  int rank,                            ///< ..
  int Nmpi)                            ///< [in]  Array with information for the other mpi nodes
{
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray() );

  const bool first_proc = (Nproc==0) || (rank==0 && Nproc==1);

  // Find number of active particles
  if (first_proc) {
    ids_active_particles.clear();
    for (int i=0; i<hydro->Nhydro; i++) {
      if (partdata[i].active) {
        ids_active_particles.push_back(i);
      }
    }
  }

  const int Nactive = ids_active_particles.size();
  const int size_export = Nactive*sizeof(ParticleType<ndim>);

  if (first_proc) {
    ptcl_export_buffer.clear();
    ptcl_export_buffer.reserve((Nmpi-1)*size_export);
    ptcl_export_buffer.resize(size_export);
  }
  else {
    ptcl_export_buffer.resize(ptcl_export_buffer.size()+size_export);
  }

//  //Copy positions of active particles inside arrays
//  int j=0;
//  for (int i=0; i< hydro->Nhydro; i++) {
//    if (partdata[i].active) {
//      for (int k=0; k<ndim; k++)
//        ptcl_export_buffer[j].r[k] = partdata[i].r[k];
//      j++;
//    }
//  }

  //Copy particles to export inside arrays
  int j = (ptcl_export_buffer.size() - size_export)/sizeof(ParticleType<ndim>);
  for (int i=0; i<hydro->Nhydro; i++) {
    if (partdata[i].active) {
      copy(&ptcl_export_buffer[j*sizeof(ParticleType<ndim>)] , &partdata[i]);
      j++;
    }
  }

  return size_export;
}


//=================================================================================================
//  BruteForceSearch::UnpackExported
/// Unpack the information exported from the other processors, contaning the positions
/// of the particles that were exported
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UnpackExported
 (vector<char > &received_array,
  vector<int> &Nbytes_from_proc,
  Hydrodynamics<ndim> *hydro)
{
  int offset = 0;
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray() );

  assert(hydro->NImportedParticles==0);


  //-----------------------------------------------------------------------------------------------
  for (int Nproc = 0; Nproc<Nbytes_from_proc.size(); Nproc++) {

    int N_received_bytes = Nbytes_from_proc[Nproc];
    int N_received_particles = N_received_bytes/sizeof(ParticleType<ndim>);

    //Ensure there is enough memory
    if (hydro->Ntot + N_received_particles > hydro->Nhydromax) {
      ExceptionHandler::getIstance().raise("Error while receiving imported particles: not enough memory!");
    }

//    //Copy particle positions inside SPH main arrays
//    for (int i=0; i<N_received_particles; i++) {
//      for (int k=0; k<ndim; k++)
//        partdata[i+hydro->Ntot].r[k] = received_array[i+offset].r[k];
//    }

    //Copy received particles inside SPH main arrays
    for (int i=0; i<N_received_particles; i++) {
      copy( &partdata[i+hydro->Ntot] , &received_array[i*sizeof(ParticleType<ndim>)+offset]);
    }

    //Update the SPH counters
    hydro->Ntot += N_received_particles;
    hydro->NImportedParticles += N_received_particles;

    //Update the offset
    offset += N_received_bytes;

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  BruteForceSearch::GetBackExportInfo
/// Return the data to transmit back to the other processors (particle acceleration etc.)
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::GetBackExportInfo
 (vector<char >& send_buffer,          ///< [inout] Buffer containing information to send
  vector<int>& Nbytes_from_proc,       ///< ..
  vector<int>& Nbytes_to_proc,         ///< ..
  Hydrodynamics<ndim> *hydro,          ///< [in] Pointer to the SPH object
  int rank)                            ///< ..
{
  int removed_particles=0;
  int Nbytes_received_exported = std::accumulate(Nbytes_from_proc.begin(), Nbytes_from_proc.end(), 0);

  send_buffer.resize(Nbytes_received_exported);

  // Loop over the processors, removing particles as we go
  //-----------------------------------------------------------------------------------------------
  for (int Nproc=0 ; Nproc < Nbytes_from_proc.size(); Nproc++ ) {

    const int N_received_particles = Nbytes_from_proc[Nproc]/sizeof(ParticleType<ndim>);

//    //Copy the accelerations and gravitational potential of the particles
//    ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetSphParticleArray() );
//    int j=0;
//    for (int i=hydro->Nhydro - N_received_particles; i<hydro->Nhydro; i++) {
//      for (int k=0; k<ndim; k++)
//        send_buffer[removed_particles+j].a[k] = partdata[i].a[k];
//      send_buffer[removed_particles+j].gpot = partdata[i].gpot;
//      j++;
//    }

    //Copy the particles inside the send buffer
    ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray() );
    int j=0;
    const int start_index = hydro->Nhydro + hydro->Nghost + removed_particles;
    for (int i=start_index; i<start_index + N_received_particles; i++) {
      copy (&send_buffer[(removed_particles+j)*sizeof(ParticleType<ndim>)],&partdata[i]);
      j++;
    }

    assert(j==N_received_particles);

    removed_particles += j;

    //Decrease the particle counter
    hydro->Ntot -= N_received_particles;
    hydro->NImportedParticles -= N_received_particles;

  }
  //-----------------------------------------------------------------------------------------------

  assert(hydro->NImportedParticles == 0);
  assert(hydro->Ntot == hydro->Nhydro + hydro->Nghost);
  assert(send_buffer.size() == removed_particles*sizeof(ParticleType<ndim>));

  return;
}



//=================================================================================================
//  BruteForceSearch::UnpackReturnedExportInfo
/// Unpack the data that was returned by the other processors, summing the accelerations to the particles
//=================================================================================================
template <int ndim, template<int> class ParticleType>
void BruteForceSearch<ndim,ParticleType>::UnpackReturnedExportInfo
 (vector<char >& received_information,   ///< ..
  vector<int>& recv_displs,              ///< ..
  Hydrodynamics<ndim> *hydro,            ///< ..
  int rank)                              ///< ..
{
  ParticleType<ndim>* partdata = static_cast<ParticleType<ndim>* > (hydro->GetParticleArray() );

//  //For each particle, sum up the accelerations returned by other processors
//  for (int i=0; i< ids_active_particles.size(); i++ ) {
//    const int j = ids_active_particles[i];
//
//    for(int Nproc=0; Nproc<recv_displs.size(); Nproc++ ) {
//
//      if (rank==Nproc)
//        continue;
//
//      for (int k=0; k<ndim; k++)
//        partdata[j].a[k] += received_information[i+recv_displs[Nproc] ].a[k];
//      partdata[j].gpot += received_information[i+recv_displs[Nproc] ].gpot;
//    }
//
//  }

  // For each particle, sum up some of the quantities returned by other processors
  //-----------------------------------------------------------------------------------------------
  for (int i=0; i<ids_active_particles.size(); i++) {
    const int j = ids_active_particles[i];

    for (int Nproc=0; Nproc<recv_displs.size(); Nproc++) {

      if (rank == Nproc) continue;

      ParticleType<ndim>* received_particle = reinterpret_cast<ParticleType<ndim>*>
        (&received_information[i * sizeof(ParticleType<ndim>) + recv_displs[Nproc] ]);


      assert(partdata[j].iorig == received_particle->iorig);

      for (int k=0; k<ndim; k++) {
        partdata[j].a[k] += received_particle->a[k];
        partdata[j].agrav[k] += received_particle->agrav[k];
      }
      partdata[j].gpot += received_particle->gpot;
      partdata[j].dudt += received_particle->dudt;
      partdata[j].div_v += received_particle->div_v;
      partdata[j].levelneib = max(partdata[j].levelneib, received_particle->levelneib);

    }

  }
  //-----------------------------------------------------------------------------------------------

}
#endif



template class BruteForceSearch<1,SphParticle>;
template class BruteForceSearch<2,SphParticle>;
template class BruteForceSearch<3,SphParticle>;
template class BruteForceSearch<1,GradhSphParticle>;
template class BruteForceSearch<2,GradhSphParticle>;
template class BruteForceSearch<3,GradhSphParticle>;
template class BruteForceSearch<1,SM2012SphParticle>;
template class BruteForceSearch<2,SM2012SphParticle>;
template class BruteForceSearch<3,SM2012SphParticle>;
template class BruteForceSearch<1,MeshlessFVParticle>;
template class BruteForceSearch<2,MeshlessFVParticle>;
template class BruteForceSearch<3,MeshlessFVParticle>;
