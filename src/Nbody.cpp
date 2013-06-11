//=============================================================================
//  Nbody.cpp
//  Contains main N-body class functions.
//=============================================================================

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Precision.h"
#include "NbodyParticle.h"
#include "StarParticle.h"
#include "SystemParticle.h"
#include "Parameters.h"
#include "Nbody.h"
#include "SphKernel.h"
#include "Debug.h"
#include "InlineFuncs.h"
using namespace std;


//template <int ndim>
//const FLOAT Nbody<ndim>::invndim;
//const int Nbody<ndim>::vdim;


//=============================================================================
//  Nbody::Nbody
/// Nbody class constructor
//=============================================================================
template <int ndim>
Nbody<ndim>::Nbody(int nbody_softening_aux, int sub_systems_aux, 
                   DOUBLE nbody_mult_aux, string KernelName, int Npec):
  nbody_softening(nbody_softening_aux),
  sub_systems(sub_systems_aux),
  nbody_mult(nbody_mult_aux),
  kerntab(TabulatedKernel<ndim>(KernelName)),
  Nstar(0),
  Nstarmax(0),
  Nsystem(0),
  Nsystemmax(0),
  Nnbody(0),
  Nnbodymax(0),
  allocated(false),
  Npec(1)
{
}



//=============================================================================
//  Nbody::AllocateMemory
/// Allocate all memory required for stars and N-body system particles.
//=============================================================================
template <int ndim>
void Nbody<ndim>::AllocateMemory(int N)
{
  debug2("[Nbody::AllocateMemory]");

  if (N > Nstarmax) {
    if (allocated) DeallocateMemory();
    Nstar = N;
    Nstarmax = N;
    Nsystem = N;
    Nsystemmax = N;
    Nnbody = N;
    Nnbodymax = Nstarmax + Nsystemmax;
    nbodydata = new struct NbodyParticle<ndim>*[Nnbodymax];
    stardata = new struct StarParticle<ndim>[Nstarmax];
    system = new struct SystemParticle<ndim>[Nsystemmax];
    allocated = true;
  }

  return;
}
 


//=============================================================================
//  Nbody::DeallocateMemory
/// Deallocate all N-body memory.
//=============================================================================
template <int ndim>
void Nbody<ndim>::DeallocateMemory(void)
{
  debug2("[Nbody::DeallocateMemory]");

  if (allocated) {
    delete[] system;
    delete[] stardata;
    delete[] nbodydata;
  }
  allocated = false;

  return;
}


//=============================================================================
//  Systems.cpp
//  Functions for integrating systems internal motions
//=============================================================================


#include <algorithm>
#include "SystemParticle.h"



//=============================================================================
//  Nbody::IntegrateInternalMotion
/// This function integrates the internal motion of a system. First integrates
/// the internal motion of its sub-systems by recursively calling their method,
/// then integrates the COM of the sub-systems.
//=============================================================================
template <int ndim>
void Nbody<ndim>::IntegrateInternalMotion(
    SystemParticle<ndim>* system, /// [inout] System that we wish to integrate the internal motion
    DOUBLE tlocal_end             /// [in]    Time to integrate the internal motion for
    )
{
  int Nchildren = system->Nchildren;
  NbodyParticle<ndim>** children = system->children;

  //-------------------------------------------------
  //First loop over children and, if they are systems, integrate
  //their internal motion
  //------------------------------------------------
  for (int ichild=0; ichild<Nchildren; ichild++) {
    if (children[ichild]->Ncomp > 1)
      //The cast is needed because the function is defined only in
      //SystemParticle, not in NbodyParticle
      //The safety of the cast relies on the correctness of the Ncomp value
      IntegrateInternalMotion(static_cast<SystemParticle<ndim>* > (children[ichild]), tlocal_end);
  }


  //-------------------------------------------------
  //Now integrate our internal motion
  //-------------------------------------------------
  DOUBLE tlocal=0.;
  DOUBLE dt;
  int nlocal_steps=0;

  NbodyParticle<ndim> children_local[Ncompmax]; //Local copy of children
  NbodyParticle<ndim>* children_local_pointers[Ncompmax]; //Array of pointers to local copy of children
  //Make local copies of children (in COM frame)
  for (int ichild=0; ichild<Nchildren; ichild++) {

    //TODO: should probably use some fancy feature of C++ (copy constructor?)
    children_local[ichild]=*(children[ichild]);

    for (int k=0; k<ndim; k++) children_local[ichild].r0[k] -=
        system->r0[k];
    for (int k=0; k<ndim; k++) children_local[ichild].v0[k] -=
        system->v0[k];
    for (int k=0; k<ndim; k++) children_local[ichild].a0[k] -=
        system->a0[k];
    for (int k=0; k<ndim; k++) children_local[ichild].adot0[k] -=
        system->adot0[k];
//    for (int k=0; k<ndim; k++) children_local[ichild].a2dot0[k] -=
//        system->a2dot0[k];
//    for (int k=0; k<ndim; k++) children_local[ichild].a3dot0[k] -=
//        system->a3dot0[k];


    children_local_pointers[ichild] = &(children_local[ichild]);

  }


  do {
    //Calculate time-step
    dt = std::min(big_number, tlocal_end-tlocal);
    for (int ichild=0; ichild<Nchildren; ichild++) {
      dt = std::min(dt, Timestep(children_local_pointers[ichild]));
    }
    tlocal += dt;
    nlocal_steps +=1;

    //Advance position and velocities
    AdvanceParticles(nlocal_steps, Nchildren, children_local_pointers, dt);

    //Iteration loop
    for (int it=0; it<Npec; it++) {

      //Zero all acceleration terms
      for (int ichild=0; ichild<Nchildren; ichild++) {
        for (int k=0; k<ndim; k++) children_local[ichild].a[k] = 0.0;
        for (int k=0; k<ndim; k++) children_local[ichild].a2dot[k] = 0.0;
        for (int k=0; k<ndim; k++) children_local[ichild].a3dot[k] = 0.0;
        children_local[ichild].gpot = 0.0;
      }
    }

    //Calculate forces, derivatives, ...
    CalculateDirectGravForces(Nchildren, children_local_pointers);

    //Apply correction terms
    CorrectionTerms(nlocal_steps, Nchildren, children_local_pointers, dt);

    //Set end-of-step variables
    EndTimestep(nlocal_steps, Nchildren, children_local_pointers);



  } while(tlocal < tlocal_end);

  //Copy systems back to our children
  for (int ichild=0; ichild<Nchildren; ichild++) {
    for (int k=0; k<ndim; k++) children[ichild]->r[k] +=
         children_local[ichild].r[k];
    for (int k=0; k<ndim; k++) children[ichild]->v[k] +=
         children_local[ichild].v[k];
    for (int k=0; k<ndim; k++) children[ichild]->a[k] +=
         children_local[ichild].a[k];
    for (int k=0; k<ndim; k++) children[ichild]->adot[k] +=
         children_local[ichild].adot[k];
    for (int k=0; k<ndim; k++) children[ichild]->a2dot[k] +=
         children_local[ichild].a2dot[k];
    for (int k=0; k<ndim; k++) children[ichild]->a3dot[k] +=
         children_local[ichild].a3dot[k];

  }

}


template class Nbody<1>;
template class Nbody<2>;
template class Nbody<3>;
