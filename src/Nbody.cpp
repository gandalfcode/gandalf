//=============================================================================
//  Nbody.cpp
//  Contains main N-body class functions.
//=============================================================================


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
				DOUBLE nbody_mult_aux, string KernelName):
  nbody_softening(nbody_softening_aux),
  sub_systems(sub_systems_aux),
  nbody_mult(nbody_mult_aux),
  kerntab(TabulatedKernel<ndim>(KernelName))
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
//  Nbody::CalculateDirectSoftenedGravForces
/// Calculate all star-star force contributions for active systems using 
/// direct summation with SPH kernel-softened gravity.
//=============================================================================
//template <int ndim>
//void Nbody<ndim>::CalculateDirectSoftenedGravForces(void)
//{
//  debug2("[Nbody::CalculateDirectSoftenedGravForces]");
//
//  return;
//}



template class Nbody<1>;
template class Nbody<2>;
template class Nbody<3>;
