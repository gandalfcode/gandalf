// ============================================================================
// SphSnapshot.cpp
// Contains all snapshot functions
// ============================================================================


#include <ctime>
#include <cstdio>
#include <iostream>
#include "Exception.h"
#include "SphSnapshot.h"
#include "Sph.h"
#include "SphParticle.h"
#include "SphSimulation.h"
#include "Debug.h"
using namespace std;

SphSnapshotBase* SphSnapshotBase::SphSnapshotFactory(string filename, SphSimulationBase* sim, int ndim) {
    if (ndim==1)
      return new SphSnapshot<1>(filename, sim);
    else if (ndim==2)
      return new SphSnapshot<2>(filename, sim);
    else if (ndim==3)
      return new SphSnapshot<3>(filename, sim);
    return NULL;
  };


// ============================================================================
// SphSnapshotBase::SphSnapshotBase
// ============================================================================
SphSnapshotBase::SphSnapshotBase(string auxfilename)
{
  allocated = false;
  nallocated = 0;
  Nsph = 0;
  t = 0.0;
  if (auxfilename != "") filename = auxfilename;
  LastUsed = time(NULL);
}

template <int ndims>
SphSnapshot<ndims>::SphSnapshot (string filename, SphSimulationBase* sim):
SphSnapshotBase(filename),
simulation(static_cast<SphSimulation<ndims>* > (sim))
{
  this->ndim = ndims;
}


// ============================================================================
// SphSnapshotBase::~SphSnapshotBase
// ============================================================================
SphSnapshotBase::~SphSnapshotBase()
{
}



// ============================================================================
// SphSnapshotBase::AllocateBufferMemory
// Allocate memory for current snapshot.  Only allocates single precision 
// to minimise memory use, even if compiled with double precision. 
// ============================================================================
void SphSnapshotBase::AllocateBufferMemory(void)
{
  debug2("[SphSnapshotBase::AllocateBufferMemory]");

  // If memory already allocated and more memory is needed for more particles,
  // deallocate now before reallocating.
  if (allocated) {
    if (Nsph > Nmax)
      DeallocateBufferMemory();
    else
      return;
  }

  // Allocate memory for all vector quantities depending on dimensionality
  if (ndim == 1) {
    x = new float[Nsph];
    vx = new float[Nsph];
    ax = new float[Nsph];
  }
  else if (ndim == 2) {
    x = new float[Nsph];
    y = new float[Nsph];
    vx = new float[Nsph];
    vy = new float[Nsph];
    ax = new float[Nsph];
    ay = new float[Nsph];
  }
  else if (ndim == 3) {
    x = new float[Nsph];
    y = new float[Nsph];
    z = new float[Nsph];
    vx = new float[Nsph];
    vy = new float[Nsph];
    vz = new float[Nsph];
    ax = new float[Nsph];
    ay = new float[Nsph];
    az = new float[Nsph];
  }
  
  // Allocate memory for other scalar quantities
  m = new float[Nsph];
  h = new float[Nsph];
  rho = new float[Nsph];
  u = new float[Nsph];
  dudt = new float[Nsph];

  // Record 3 vectors of size ndim (r,v,a) and 5 scalars (m,h,rho,u,dudt)
  nallocated = 3*ndim + 5;
  allocated = true;
  Nmax = Nsph;

  return;
}



// ============================================================================
// SphSnapshotBase::DeallocateBufferMemory
// Deallocate memory for current snapshot.
// ============================================================================
void SphSnapshotBase::DeallocateBufferMemory(void)
{
  debug2("[SphSnapshotBase::DeallocateBufferMemory]");

  // Deallocate scalar array memory
  delete[] dudt;
  delete[] u;
  delete[] rho;
  delete[] h;
  delete[] m;

  // Deallocate vector array memory
  if (ndim == 1) {
    delete[] ax;
    delete[] vx;
    delete[] x;
  }
  else if (ndim == 2) {
    delete[] ay;
    delete[] ax;
    delete[] vy;
    delete[] vx;
    delete[] y;
    delete[] x;
  }
  else if (ndim == 3) {
    delete[] az;
    delete[] ay;
    delete[] ax;
    delete[] vz;
    delete[] vy;
    delete[] vx;
    delete[] z;
    delete[] y;
    delete[] x;
  }
  
  allocated = false;
  nallocated = 0;

  return;
}



// ============================================================================
// SphSnapshotBase::CalculateMemoryUsage
// Returns no. of bytes required for current snapshot
// ============================================================================
int SphSnapshotBase::CalculateMemoryUsage(void)
{
  return Nsph*nallocated*sizeof(float);
}



// ============================================================================
// SphSnapshotBase::CopyDataFromSimulation
// Copy particle data from main memory to current snapshot arrays.
// ============================================================================
template <int ndims>
void SphSnapshot<ndims>::CopyDataFromSimulation()
{
  debug2("[SphSnapshotBase::CopyDataFromSimulation]");

  SphParticle<ndims>* sphaux = simulation->sph->sphdata;

  Nsph = simulation->sph->Nsph;

  AllocateBufferMemory();

  // Loop over all SPH particles and record particle data
  for (int i=0; i<Nsph; i++) {

    if (ndim == 1) {
      x[i] = (float) sphaux[i].r[0];
      vx[i] = (float) sphaux[i].v[0];
      ax[i] = (float) sphaux[i].a[0];
    }
    else if (ndim == 2) {
      x[i] = (float) sphaux[i].r[0];
      y[i] = (float) sphaux[i].r[1];
      vx[i] = (float) sphaux[i].v[0];
      vy[i] = (float) sphaux[i].v[1];
      ax[i] = (float) sphaux[i].a[0];
      ay[i] = (float) sphaux[i].a[1];
    }
    else if (ndim == 3) {
      x[i] = (float) sphaux[i].r[0];
      y[i] = (float) sphaux[i].r[1];
      z[i] = (float) sphaux[i].r[2];
      vx[i] = (float) sphaux[i].v[0];
      vy[i] = (float) sphaux[i].v[1];
      vz[i] = (float) sphaux[i].v[2];
      ax[i] = (float) sphaux[i].a[0];
      ay[i] = (float) sphaux[i].a[1];
      az[i] = (float) sphaux[i].a[2];
    }

    m[i] = (float) sphaux[i].press/powf(sphaux[i].rho,1.4); //sphaux[i].m;
    h[i] = (float) sphaux[i].h;
    rho[i] = (float) sphaux[i].rho;
    u[i] = (float) sphaux[i].u;
    dudt[i] = (float) sphaux[i].dudt;

  }

  LastUsed = time(NULL);
  return;
}



// ============================================================================
// SphSnapshotBase::ExtractArray
// Returns pointer to required array stored in snapshot buffer memory.
// Currently also returns scaling factors for that array.
// ============================================================================
void SphSnapshotBase::ExtractArray(string name, float** out_array, int* size_array,
                               float& scaling_factor, string RequestedUnit)
{

  //Check that the memory is allocated. If not, fails very rumorously
  if (!allocated){
    cout << "Error: requested a snapshot that it's not allocated!!!!" << endl;
    cout << "This means there's a bug in the memory management: please inform the authors" << endl;
    exit(-2);
  }

  SimUnit* unit;                            // Unit pointer

  LastUsed = time(NULL);

  // If array name is valid, pass pointer to array and also set unit
  if (name == "x") {
    *out_array = x;
    unit = &(units->r);
  }
  else if (name == "y") {
    *out_array = y;
    unit = &(units->r);
  }
  else if (name == "z") {
    *out_array = z;
    unit = &(units->r);
  }
  else if (name == "vx") {
    *out_array = vx;
    unit = &(units->v);
  }
  else if (name == "vy") {
    *out_array = vy;
    unit = &(units->v);
  }
  else if (name == "vz") {
    *out_array = vz;
    unit = &(units->v);
  }
  else if (name == "ax") {
    *out_array = ax;
    unit = &(units->a);
  }
  else if (name == "ay") {
    *out_array = ay;
    unit = &(units->a);
  }
  else if (name == "az") {
    *out_array = az;
    unit = &(units->a);
  }
  else if (name == "m") {
    *out_array = m;
    unit = &(units->m);
  }
  else if (name == "h") {
    *out_array = h;
    unit = &(units->r);
  }
  else if (name == "rho") {
    *out_array = rho;
    unit = &(units->rho);
  }
  else if (name == "u") {
    *out_array = u;
    unit = &(units->u);
  }
  else if (name == "dudt") {
    *out_array = dudt;
    unit = &(units->dudt);
  }
  else {
    string message = "Warning: the selected array: " + name + 
      " has not been recognized";
    ExceptionHandler::getIstance().raise(message);
    *size_array = 0;
  }

  //set the size now that we have the array
  *size_array = Nsph;

  // If no new unit is requested, pass the default scaling values.
  // Otherwise, calculate new scaling factor plus latex label.
  if (RequestedUnit == "default") {
    unitname = unit->outunit;
    RequestedUnit=unitname;
  }
  else {
    unitname=RequestedUnit;
  }
  label = unit->LatexLabel(RequestedUnit);
  scaling_factor = unit->OutputScale(RequestedUnit);

  return;
}


// ============================================================================
// SphSnapshotBase::ReadSnapshot
// Read snapshot into main memory and then copy into snapshot buffer.
// ============================================================================
template <int ndims>
void SphSnapshot<ndims>::ReadSnapshot(string format)
{
  debug2("[SphSnapshotBase::ReadSnapshot]");

  // Set pointer to units object
  units = &(simulation->simunits);

  // Read simulation into main memory
  simulation->ReadSnapshotFile(filename, format);

  // Now copy from main memory to current snapshot
  CopyDataFromSimulation();

  // Record simulation snapshot time
  t = simulation->t;

}
