//=============================================================================
//  SphSnapshot.cpp
//  Contains all snapshot functions
//=============================================================================


#include <ctime>
#include <cstdio>
#include <iostream>
#include "Exception.h"
#include "SphSnapshot.h"
#include "Sph.h"
#include "SphParticle.h"
#include "Simulation.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "UnitInfo.h"
#include "HeaderInfo.h"
using namespace std;



//=============================================================================
//  SphSnapshotBase::SphSnapshotFactory
/// Creates and returns a new snapshot object based on dimensionality of sim.
//=============================================================================
SphSnapshotBase* SphSnapshotBase::SphSnapshotFactory
(string filename,                   ///< Filename containing snapshot data
 SimulationBase* sim,               ///< Simulation object pointer
 int ndim)                          ///< Dimensionality of simulation
{
  if (ndim==1)
    return new SphSnapshot<1>(filename, sim);
  else if (ndim==2)
    return new SphSnapshot<2>(filename, sim);
  else if (ndim==3)
    return new SphSnapshot<3>(filename, sim);
  return NULL;
};



//=============================================================================
//  SphSnapshotBase::SphSnapshotBase
/// Constructor for SphSnapshotBase class.
//=============================================================================
SphSnapshotBase::SphSnapshotBase(SimUnits* _units, string auxfilename):
    units(_units)
{
  allocated = false;
  allocatedstar = false;
  allocatedsph = false;
  nallocatedsph = 0;
  nallocatedstar = 0;
  Nsph = 0;
  Nsphmax = 0;
  Nstar = 0;
  Nstarmax = 0;
  t = 0.0;
  LastUsed = time(NULL);
  if (auxfilename != "") filename = auxfilename;
}



//=============================================================================
//  SphSnapshot::SphSnapshot
/// Constructor for SphSnapshot class.
//=============================================================================
template <int ndims>
SphSnapshot<ndims>::SphSnapshot
(string filename,                   ///< Filename containing snapshot data
 SimulationBase* sim):              ///< Simulation object pointer
SphSnapshotBase(&(sim->simunits), filename),
simulation(static_cast<Simulation<ndims>* > (sim))
{
  this->ndim = ndims;

  this->fileform = sim->GetParam("in_file_form");

  // Computes how numbers we need to store for each sph/star particle
  nneededsph = 3*ndims + 5;
  nneededstar = 3*ndims + 2;

  if (filename != "") {
    HeaderInfo info;
    info = sim->ReadHeaderSnapshotFile(filename, this->fileform);
    this->t = info.t;
    this->Nstar = info.Nstar;
    this->Nsph = info.Nsph;
  }
}


//=============================================================================
//  SphSnapshotBase::~SphSnapshotBase
/// Deallocate any heap arrays to avoid leaking memory.
//=============================================================================
SphSnapshotBase::~SphSnapshotBase()
{
  DeallocateBufferMemory();
}



//=============================================================================
//  SphSnapshotBase::AllocateBufferMemory
/// Allocate memory for current snapshot.  Only allocates single precision 
/// to minimise memory use, even if compiled with double precision.
/// Wrapper around AllocateBufferMemoryStar and AllocateBufferMemorySph
//=============================================================================
void SphSnapshotBase::AllocateBufferMemory(void)
{
  debug2("[SphSnapshotBase::AllocateBufferMemory]");

  AllocateBufferMemoryStar();
  AllocateBufferMemorySph();

  allocated = true;

  return;

}

//=============================================================================
//  SphSnapshotBase::AllocateBufferMemoryStar
/// Allocate memory for stars in current snapshot.
//=============================================================================
void SphSnapshotBase::AllocateBufferMemoryStar(void) 
{
  // If memory is already allocated and more memory is needed for more 
  // particles, deallocate now before reallocating.
  if (allocatedstar) {
    if (Nstar > Nstarmax)
      DeallocateBufferMemoryStar();
    else
      return;
  }

  // Allocate memory for all vector quantities depending on dimensionality
  if (ndim == 1) {
    xstar = new float[Nstar];
    vxstar = new float[Nstar];
    axstar = new float[Nstar];
  }
  else if (ndim == 2) {
    xstar = new float[Nstar];
    ystar = new float[Nstar];
    vxstar = new float[Nstar];
    vystar = new float[Nstar];
    axstar = new float[Nstar];
    aystar = new float[Nstar];

  }
  else if (ndim == 3) {
    xstar = new float[Nstar];
    ystar = new float[Nstar];
    zstar = new float[Nstar];
    vxstar = new float[Nstar];
    vystar = new float[Nstar];
    vzstar = new float[Nstar];
    axstar = new float[Nstar];
    aystar = new float[Nstar];
    azstar = new float[Nstar];

  }

  // Stars scalar quantities
  mstar = new float[Nstar];
  hstar = new float[Nstar];

  // Record 3 vectors of size ndim (r,v,a) and 2 scalars (m,h)
  nallocatedstar = 3*ndim + 2;
  allocatedstar = true;
  Nstarmax = Nstar;

  return;
}



//=============================================================================
//  SphSnapshotBase::AllocateBufferMemorySph
/// Allocate memory for sph particles in current snapshot.
//=============================================================================
void SphSnapshotBase::AllocateBufferMemorySph(void) 
{
  // If memory already allocated and more memory is needed for more particles,
  // deallocate now before reallocating.
  if (allocatedsph) {
    if (Nsph > Nsphmax)
      DeallocateBufferMemorySph();
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
  nallocatedsph = 3*ndim + 5;
  allocatedsph = true;
  Nsphmax = Nsph;

  return;
}



//=============================================================================
//  SphSnapshotBase::DeallocateBufferMemory
/// Deallocate memory for current snapshot.
//=============================================================================
void SphSnapshotBase::DeallocateBufferMemory(void)
{
  debug2("[SphSnapshotBase::DeallocateBufferMemory]");

  DeallocateBufferMemorySph();
  DeallocateBufferMemoryStar();

  allocated=false;
  return;
}



//=============================================================================
//  SphSnapshotBase::DeallocateBufferMemorySph
/// Deallocate sph particles memory for current snapshot.
//=============================================================================
void SphSnapshotBase::DeallocateBufferMemorySph(void)
{
  // If we are not allocated, return immediately,
  // to avoid double deleting a pointer
  if (!allocatedsph)
    return;

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
  
  allocatedsph = false;
  nallocatedsph = 0;
  
  return;
}



//=============================================================================
//  SphSnapshotBase::DeallocateBufferMemoryStar
/// Deallocate star particles memory for current snapshot.
//=============================================================================
void SphSnapshotBase::DeallocateBufferMemoryStar(void)
{
  // If we are not allocated, return immediately,
  // to avoid double deleting a pointer
  if (!allocatedstar)
    return;

  // Deallocate scalar array memory
  delete[] hstar;
  delete[] mstar;

  // Deallocate vector array memory
  if (ndim == 1) {
    delete[] axstar;
    delete[] vxstar;
    delete[] xstar;
  }
  else if (ndim == 2) {
    delete[] aystar;
    delete[] axstar;
    delete[] vystar;
    delete[] vxstar;
    delete[] ystar;
    delete[] xstar;
  }
  else if (ndim == 3) {
    delete[] azstar;
    delete[] aystar;
    delete[] axstar;
    delete[] vzstar;
    delete[] vystar;
    delete[] vxstar;
    delete[] zstar;
    delete[] ystar;
    delete[] xstar;
  }
  
  allocatedstar = false;
  nallocatedstar = 0;
  
  return;
}



//=============================================================================
//  SphSnapshotBase::CalculateMemoryUsage
/// Returns no. of bytes allocated for current snapshot
//=============================================================================
int SphSnapshotBase::CalculateMemoryUsage(void)
{
  return Nsph*nallocatedsph*sizeof(float) + 
    Nstar*nallocatedstar*sizeof(float);
}



//=============================================================================
//  SphSnapshotBase::CalculatePredictedMemoryUsage
/// Returns no. of bytes that the current snapshot would use, if allocated
//=============================================================================
int SphSnapshotBase::CalculatePredictedMemoryUsage(void)
{
  return Nsph*nneededsph*sizeof(float) + Nstar*nneededstar*sizeof(float);
}



//=============================================================================
//  SphSnapshot::CopyDataFromSimulation
/// Copy particle data from main memory to current snapshot arrays.
//=============================================================================
template <int ndims>
void SphSnapshot<ndims>::CopyDataFromSimulation()
{
  debug2("[SphSnapshotBase::CopyDataFromSimulation]");
  SphParticle<ndims>* sphaux;
  StarParticle<ndims>* staraux;

  // Reset the species
  _species.clear();

  // Read which species are there
  if (simulation->sph != NULL && simulation->sph->sphdata != NULL) {
    sphaux = simulation->sph->sphdata;
    Nsph = simulation->sph->Nsph;
    if (Nsph != 0) {
      _species.push_back("sph");
    }
  }
  if (simulation->nbody != NULL && simulation->nbody->stardata != NULL) {
    staraux = simulation->nbody->stardata;
    Nstar = simulation->nbody->Nstar;
    if (Nstar != 0) {
     _species.push_back("star");
    }
  }

  AllocateBufferMemory();

  // Loop over all SPH particles and record particle data
  for (int i=0; i<Nsph; i++) {

    if (ndim == 1) {
      x[i] = (float) sphaux[i].r[0];
      vx[i] = (float) sphaux[i].v[0];
      ax[i] = (float) pow(2,simulation->level_step - sphaux[i].level)*
	simulation->timestep; //(float) sphaux[i].a[0];
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
    m[i] = (float) sphaux[i].m;
    h[i] = (float) sphaux[i].h;
    rho[i] = (float) sphaux[i].rho;
    u[i] = (float) sphaux[i].u;
    dudt[i] = sphaux[i].dudt;

  }

  // Loop over star particles and record particle data
  for (int i=0; i<Nstar; i++) {
    if (ndim == 1) {
      xstar[i] = (float) staraux[i].r[0];
      vxstar[i] = (float) staraux[i].v[0];
      axstar[i] = (float) staraux[i].a[0];
    }
    else if (ndim == 2) {
      xstar[i] = (float) staraux[i].r[0];
      ystar[i] = (float) staraux[i].r[1];
      vxstar[i] = (float) staraux[i].v[0];
      vystar[i] = (float) staraux[i].v[1];
      axstar[i] = (float) staraux[i].a[0];
      aystar[i] = (float) staraux[i].a[1];
    }
    else if (ndim == 3) {
      xstar[i] = (float) staraux[i].r[0];
      ystar[i] = (float) staraux[i].r[1];
      zstar[i] = (float) staraux[i].r[2];
      vxstar[i] = (float) staraux[i].v[0];
      vystar[i] = (float) staraux[i].v[1];
      vzstar[i] = (float) staraux[i].v[2];
      axstar[i] = (float) staraux[i].a[0];
      aystar[i] = (float) staraux[i].a[1];
      azstar[i] = (float) staraux[i].a[2];
    }

    mstar[i] = (float) staraux[i].m;
    hstar[i] = (float) staraux[i].h;
  }

  LastUsed = time(NULL);
  return;
}



//=============================================================================
//  SphSnapshotBase::GetRealType
/// Convert the given type into the 'real' type, meaning, if we pass default,
/// gives back the true underlying type
//=============================================================================
string SphSnapshotBase::GetRealType(string type) 
{
  // Default is: if there is only one species, then we use that one
  // Otherwise, we return sph particles
  if (type=="default") {
    if (GetNTypes() == 0){
      string message = "Error: the requested simulation has no species!!!";
      ExceptionHandler::getIstance().raise(message);
    }
    else if (GetNTypes() == 1)
      type = GetSpecies(0);
    else
      type = "sph";
  }

  return type;
}



//=============================================================================
//  SphSnapshotBase::ExtractArray
/// Returns pointer to required array stored in snapshot buffer memory.
/// Currently also returns scaling factors for that array.
//=============================================================================
UnitInfo SphSnapshotBase::ExtractArray
(string name,                       ///< ..
 string type,                       ///< ..
 float** out_array,                 ///< ..
 int* size_array,                   ///< ..
 float& scaling_factor,             ///< ..
 string RequestedUnit)              ///< ..
{
  string unitname;                  // ..
  UnitInfo unitinfo;                // ..
  SimUnit* unit;                    // Unit pointer

  // Zero initial array pointers and size
  *out_array = NULL;
  *size_array = 0;


  // Check that the memory is allocated. If not, fails very rumorously
  if (!allocated){
    cout << "Error: requested a snapshot that is not allocated!!!!" << endl;
    cout << "This means there's a bug in the memory management: please inform the authors" << endl;
    exit(-2);
  }

  // Set last time used
  LastUsed = time(NULL);

  // Get the real underlying type, in case we are passing "default"
  type = GetRealType(type);

  // Check type
  if (type != "sph" && type != "star") {
    string message = "Error: the type " + type + " was not recognized!";
    ExceptionHandler::getIstance().raise(message);
  }


  // If array type and name is valid, pass pointer to array and also set unit
  if (name == "x") {
    if (type == "sph") 
      *out_array = x;
    else if (type == "star") 
      *out_array = xstar;
    unit = &(units->r);
  }
  else if (name == "y") {
    if (type == "sph")
      *out_array = y;
    else if (type == "star")
      *out_array = ystar;
    unit = &(units->r);
  }
  else if (name == "z") {
    if (type == "sph")
      *out_array = z;
    else if (type == "star")
      *out_array = zstar;
    unit = &(units->r);
  }
  else if (name == "vx") {
    if (type == "sph")
      *out_array = vx;
    else if (type == "star")
      *out_array= vxstar;
    unit = &(units->v);
  }
  else if (name == "vy") {
    if (type == "sph")
      *out_array = vy;
    else if (type == "star")
      *out_array = vystar;
    unit = &(units->v);
  }
  else if (name == "vz") {
    if (type == "sph")
      *out_array = vz;
    else if (type == "star")
      *out_array = vzstar;
    unit = &(units->v);
  }
  else if (name == "ax") {
    if (type == "sph")
      *out_array = ax;
    else if (type == "star")
      *out_array = axstar;
    unit = &(units->a);
  }
  else if (name == "ay") {
    if (type == "sph")
      *out_array = ay;
    else if (type == "star")
      *out_array = aystar;
    unit = &(units->a);
  }
  else if (name == "az") {
    if (type == "sph")
      *out_array = az;
    else if (type == "star")
      *out_array = azstar;
    unit = &(units->a);
  }
  else if (name == "m") {
    if (type == "sph")
      *out_array = m;
    else if (type == "star")
      *out_array = mstar;
    unit = &(units->m);
  }
  else if (name == "h") {
    if (type == "sph")
      *out_array = h;
    else if (type == "star")
      *out_array = hstar;
    unit = &(units->r);
  }
  else if (name == "rho") {
    if (type == "sph") {
      *out_array = rho;
    }
    unit = &(units->rho);
  }
  else if (name == "u") {
    if (type == "sph")
      *out_array = u;
    unit = &(units->u);
  }
  else if (name == "dudt") {
    if (type == "sph")
      *out_array = dudt;
    unit = &(units->dudt);
  }
  else {
    string message = "Warning: the selected array: " + name + 
      " has not been recognized";
    ExceptionHandler::getIstance().raise(message);
    *size_array = 0;
  }


  // Check that we did not get a NULL
  if (out_array == NULL) {
    string message;
    if (type == "star" && (name == "rho" || name == "u" || name == "dudt"))
      message = "Error: for stars, you cannot request the array " + name;
    else
      message = "Error: the requested array: " + name + " is not allocated! "
        "Probably a dimensionality problem";
    ExceptionHandler::getIstance().raise(message);
  }

  // Set the size now that we have the array
  if (type == "sph")
    *size_array = Nsph;
  else if (type == "star")
    *size_array = Nstar;

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

  unitinfo.name = unitname;
  unitinfo.label = label;

  return unitinfo;
}



//=============================================================================
//  SphSnapshot::ReadSnapshot
/// Read snapshot into main memory and then copy into snapshot buffer.
//=============================================================================
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

  return;
}
