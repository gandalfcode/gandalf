//=================================================================================================
//  SphSnapshot.cpp
//  Contains all functions for managing SPH snapshot objects.
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


#include <ctime>
#include <cstdio>
#include <iostream>
#include "Exception.h"
#include "SphSnapshot.h"
#include "Sph.h"
#include "Particle.h"
#include "Simulation.h"
#include "Debug.h"
#include "InlineFuncs.h"
#include "UnitInfo.h"
#include "HeaderInfo.h"
#include "Precision.h"
using namespace std;



//=================================================================================================
//  SphSnapshotBase::SphSnapshotFactory
/// Creates and returns a new snapshot object based on dimensionality of sim.
//=================================================================================================
SphSnapshotBase* SphSnapshotBase::SphSnapshotFactory
 (string filename,                   ///< Filename containing snapshot data
  SimulationBase* sim,               ///< Simulation object pointer
  int ndim)                          ///< Dimensionality of simulation
{
  if (ndim == 1) {
    return new SphSnapshot<1>(filename, sim);
  }
  else if (ndim == 2) {
    return new SphSnapshot<2>(filename, sim);
  }
  else if (ndim == 3) {
    return new SphSnapshot<3>(filename, sim);
  }
  return NULL;
};



//=================================================================================================
//  SphSnapshotBase::SphSnapshotBase
/// Constructor for SphSnapshotBase class.
//=================================================================================================
SphSnapshotBase::SphSnapshotBase(SimUnits* _units, string auxfilename): units(_units)
{
  allocated        = false;
  t                = 0.0;
  LastUsed         = time(NULL);
  if (auxfilename != "") filename = auxfilename;
}



//=================================================================================================
//  SphSnapshot::SphSnapshot
/// Constructor for SphSnapshot class.
//=================================================================================================
template <int ndims>
SphSnapshot<ndims>::SphSnapshot
 (string filename,                     ///< Filename containing snapshot data
  SimulationBase* sim):                ///< Simulation object pointer
  SphSnapshotBase(&(sim->simunits), filename),
  simulation(static_cast<Simulation<ndims>* > (sim))
{
  this->ndim = ndims;
  this->fileform = sim->GetParam("out_file_form");

  // Computes how numbers we need to store for each sph/star particle
//  nneededsph = 3*ndims + 5;
//  nneededstar = 3*ndims + 2;
//  nneededbinary = 5;

  if (filename != "") {
    HeaderInfo info;
    info = sim->ReadHeaderSnapshotFile(filename, this->fileform);
    this->t = info.t;
    int Nhydro = info.Nhydro;
    if (Nhydro>0) {
      this->data["sph"]=Species(Nhydro,"sph");
      Species::maptype& sph_values=data["sph"].values;
      sph_values["iorig"]=vector<SNAPFLOAT>();
      sph_values["x"]=vector<SNAPFLOAT>();
      sph_values["vx"]=vector<SNAPFLOAT>();
      sph_values["ax"]=vector<SNAPFLOAT>();
      sph_values["gpot"]=vector<SNAPFLOAT>();
      sph_values["m"]=vector<SNAPFLOAT>();
      sph_values["h"]=vector<SNAPFLOAT>();
      sph_values["rho"]=vector<SNAPFLOAT>();
      sph_values["u"]=vector<SNAPFLOAT>();
      sph_values["dudt"]=vector<SNAPFLOAT>();
      if (ndim>1) {
        sph_values["y"]=vector<SNAPFLOAT>();
        sph_values["vy"]=vector<SNAPFLOAT>();
        sph_values["ay"]=vector<SNAPFLOAT>();
      }
      if (ndim>2) {
        sph_values["z"]=vector<SNAPFLOAT>();
        sph_values["vz"]=vector<SNAPFLOAT>();
        sph_values["az"]=vector<SNAPFLOAT>();
      }
    }
    int Ndust = info.Ndust;
    if (Ndust>0) {
      this->data["dust"]=Species(Ndust,"dust");
      Species::maptype& dust_values=data["dust"].values;
      dust_values["iorig"]=vector<SNAPFLOAT>();
      dust_values["x"]=vector<SNAPFLOAT>();
      dust_values["vx"]=vector<SNAPFLOAT>();
      dust_values["ax"]=vector<SNAPFLOAT>();
      dust_values["gpot"]=vector<SNAPFLOAT>();
      dust_values["m"]=vector<SNAPFLOAT>();
      dust_values["h"]=vector<SNAPFLOAT>();
      dust_values["rho"]=vector<SNAPFLOAT>();
      if (ndim>1) {
        dust_values["y"]=vector<SNAPFLOAT>();
        dust_values["vy"]=vector<SNAPFLOAT>();
        dust_values["ay"]=vector<SNAPFLOAT>();
      }
      if (ndim>2) {
        dust_values["z"]=vector<SNAPFLOAT>();
        dust_values["vz"]=vector<SNAPFLOAT>();
        dust_values["az"]=vector<SNAPFLOAT>();
      }
    }
    int Nstar = info.Nstar;
    if (Nstar>0) {
      this->data["star"]=Species(Nstar,"star");
      Species::maptype& star_values = data["star"].values;
      star_values["x"]=vector<SNAPFLOAT>(Nstar);
      star_values["vx"]=vector<SNAPFLOAT>(Nstar);
      star_values["ax"]=vector<SNAPFLOAT>(Nstar);
      star_values["gpot"]=vector<SNAPFLOAT>(Nstar);
      star_values["m"]=vector<SNAPFLOAT>(Nstar);
      star_values["h"]=vector<SNAPFLOAT>(Nstar);
      if (ndim>1) {
        star_values["y"]=vector<SNAPFLOAT>(Nstar);
        star_values["vy"]=vector<SNAPFLOAT>(Nstar);
        star_values["ay"]=vector<SNAPFLOAT>(Nstar);
      }
      if (ndim>2) {
        star_values["z"]=vector<SNAPFLOAT>(Nstar);
        star_values["vz"]=vector<SNAPFLOAT>(Nstar);
        star_values["az"]=vector<SNAPFLOAT>(Nstar);
      }
    }
  }
}



//=================================================================================================
//  SphSnapshotBase::DeallocateBufferMemory
/// Deallocate memory for current snapshot.
//=================================================================================================
void SphSnapshotBase::DeallocateBufferMemory(void)
{
  debug2("[SphSnapshotBase::DeallocateBufferMemory]");

  for (DataIterator it=data.begin(); it != data.end(); it++) {
    it->second.DeallocateMemory();
  }

  allocated = false;
  return;
}


//=================================================================================================
//  SphSnapshotBase::CalculateMemoryUsage
/// Returns no. of bytes allocated for current snapshot.
//=================================================================================================
int SphSnapshotBase::CalculateMemoryUsage(void)
{
  int result=0;
  for (DataIterator it=data.begin(); it != data.end(); it++) {
    result += it->second.CalculateMemoryUsage();
  }
  return result;
}



//=================================================================================================
//  SphSnapshotBase::CalculatePredictedMemoryUsage
/// Returns no. of bytes that the current snapshot would use, if allocated.
//=================================================================================================
int SphSnapshotBase::CalculatePredictedMemoryUsage(void)
{
  int result=0;
  for (DataIterator it=data.begin(); it != data.end(); it++) {
    result += it->second.CalculatePredictedMemoryUsage();
  }
  return result;
}



//=================================================================================================
//  SphSnapshot::CopyDataFromSimulation
/// Copy particle data from main memory to current snapshot arrays.
//=================================================================================================
template <int ndims>
void SphSnapshot<ndims>::CopyDataFromSimulation()
{
  StarParticle<ndims>* staraux = 0;    ///< ..
  //BinaryOrbit *orbitaux = 0;           ///< ..

  int Nhydro=0;
  int Nstar=0;
  int Ndust=0;

  debug2("[SphSnapshotBase::CopyDataFromSimulation]");

  // Reset the species
  _species.clear();
  data.clear();

  // Read which species are there
  // if (simulation->hydro != NULL && simulation->hydro->GetParticleArray() != NULL) {
  if (simulation->hydro != NULL) {
    // Compute Ndust - note that here Nhydro is Ngas+Ndust
    Nhydro = simulation->hydro->Nhydro;
    for (int n=0; n < Nhydro; n++){
      int ptype = simulation->hydro->GetParticlePointer(n).ptype;
      if (ptype == dust_type) Ndust++;
    }
    // And now Nhydro becomes Ngas
    Nhydro -= Ndust;
    if (Nhydro != 0) {
      _species.push_back("sph");
      data["sph"]=Species(Nhydro,"sph");

      Species& sph = data["sph"];
      Species::maptype& sph_values = data["sph"].values;
      sph_values["iorig"]=vector<SNAPFLOAT>(sph.N);
      sph_values["x"]=vector<SNAPFLOAT>(sph.N);
      sph_values["vx"]=vector<SNAPFLOAT>(sph.N);
      sph_values["ax"]=vector<SNAPFLOAT>(sph.N);
      sph_values["gpot"]=vector<SNAPFLOAT>(sph.N);
      sph_values["m"]=vector<SNAPFLOAT>(sph.N);
      sph_values["h"]=vector<SNAPFLOAT>(sph.N);
      sph_values["rho"]=vector<SNAPFLOAT>(sph.N);
      sph_values["u"]=vector<SNAPFLOAT>(sph.N);
      sph_values["dudt"]=vector<SNAPFLOAT>(sph.N);
      if (ndim>1) {
        sph_values["y"]=vector<SNAPFLOAT>(sph.N);
        sph_values["vy"]=vector<SNAPFLOAT>(sph.N);
        sph_values["ay"]=vector<SNAPFLOAT>(sph.N);
      }
      if (ndim>2) {
        sph_values["z"]=vector<SNAPFLOAT>(sph.N);
        sph_values["vz"]=vector<SNAPFLOAT>(sph.N);
        sph_values["az"]=vector<SNAPFLOAT>(sph.N);
      }
    }
    if (Ndust != 0) {
      _species.push_back("dust");
      data["dust"]=Species(Ndust,"dust");

      Species::maptype& dust_values = data["dust"].values;
      dust_values["iorig"]=vector<SNAPFLOAT>(Ndust);
      dust_values["x"]=vector<SNAPFLOAT>(Ndust);
      dust_values["vx"]=vector<SNAPFLOAT>(Ndust);
      dust_values["ax"]=vector<SNAPFLOAT>(Ndust);
      dust_values["gpot"]=vector<SNAPFLOAT>(Ndust);
      dust_values["m"]=vector<SNAPFLOAT>(Ndust);
      dust_values["h"]=vector<SNAPFLOAT>(Ndust);
      dust_values["rho"]=vector<SNAPFLOAT>(Ndust);
      if (ndim>1) {
        dust_values["y"]=vector<SNAPFLOAT>(Ndust);
        dust_values["vy"]=vector<SNAPFLOAT>(Ndust);
        dust_values["ay"]=vector<SNAPFLOAT>(Ndust);
      }
      if (ndim>2) {
        dust_values["z"]=vector<SNAPFLOAT>(Ndust);
        dust_values["vz"]=vector<SNAPFLOAT>(Ndust);
        dust_values["az"]=vector<SNAPFLOAT>(Ndust);
      }
    }
  }
  if (simulation->nbody != NULL && simulation->nbody->stardata != NULL) {
    staraux  = simulation->nbody->stardata;
//    orbitaux = simulation->nbodytree.orbit;
    Nstar    = simulation->nbody->Nstar;
//    int Norbit   = simulation->nbodytree.Norbit;
    if (Nstar != 0) {
      _species.push_back("star");
      data["star"]=Species(Nstar,"star");

      Species::maptype& star_values = data["star"].values;
      star_values["x"]=vector<SNAPFLOAT>(Nstar);
      star_values["vx"]=vector<SNAPFLOAT>(Nstar);
      star_values["ax"]=vector<SNAPFLOAT>(Nstar);
      star_values["gpot"]=vector<SNAPFLOAT>(Nstar);
      star_values["m"]=vector<SNAPFLOAT>(Nstar);
      star_values["h"]=vector<SNAPFLOAT>(Nstar);
      if (ndim>1) {
        star_values["y"]=vector<SNAPFLOAT>(Nstar);
        star_values["vy"]=vector<SNAPFLOAT>(Nstar);
        star_values["ay"]=vector<SNAPFLOAT>(Nstar);
      }
      if (ndim>2) {
        star_values["z"]=vector<SNAPFLOAT>(Nstar);
        star_values["vz"]=vector<SNAPFLOAT>(Nstar);
        star_values["az"]=vector<SNAPFLOAT>(Nstar);
      }

    }
//    if (Norbit != 0) {
//      _species.push_back("binary");
//      data["binary"]=Species(Norbit,"binary");
//    }
  }



  // Loop over all SPH particles and record particle data
  //-----------------------------------------------------------------------------------------------
  int igas=0; int idust=0;
  for (int i=0; i<simulation->hydro->Nhydro; i++) {
    Particle<ndims>& part = simulation->hydro->GetParticlePointer(i);
    int ptype = part.ptype;

    if (ptype == gas_type) {
      Species::maptype& sph_values = data["sph"].values;
      sph_values["iorig"][igas] = reinterpret_cast<SNAPFLOAT&>( part.iorig);
      sph_values["x"][igas] = (SNAPFLOAT) part.r[0];
      sph_values["vx"][igas] = (SNAPFLOAT) part.v[0];
      sph_values["ax"][igas] = (SNAPFLOAT) part.a[0];
      sph_values["gpot"][igas] = (SNAPFLOAT) part.gpot;
      if (ndims > 1) {
        sph_values["y"][igas] = (SNAPFLOAT) part.r[1];
        sph_values["vy"][igas] = (SNAPFLOAT) part.v[1];
        sph_values["ay"][igas] = (SNAPFLOAT) part.a[1];
      }
      if (ndims > 2) {
        sph_values["z"][igas] = (SNAPFLOAT) part.r[2];
        sph_values["vz"][igas] = (SNAPFLOAT) part.v[2];
        sph_values["az"][igas] = (SNAPFLOAT) part.a[2];
      }
      sph_values["m"][igas]    = (SNAPFLOAT) part.m;
      sph_values["h"][igas]    = (SNAPFLOAT) part.h;
      sph_values["rho"][igas]  = (SNAPFLOAT) part.rho;
      sph_values["u"][igas]    = (SNAPFLOAT) part.u;
      sph_values["dudt"][igas] = (SNAPFLOAT) part.dudt;
      igas++;
    }
    else if (ptype == dust_type) {
      Species::maptype& dust_values = data["dust"].values;
      dust_values["iorig"][idust] = reinterpret_cast<SNAPFLOAT&>( part.iorig);
      dust_values["x"][idust] = (SNAPFLOAT) part.r[0];
      dust_values["vx"][idust] = (SNAPFLOAT) part.v[0];
      dust_values["ax"][idust] = (SNAPFLOAT) part.a[0];
      dust_values["gpot"][idust] = (SNAPFLOAT) part.gpot;
      if (ndims > 1) {
        dust_values["y"][idust] = (SNAPFLOAT) part.r[1];
        dust_values["vy"][idust] = (SNAPFLOAT) part.v[1];
        dust_values["ay"][idust] = (SNAPFLOAT) part.a[1];
      }
      if (ndims > 2) {
        dust_values["z"][idust] = (SNAPFLOAT) part.r[2];
        dust_values["vz"][idust] = (SNAPFLOAT) part.v[2];
        dust_values["az"][idust] = (SNAPFLOAT) part.a[2];
      }
      dust_values["m"][idust]    = (SNAPFLOAT) part.m;
      dust_values["h"][idust]    = (SNAPFLOAT) part.h;
      dust_values["rho"][idust]  = (SNAPFLOAT) part.rho;
      idust++;
    }

  }

  // Loop over star particles and record particle data
  for (int i=0; i<Nstar; i++) {
    Species::maptype& star_values = data["star"].values;
    star_values["x"][i] = (SNAPFLOAT) staraux[i].r[0];
    star_values["vx"][i] = (SNAPFLOAT) staraux[i].v[0];
    star_values["ax"][i] = (SNAPFLOAT) staraux[i].a[0];
    star_values["gpot"][i] = (SNAPFLOAT) staraux[i].gpot;
    if (ndims > 1) {
      star_values["y"][i] = (SNAPFLOAT) staraux[i].r[1];
      star_values["vy"][i] = (SNAPFLOAT) staraux[i].v[1];
      star_values["ay"][i] = (SNAPFLOAT) staraux[i].a[1];
    }
    if (ndims > 2) {
      star_values["z"][i] = (SNAPFLOAT) staraux[i].r[2];
      star_values["vz"][i] = (SNAPFLOAT) staraux[i].v[2];
      star_values["az"][i] = (SNAPFLOAT) staraux[i].a[2];
    }

    star_values["m"][i] = (SNAPFLOAT) staraux[i].m;
    star_values["h"][i] = (SNAPFLOAT) staraux[i].h;
  }

  // Loop over all binary orbits and record data
//  for (int i=0; i<Norbit; i++) {
//    ecc[i]    = (float) orbitaux[i].ecc;
//    mbin[i]   = (float) orbitaux[i].m;
//    period[i] = (float) orbitaux[i].period;
//    qbin[i]   = (float) orbitaux[i].q;
//    sma[i]    = (float) orbitaux[i].sma;
//  }

  LastUsed = time(NULL);

  allocated=true;

  return;
}



//=================================================================================================
//  SphSnapshotBase::GetRealType
/// Convert the given type into the 'real' type, meaning, if we pass default,
/// gives back the true underlying type
//=================================================================================================
string SphSnapshotBase::GetRealType(string type)
{
  // Default is: if there is only one species, then we use that one
  // Otherwise, we return sph particles
  if (type == "default") {
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



//=================================================================================================
///  SphSnapshotBase::GetNparticlesType
///  Get the number of particles for the given type
//=================================================================================================
int SphSnapshotBase::GetNparticlesType(string type)
{
  type = GetRealType(type);
  return data[type].N;
}



//=================================================================================================
//  SphSnapshotBase::ExtractArray
/// Returns pointer to required array stored in snapshot buffer memory.
/// Currently also returns scaling factors for that array.
//=================================================================================================
#ifdef GANDALF_SNAPSHOT_SINGLE_PRECISION
UnitInfo SphSnapshotBase::ExtractArray
 (string name,                         ///< Name of variable to extract
  string type,                         ///< Particle type
  float** out_array,                   ///< Outputted array
  int* size_array,                     ///< No. of elements in outputted array
  float& scaling_factor,               ///< Scaling factor for outputted variable
  string RequestedUnit)                ///< Requested unit for outputted variable
#else
UnitInfo SphSnapshotBase::ExtractArray
 (string name,                         ///< Name of variable to extract
  string type,                         ///< Particle type
  double** out_array,                   ///< Outputted array
  int* size_array,                     ///< No. of elements in outputted array
  double& scaling_factor,               ///< Scaling factor for outputted variable
  string RequestedUnit)                ///< Requested unit for outputted variable
#endif
{
  string unitname;                     // Name of unit
  UnitInfo unitinfo;                   // All data for units and scaling
  SimUnit* unit;                       // Unit object pointer

  // Zero initial array pointers and size
  *out_array = NULL;
  *size_array = 0;
  unit = 0;

  // Check that the memory is allocated. If not, fails very rumorously
  if (!allocated){
    cout << "Error: requested a snapshot that is not allocated!!!!" << endl;
    cout << "This means there's a bug in the memory management: "
      "please inform the authors" << endl;
    exit(-2);
  }

  // Set last time used
  LastUsed = time(NULL);

  // Get the real underlying type, in case we are passing "default"
  type = GetRealType(type);

  // Check type
  if (type != "sph" && type != "star" && type != "binary" && type != "dust") {
    string message = "Error: the type " + type + " was not recognized!";
    ExceptionHandler::getIstance().raise(message);
  }

  *out_array=&(data[type].values[name][0]);

  // If array type and name is valid, pass pointer to array and also set unit
  if (name == "x") {
    unit = &(units->r);
  }
  else if (name == "iorig") {
	  unit = &(units->nounits);
  }
  else if (name == "y") {
    unit = &(units->r);
  }
  else if (name == "z") {
    unit = &(units->r);
  }
  else if (name == "vx") {
    unit = &(units->v);
  }
  else if (name == "vy") {
    unit = &(units->v);
  }
  else if (name == "vz") {
    unit = &(units->v);
  }
  else if (name == "ax") {
    unit = &(units->a);
  }
  else if (name == "ay") {
    unit = &(units->a);
  }
  else if (name == "az") {
    unit = &(units->a);
  }
  else if (name == "m") {
    unit = &(units->m);
  }
  else if (name == "h") {
    unit = &(units->r);
  }
  else if (name == "rho") {
    unit = &(units->rho);
  }
  else if (name == "u") {
    unit = &(units->u);
  }
  else if (name == "dudt") {
    unit = &(units->dudt);
  }
  else if (name == "ecc") {
    unit = &(units->nounits);
  }
  else if (name == "mbin") {
    unit = &(units->m);
  }
  else if (name == "period") {
    unit = &(units->t);
  }
  else if (name == "qbin") {
    unit = &(units->nounits);
  }
  else if (name == "sma") {
    unit = &(units->r);
  }
  else if (name == "gpot") {
    unit = &(units->E);
  }
  else {
    string message = "Warning: the selected array: " + name + " has not been recognized";
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
  *size_array = data[type].N;

  // If no new unit is requested, pass the default scaling values.
  // Otherwise, calculate new scaling factor plus latex label.
  if (RequestedUnit == "default") {
    unitname = unit->outunit;
    RequestedUnit = unitname;
  }
  else {
    unitname = RequestedUnit;
  }
  label = unit->LatexLabel(RequestedUnit);
  scaling_factor = unit->OutputScale(RequestedUnit);

  unitinfo.name = unitname;
  unitinfo.label = label;

  return unitinfo;
}



//=================================================================================================
//  SphSnapshot::ReadSnapshot
/// Read snapshot into main memory and then copy into snapshot buffer.
//=================================================================================================
template <int ndims>
void SphSnapshot<ndims>::ReadSnapshot(string format)
{
  debug2("[SphSnapshotBase::ReadSnapshot]");

  // Set pointer to units object
  units = &(simulation->simunits);

  // Read simulation into main memory
  simulation->ReadSnapshotFile(filename, format);

  // Recalculate input units if required
  units->SetupUnits(simulation->simparams);

  // Scale particle data to dimensionless code units
  simulation->ConvertToCodeUnits();

  // Now copy from main memory to current snapshot
  CopyDataFromSimulation();

  // Record simulation snapshot time
  t = simulation->t;

  return;
}
