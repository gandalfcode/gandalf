//=================================================================================================
//  SphSnapshot.h
//  Contains definitions for SphSnapshot class
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


#ifndef _SPH_SNAPSHOT_H_
#define _SPH_SNAPSHOT_H_


#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "Precision.h"
#include "Sph.h"
#include "Simulation.h"
#include "UnitInfo.h"
#include "BinaryOrbit.h"
using namespace std;

class Species {
public:
  typedef map<string,vector<SNAPFLOAT> > maptype;
  map<string,vector<SNAPFLOAT> > values;
  int N;
  string name;

  Species(): N(0) {};

  Species(int _N, string _name): N(_N), name(_name) {};

  void DeallocateMemory() {
    for (maptype::iterator it=values.begin(); it != values.end(); it++) {
      it->second.clear();
    }
  }

  bool IsAllocated() {
    bool result=false;
    if (values.size() > 0) {
      maptype::iterator it=values.begin();
      if (it->second.size()>0)
        result=true;
    }
    return result;
  }

  int CalculateMemoryUsage() {
    if (IsAllocated() )
      return CalculatePredictedMemoryUsage();
    else
      return 0;
  }

  int CalculatePredictedMemoryUsage() {
    return N*sizeof(SNAPFLOAT)*values.size();
  }
};



//=================================================================================================
/// Class SphSnapshotBase
/// \brief   Definition for Sph Snapshot (Base) for recording simulation data.
/// \details Snapshot class contains a copy of the simulation data,
///          either from main memory or from a loaded file.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
class SphSnapshotBase
{
 private:
//  void AllocateBufferMemoryBinary();
//  void AllocateBufferMemorySph();
//  void AllocateBufferMemoryStar();
//  void DeallocateBufferMemoryBinary();
//  void DeallocateBufferMemorySph();
//  void DeallocateBufferMemoryStar();

 protected:
  //int nneededbinary;                   ///< No. of variables needed to store binary orbit
  //int nneededsph;                      ///< No. of variables needed to store for sph ptcl
  //int nneededstar;                     ///< No. of variables needed to store for star ptcl
  vector<string> _species;             ///< ..
  typedef map<string, Species> MapData;
  typedef MapData::iterator DataIterator;
  map<string, Species> data;

 public:

  static SphSnapshotBase* SphSnapshotFactory(string filename,
                                             SimulationBase* sim, int ndim);

  SphSnapshotBase(SimUnits*, string="");
  //virtual ~SphSnapshotBase();


  // Snapshot function prototypes
  //-----------------------------------------------------------------------------------------------
  //void AllocateBufferMemory(void);
  void DeallocateBufferMemory(void);
  int CalculateMemoryUsage(void);
  int CalculatePredictedMemoryUsage(void);
  virtual void CopyDataFromSimulation()=0;
#if defined(GANDALF_SNAPSHOT_SINGLE_PRECISION)
  UnitInfo ExtractArray(string, string, float** out_array, int* size_array,
                        float& scaling_factor, string RequestedUnit);
#else
  UnitInfo ExtractArray(string, string, double** out_array, int* size_array,
                        double& scaling_factor, string RequestedUnit);
#endif
  virtual void ReadSnapshot(string)=0;
  int GetNTypes() {return _species.size(); };
  string GetSpecies(int ispecies) { return _species.at(ispecies); };
  string GetRealType(string);
  int GetNparticlesType(string species);

  // All variables
  //-----------------------------------------------------------------------------------------------
  bool allocated;                   ///< Is snapshot memory allocated?
  //bool allocatedbinary;             ///< Is SPH particle memory allocated?
  //bool allocatedsph;                ///< Is SPH particle memory allocated?
  //bool allocatedstar;               ///< Is star particle memory allocated?
  //bool computedbinary;              ///< Are binary properties computed?
  //bool computedsph;                 ///< Are additional SPH values computed?
  //bool computednbody;               ///< Are additional star values computed?
  int LastUsed;                     ///< ??
  //int nallocatedbinary;             ///< No. of floats allocated for SPH
  //int nallocatedsph;                ///< No. of floats allocated for SPH
  //int nallocatedstar;               ///< No. of floats allocated for stars
  int ndim;                         ///< Local copy of ndim
  //int Nbinary;                      ///< No. of binary stars
  //int Norbit;                       ///< No. of orbits in memory
  //int Norbitmax;                    ///< Max. no. of orbits
  //int Nquadruple;                   ///< No. of quadruple systems
  //int Nhydro;                         ///< No. of SPH particles
  //int Nhydromax;                      ///< Max. no. of SPH particles
  //int Nstar;                        ///< No. of star particles
  //int Nstarmax;                     ///< Max. no. of star particles
  //int Ntriple;                      ///< No. of triple systems
  SNAPFLOAT t;                         ///< Simulation time of snapshot

  string filename;                  ///< Filename of snapshot
  string fileform;                  ///< File format of snapshot
  //string unitname;                  ///< Aux. unit string
  string label;                     ///< Aux. latex label

  SimUnits* units;                  ///< Pointer to units object


};



//=================================================================================================
/// Class SphSnapshot
/// \brief   Definition for Sph Snapshot for recording simulation data.
/// \details Snapshot class contains a copy of the simulation data,
///          either from main memory or from a loaded file.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndims>
class SphSnapshot : public SphSnapshotBase {
public:
  SphSnapshot (string, SimulationBase* );
  //~SphSnapshot() {};
  void CopyDataFromSimulation();
  void ReadSnapshot(string);

  Simulation<ndims>* simulation;
};
#endif
