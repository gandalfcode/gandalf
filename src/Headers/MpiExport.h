//=============================================================================
//  MpiExport.h
//  Header file containing definitions of structures needed to export particles
//  and tree cells
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

#ifndef MPIEXPORT_H_
#define MPIEXPORT_H_

#include <cstring>
#include <vector>
#include "Precision.h"



//=============================================================================
//  Struct ExportParticleInfo
/// Information needed to export a particle that needs only gravity calculation
//=============================================================================
template <int ndim>
struct ExportParticleInfo {
  FLOAT r[ndim];

  static MPI_Datatype CreateMpiDataType () {
    MPI_Datatype export_particle_type;
    MPI_Datatype types[1] = {MPI_BYTE};
    MPI_Aint offsets[1] = {0};
    int blocklen[1] = {sizeof(ExportParticleInfo<ndim>)};

    MPI_Type_create_struct(1,blocklen,offsets,types,&export_particle_type);

    return export_particle_type;
  }

};

//=============================================================================
//  Struct ExportBackParticleInfo
/// Information needed to export back from a particle that needs only gravity calculation
//=============================================================================
template <int ndim>
struct ExportBackParticleInfo {
  FLOAT a[ndim];
  FLOAT gpot;

  static MPI_Datatype CreateMpiDataType () {
    MPI_Datatype export_back_particle_type;
    MPI_Datatype types[1] = {MPI_BYTE};
    MPI_Aint offsets[1] = {0};
    int blocklen[1] = {sizeof(ExportBackParticleInfo<ndim>)};

    MPI_Type_create_struct(1,blocklen,offsets,types,&export_back_particle_type);

    return export_back_particle_type;
  }


};


//=============================================================================
//  Struct ExportCellInfo
/// Information needed to export a cell to another processor
//=============================================================================
template <int ndim>
struct ExportCellInfo {

};


//=============================================================================
//  Struct BytesChunk
/// Generic information needed to identify a chunk of bytes: pointer to location
/// in memory and size of the chunk
//=============================================================================
struct BytesChunk {
  void* data;
  int size;
};

template <class T>
inline void copy (char* pointer, T* element) {
  void*  element_unsafe = reinterpret_cast<void*> (element);
  void* pointer_unsafe = reinterpret_cast<void*> (pointer);
  memcpy(pointer_unsafe,element_unsafe,sizeof(T));
}

template <class T>
inline void copy (T* element, char* pointer) {
  void* element_unsafe = reinterpret_cast<void*> (element);
  void* pointer_unsafe = reinterpret_cast<void*> (pointer);
  memcpy(element_unsafe,pointer_unsafe,sizeof(T));
}



#endif /* MPIEXPORT_H_ */
