//=================================================================================================
//  SimulationIO.hpp
//  Contains all functions for reading and writing simulation data to/from snapshot files.
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


#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstring>
#include "Simulation.h"
#include "Parameters.h"
#include "Debug.h"
#include "HeaderInfo.h"
#include "formatted_output.h"
#include "BinaryIO.h"
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif
using namespace std;


static const string binary_tag("SERENBINARYDUMPV3");
static const string ascii_tag("SERENASCIIDUMPV2");
static const int string_length = 20;



//=================================================================================================
//  SimulationBase::ReadSnapshotFile
/// Read snapshot file with specified filename and format.
//=================================================================================================
bool SimulationBase::ReadSnapshotFile
 (string filename,                     ///< [in] Name of input snapshot file
  string fileform)                     ///< [in] Format of input snapshot file
{
  debug2("[Simulation::ReadSnapshotFile]");

  cout << "Reading snapshot : " << filename << "   format : " << fileform << endl;

  // Read in snapshot file to main memory
  if (fileform == "column") {
    return ReadColumnSnapshotFile(filename);
  }
  else if (fileform == "sf" || fileform == "seren_form") {
    return ReadSerenFormSnapshotFile(filename);
  }
  else if (fileform == "su" || fileform == "seren_unform") {
    return ReadSerenUnformSnapshotFile(filename);
  }
  else {
    cout << "Unrecognised file format" << endl;
    return false;
  }

  // Set flag to ensure variables are converted to code units
  rescale_particle_data = true;

}



//=================================================================================================
//  SimulationBase::WriteSnapshotFile
/// Write snapshot file with specified filename and format.
//=================================================================================================
bool SimulationBase::WriteSnapshotFile
 (string filename,                     ///< [in] Name of output snapshot file
  string fileform)                     ///< [in] Format of output snapshot file
{
  debug2("[Simulation::WriteSnapshotFile]");

  if (fileform == "column") {
    return WriteColumnSnapshotFile(filename);
  }
  else if (fileform == "sf" || fileform == "seren_form") {
    return WriteSerenFormSnapshotFile(filename);
  }
  else if (fileform == "su" || fileform == "seren_unform") {
    return WriteSerenUnformSnapshotFile(filename);
  }
  else if (fileform == "slite" || fileform == "seren_lite") {
    return WriteSerenLiteSnapshotFile(filename);
  }
  else {
    cout << "Unrecognised file format" << endl;
    return false;
  }
}



//=================================================================================================
//  SimulationBase::ReadHeaderSnapshotFile
/// Read the header of a snapshot file, given the filename and the format.  Return information
/// in a HeaderInfo struct. This is designed to be called from Python.
//=================================================================================================
HeaderInfo SimulationBase::ReadHeaderSnapshotFile
 (string filename,                     ///< [in] Name of output snapshot file
  string fileform)                     ///< [in] Format of output snapshot file
{
  HeaderInfo info;                     // Important information contained in header file
  ifstream infile;                     // File stream object

  debug2("[Simulation::ReadHeaderSnapshotFile]");

  infile.open(filename.c_str());

  if (fileform == "column") {
    ReadColumnHeaderFile(infile, info);
  }
  else if (fileform == "sf" || fileform == "seren_form") {
    ReadSerenFormHeaderFile(infile, info);
  }
  else if (fileform == "su" || fileform == "seren_unform") {
    ReadSerenUnformHeaderFile(infile, info);
  }
  else {
    ExceptionHandler::getIstance().raise("Unrecognised file format");
  }

  return info;
}



//=================================================================================================
//  Simulation::ReadColumnHeaderFile
/// Function for reading the header file of a snapshot. Does not modify the variables of the
/// Simulation class, but rather returns information in a HeaderInfo struct.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ReadColumnHeaderFile
 (ifstream& infile,                    ///< [in] Reference to input file stream object
  HeaderInfo& info)                    ///< [out] Header data structure
{
  debug2("[Simulation::ReadColumnHeaderFile]");

  // Open file and read header information
  infile >> info.Nhydro;
  infile >> info.Nstar;
  infile >> info.ndim;
  infile >> info.t;
  info.t /= simunits.t.inscale;

  // Check dimensionality matches if using fixed dimensions
  if (info.ndim != ndim) {
    std::ostringstream stream;
    stream << "Incorrect no. of dimensions in file : "
     << info.ndim << "  [ndim : " << ndim << "]" << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }

  return;
}



//=================================================================================================
//  Simulation::ReadColumnSnapshotFile
/// Reads a column format data snapshot of given filename.
//=================================================================================================
template <int ndim>
bool Simulation<ndim>::ReadColumnSnapshotFile
 (string filename)                     ///< Filename of column data snapshot file
{
  int i;                               // Particle counter
  FLOAT raux;                          // ..
  ifstream infile;                     // ..
  HeaderInfo info;                     // ..

  debug2("[Simulation::ReadColumnSnapshotFile]");

  infile.open(filename.c_str());

  ReadColumnHeaderFile(infile, info);
  t = info.t;
  tsnaplast = t;

  hydro->Nhydro = info.Nhydro;
  AllocateParticleMemory();
  i = 0;

  // Read in data depending on dimensionality
  //-----------------------------------------------------------------------------------------------
  while (infile.good() && i < hydro->Nhydro) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    if (ndim == 1) {
      infile >> part.r[0] >> part.v[0] >> part.m >> part.h >> part.rho >> part.u;
    }
    else if (ndim == 2) {
      infile >> part.r[0] >> part.r[1] >> part.v[0] >> part.v[1]
             >> part.m >> part.h >> part.rho >> part.u;
    }
    else if (ndim == 3) {
      infile >> part.r[0] >> part.r[1] >> part.r[2] >> part.v[0] >> part.v[1] >> part.v[2]
             >> part.m >> part.h >> part.rho >> part.u;
    }
    i++;
  }


  nbody->Nstar = info.Nstar;
  sinks->Nsink = info.Nstar;
  nbody->AllocateMemory(nbody->Nstar);
  sinks->AllocateMemory(nbody->Nstar);
  i = 0;

  // Read in data depending on dimensionality
  //-----------------------------------------------------------------------------------------------
  while (infile.good() && i < nbody->Nstar) {
    if (ndim == 1) {
      infile >> nbody->stardata[i].r[0] >> nbody->stardata[i].v[0] >> nbody->stardata[i].m
             >> nbody->stardata[i].h >> raux >> raux;
    }
    else if (ndim == 2) {
      infile >> nbody->stardata[i].r[0] >> nbody->stardata[i].r[1] >> nbody->stardata[i].v[0]
             >> nbody->stardata[i].v[1] >> nbody->stardata[i].m >> nbody->stardata[i].h
             >> raux >> raux;
    }
    else if (ndim == 3) {
      infile >> nbody->stardata[i].r[0] >> nbody->stardata[i].r[1] >> nbody->stardata[i].r[2]
             >> nbody->stardata[i].v[0] >> nbody->stardata[i].v[1] >> nbody->stardata[i].v[2]
             >> nbody->stardata[i].m >> nbody->stardata[i].h >> raux >> raux;
    }
    sinks->sink[i].radius = nbody->kernp->kernrange*nbody->stardata[i].h;
    sinks->sink[i].star = &(nbody->stardata[i]);
    i++;
  }


  infile.close();

  return true;
}



#ifdef MPI_PARALLEL
//=================================================================================================
//  Simulation::WriteColumnSnapshotFile
/// Write hydro and N-body particle data to column data snapshot file (MPI parallelised version)
//=================================================================================================
template <int ndim>
bool Simulation<ndim>::WriteColumnSnapshotFile(string filename)
{
  int i;
  ostringstream outfile;

  debug2("[Simulation::WriteColumnSnapshotFileMPI]");

  if (rank == 0) cout << "Writing current data to snapshot file : " << filename << endl;

  // Open file
  MPI_File file;
  char* filename_str = new char[strlen(filename.c_str())+1];
  strcpy(filename_str,filename.c_str());
  MPI_File_open(MPI_COMM_WORLD, filename_str, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
  MPI_File_set_size(file,0);
  delete[] filename_str;

  // Collect total number of particles
  int Ntotsph;
  MPI_Allreduce(&hydro->Nhydro, &Ntotsph, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  int Ntotstar;
  MPI_Allreduce(&nbody->Nstar, &Ntotstar, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Root node writes header
  if (rank == 0) {
    outfile << Ntotsph << endl;
    outfile << Ntotstar << endl;
    outfile << ndim << endl;
    outfile << t*simunits.t.outscale << endl;
  }

  // Write data for hydro particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    if (ndim == 1)
      outfile << part.r[0]*simunits.r.outscale << "   "
              << part.v[0]*simunits.v.outscale << "   "
              << part.m*simunits.m.outscale << "   "
              << part.h*simunits.r.outscale << "   "
              << part.rho*simunits.rho.outscale << "   "
              << part.u*simunits.u.outscale
              << endl;
    else if (ndim == 2)
      outfile << part.r[0]*simunits.r.outscale << "   "
              << part.r[1]*simunits.r.outscale << "   "
              << part.v[0]*simunits.v.outscale << "   "
              << part.v[1]*simunits.v.outscale << "   "
              << part.m*simunits.m.outscale << "   "
              << part.h*simunits.r.outscale << "   "
              << part.rho*simunits.rho.outscale << "   "
              << part.u*simunits.u.outscale
              << endl;
    else if (ndim == 3)
      outfile << part.r[0]*simunits.r.outscale << "   "
              << part.r[1]*simunits.r.outscale << "   "
              << part.r[2]*simunits.r.outscale << "   "
              << part.v[0]*simunits.v.outscale << "   "
              << part.v[1]*simunits.v.outscale << "   "
              << part.v[2]*simunits.v.outscale << "   "
              << part.m*simunits.m.outscale << "   "
              << part.h*simunits.r.outscale << "   "
              << part.rho*simunits.rho.outscale << "   "
              << part.u*simunits.u.outscale
              << endl;
  }

  //Now all nodes write to the file their portion
  //To do that, we need to know the offset of each node, summin up the length of each bit
  {
    std::string content = outfile.str();
    int offset;
    int length_char = content.length();
    MPI_Exscan(&length_char, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (rank==0) {
      offset = 0;
    }
    MPI_Offset offset_mpi = offset;
    MPI_File_seek(file,offset_mpi,MPI_SEEK_SET);
    //Now we can do the actual writing. Extract information from the string and pass it to MPI
    char* data = new char[strlen(content.c_str())+1];
    strcpy(data, content.c_str());
    MPI_Status status;
    MPI_File_write(file, data, length_char, MPI_CHAR, &status);
    delete[] data;
  }

  //Now clear the stream object
  outfile.clear();
  outfile.str("");

  //We need to know where the last process got to seek to that point
  MPI_Offset end_sph_mpi; int end_sph;
  if (rank==Nmpi-1) {
    MPI_File_get_position(file, &end_sph_mpi);
    end_sph = end_sph_mpi;
  }
  MPI_Bcast(&end_sph, 1, MPI_INT, Nmpi-1, MPI_COMM_WORLD);


  // Write data for Nbody particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<nbody->Nstar; i++) {
    if (ndim == 1)
      outfile << nbody->stardata[i].r[0]*simunits.r.outscale << "   "
              << nbody->stardata[i].v[0]*simunits.v.outscale << "   "
              << nbody->stardata[i].m*simunits.m.outscale << "   "
              << nbody->stardata[i].h*simunits.r.outscale << "   "
              << 0.0 << "   "
              << 0.0
              << endl;
    else if (ndim == 2)
      outfile << nbody->stardata[i].r[0]*simunits.r.outscale << "   "
              << nbody->stardata[i].r[1]*simunits.r.outscale << "   "
              << nbody->stardata[i].v[0]*simunits.v.outscale << "   "
              << nbody->stardata[i].v[1]*simunits.v.outscale << "   "
              << nbody->stardata[i].m*simunits.m.outscale << "   "
              << nbody->stardata[i].h*simunits.r.outscale << "   "
              << 0.0 << "   "
              << 0.0
              << endl;
    else if (ndim == 3)
      outfile << nbody->stardata[i].r[0]*simunits.r.outscale << "   "
              << nbody->stardata[i].r[1]*simunits.r.outscale << "   "
              << nbody->stardata[i].r[2]*simunits.r.outscale << "   "
              << nbody->stardata[i].v[0]*simunits.v.outscale << "   "
              << nbody->stardata[i].v[1]*simunits.v.outscale << "   "
              << nbody->stardata[i].v[2]*simunits.v.outscale << "   "
              << nbody->stardata[i].m*simunits.m.outscale << "   "
              << nbody->stardata[i].h*simunits.r.outscale << "   "
              << 0.0 << "   "
              << 0.0
              << endl;
  }


  // Now all nodes write to the file their portion
  // To do that, we need to know the offset of each node, summing up the length of each bit
  {
    std::string content = outfile.str();
    int offset;
    int length_char = content.length();
    MPI_Exscan(&length_char, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == 0) offset = 0;
    MPI_Offset offset_mpi = offset;

    //Offset the position by the end of the hydro information
    offset_mpi += end_sph_mpi;
    MPI_File_seek(file, offset_mpi, MPI_SEEK_SET);

    // Now we can do the actual writing. Extract information from the string and pass it to MPI
    char* data = new char[strlen(content.c_str())+1];
    strcpy(data, content.c_str());
    delete[] data;
  }

  MPI_File_close(&file);

  return true;
}



#else
//=================================================================================================
//  Simulation::WriteColumnSnapshotFile
/// Write hydro and N-body particle data to column data snapshot file.
//=================================================================================================
template <int ndim>
bool Simulation<ndim>::WriteColumnSnapshotFile(string filename)
{
  int i;
  ofstream outfile;

  debug2("[Simulation::WriteColumnSnapshotFile]");

  cout << "Writing current data to snapshot file : " << filename << endl;

  // Open file and read header information
  outfile.open(filename.c_str());
  outfile << hydro->Nhydro << endl;
  outfile << nbody->Nstar << endl;
  outfile << ndim << endl;
  outfile << t*simunits.t.outscale << endl;

  // Write data for hydro particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    if (ndim == 1) {
      outfile << part.r[0]*simunits.r.outscale << "   "
              << part.v[0]*simunits.v.outscale << "   "
              << part.m*simunits.m.outscale << "   "
              << part.h*simunits.r.outscale << "   "
              << part.rho*simunits.rho.outscale << "   "
              << part.u*simunits.u.outscale << "   "
              << endl;
    }
    else if (ndim == 2) {
      outfile << part.r[0]*simunits.r.outscale << "   "
              << part.r[1]*simunits.r.outscale << "   "
              << part.v[0]*simunits.v.outscale << "   "
              << part.v[1]*simunits.v.outscale << "   "
              << part.m*simunits.m.outscale << "   "
              << part.h*simunits.r.outscale << "   "
              << part.rho*simunits.rho.outscale << "   "
              << part.u*simunits.u.outscale
              << endl;
    }
    else if (ndim == 3) {
      outfile << part.r[0]*simunits.r.outscale << "   "
              << part.r[1]*simunits.r.outscale << "   "
              << part.r[2]*simunits.r.outscale << "   "
              << part.v[0]*simunits.v.outscale << "   "
              << part.v[1]*simunits.v.outscale << "   "
              << part.v[2]*simunits.v.outscale << "   "
              << part.m*simunits.m.outscale << "   "
              << part.h*simunits.r.outscale << "   "
              << part.rho*simunits.rho.outscale << "   "
              << part.u*simunits.u.outscale
              << endl;
    }
  }

  // Write data for N-body particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<nbody->Nstar; i++) {
    if (ndim == 1) {
      outfile << nbody->stardata[i].r[0]*simunits.r.outscale << "   "
              << nbody->stardata[i].v[0]*simunits.v.outscale << "   "
              << nbody->stardata[i].m*simunits.m.outscale << "   "
              << nbody->stardata[i].h*simunits.r.outscale << "   "
              << 0.0 << "   "
              << 0.0
              << endl;
    }
    else if (ndim == 2) {
      outfile << nbody->stardata[i].r[0]*simunits.r.outscale << "   "
              << nbody->stardata[i].r[1]*simunits.r.outscale << "   "
              << nbody->stardata[i].v[0]*simunits.v.outscale << "   "
              << nbody->stardata[i].v[1]*simunits.v.outscale << "   "
              << nbody->stardata[i].m*simunits.m.outscale << "   "
              << nbody->stardata[i].h*simunits.r.outscale << "   "
              << 0.0 << "   "
              << 0.0
              << endl;
    }
    else if (ndim == 3) {
      outfile << nbody->stardata[i].r[0]*simunits.r.outscale << "   "
              << nbody->stardata[i].r[1]*simunits.r.outscale << "   "
              << nbody->stardata[i].r[2]*simunits.r.outscale << "   "
              << nbody->stardata[i].v[0]*simunits.v.outscale << "   "
              << nbody->stardata[i].v[1]*simunits.v.outscale << "   "
              << nbody->stardata[i].v[2]*simunits.v.outscale << "   "
              << nbody->stardata[i].m*simunits.m.outscale << "   "
              << nbody->stardata[i].h*simunits.r.outscale << "   "
              << 0.0 << "   "
              << 0.0
              << endl;
    }
  }

  outfile.close();

  return true;
}
#endif



//=================================================================================================
//  Simulation::ReadSerenFormHeaderFile
/// Function for reading the header file of a snapshot. Does not modify the variables of the
/// Simulation class, but rather returns information in a HeaderInfo struct.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ReadSerenFormHeaderFile
 (ifstream& infile,                    ///< [in] Input file stream
  HeaderInfo& info)                    ///< [out] Header info data structure
{
  int i;                               // Aux. counter
  int idata[50];                       // Integer data array
  int ilpdata[50];                     // Long integer data array
  FLOAT rdata[50];                     // Real data array
  DOUBLE ddata[50];                    // Double precision data array
  string dummy;                        // Dummy string variable

  // Read information identifying format and precision of file.
  // Then check if each value corresponds to the current values.
  //-----------------------------------------------------------------------------------------------
  infile >> dummy;
  infile >> dummy;
  infile >> info.ndim;
  infile >> dummy;
  infile >> dummy;

  // Read infile header integer data
  for (i=0; i<50; i++) infile >> idata[i];
  for (i=0; i<50; i++) infile >> ilpdata[i];
  for (i=0; i<50; i++) infile >> rdata[i];
  for (i=0; i<50; i++) infile >> ddata[i];

  info.Nhydro = idata[0];
  info.Nstar  = idata[1];
  info.t      = ddata[0];
  info.t      /= simunits.t.inscale;

  // Check dimensionality matches if using fixed dimensions
  if (info.ndim != ndim) {
    std::ostringstream stream;
    stream << "No. of dimensions in file : "
     << info.ndim << "  [ndim : " << ndim << "]" << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }

  return;
}



//=================================================================================================
//  Simulation::ReadSerenFormSnapshotFile
/// Read data from Seren format snapshot file into main arrays.
//=================================================================================================
template <int ndim>
bool Simulation<ndim>::ReadSerenFormSnapshotFile(string filename)
{
  int dim_check;               // Dimension check
  //int dmdt_range_aux;          // Accretion history array size
  int i;                       // Aux. counter
  int ifirst;                  // i.d. of first particle
  int ilast;                   // i.d. of last particle
  int j;                       // Aux. counter
  int k;                       // Dimension counter
  int ndata;                   // No. of data arrays written
  int nunit;                   // No. of unit strings
  int pr_check;                // Precision check
  int sink_data_length;        // Length of float array for sink data
  int idata[50];               // Integer data array
  int ilpdata[50];             // Long data array
  int typedata[50][5];         // Hydro Particle data array information
  FLOAT rdata[50];             // Real data array
  FLOAT rtemp;                 // Dummy float variable
  DOUBLE ddata[50];            // Double float data array
  string format_id;            // File format (for verification)
  string data_id[50];          // String ids of arrays written
  string unit_data[50];        // String ids of units written
  string dummystring;          // Dummy string variable
  ifstream infile;             // Stream of input file

  debug2("[Simulation::ReadSerenFormSnapshotFile]");

  infile.open(filename.c_str());


  // Read information identifying format and precision of file.
  // Then check if each value corresponds to the current values.
  //----------------------------------------------------------------------------------------------
  infile >> format_id;
  simparams->TrimWhiteSpace(format_id);
  if (format_id != ascii_tag && format_id != "SERENASCIIDUMPV3") {
    std::ostringstream stream;
    stream << "Incorrect format of snapshot file " << filename << ": " << format_id << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }
  infile >> pr_check;
  infile >> dim_check;
  if (dim_check != ndim) {
    std::ostringstream stream;
    stream << "Incorrect NDIM in file " << filename << " : got "
           << dim_check << ", expected " << ndim << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }
  infile >> dim_check;
    if (dim_check != ndim) {
      std::ostringstream stream;
      stream << "Incorrect NDIM in file " << filename << " : got "
             << dim_check << ", expected " << ndim << endl;
      ExceptionHandler::getIstance().raise(stream.str());
  }
  infile >> dim_check;
  if (dim_check != ndim) {
    std::ostringstream stream;
    stream << "Incorrect NDIM in file " << filename << " : got "
           << dim_check << ", expected " << ndim << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }

  // Read infile header integer data
  for (i=0; i<50; i++) infile >> idata[i];
  for (i=0; i<50; i++) infile >> ilpdata[i];
  for (i=0; i<50; i++) infile >> rdata[i];
  for (i=0; i<50; i++) infile >> ddata[i];

  hydro->Nhydro      = idata[0];
  nbody->Nstar   = idata[1];
  sinks->Nsink    = idata[1];
  //dmdt_range_aux = idata[29];
  nunit          = idata[19];
  ndata          = idata[20];

  // Variables that should be remembered for restarts
  if (restart) {
    Noutsnap      = ilpdata[0];
    Nsteps        = ilpdata[1];
    Noutlitesnap  = ilpdata[10];
    t             = ddata[0];
    tsnaplast     = ddata[1];
    hydro->mmean  = ddata[2];
    tlitesnaplast = ddata[10];
  }

  //hydro->h_fac = rdata[0];
  //hydroptr->tlastsnap = ddata[1];
  //hydro->mgas_orig = ddata[2];

  // Read in unit data
  if (nunit > 0) {
    for (i=0; i<nunit; i++) infile >> unit_data[i];
  }

  // Read in unit variables here
  simunits.r.inunit      = unit_data[0];
  simunits.m.inunit      = unit_data[1];
  simunits.t.inunit      = unit_data[2];
  simunits.v.inunit      = unit_data[3];
  simunits.a.inunit      = unit_data[4];
  simunits.rho.inunit    = unit_data[5];
  simunits.sigma.inunit  = unit_data[6];
  simunits.press.inunit  = unit_data[7];
  simunits.f.inunit      = unit_data[8];
  simunits.E.inunit      = unit_data[9];
  simunits.mom.inunit    = unit_data[10];
  simunits.angmom.inunit = unit_data[11];
  simunits.angvel.inunit = unit_data[12];
  simunits.dmdt.inunit   = unit_data[13];
  simunits.L.inunit      = unit_data[14];
  simunits.kappa.inunit  = unit_data[15];
  simunits.B.inunit      = unit_data[16];
  simunits.Q.inunit      = unit_data[17];
  simunits.Jcur.inunit   = unit_data[18];
  simunits.u.inunit      = unit_data[19];
  simunits.temp.inunit   = unit_data[20];
  if (nunit > 21) simunits.dudt.inunit = unit_data[21];

  // Read infile ids of arrays contained in file
  if (ndata > 0) {
    for (i=0; i<ndata; i++) infile >> data_id[i];
  }

  // Read infile ids of arrays contained in file
  if (ndata > 0) {
    for (i=0; i<ndata; i++) {
      infile >> typedata[i][0] >> typedata[i][1]
             >> typedata[i][2] >> typedata[i][3] >> typedata[i][4];
    }
  }

  // Allocate memory now we know particle numbers
  AllocateParticleMemory();


  // Loop through array ids and read each array in turn
  //===============================================================================================
  for (j=0; j<ndata; j++) {

    ifirst = typedata[j][1];
    ilast = typedata[j][2];

    // porig
    //---------------------------------------------------------------------------------------------
    if (data_id[j] == "porig") {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        infile >> part.iorig;
      }
    }

    // Positions
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "r") {
      if (ndim == 1) {
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          infile >> part.r[0];
        }
      }
      else if (ndim == 2) {
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          infile >> part.r[0] >> part.r[1];
        }
      }
      else if (ndim == 3) {
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          infile >> part.r[0] >> part.r[1] >> part.r[2];
        }
      }
    }

    // Masses
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "m") {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        infile >> part.m;
        assert(part.m > (FLOAT) 0.0);
      }
    }

    // Smoothing lengths
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "h") {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        infile >> part.h;
      }
      initial_h_provided = true;
    }

    // Velocities
    //--------------------------------------------------------------------------------------------
    else if (data_id[j] == "v") {
      if (ndim == 1) {
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          infile >> part.v[0];
        }
      }
      else if (ndim == 2) {
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          infile >> part.v[0] >> part.v[1];
        }
      }
      else if (ndim == 3) {
        for (i=0; i<hydro->Nhydro; i++) {
          Particle<ndim>& part = hydro->GetParticlePointer(i);
          infile >> part.v[0] >> part.v[1] >> part.v[2];
        }
      }
    }

    // Other 1-D redundant information
    //---------------------------------------------------------------------------------------------
    else if(data_id[j] == "temp") {
      for (i=0; i<hydro->Nhydro; i++) infile >> rtemp;
    }

    // Densities
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "rho") {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        infile >> part.rho;
      }
    }

    // Specific internal energies
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "u") {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        infile >> part.u;
      }
    }

    // Sinks/stars
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "sink_v1") {
      sink_data_length = 12 + 2*ndim;  //+ 2*dmdt_range_aux;
      int ii;
      FLOAT sdata[sink_data_length];
      for (ii=0; ii<6; ii++) infile >> idata[ii];
      if (nbody->Nstar > 0) {
        for (i=0; i<nbody->Nstar; i++) {
          for (ii=0; ii<2; ii++) {
            char boolean;
            infile >> boolean;
          }
          for (ii=0; ii<2; ii++) infile >> idata[ii];
          for (ii=0; ii<sink_data_length; ii++) infile >> sdata[ii];
          for (k=0; k<ndim; k++) nbody->stardata[i].r[k] = sdata[k+1];
          for (k=0; k<ndim; k++) nbody->stardata[i].v[k] = sdata[k+1+ndim];
          nbody->stardata[i].m = sdata[1+2*ndim];
          nbody->stardata[i].h = sdata[2+2*ndim];
          nbody->stardata[i].radius = sdata[3+2*ndim];
          nbody->nbodydata[i] = &(nbody->stardata[i]);
          sinks->sink[i].star = &(nbody->stardata[i]);
          sinks->sink[i].radius = sdata[3+2*ndim];
        }
      }
    }

    // Skip through arbitrary 1D or 2D array
    //---------------------------------------------------------------------------------------------
    else if (typedata[j][0] >= 1) {
      int kk = typedata[j][0];
      for (i=ifirst-1; i<ilast; i++) {
        for (k=0; k<kk; k++) infile >> dummystring;
      }
    }

  }
  //===============================================================================================

  // Close file
  infile.close();

  return true;
}



//=================================================================================================
//  Simulation::WriteSerenFormSnapshotFile
/// Write hydro and N-body particle data to snapshot file in Seren format.
//=================================================================================================
template <int ndim>
bool Simulation<ndim>::WriteSerenFormSnapshotFile
 (string filename)                     ///< [in] Filename to write new snapshot to
{
  //int dmdt_range_aux;          // Accretion history array size
  int i;                       // Aux. counter
  //int ifirst;                  // i.d. of first particle
  //int ilast;                   // i.d. of last particle
  //int j;                       // Aux. counter
  int k;                       // Dimension counter
  int ndata;                   // No. of data arrays written
  int nunit;                   // No. of unit strings
  int idata[50];               // Integer data array
  int ilpdata[50];             // Long integer data array
  int typedata[50][5];         // Hydro particle data array information
  FLOAT rdata[50];             // Real data array
  DOUBLE ddata[50];            // Double float data array
  string format_id;            // File format (for verification)
  string data_id[50];          // String ids of arrays written
  string unit_data[50];        // String ids of units written
  int sink_data_length;        // Length of sink float array
  string dummystring;          // Dummy string variable
  ofstream outfile;            // Output file stream

  debug2("[Simulation::WriteSerenFormSnapshotFile]");

  cout << "Writing snapshot file : " << filename << endl;

  outfile.open(filename.c_str());
  outfile.setf (std::ios::scientific, std::ios::floatfield);
  formatted_output outfile_format(outfile, 18, 2, 10);

  for (i=0; i<50; i++) idata[i] = 0;
  for (i=0; i<50; i++) ilpdata[i] = 0;
  for (i=0; i<50; i++) rdata[i] = (FLOAT) 0.0;
  for (i=0; i<50; i++) ddata[i] = 0.0;
  nunit = 0;
  ndata = 0;

  // Set units
  if (!simunits.dimensionless) {
    unit_data[0]  = simunits.r.outunit;
    unit_data[1]  = simunits.m.outunit;
    unit_data[2]  = simunits.t.outunit;
    unit_data[3]  = simunits.v.outunit;
    unit_data[4]  = simunits.a.outunit;
    unit_data[5]  = simunits.rho.outunit;
    unit_data[6]  = simunits.sigma.outunit;
    unit_data[7]  = simunits.press.outunit;
    unit_data[8]  = simunits.f.outunit;
    unit_data[9]  = simunits.E.outunit;
    unit_data[10] = simunits.mom.outunit;
    unit_data[11] = simunits.angmom.outunit;
    unit_data[12] = simunits.angvel.outunit;
    unit_data[13] = simunits.dmdt.outunit;
    unit_data[14] = simunits.L.outunit;
    unit_data[15] = simunits.kappa.outunit;
    unit_data[16] = simunits.B.outunit;
    unit_data[17] = simunits.Q.outunit;
    unit_data[18] = simunits.Jcur.outunit;
    unit_data[19] = simunits.u.outunit;
    unit_data[20] = simunits.temp.outunit;
    nunit = 21;
  }


  // Set array ids and array information data if there are any hydro particles
  //-----------------------------------------------------------------------------------------------
  if (hydro->Nhydro > 0) {
    data_id[ndata] = "porig";
    typedata[ndata][0] = 1;              typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro;  typedata[ndata][3] = 2;
    typedata[ndata][4] = 0;              ndata++;

    data_id[ndata] = "r";
    typedata[ndata][0] = ndim;           typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro;  typedata[ndata][3] = 4;
    typedata[ndata][4] = 1;              ndata++;

    data_id[ndata] = "m";
    typedata[ndata][0] = 1;              typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro;  typedata[ndata][3] = 4;
    typedata[ndata][4] = 2;              ndata++;

    data_id[ndata] = "h";
    typedata[ndata][0] = 1;              typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro;  typedata[ndata][3] = 4;
    typedata[ndata][4] = 1;              ndata++;

    data_id[ndata] = "v";
    typedata[ndata][0] = ndim;           typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro;  typedata[ndata][3] = 4;
    typedata[ndata][4] = 4;              ndata++;

    data_id[ndata] = "rho";
    typedata[ndata][0] = 1;              typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro;  typedata[ndata][3] = 4;
    typedata[ndata][4] = 6;              ndata++;

    data_id[ndata] = "u";
    typedata[ndata][0] = 1;              typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro;  typedata[ndata][3] = 4;
    typedata[ndata][4] = 20;             ndata++;
  }

  if (nbody->Nstar > 0) {
    data_id[ndata] = "sink_v1";
    typedata[ndata][0] = 1;             typedata[ndata][1] = 1;
    typedata[ndata][2] = nbody->Nstar;  typedata[ndata][3] = 7;
    typedata[ndata][4] = 0;             ndata++;
  }

  // Set important header information
  idata[0]    = hydro->Nhydro;
  idata[1]    = nbody->Nstar;
  idata[4]    = hydro->Nhydro;
  idata[19]   = nunit;
  idata[20]   = ndata;
  ilpdata[0]  = Noutsnap;
  ilpdata[1]  = Nsteps;
  ilpdata[10] = Noutlitesnap;
  rdata[0]    = hydro->h_fac;
  rdata[1]    = (FLOAT) 0.0;
  ddata[0]    = t*simunits.t.outscale;
  ddata[1]    = tsnaplast*simunits.t.outscale;
  ddata[2]    = hydro->mmean*simunits.m.outscale;
  ddata[10]   = tlitesnaplast*simunits.t.outscale;

  // Write header information to file
  //-----------------------------------------------------------------------------------------------
  outfile_format << ascii_tag << endl;;
  outfile_format << 4 << endl;
  outfile_format << ndim << endl;
  outfile_format << ndim << endl;
  outfile_format << ndim << endl;
  outfile_format.set_width_integer(10);
  for (i=0; i<50; i++) outfile_format << idata[i] << endl;
  for (i=0; i<50; i++) outfile_format << ilpdata[i] << endl;
  for (i=0; i<50; i++) outfile_format << rdata[i] << endl;
  for (i=0; i<50; i++) outfile_format << ddata[i] << endl;
  if (nunit > 0) {
    for (i=0; i<nunit; i++) outfile_format << unit_data[i] << endl;
  }
  if (ndata > 0) {
    for (i=0; i<ndata; i++) outfile_format << data_id[i] << endl;
  }
  if (ndata > 0) {
    for (i=0; i<ndata; i++) {
      outfile_format << typedata[i][0] <<  typedata[i][1] << typedata[i][2]
                     << typedata[i][3] << typedata[i][4] << endl;
    }
  }

  // Write arrays for hydro particles
  //-----------------------------------------------------------------------------------------------
  if (hydro->Nhydro > 0) {

    // porig
    //---------------------------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      outfile_format << part.iorig << endl;
    }

    // Positions
    //---------------------------------------------------------------------------------------------
    if (ndim == 1) {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        outfile_format << part.r[0]*simunits.r.outscale << endl;
      }
    }
    else if (ndim == 2) {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        outfile_format << part.r[0]*simunits.r.outscale << part.r[1]*simunits.r.outscale << endl;
      }
    }
    else if (ndim == 3) {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        outfile_format << part.r[0]*simunits.r.outscale << part.r[1]*simunits.r.outscale
                       << part.r[2]*simunits.r.outscale << endl;
      }
    }

    // Masses
    //---------------------------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      outfile_format << part.m*simunits.m.outscale << endl;
    }

    // Smoothing lengths
    //---------------------------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      outfile_format << part.h*simunits.r.outscale << endl;
    }

    // Velocities
    //---------------------------------------------------------------------------------------------
    if (ndim == 1) {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        outfile_format << part.v[0]*simunits.v.outscale << endl;
      }
    }
    else if (ndim == 2) {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        outfile_format << part.v[0]*simunits.v.outscale
                << part.v[1]*simunits.v.outscale << endl;
      }
    }
    else if (ndim == 3) {
      for (i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        outfile_format << part.v[0]*simunits.v.outscale
                << part.v[1]*simunits.v.outscale
                << part.v[2]*simunits.v.outscale << endl;
      }
    }

    // Densities
    //---------------------------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      outfile_format << part.rho*simunits.rho.outscale << endl;;
    }

    // Specific internal energies
    //---------------------------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      outfile_format << part.u*simunits.u.outscale << endl;
    }

  }
  //-----------------------------------------------------------------------------------------------


  // Sinks/stars
  //-----------------------------------------------------------------------------------------------
  if (nbody->Nstar > 0) {
    sink_data_length = 12 + 2*ndim; //+ 2*dmdt_range_aux;
    int ii;
    FLOAT sdata[sink_data_length];
    for (k=0; k<sink_data_length; k++) sdata[k] = 0.0;
    outfile_format << 2  << 2  << 0 << sink_data_length  << 0 << 0 << endl;
    outfile_format.set_width_integer(8);
    for (i=0; i<nbody->Nstar; i++) {
      outfile_format << true << true << endl;
      outfile_format << i+1  << 0 << endl;
      for (k=0; k<ndim; k++) sdata[k+1] = nbody->stardata[i].r[k]*simunits.r.outscale;
      for (k=0; k<ndim; k++) sdata[k+1+ndim] = nbody->stardata[i].v[k]*simunits.v.outscale;
      sdata[1+2*ndim] = nbody->stardata[i].m*simunits.m.outscale;
      sdata[2+2*ndim] = nbody->stardata[i].h*simunits.r.outscale;
      sdata[3+2*ndim] = nbody->stardata[i].radius*simunits.r.outscale;
      for (ii=0; ii<sink_data_length; ii++) outfile_format << sdata[ii] ;
      outfile << endl;
    }
  }
  //-----------------------------------------------------------------------------------------------

  // Close file
  outfile.close();

  return true;
}



//=================================================================================================
//  Simulation::ReadSerenUnformHeaderFile
/// Function for reading the header file of a snapshot. Does not modify the variables
/// of the Simulation class, but rather returns information in a HeaderInfo struct.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ReadSerenUnformHeaderFile
 (ifstream& infile,                    ///< Input file stream
  HeaderInfo& info)                    ///< Header info data structure
{
  BinaryReader reader(infile);         // ..

  debug2("[Simulation::ReadSerenUnformHeaderFile]");

  // Skip the first bits (tag + the 4)
  int id_length = 20; //binary_tag.size();
  infile.seekg(id_length);
  infile.seekg(4, ios_base::cur);

  // Read number of dimensions
  reader.read_value(info.ndim);

  // Skip the following two integers...
  infile.seekg(8, ios_base::cur);

  // Read number of hydro particles
  reader.read_value(info.Nhydro);

  // Read number of star particles
  reader.read_value(info.Nstar);

  // Skip the remaining 48 idata, ilpdata and rdata
  infile.seekg(48*sizeof(int) + 50*sizeof(long) + 50*sizeof(FLOAT), ios_base::cur);

  // Read time and convert to code units
  reader.read_value(info.t);
  info.t /= simunits.t.inscale;

  // Check dimensionality matches if using fixed dimensions
  if (info.ndim != ndim) {
    std::ostringstream stream;
    stream << "Incorrect no. of dimensions in file : "
           << info.ndim << "  [ndim : " << ndim << "]" << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }

  return;
}


//=================================================================================================
//  Simulation::ReadSerenUnformSnapshotFile
/// Read data from Seren binary format snapshot file into main arrays.
//=================================================================================================
template <int ndim>
bool Simulation<ndim>::ReadSerenUnformSnapshotFile(string filename)
{
  int dummy;                        // ..
  //int dmdt_range_aux;               // Accretion history array size
  int ndata;                        // No. of data arrays written
  int nunit;                        // No. of unit strings
  int typedata[50][5];              // ..
  int idata[50];                    // ..
  long ilpdata[50];                 // ..
  FLOAT rdata[50];                  // ..
  DOUBLE ddata[50];                 // ..
  string data_id[50];               // String ids of arrays written
  string unit_data[50];             // String ids of units written

  debug2("[Simulation::ReadSerenUnformSnapshotFile]");

  ifstream infile(filename.c_str());
  BinaryReader reader(infile);

  // Read information identifying format and precision of file.
  // Then check if each value corresponds to the current values.
  //-----------------------------------------------------------------------------------------------

  //Tag
  std::vector<char> file_tag(string_length);
  infile.read(&file_tag[0], string_length);
  string file_tag_string (&file_tag[0],string_length);
  file_tag_string = simparams->TrimWhiteSpace(file_tag_string);
  if (file_tag_string != binary_tag) {
    std::ostringstream stream;
    stream << "Incorrect format of snapshot file " << filename << ": " << file_tag_string <<
        " instead of " << binary_tag << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }

  // Precision
  reader.read_value(dummy);
#if defined GANDALF_DOUBLE_PRECISION
  if (dummy != 8) {
#else
  if (dummy !=4) {
#endif
    std::ostringstream stream;
    stream << "Incorrect precision of snapshot file " << filename << ": " << dummy <<
        " but the code was compiled with a FLOAT size of " << sizeof(FLOAT) << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }

  // Dimensions
  int ndim_file;
  reader.read_value(ndim_file);
  if (ndim_file != ndim) {
    std::ostringstream stream;
    stream << "Incorrect NDIM in file " << filename << " : got " << ndim_file << ", expected " << ndim << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }
  reader.read_value(ndim_file);
  if (ndim_file != ndim) {
    std::ostringstream stream;
    stream << "Incorrect NDIM in file " << filename << " : got " << ndim_file << ", expected " << ndim << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }
  reader.read_value(ndim_file);
  if (ndim_file != ndim) {
    std::ostringstream stream;
    stream << "Incorrect NDIM in file " << filename << " : got " << ndim_file << ", expected " << ndim << endl;
    ExceptionHandler::getIstance().raise(stream.str());
  }

  // Read infile header integer data
  {
    for (int i=0; i<50; i++) reader.read_value(idata[i]);
    hydro->Nhydro  = idata[0];
    nbody->Nstar   = idata[1];
    sinks->Nsink    = idata[1];
    //dmdt_range_aux = idata[29];
    nunit          = idata[19];
    ndata          = idata[20];
  }
  {
    for (int i=0; i<50; i++) reader.read_value(ilpdata[i]);
    Nsteps         = ilpdata[1];
  }

  // Read infile header floating point data
  {
    for (int i=0; i<50; i++) reader.read_value(rdata[i]);
    for (int i=0; i<50; i++) reader.read_value(ddata[i]);
    t = ddata[0];
  }

  // Variables that should be remembered for restarts
  if (restart) {
    Noutsnap      = ilpdata[0];
    Nsteps        = ilpdata[1];
    Noutlitesnap  = ilpdata[10];
    t             = ddata[0];
    tsnaplast     = ddata[1];
    hydro->mmean  = ddata[2];
    tlitesnaplast = ddata[10];
  }

  // Read unit_data
  for (int i=0; i<nunit; i++) {
    char buffer[string_length];
    infile.read(buffer,string_length);
    unit_data[i] = std::string(buffer,string_length);
    unit_data[i] = simparams->TrimWhiteSpace(unit_data[i]);
  }

  // Read data_id
  for (int i=0; i< ndata; i++) {
    char buffer[string_length];
    infile.read(buffer,string_length);
    data_id[i] = std::string(buffer,string_length);
    data_id[i] = simparams->TrimWhiteSpace(data_id[i]);
  }

  // Read infile ids of arrays contained in file
  for (int i=0; i<ndata; i++) {
    for (int j=0; j<5; j++) reader.read_value(typedata[i][j]);
  }

  // Allocate memory now we know particle numbers
  AllocateParticleMemory();


  // Loop through array ids and read each array in turn
  //===============================================================================================
  for (int j=0; j<ndata; j++) {

    // Comment these out for now (may need in future for different particle types)
    //int ifirst = typedata[j][1];
    //int ilast = typedata[j][2];

    // porig
    //---------------------------------------------------------------------------------------------
    if (data_id[j] == "porig") {
      for (int i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        reader.read_value(part.iorig);
      }
    }

    // Positions
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "r") {
      for (int i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (int k=0; k<ndim; k++)
          reader.read_value(part.r[k]);
      }
    }

    // Masses
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "m") {
      for (int i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        reader.read_value(part.m);
        assert(part.m > (FLOAT) 0.0);
      }
    }

    // Smoothing lengths
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "h") {
      for (int i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        reader.read_value(part.h);
      }
    }

    // Velocities
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "v") {
      for (int i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        for (int k=0; k< ndim; k++)
          reader.read_value(part.v[k]);
      }
    }

    // Densities
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "rho") {
      for (int i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        reader.read_value(part.rho);
      }
    }

    // Specific internal energies
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "u") {
      for (int i=0; i<hydro->Nhydro; i++) {
        Particle<ndim>& part = hydro->GetParticlePointer(i);
        reader.read_value(part.u);
      }
    }

    // Skip 1-D redundant information
    //---------------------------------------------------------------------------------------------
    else if(data_id[j] == "temp") {
      infile.seekg(sizeof(FLOAT)*hydro->Nhydro, ios_base::cur);
    }

    // Sinks/stars
    //---------------------------------------------------------------------------------------------
    else if (data_id[j] == "sink_v1") {
      int sink_data_length = 12 + 2*ndim;  //+ 2*dmdt_range_aux;
      int ii;
      FLOAT sdata[sink_data_length];
      int idata[50];
      for (ii=0; ii<6; ii++) reader.read_value(idata[ii]);
      if (nbody->Nstar > 0) {
        for (int i=0; i<nbody->Nstar; i++) {
          //Logical flags
          for (ii=0; ii<2; ii++) {
            int boolean;
            reader.read_value(boolean);
          }
          for (ii=0; ii<2; ii++) reader.read_value(idata[ii]);
          for (ii=0; ii<sink_data_length; ii++) reader.read_value(sdata[ii]);
          for (int k=0; k<ndim; k++) nbody->stardata[i].r[k] = sdata[k+1];
          for (int k=0; k<ndim; k++) nbody->stardata[i].v[k] = sdata[k+1+ndim];
          nbody->stardata[i].m = sdata[1+2*ndim];
          nbody->stardata[i].h = sdata[2+2*ndim];
          nbody->stardata[i].radius = sdata[3+2*ndim];
          sinks->sink[i].radius = sdata[3+2*ndim];
          sinks->sink[i].star = &(nbody->stardata[i]);
          nbody->nbodydata[i] = &(nbody->stardata[i]);
          assert(nbody->stardata[i].m > 0.0);
        }
      }
    }

    // Skip through arbitrary 1D or 2D array
    //---------------------------------------------------------------------------------------------
    else if (typedata[j][0] >= 1) {
      ExceptionHandler::getIstance().raise("Arbitrary data reading not implemented!");
    }

  }
  //===============================================================================================


  infile.close();

  return true;

}



//=================================================================================================
//  Simulation::WriteSerenUnformSnapshotFile
/// Write hydro and N-body particle data to snapshot file in Seren binary format.
//=================================================================================================
#ifdef MPI_PARALLEL
template <int ndim>
bool Simulation<ndim>::WriteSerenUnformSnapshotFile(string filename)
{
  int i;                            // Aux. counter
  int idata[50];                    // Integer data array
  int ii;                           // Aux. counter
  int k;                            // Aux. loop counter
  int typedata[50][5];              // Hydro particle data array information
  int ndata;                        // No. of data arrays written
  int nunit;                        // No. of unit strings
  int sink_data_length = 12+2*ndim; // (+ 2*dmdt_range_aux);
  long ilpdata[50];                 // Long integer data array
  FLOAT rdata[50];                  // Real data array
  FLOAT sdata[sink_data_length];    // Sink data packet
  DOUBLE ddata[50];                 // Double float data array
  string unit_data[50];             // String ids of units written
  string data_id[50];               // String ids of arrays written

  debug2("[Simulation::WriteSerenUnformSnapshotFile]");

  if (rank==0)
    cout << "Writing snapshot file : " << filename << endl;

  // Total number of particles
  int Ntot_hydro;
  MPI_Allreduce(&hydro->Nhydro,&Ntot_hydro,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  // Zero arrays
  for (i=0; i<50; i++) idata[i] = 0;
  for (i=0; i<50; i++) ilpdata[i] = 0;
  for (i=0; i<50; i++) rdata[i] = 0.0;
  for (i=0; i<50; i++) ddata[i] = 0.0;
  nunit = 0;
  ndata = 0;

  // Set units
  if (!simunits.dimensionless) {
    unit_data[0] = simunits.r.outunit;
    unit_data[1] = simunits.m.outunit;
    unit_data[2] = simunits.t.outunit;
    unit_data[3] = simunits.v.outunit;
    unit_data[4] = simunits.a.outunit;
    unit_data[5] = simunits.rho.outunit;
    unit_data[6] = simunits.sigma.outunit;
    unit_data[7] = simunits.press.outunit;
    unit_data[8] = simunits.f.outunit;
    unit_data[9] = simunits.E.outunit;
    unit_data[10] = simunits.mom.outunit;
    unit_data[11] = simunits.angmom.outunit;
    unit_data[12] = simunits.angvel.outunit;
    unit_data[13] = simunits.dmdt.outunit;
    unit_data[14] = simunits.L.outunit;
    unit_data[15] = simunits.kappa.outunit;
    unit_data[16] = simunits.B.outunit;
    unit_data[17] = simunits.Q.outunit;
    unit_data[18] = simunits.Jcur.outunit;
    unit_data[19] = simunits.u.outunit;
    unit_data[20] = simunits.temp.outunit;
    nunit = 21;
  }

  // Set array ids and array information data if there are any hydro particles
  //---------------------------------------------------------------------------
  if (Ntot_hydro > 0) {
    data_id[ndata] = "porig";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = Ntot_hydro; typedata[ndata][3] = 2;
    typedata[ndata][4] = 0; ndata++;

    data_id[ndata] = "r";
    typedata[ndata][0] = ndim; typedata[ndata][1] = 1;
    typedata[ndata][2] = Ntot_hydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 1; ndata++;

    data_id[ndata] = "m";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = Ntot_hydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 2; ndata++;

    data_id[ndata] = "h";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = Ntot_hydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 1; ndata++;

    data_id[ndata] = "v";
    typedata[ndata][0] = ndim; typedata[ndata][1] = 1;
    typedata[ndata][2] = Ntot_hydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 4; ndata++;

    data_id[ndata] = "rho";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = Ntot_hydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 6; ndata++;

    data_id[ndata] = "u";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = Ntot_hydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 20; ndata++;
  }

  if (nbody->Nstar > 0) {
    data_id[ndata] = "sink_v1";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = nbody->Nstar; typedata[ndata][3] = 7;
    typedata[ndata][4] = 0; ndata++;
  }

  // Set important header information
  idata[0]    = Ntot_hydro;
  idata[1]    = nbody->Nstar;
  idata[4]    = Ntot_hydro;
  idata[19]   = nunit;
  idata[20]   = ndata;
  ilpdata[0]  = Noutsnap;
  ilpdata[1]  = Nsteps;
  ilpdata[10] = Noutlitesnap;
  rdata[0]    = hydro->h_fac;
  rdata[1]    = 0.0;
  ddata[0]    = t*simunits.t.outscale;
  ddata[1]    = tsnaplast*simunits.t.outscale;
  ddata[2]    = hydro->mmean*simunits.m.outscale;
  ddata[10]   = tlitesnaplast*simunits.t.outscale;


  // Write header information to file
  // Only cpu 0 does it
  //---------------------------------------------------------------------------

  if (rank==0) {
    ofstream outfile(filename.c_str(),ios::binary);
    BinaryWriter writer(outfile);
    std::ostringstream stream;
    stream << std::left << std::setw(string_length) << std::setfill(' ')
                      << binary_tag;
    outfile << stream.str();

#if defined GANDALF_DOUBLE_PRECISION
    writer.write_value(8);
#else
    writer.write_value(4);
#endif
    writer.write_value(ndim);
    writer.write_value(ndim);
    writer.write_value(ndim);
    for (i=0; i<50; i++) writer.write_value(idata[i]);
    for (i=0; i<50; i++) writer.write_value(ilpdata[i]);
    for (i=0; i<50; i++) writer.write_value(rdata[i]);
    for (i=0; i<50; i++) writer.write_value(ddata[i]);
    for (i=0; i<nunit; i++) {
      std::ostringstream stream;
      stream << std::left << std::setw(string_length) << std::setfill(' ')
                    << unit_data[i];
      outfile << stream.str();
    }
    for (i=0; i<ndata; i++) {
      std::ostringstream stream;
      stream << std::left << std::setw(string_length) << std::setfill(' ')
                << data_id[i];
      outfile << stream.str();
    }
    for (i=0; i<ndata; i++) {
      for (int j=0; j< 5; j++) writer.write_value(typedata[i][j]);
    }
  }


  // Write arrays for hydro particles
  //---------------------------------------------------------------------------
  if (Ntot_hydro > 0) {

    // Let's make sure root has finished writing the header before we start writing
    MPI_Barrier(MPI_COMM_WORLD);

    // Now everyone opens the output file
    MPI_File file;
    char* filename_str = new char[strlen(filename.c_str())+1];
    strcpy(filename_str,filename.c_str());
    MPI_File_open(MPI_COMM_WORLD, filename_str, MPI_MODE_RDWR,MPI_INFO_NULL, &file);
    delete[] filename_str;
    MPI_Offset end_header;
    MPI_File_get_size(file,&end_header);
    MPI_File_seek(file,end_header,MPI_SEEK_SET);

    // We need to know how much to seek to start writing -
    // need to know the cumulative number of particles
    int Nhydro_before;
    MPI_Exscan(&hydro->Nhydro,&Nhydro_before,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (rank==0) {
      Nhydro_before = 0;
    }

    // Starting point - let's remember it
    MPI_Offset end_previous_write;

    // Buffer
    void* buffer = malloc(sizeof(FLOAT)*ndim*hydro->Nhydro);

    // porig
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      int* buffer_int = (int*) buffer;
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      buffer_int[i] = part.iorig;
    }
    // Seek at the right position in the file
    MPI_File_seek(file, sizeof(int)*Nhydro_before  ,MPI_SEEK_CUR);
    // Write date
    MPI_Status status;
    MPI_File_write_all (file, buffer, hydro->Nhydro, MPI_INT, &status);
    // Seek at the end of the porig section
    end_previous_write = end_header+sizeof(int)*Ntot_hydro;
    MPI_File_seek(file, end_previous_write, MPI_SEEK_SET);

    // Positions
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      FLOAT* buffer_float = (FLOAT*) buffer;
      for (int k=0; k<ndim; k++) buffer_float[ndim*i+k] = part.r[k]*simunits.r.outscale;
    }
    MPI_File_seek(file, sizeof(FLOAT)*ndim*Nhydro_before,MPI_SEEK_CUR);
    MPI_File_write_all (file, buffer, ndim*hydro->Nhydro, GANDALF_MPI_FLOAT, &status);
    end_previous_write += sizeof(FLOAT)*ndim*Ntot_hydro;
    MPI_File_seek(file, end_previous_write, MPI_SEEK_SET);


    // Masses
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      FLOAT* buffer_float = (FLOAT*) buffer;
      buffer_float[i] = part.m*simunits.m.outscale;
      assert(part.m > 0.0);
    }
    MPI_File_seek(file, sizeof(FLOAT)*Nhydro_before,MPI_SEEK_CUR);
    MPI_File_write_all (file, buffer, hydro->Nhydro, GANDALF_MPI_FLOAT, &status);
    end_previous_write += sizeof(FLOAT)*Ntot_hydro;
    MPI_File_seek(file, end_previous_write, MPI_SEEK_SET);

    // Smoothing lengths
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      FLOAT* buffer_float = (FLOAT*) buffer;
      buffer_float[i] = part.h*simunits.r.outscale;
    }
    MPI_File_seek(file, sizeof(FLOAT)*Nhydro_before,MPI_SEEK_CUR);
    MPI_File_write_all (file, buffer, hydro->Nhydro, GANDALF_MPI_FLOAT, &status);
    end_previous_write += sizeof(FLOAT)*Ntot_hydro;
    MPI_File_seek(file, end_previous_write, MPI_SEEK_SET);

    // Velocities
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      FLOAT* buffer_float = (FLOAT*) buffer;
      for (int k=0; k<ndim; k++)
        buffer_float[ndim*i+k] = part.v[k]*simunits.v.outscale;
    }
    MPI_File_seek(file, sizeof(FLOAT)*ndim*Nhydro_before,MPI_SEEK_CUR);
    MPI_File_write_all (file, buffer, ndim*hydro->Nhydro, GANDALF_MPI_FLOAT, &status);
    end_previous_write += sizeof(FLOAT)*ndim*Ntot_hydro;
    MPI_File_seek(file, end_previous_write, MPI_SEEK_SET);

    // Densities
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      FLOAT* buffer_float = (FLOAT*) buffer;
      buffer_float[i] = part.rho*simunits.rho.outscale;
    }
    MPI_File_seek(file, sizeof(FLOAT)*Nhydro_before,MPI_SEEK_CUR);
    MPI_File_write_all (file, buffer, hydro->Nhydro, GANDALF_MPI_FLOAT, &status);
    end_previous_write += sizeof(FLOAT)*Ntot_hydro;
    MPI_File_seek(file, end_previous_write, MPI_SEEK_SET);

    // Specific internal energies
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      FLOAT* buffer_float = (FLOAT*) buffer;
      buffer_float[i] = part.u*simunits.u.outscale;
    }
    MPI_File_seek(file, sizeof(FLOAT)*Nhydro_before,MPI_SEEK_CUR);
    MPI_File_write_all (file, buffer, hydro->Nhydro, GANDALF_MPI_FLOAT, &status);
    end_previous_write += sizeof(FLOAT)*Ntot_hydro;
    MPI_File_seek(file, end_previous_write, MPI_SEEK_SET);

    free(buffer);

    MPI_Offset end_file;
    MPI_File_get_position(file, &end_file);
    assert(end_file == end_header+(sizeof(FLOAT)*ndim*2+sizeof(FLOAT)*4+sizeof(int)*1 ) *Ntot_hydro );

    MPI_File_close(&file);


  }


  // Sinks/stars
  //---------------------------------------------------------------------------
  if (rank==0 && nbody->Nstar > 0) {

    ofstream outfile(filename.c_str(),ios::binary|ios::app);
    BinaryWriter writer(outfile);

    for (k=0; k<sink_data_length; k++) sdata[k] = 0.0;
    int values[6] = {2,2,0,sink_data_length,0,0};
    for (i=0; i<6;i++)
      writer.write_value(values[i]);
    for (i=0; i<nbody->Nstar; i++) {
      writer.write_value(true); writer.write_value(true);
      writer.write_value(i+1); writer.write_value(0);
      for (k=0; k<ndim; k++)
        sdata[k+1] = nbody->stardata[i].r[k]*simunits.r.outscale;
      for (k=0; k<ndim; k++)
        sdata[k+1+ndim] = nbody->stardata[i].v[k]*simunits.v.outscale;
      sdata[1+2*ndim] = nbody->stardata[i].m*simunits.m.outscale;
      sdata[2+2*ndim] = nbody->stardata[i].h*simunits.r.outscale;
      sdata[3+2*ndim] = nbody->stardata[i].radius*simunits.r.outscale;
      for (ii=0; ii<sink_data_length; ii++) {
        writer.write_value(sdata[ii]);
      }
    }
  }
  //---------------------------------------------------------------------------



  return true;
}
#else
template <int ndim>
bool Simulation<ndim>::WriteSerenUnformSnapshotFile(string filename)
{
  int i;                            // Aux. counter
  int idata[50];                    // Integer data array
  int ii;                           // Aux. counter
  int k;                            // Aux. loop counter
  int typedata[50][5];              // Hydro particle data array information
  int ndata;                        // No. of data arrays written
  int nunit;                        // No. of unit strings
  int sink_data_length = 12+2*ndim; // (+ 2*dmdt_range_aux);
  long ilpdata[50];                 // Long integer data array
  FLOAT rdata[50];                  // Real data array
  FLOAT sdata[sink_data_length];    // Sink data packet
  DOUBLE ddata[50];                 // Double float data array
  string unit_data[50];             // String ids of units written
  string data_id[50];               // String ids of arrays written

  debug2("[Simulation::WriteSerenUnformSnapshotFile]");

  cout << "Writing snapshot file : " << filename << endl;

  ofstream outfile(filename.c_str(),ios::binary);
  BinaryWriter writer(outfile);

  // Zero arrays
  for (i=0; i<50; i++) idata[i] = 0;
  for (i=0; i<50; i++) ilpdata[i] = 0;
  for (i=0; i<50; i++) rdata[i] = 0.0;
  for (i=0; i<50; i++) ddata[i] = 0.0;
  nunit = 0;
  ndata = 0;

  // Set units
  if (!simunits.dimensionless) {
    unit_data[0] = simunits.r.outunit;
    unit_data[1] = simunits.m.outunit;
    unit_data[2] = simunits.t.outunit;
    unit_data[3] = simunits.v.outunit;
    unit_data[4] = simunits.a.outunit;
    unit_data[5] = simunits.rho.outunit;
    unit_data[6] = simunits.sigma.outunit;
    unit_data[7] = simunits.press.outunit;
    unit_data[8] = simunits.f.outunit;
    unit_data[9] = simunits.E.outunit;
    unit_data[10] = simunits.mom.outunit;
    unit_data[11] = simunits.angmom.outunit;
    unit_data[12] = simunits.angvel.outunit;
    unit_data[13] = simunits.dmdt.outunit;
    unit_data[14] = simunits.L.outunit;
    unit_data[15] = simunits.kappa.outunit;
    unit_data[16] = simunits.B.outunit;
    unit_data[17] = simunits.Q.outunit;
    unit_data[18] = simunits.Jcur.outunit;
    unit_data[19] = simunits.u.outunit;
    unit_data[20] = simunits.temp.outunit;
    nunit = 21;
  }

  // Set array ids and array information data if there are any hydro particles
  //---------------------------------------------------------------------------
  if (hydro->Nhydro > 0) {
    data_id[ndata] = "porig";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 2;
    typedata[ndata][4] = 0; ndata++;

    data_id[ndata] = "r";
    typedata[ndata][0] = ndim; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 1; ndata++;

    data_id[ndata] = "m";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 2; ndata++;

    data_id[ndata] = "h";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 1; ndata++;

    data_id[ndata] = "v";
    typedata[ndata][0] = ndim; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 4; ndata++;

    data_id[ndata] = "rho";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 6; ndata++;

    data_id[ndata] = "u";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 20; ndata++;
  }

  if (nbody->Nstar > 0) {
    data_id[ndata] = "sink_v1";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = nbody->Nstar; typedata[ndata][3] = 7;
    typedata[ndata][4] = 0; ndata++;
  }

  // Set important header information
  idata[0]    = hydro->Nhydro;
  idata[1]    = nbody->Nstar;
  idata[4]    = hydro->Nhydro;
  idata[19]   = nunit;
  idata[20]   = ndata;
  ilpdata[0]  = Noutsnap;
  ilpdata[1]  = Nsteps;
  ilpdata[10] = Noutlitesnap;
  rdata[0]    = hydro->h_fac;
  rdata[1]    = 0.0;
  ddata[0]    = t*simunits.t.outscale;
  ddata[1]    = tsnaplast*simunits.t.outscale;
  ddata[2]    = hydro->mmean*simunits.m.outscale;
  ddata[10]   = tlitesnaplast*simunits.t.outscale;


  // Write header information to file
  //---------------------------------------------------------------------------
  {
    std::ostringstream stream;
    stream << std::left << std::setw(string_length) << std::setfill(' ')
                      << binary_tag;
    outfile << stream.str();
  }


#if defined GANDALF_DOUBLE_PRECISION
  writer.write_value(8);
#else
  writer.write_value(4);
#endif
  writer.write_value(ndim);
  writer.write_value(ndim);
  writer.write_value(ndim);
  for (i=0; i<50; i++) writer.write_value(idata[i]);
  for (i=0; i<50; i++) writer.write_value(ilpdata[i]);
  for (i=0; i<50; i++) writer.write_value(rdata[i]);
  for (i=0; i<50; i++) writer.write_value(ddata[i]);
  for (i=0; i<nunit; i++) {
    std::ostringstream stream;
    stream << std::left << std::setw(string_length) << std::setfill(' ')
                  << unit_data[i];
    outfile << stream.str();
  }
  for (i=0; i<ndata; i++) {
    std::ostringstream stream;
    stream << std::left << std::setw(string_length) << std::setfill(' ')
              << data_id[i];
    outfile << stream.str();
  }
  for (i=0; i<ndata; i++) {
    for (int j=0; j< 5; j++) writer.write_value(typedata[i][j]);
  }

  // Write arrays for hydro particles
  //---------------------------------------------------------------------------
  if (hydro->Nhydro > 0) {

    // porig
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      writer.write_value(part.iorig);
    }

    // Positions
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (int k=0; k<ndim; k++) writer.write_value((FLOAT)(part.r[k]*simunits.r.outscale));
    }

    // Masses
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      writer.write_value((FLOAT)(part.m*simunits.m.outscale));
      assert(part.m > 0.0);
    }

    // Smoothing lengths
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      writer.write_value((FLOAT)(part.h*simunits.r.outscale));
    }

    // Velocities
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (int k=0; k<ndim; k++)
        writer.write_value((FLOAT)(part.v[k]*simunits.v.outscale));
    }

    // Densities
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      writer.write_value((FLOAT)(part.rho*simunits.rho.outscale));
    }

    // Specific internal energies
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      writer.write_value((FLOAT)(part.u*simunits.u.outscale));
    }

  }


  // Sinks/stars
  //---------------------------------------------------------------------------
  if (nbody->Nstar > 0) {
    for (k=0; k<sink_data_length; k++) sdata[k] = 0.0;
    int values[6] = {2,2,0,sink_data_length,0,0};
    for (i=0; i<6;i++)
      writer.write_value(values[i]);
    for (i=0; i<nbody->Nstar; i++) {
      writer.write_value(true); writer.write_value(true);
      writer.write_value(i+1); writer.write_value(0);
      for (k=0; k<ndim; k++)
        sdata[k+1] = nbody->stardata[i].r[k]*simunits.r.outscale;
      for (k=0; k<ndim; k++)
        sdata[k+1+ndim] = nbody->stardata[i].v[k]*simunits.v.outscale;
      sdata[1+2*ndim] = nbody->stardata[i].m*simunits.m.outscale;
      sdata[2+2*ndim] = nbody->stardata[i].h*simunits.r.outscale;
      sdata[3+2*ndim] = nbody->stardata[i].radius*simunits.r.outscale;
      for (ii=0; ii<sink_data_length; ii++) {
        writer.write_value(sdata[ii]);
      }
    }
  }
  //---------------------------------------------------------------------------

  outfile.close();

  return true;
}
#endif


//=================================================================================================
//  Simulation::WriteSerenLiteSnapshotFile
/// Write hydro and N-body particle data to snapshot file in Seren 'lite' format
/// (i.e. stripped down, low-memory data for basic visualisation and movies).
//=================================================================================================
template <int ndim>
bool Simulation<ndim>::WriteSerenLiteSnapshotFile(string filename)
{
  int i;                            // Aux. counter
  int idata[50];                    // Integer data array
  int ii;                           // Aux. counter
  int k;                            // Aux. loop counter
  int typedata[50][5];              // Hydro particle data array information
  int ndata;                        // No. of data arrays written
  int nunit;                        // No. of unit strings
  int sink_data_length = 12+2*ndim; // (+ 2*dmdt_range_aux);
  long ilpdata[50];                 // Long integer data array
  float rdata[50];                  // Real data array
  float sdata[sink_data_length];    // Sink data packet
  double ddata[50];                 // Double float data array
  string unit_data[50];             // String ids of units written
  string data_id[50];               // String ids of arrays written


  debug2("[Simulation::WriteSerenLiteSnapshotFile]");

  cout << "Writing snapshot file : " << filename << endl;

  ofstream outfile(filename.c_str(),ios::binary);
  BinaryWriter writer(outfile);

  // Zero arrays
  for (i=0; i<50; i++) idata[i] = 0;
  for (i=0; i<50; i++) ilpdata[i] = 0;
  for (i=0; i<50; i++) rdata[i] = 0.0;
  for (i=0; i<50; i++) ddata[i] = 0.0;
  nunit = 0;
  ndata = 0;

  // Set units
  if (!simunits.dimensionless) {
    unit_data[0] = simunits.r.outunit;
    unit_data[1] = simunits.m.outunit;
    unit_data[2] = simunits.t.outunit;
    unit_data[3] = simunits.v.outunit;
    unit_data[4] = simunits.a.outunit;
    unit_data[5] = simunits.rho.outunit;
    unit_data[6] = simunits.sigma.outunit;
    unit_data[7] = simunits.press.outunit;
    unit_data[8] = simunits.f.outunit;
    unit_data[9] = simunits.E.outunit;
    unit_data[10] = simunits.mom.outunit;
    unit_data[11] = simunits.angmom.outunit;
    unit_data[12] = simunits.angvel.outunit;
    unit_data[13] = simunits.dmdt.outunit;
    unit_data[14] = simunits.L.outunit;
    unit_data[15] = simunits.kappa.outunit;
    unit_data[16] = simunits.B.outunit;
    unit_data[17] = simunits.Q.outunit;
    unit_data[18] = simunits.Jcur.outunit;
    unit_data[19] = simunits.u.outunit;
    unit_data[20] = simunits.temp.outunit;
    nunit = 21;
  }

  // Set array ids and array information data if there are any hydro particles
  //---------------------------------------------------------------------------
  if (hydro->Nhydro > 0) {

    data_id[ndata] = "r";
    typedata[ndata][0] = ndim; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 1; ndata++;

    data_id[ndata] = "m";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 2; ndata++;

    data_id[ndata] = "h";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 1; ndata++;

    data_id[ndata] = "rho";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 6; ndata++;

    data_id[ndata] = "u";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = hydro->Nhydro; typedata[ndata][3] = 4;
    typedata[ndata][4] = 20; ndata++;
  }

  if (nbody->Nstar > 0) {
    data_id[ndata] = "sink_v1";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = nbody->Nstar; typedata[ndata][3] = 7;
    typedata[ndata][4] = 0; ndata++;
  }

  // Set important header information
  idata[0]   = hydro->Nhydro;
  idata[1]   = nbody->Nstar;
  idata[4]   = hydro->Nhydro;
  idata[19]  = nunit;
  idata[20]  = ndata;
  ilpdata[0] = Noutsnap;
  ilpdata[1] = Nsteps;
  rdata[0]   = hydro->h_fac;
  rdata[1]   = 0.0;
  ddata[0]   = t*simunits.t.outscale;
  ddata[1]   = tsnaplast*simunits.t.outscale;
  ddata[2]   = hydro->mmean*simunits.m.outscale;


  // Write header information to file
  //---------------------------------------------------------------------------
  {
    std::ostringstream stream;
    stream << std::left << std::setw(string_length) << std::setfill(' ') << binary_tag;
    outfile << stream.str();
  }

  // Hard-wired to single precision for low-memory usage
  writer.write_value(4);
  writer.write_value(ndim);
  writer.write_value(ndim);
  writer.write_value(ndim);
  for (i=0; i<50; i++) writer.write_value(idata[i]);
  for (i=0; i<50; i++) writer.write_value(ilpdata[i]);
  for (i=0; i<50; i++) writer.write_value(rdata[i]);
  for (i=0; i<50; i++) writer.write_value(ddata[i]);
  for (i=0; i<nunit; i++) {
    std::ostringstream stream;
    stream << std::left << std::setw(string_length) << std::setfill(' ') << unit_data[i];
    outfile << stream.str();
  }
  for (i=0; i<ndata; i++) {
    std::ostringstream stream;
    stream << std::left << std::setw(string_length) << std::setfill(' ') << data_id[i];
    outfile << stream.str();
  }
  for (i=0; i<ndata; i++) {
    for (int j=0; j< 5; j++) writer.write_value(typedata[i][j]);
  }


  // Write arrays for hydro particles
  //---------------------------------------------------------------------------
  if (hydro->Nhydro > 0) {

    // Positions
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      for (int k=0; k<ndim; k++) writer.write_value((float) (part.r[k]*simunits.r.outscale));
    }

    // Masses
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      writer.write_value((float) (part.m*simunits.m.outscale));
    }

    // Smoothing lengths
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      writer.write_value((float) (part.h*simunits.r.outscale));
    }

    // Densities
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      writer.write_value((float) (part.rho*simunits.rho.outscale));
    }

    // Specific internal energies
    //-------------------------------------------------------------------------
    for (i=0; i<hydro->Nhydro; i++) {
      Particle<ndim>& part = hydro->GetParticlePointer(i);
      writer.write_value((float) (part.u*simunits.u.outscale));
    }

  }


  // Sinks/stars
  //---------------------------------------------------------------------------
  if (nbody->Nstar > 0) {
    for (k=0; k<sink_data_length; k++) sdata[k] = 0.0;
    int values[6] = {2,2,0,sink_data_length,0,0};
    for (i=0; i<6;i++)
      writer.write_value(values[i]);
    for (i=0; i<nbody->Nstar; i++) {
      writer.write_value(true); writer.write_value(true);
      writer.write_value(i+1); writer.write_value(0);
      for (k=0; k<ndim; k++) sdata[k+1] = (float) (nbody->stardata[i].r[k]*simunits.r.outscale);
      for (k=0; k<ndim; k++) {
        sdata[k+1+ndim] = (float) (nbody->stardata[i].v[k]*simunits.v.outscale);
      }
      sdata[1+2*ndim] = (float) (nbody->stardata[i].m*simunits.m.outscale);
      sdata[2+2*ndim] = (float) (nbody->stardata[i].h*simunits.r.outscale);
      sdata[3+2*ndim] = (float)(nbody->stardata[i].radius*simunits.r.outscale);
      for (ii=0; ii<sink_data_length; ii++) {
        writer.write_value(sdata[ii]);
      }
    }
  }
  //---------------------------------------------------------------------------

  outfile.close();

  return true;
}



//=================================================================================================
//  Simulation::ConvertToCodeUnits
/// For any simulations loaded into memory via a snapshot file, all particle
/// variables are converted into dimensionless code units here.
//=================================================================================================
template <int ndim>
void Simulation<ndim>::ConvertToCodeUnits(void)
{
  int i;                            ///< Particle counter
  int k;                            ///< Dimension counter

  debug2("[Simulation::ConvertToCodeUnits]");

  // Rescale all hydro particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<hydro->Nhydro; i++) {
    Particle<ndim>& part = hydro->GetParticlePointer(i);
    for (k=0; k<ndim; k++) part.r[k] /= simunits.r.inscale;
    for (k=0; k<ndim; k++) part.v[k] /= simunits.v.inscale;
    part.m    /= simunits.m.inscale;
    part.h    /= simunits.r.inscale;
    part.u    /= simunits.u.inscale;
    part.rho  /= simunits.rho.inscale;
    part.dudt /= simunits.dudt.inscale;
  }


  // Rescale all N-body particles
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<nbody->Nstar; i++) {
    for (k=0; k<ndim; k++) nbody->stardata[i].r[k] /= simunits.r.inscale;
    for (k=0; k<ndim; k++) nbody->stardata[i].v[k] /= simunits.v.inscale;
    nbody->stardata[i].m      /= simunits.m.inscale;
    nbody->stardata[i].h      /= simunits.r.inscale;
    nbody->stardata[i].radius /= simunits.r.inscale;
    nbody->stardata[i].invh   = 1.0/nbody->stardata[i].h;
  }


  // Rescale all sink information
  //-----------------------------------------------------------------------------------------------
  for (i=0; i<sinks->Nsink; i++) {
    sinks->sink[i].radius /= simunits.r.inscale;
  }


  // Rescale other variables
  t             /= simunits.t.inscale;
  tsnaplast     /= simunits.t.inscale;
  tlitesnaplast /= simunits.t.inscale;
  hydro->mmean  /= simunits.m.inscale;
  if (restart) {
    tsnapnext = tsnaplast + dt_snap;
    tlitesnapnext = tlitesnaplast + dt_litesnap;
  }

  return;
}
