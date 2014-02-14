//=============================================================================
//  SimulationIO.hpp
//  Contains all functions for reading and writing simulation data
//  to/from snapshot files.
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
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif
using namespace std;



//=============================================================================
//  SimulationBase::ReadSnapshotFile
/// Read snapshot file with specified filename and format.
//=============================================================================
bool SimulationBase::ReadSnapshotFile
(string filename,                   ///< [in] Name of input snapshot file
 string fileform)                   ///< [in] Format of input snapshot file
{
  debug2("[Simulation::ReadSnapshotFile]");

  cout << "Reading snapshot : " << filename << "   format : " << fileform << endl;

  // Read in snapshot file to main memory
  if (fileform == "column")
    return ReadColumnSnapshotFile(filename);
  else if (fileform == "sf" || fileform == "seren_form")
    return ReadSerenFormSnapshotFile(filename);
  else {
    cout << "Unrecognised file format" << endl;
    return false;
  }

  // Set flag to ensure variables are converted to code units
  rescale_particle_data = true;

}



//=============================================================================
//  SimulationBase::WriteSnapshotFile
/// Write snapshot file with specified filename and format.
//=============================================================================
bool SimulationBase::WriteSnapshotFile
(string filename,                   ///< [in] Name of output snapshot file
 string fileform)                   ///< [in] Format of output snapshot file
{
  debug2("[Simulation::WriteSnapshotFile]");

  if (fileform == "column")
    return WriteColumnSnapshotFile(filename);
  else if (fileform == "sf" || fileform == "seren_form")
    return WriteSerenFormSnapshotFile(filename);
  else {
    cout << "Unrecognised file format" << endl;
    return false;
  }
}



//=============================================================================
//  SimulationBase::ReadHeaderSnapshotFile
/// Read the header of a snapshot file, given the filename and the format.
/// Return information in a HeaderInfo struct. This is designed to be called
/// from Python.
//=============================================================================
HeaderInfo SimulationBase::ReadHeaderSnapshotFile
(string filename,                   ///< [in] Name of output snapshot file
 string fileform)                   ///< [in] Format of output snapshot file
{
  HeaderInfo info;                  // ..
  ifstream infile;                  // ..

  debug2("[Simulation::ReadHeaderSnapshotFile]");

  infile.open(filename.c_str());

  if (fileform == "column")
    ReadColumnHeaderFile(infile, info);
  else if (fileform == "sf" || fileform == "seren_form")
    ReadSerenFormHeaderFile(infile, info);
  else
    ExceptionHandler::getIstance().raise("Unrecognised file format");

  return info;
}



//=============================================================================
//  Simulation::ReadColumnHeaderFile
/// Function for reading the header file of a snapshot. Does not modify the
/// variables of the Simulation class, but rather returns information in a 
/// HeaderInfo struct.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ReadColumnHeaderFile
(ifstream& infile,                  ///< ..
 HeaderInfo& info)                  ///< ..
{
  debug2("[Simulation::ReadColumnHeaderFile]");

  // Open file and read header information
  infile >> info.Nsph;
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



//=============================================================================
//  Simulation::ReadColumnSnapshotFile
/// Reads a column format data snapshot of given filename.
//=============================================================================
template <int ndim>
bool Simulation<ndim>::ReadColumnSnapshotFile
(string filename)                  ///< Filename of column data snapshot file
{
  int i;
  ifstream infile;
  FLOAT raux;
  HeaderInfo info;

  debug2("[Simulation::ReadColumnSnapshotFile]");

  infile.open(filename.c_str());

  ReadColumnHeaderFile(infile, info);
  t = info.t;

  sph->Nsph = info.Nsph;
  sph->AllocateMemory(sph->Nsph);
  i = 0;

  // Read in data depending on dimensionality
  //---------------------------------------------------------------------------
  while (infile.good() && i < sph->Nsph) {
    SphParticle<ndim>* part = sph->GetParticleIPointer(i);
    if (ndim == 1) 
      infile >> part->r[0] >> part->v[0] >> part->m >> part->h
	     >> part->rho >> part->u;
    else if (ndim == 2) 
      infile >> part->r[0] >> part->r[1] >> part->v[0] >> part->v[1]
             >> part->m >> part->h >> part->rho >> part->u;
    else if (ndim == 3) 
      infile >> part->r[0] >> part->r[1] >> part->r[2] >> part->v[0]
             >> part->v[1] >> part->v[2] >> part->m >> part->h
	     >> part->rho >> part->u;
    i++;
  }

  nbody->Nstar = info.Nstar;
  nbody->AllocateMemory(nbody->Nstar);
  i = 0;

  // Read in data depending on dimensionality
  //---------------------------------------------------------------------------
  while (infile.good() && i < nbody->Nstar) {
    if (ndim == 1) 
      infile >> nbody->stardata[i].r[0] >> nbody->stardata[i].v[0] 
	     >> nbody->stardata[i].m >> nbody->stardata[i].h
	     >> raux >> raux;
    else if (ndim == 2) 
      infile >> nbody->stardata[i].r[0] >> nbody->stardata[i].r[1] 
	     >> nbody->stardata[i].v[0] >> nbody->stardata[i].v[1] 
	     >> nbody->stardata[i].m >> nbody->stardata[i].h
	     >> raux >> raux;
    else if (ndim == 3) 
      infile >> nbody->stardata[i].r[0] >> nbody->stardata[i].r[1] 
	     >> nbody->stardata[i].r[2] >> nbody->stardata[i].v[0] 
	     >> nbody->stardata[i].v[1] >> nbody->stardata[i].v[2] 
	     >> nbody->stardata[i].m >> nbody->stardata[i].h
	     >> raux >> raux;
    i++;
  }


  infile.close();

  return true;
}



//=============================================================================
//  Simulation::WriteColumnSnapshotFile
/// Write SPH and N-body particle data to column data snapshot file.
//=============================================================================
#ifdef MPI_PARALLEL
template <int ndim>
bool Simulation<ndim>::WriteColumnSnapshotFile(string filename)
{
  int i;
  ostringstream outfile;

  debug2("[Simulation::WriteColumnSnapshotFileMPI]");

  if (rank == 0)
    cout << "Writing current data to snapshot file : " << filename << endl;

  // Open file
  MPI_File file;
  char* filename_str = new char[strlen(filename.c_str())+1];
  strcpy(filename_str,filename.c_str());
  MPI_File_open(MPI_COMM_WORLD, filename_str, MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL, &file);
  delete[] filename_str;
  //Collect total number of particles
  int Ntotsph;
  MPI_Allreduce(&sph->Nsph,&Ntotsph,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  int Ntotstar;
  MPI_Allreduce(&nbody->Nstar,&Ntotstar,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  //Root node writes header
  if (rank==0) {
    outfile << Ntotsph << endl;
    outfile << Ntotstar << endl;
    outfile << ndim << endl;
    outfile << t*simunits.t.outscale << endl;
  }

  // Write data for SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    SphParticle<ndim>* part = sph->GetParticleIPointer(i);
    if (ndim == 1) 
      outfile << part->r[0]*simunits.r.outscale << "   "
              << part->v[0]*simunits.v.outscale << "   "
              << part->m*simunits.m.outscale << "   "
              << part->h*simunits.r.outscale << "   "
              << part->rho*simunits.rho.outscale << "   "
              << part->u*simunits.u.outscale
              << endl;
    else if (ndim == 2) 
      outfile << part->r[0]*simunits.r.outscale << "   "
              << part->r[1]*simunits.r.outscale << "   "
              << part->v[0]*simunits.v.outscale << "   "
              << part->v[1]*simunits.v.outscale << "   "
              << part->m*simunits.m.outscale << "   "
              << part->h*simunits.r.outscale << "   "
              << part->rho*simunits.rho.outscale << "   "
              << part->u*simunits.u.outscale
              << endl;
    else if (ndim == 3) 
      outfile << part->r[0]*simunits.r.outscale << "   "
              << part->r[1]*simunits.r.outscale << "   "
              << part->r[2]*simunits.r.outscale << "   "
              << part->v[0]*simunits.v.outscale << "   "
              << part->v[1]*simunits.v.outscale << "   "
              << part->v[2]*simunits.v.outscale << "   "
              << part->m*simunits.m.outscale << "   "
              << part->h*simunits.r.outscale << "   "
              << part->rho*simunits.rho.outscale << "   "
              << part->u*simunits.u.outscale
              << endl;
  }

  //Now all nodes write to the file their portion
  //To do that, we need to know the offset of each node, summin up the length of each bit
  {
    std::string content = outfile.str();
    int offset;
    int length_char = content.length();
    MPI_Exscan(&length_char,&offset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
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
    MPI_File_get_position(file,&end_sph_mpi);
    end_sph = end_sph_mpi;
  }
  MPI_Bcast(&end_sph,1,MPI_INT,Nmpi-1,MPI_COMM_WORLD);

  // Write data for Nbody particles
  //---------------------------------------------------------------------------
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


  //Now all nodes write to the file their portion
  //To do that, we need to know the offset of each node, summing up the length of each bit
  {
    std::string content = outfile.str();
    int offset;
    int length_char = content.length();
    MPI_Exscan(&length_char,&offset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (rank==0) {
      offset = 0;
    }
    MPI_Offset offset_mpi = offset;
    //Offset the position by the end of the sph information
    offset_mpi += end_sph_mpi;
    MPI_File_seek(file,offset_mpi,MPI_SEEK_SET);
    //Now we can do the actual writing. Extract information from the string and pass it to MPI
    char* data = new char[strlen(content.c_str())+1];
    strcpy(data, content.c_str());
    MPI_Status status;
//    MPI_File_write(file, data, length_char, MPI_CHAR, &status);
    delete[] data;
  }

  MPI_File_close(&file);

  return true;
}
#else
template <int ndim>
bool Simulation<ndim>::WriteColumnSnapshotFile(string filename)
{
  int i;
  ofstream outfile;

  debug2("[Simulation::WriteColumnSnapshotFile]");

  cout << "Writing current data to snapshot file : " << filename << endl;

  // Open file and read header information
  outfile.open(filename.c_str());
  outfile << sph->Nsph << endl;
  outfile << nbody->Nstar << endl;
  outfile << ndim << endl;
  outfile << t*simunits.t.outscale << endl;

  // Write data for SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    SphParticle<ndim>* part = sph->GetParticleIPointer(i);
    if (ndim == 1) 
      outfile << part->r[0]*simunits.r.outscale << "   "
	      << part->v[0]*simunits.v.outscale << "   "
	      << part->m*simunits.m.outscale << "   "
	      << part->h*simunits.r.outscale << "   "
	      << part->rho*simunits.rho.outscale << "   "
	      << part->u*simunits.u.outscale
	      << endl;
    else if (ndim == 2) 
      outfile << part->r[0]*simunits.r.outscale << "   "
	      << part->r[1]*simunits.r.outscale << "   "
	      << part->v[0]*simunits.v.outscale << "   "
	      << part->v[1]*simunits.v.outscale << "   "
	      << part->m*simunits.m.outscale << "   "
	      << part->h*simunits.r.outscale << "   "
	      << part->rho*simunits.rho.outscale << "   "
	      << part->u*simunits.u.outscale
	      << endl;
    else if (ndim == 3) 
      outfile << part->r[0]*simunits.r.outscale << "   "
	      << part->r[1]*simunits.r.outscale << "   "
	      << part->r[2]*simunits.r.outscale << "   "
	      << part->v[0]*simunits.v.outscale << "   "
	      << part->v[1]*simunits.v.outscale << "   "
	      << part->v[2]*simunits.v.outscale << "   "
	      << part->m*simunits.m.outscale << "   "
	      << part->h*simunits.r.outscale << "   "
	      << part->rho*simunits.rho.outscale << "   "
	      << part->u*simunits.u.outscale
	      << endl;
  }

  // Write data for SPH particles
  //---------------------------------------------------------------------------
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

  outfile.close();

  return true;
}
#endif


//=============================================================================
//  Simulation::ReadSerenFormHeaderFile
/// Function for reading the header file of a snapshot. Does not modify the
/// variables of the Simulation class, but rather returns information in a 
/// HeaderInfo struct.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ReadSerenFormHeaderFile
(ifstream& infile,                 ///< Input file stream
 HeaderInfo& info)                 ///< Header info data structure
{
  int i;                           // Aux. counter
  int idata[50];                   // Integer data array
  int ilpdata[50];                 // Long integer data array
  FLOAT rdata[50];                 // Real data array
  DOUBLE ddata[50];                // Double precision data array
  string dummy;                    // Dummy string variable

  // Read information identifying format and precision of file.
  // Then check if each value corresponds to the current values.
  //---------------------------------------------------------------------------
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

  info.Nsph  = idata[0];
  info.Nstar = idata[1];
  info.t     = ddata[0];
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



//=============================================================================
//  Simulation::ReadSerenFormSnapshotFile
/// Read data from Seren format snapshot file into main arrays.
//=============================================================================
template <int ndim>
bool Simulation<ndim>::ReadSerenFormSnapshotFile(string filename)
{
  int dim_check;               // Dimension check
  int dmdt_range_aux;          // Accretion history array size
  int i;                       // Aux. counter
  int ifirst;                  // i.d. of first particle
  int ilast;                   // i.d. of last particle
  int j;                       // Aux. counter
  int k;                       // Dimension counter
  int ndata;                   // No. of data arrays written
  int nunit;                   // No. of unit strings
  int pr_check;                // Precision check
  int idata[50];               // Integer data array
  int ilpdata[50];             // Long data array
  int typedata[50][5];         // SPH Particle data array information
  FLOAT rdata[50];             // Real data array
  FLOAT rtemp;                 // Dummy float variable
  DOUBLE ddata[50];            // Double float data array
  string format_id;            // File format (for verification)
  string data_id[50];          // String ids of arrays written
  string unit_data[50];        // String ids of units written
  ifstream infile;             // Stream of input file
  int sink_data_length;        // Length of float array for sink data
  string dummystring;          // Dummy string variable
  bool booldummy;              // Dummy boolean variable

  debug2("[Simulation::ReadSerenFormSnapshotFile]");

  infile.open(filename.c_str());


  // Read information identifying format and precision of file.
  // Then check if each value corresponds to the current values.
  //---------------------------------------------------------------------------
  infile >> format_id;
  simparams->TrimWhiteSpace(format_id);
  if (format_id != "SERENASCIIDUMPV2" && format_id != "SERENASCIIDUMPV3") {
    cout << "Incorrect format of IC file : " << format_id << endl;
    exit(0);
  }
  infile >> pr_check;
  infile >> dim_check;
  if (dim_check != ndim) {
    cout << "Incorrect NDIM of IC file" << endl;
    exit(0);
  }
  infile >> dim_check;
    if (dim_check != ndim) {
    cout << "Incorrect BDIM of IC file" << endl;
    exit(0);
  }
  infile >> dim_check;
  if (dim_check != ndim) {
    cout << "Incorrect BDIM of IC file" << endl;
    exit(0);
  }

  // Read infile header integer data
  for (i=0; i<50; i++) infile >> idata[i];
  for (i=0; i<50; i++) infile >> ilpdata[i];
  sph->Nsph      = idata[0];
  nbody->Nstar   = idata[1];
  sinks.Nsink    = idata[1];
  dmdt_range_aux = idata[29];
  nunit          = idata[19];
  ndata          = idata[20];
  //Nsnap          = ilpdata[0];
  Nsteps         = ilpdata[1];

  // Read infile head real data
  for (i=0; i<50; i++) infile >> rdata[i];
  for (i=0; i<50; i++) infile >> ddata[i];
  //sph->h_fac = rdata[0];
  t = ddata[0];
  //sphptr->tlastsnap = ddata[1];
  //sph->mgas_orig = ddata[2];

  // Read in unit data
  if (nunit > 0) 
    for (i=0; i<nunit; i++) infile >> unit_data[i];

  // Read in unit variables here
  simunits.r.inunit = unit_data[0];
  simunits.m.inunit = unit_data[1];
  simunits.t.inunit = unit_data[2];
  simunits.v.inunit = unit_data[3];
  simunits.a.inunit = unit_data[4];
  simunits.rho.inunit = unit_data[5];
  simunits.sigma.inunit = unit_data[6];
  simunits.press.inunit = unit_data[7];
  simunits.f.inunit = unit_data[8];
  simunits.E.inunit = unit_data[9];
  simunits.mom.inunit = unit_data[10];
  simunits.angmom.inunit = unit_data[11];
  simunits.angvel.inunit = unit_data[12];
  simunits.dmdt.inunit = unit_data[13];
  simunits.L.inunit = unit_data[14];
  simunits.kappa.inunit = unit_data[15];
  simunits.B.inunit = unit_data[16];
  simunits.Q.inunit = unit_data[17];
  simunits.Jcur.inunit = unit_data[18];
  simunits.u.inunit = unit_data[19];
  simunits.temp.inunit = unit_data[20];
  if (nunit > 21) simunits.dudt.inunit = unit_data[21];

  // Read infile ids of arrays contained in file
  if (ndata > 0)
    for (i=0; i<ndata; i++) infile >> data_id[i];

  // Read infile ids of arrays contained in file
  if (ndata > 0)
    for (i=0; i<ndata; i++) 
      infile >> typedata[i][0] >> typedata[i][1] 
	     >> typedata[i][2] >> typedata[i][3] >> typedata[i][4];

  // Allocate memory now we know particle numbers
  AllocateParticleMemory();


  // Loop through array ids and read each array in turn
  //===========================================================================
  for (j=0; j<ndata; j++) {

    ifirst = typedata[j][1];
    ilast = typedata[j][2];

    // porig
    //-------------------------------------------------------------------------
    if (data_id[j] == "porig") {
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        infile >> part->iorig;
      }
    }    

    // Positions
    //-------------------------------------------------------------------------
    else if (data_id[j] == "r") {
      if (ndim == 1)
	for (i=0; i<sph->Nsph; i++) {
	  SphParticle<ndim>* part = sph->GetParticleIPointer(i);
	  infile >> part->r[0];
	}
      else if (ndim == 2)
	for (i=0; i<sph->Nsph; i++) {
	  SphParticle<ndim>* part = sph->GetParticleIPointer(i);
	  infile >> part->r[0] >> part->r[1];
	}
      else if (ndim == 3)
	for (i=0; i<sph->Nsph; i++) {
	  SphParticle<ndim>* part = sph->GetParticleIPointer(i);
	  infile >> part->r[0] >> part->r[1]
	           >> part->r[2];
	}
    }

    // Masses
    //-------------------------------------------------------------------------
    else if (data_id[j] == "m") {
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        infile >> part->m;
      }
    }

    // Smoothing lengths
    //-------------------------------------------------------------------------
    else if (data_id[j] == "h") {
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        infile >> part->h;
      }
    }

    // Velocities
    //-------------------------------------------------------------------------
    else if (data_id[j] == "v") {
      if (ndim == 1)
	for (i=0; i<sph->Nsph; i++) {
	  SphParticle<ndim>* part = sph->GetParticleIPointer(i);
	  infile >> part->v[0];
	}
      else if (ndim == 2)
	for (i=0; i<sph->Nsph; i++) {
	  SphParticle<ndim>* part = sph->GetParticleIPointer(i);
	  infile >> part->v[0] >> part->v[1];
	}

      else if (ndim == 3)
	for (i=0; i<sph->Nsph; i++) {
	  SphParticle<ndim>* part = sph->GetParticleIPointer(i);
	  infile >> part->v[0] >> part->v[1] >> part->v[2];
	}

    }

    // Other 1-D redundant information
    //-------------------------------------------------------------------------
    else if(data_id[j] == "temp") {
      for (i=0; i<sph->Nsph; i++) infile >> rtemp;
    }

    // Densities
    //-------------------------------------------------------------------------
    else if (data_id[j] == "rho") {
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        infile >> part->rho;
      }
    }

    // Specific internal energies
    //-------------------------------------------------------------------------
    else if (data_id[j] == "u") {
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        infile >> part->u;
      }
    }

    // Sinks/stars
    //-------------------------------------------------------------------------
    else if (data_id[j] == "sink_v1") {
      sink_data_length = 12 + 2*ndim;  //+ 2*dmdt_range_aux;
      int ii;
      FLOAT sdata[sink_data_length];
      for (ii=0; ii<6; ii++) infile >> idata[ii];
      if (nbody->Nstar > 0) {
        for (i=0; i<nbody->Nstar; i++) {
          for (ii=0; ii<2; ii++) infile >> idata[ii];
          for (ii=0; ii<2; ii++) infile >> idata[ii];
          for (ii=0; ii<sink_data_length; ii++) infile >> sdata[ii];
          for (k=0; k<ndim; k++) nbody->stardata[i].r[k] = sdata[k+1];
          for (k=0; k<ndim; k++) nbody->stardata[i].v[k] = sdata[k+1+ndim];
          nbody->stardata[i].m = sdata[1+2*ndim];
          nbody->stardata[i].h = sdata[2+2*ndim];
          nbody->stardata[i].radius = sdata[3+2*ndim];
        }
      }
    }

    // Skip through arbitrary 1D or 2D array
    //-------------------------------------------------------------------------
    else if (typedata[i][0] >= 1) {
      int kk = typedata[i][0];
      for (i=ifirst-1; i<ilast; i++)
        for (k=0; k<kk; k++) infile >> dummystring;
    }

  }
  //===========================================================================

  // Close file
  infile.close();
  
  return true;
}



//=============================================================================
//  Simulation::WriteSerenFormSnapshotFile
/// Write SPH and N-body particle data to snapshot file in Seren format.
//=============================================================================
template <int ndim>
bool Simulation<ndim>::WriteSerenFormSnapshotFile(string filename)
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
  int typedata[50][5];         // SPH Particle data array information
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

  for (i=0; i<50; i++) idata[i] = 0;
  for (i=0; i<50; i++) ilpdata[i] = 0;
  for (i=0; i<50; i++) rdata[i] = 0.0;
  for (i=0; i<50; i++) ddata[i] = 0.0;
  for (i=0; i<50; i++) unit_data[i] = 'null';
  //for (i=0; i<50; i++) data_id[i] = '';
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


  // Set array ids and array information data if there are any SPH particles
  //---------------------------------------------------------------------------
  if (sph->Nsph > 0) {
    data_id[ndata] = "porig";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = sph->Nsph; typedata[ndata][3] = 2;
    typedata[ndata][4] = 0; ndata++;

    data_id[ndata] = "r";
    typedata[ndata][0] = ndim; typedata[ndata][1] = 1;
    typedata[ndata][2] = sph->Nsph; typedata[ndata][3] = 4;
    typedata[ndata][4] = 1; ndata++;

    data_id[ndata] = "m";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = sph->Nsph; typedata[ndata][3] = 4;
    typedata[ndata][4] = 2; ndata++;

    data_id[ndata] = "h";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = sph->Nsph; typedata[ndata][3] = 4;
    typedata[ndata][4] = 1; ndata++;

    data_id[ndata] = "v";
    typedata[ndata][0] = ndim; typedata[ndata][1] = 1;
    typedata[ndata][2] = sph->Nsph; typedata[ndata][3] = 4;
    typedata[ndata][4] = 4; ndata++;

    data_id[ndata] = "rho";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = sph->Nsph; typedata[ndata][3] = 4;
    typedata[ndata][4] = 6; ndata++;

    data_id[ndata] = "u";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = sph->Nsph; typedata[ndata][3] = 4;
    typedata[ndata][4] = 20; ndata++;
  }

  if (nbody->Nstar > 0) {
    data_id[ndata] = "sink_v1";
    typedata[ndata][0] = 1; typedata[ndata][1] = 1;
    typedata[ndata][2] = nbody->Nstar; typedata[ndata][3] = 7;
    typedata[ndata][4] = 0; ndata++;
  }

  // Set important header information
  idata[0] = sph->Nsph;
  idata[1] = nbody->Nstar;
  idata[4] = sph->Nsph;
  idata[19] = nunit;
  idata[20] = ndata;
  ilpdata[0] = 0;
  ilpdata[1] = Nsteps;
  rdata[0] = sph->h_fac;
  rdata[1] = 0.0;
  ddata[0] = t*simunits.t.outscale;


  // Write header information to file
  //---------------------------------------------------------------------------
  outfile << "SERENASCIIDUMPV2" << endl;;
  outfile << 4 << endl;
  outfile << ndim << endl;
  outfile << ndim << endl;
  outfile << ndim << endl;
  for (i=0; i<50; i++) outfile << idata[i] << endl;
  for (i=0; i<50; i++) outfile << ilpdata[i] << endl;
  for (i=0; i<50; i++) outfile << rdata[i] << endl;
  for (i=0; i<50; i++) outfile << ddata[i] << endl;
  if (nunit > 0)
    for (i=0; i<nunit; i++) outfile << unit_data[i] << endl;
  if (ndata > 0)
    for (i=0; i<ndata; i++) outfile << data_id[i] << endl;
  if (ndata > 0)
    for (i=0; i<ndata; i++)
      outfile << typedata[i][0] << "    " << typedata[i][1] << "    "
              << typedata[i][2] << "    " << typedata[i][3] << "    "
              << typedata[i][4] << endl;

  // Write arrays for SPH particles
  //---------------------------------------------------------------------------
  if (sph->Nsph > 0) {

    // porig
    //-------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>* part = sph->GetParticleIPointer(i);
      outfile << part->iorig << endl;
    }

    // Positions
    //-------------------------------------------------------------------------
    if (ndim == 1)
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        outfile << part->r[0]*simunits.r.outscale << endl;
      }

    else if (ndim == 2)
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        outfile << part->r[0]*simunits.r.outscale << "    "
                << part->r[1]*simunits.r.outscale << endl;
      }

    else if (ndim == 3)
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        outfile << part->r[0]*simunits.r.outscale << "    "
                << part->r[1]*simunits.r.outscale << "    "
                << part->r[2]*simunits.r.outscale << endl;

      }

    // Masses
    //-------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>* part = sph->GetParticleIPointer(i);
      outfile << part->m*simunits.m.outscale << endl;
    }


    // Smoothing lengths
    //-------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>* part = sph->GetParticleIPointer(i);
      outfile << part->h*simunits.r.outscale << endl;
    }


    // Velocities
    //-------------------------------------------------------------------------
    if (ndim == 1)
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        outfile << part->v[0]*simunits.v.outscale << endl;
      }

    else if (ndim == 2)
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        outfile << part->v[0]*simunits.v.outscale << "    "
                << part->v[1]*simunits.v.outscale << endl;
      }

    else if (ndim == 3)
      for (i=0; i<sph->Nsph; i++) {
        SphParticle<ndim>* part = sph->GetParticleIPointer(i);
        outfile << part->v[0]*simunits.v.outscale << "    "
                << part->v[1]*simunits.v.outscale << "    "
                << part->v[2]*simunits.v.outscale << endl;
      }


    // Densities
    //-------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>* part = sph->GetParticleIPointer(i);
      outfile << part->rho*simunits.rho.outscale << endl;;
    }


    // Specific internal energies
    //-------------------------------------------------------------------------
    for (i=0; i<sph->Nsph; i++) {
      SphParticle<ndim>* part = sph->GetParticleIPointer(i);
      outfile << part->u*simunits.u.outscale << endl;
    }


  }

  // Sinks/stars
  //---------------------------------------------------------------------------
  if (nbody->Nstar > 0) {
    sink_data_length = 12 + 2*ndim; //+ 2*dmdt_range_aux;
    int ii;
    FLOAT sdata[sink_data_length];
    for (k=0; k<sink_data_length; k++) sdata[k] = 0.0;
    outfile << 2 << "    " << 2 << "    " << 0 << "    "
	    << sink_data_length << "    " << 0 << "    "
	    << 0 << endl;
    for (i=0; i<nbody->Nstar; i++) {
      outfile << true << "   " << true << endl;
      outfile << i+1 << "    " << 0 << endl;
      for (k=0; k<ndim; k++) 
        sdata[k+1] = nbody->stardata[i].r[k]*simunits.r.outscale;
      for (k=0; k<ndim; k++)
        sdata[k+1+ndim] = nbody->stardata[i].v[k]*simunits.v.outscale;
      sdata[1+2*ndim] = nbody->stardata[i].m*simunits.m.outscale;
      sdata[2+2*ndim] = nbody->stardata[i].h*simunits.r.outscale;
      sdata[3+2*ndim] = nbody->stardata[i].radius*simunits.r.outscale;
      for (ii=0; ii<sink_data_length; ii++) outfile << sdata[ii] << "    ";
      outfile << endl;
    }
  }
  //---------------------------------------------------------------------------

  // Close file
  outfile.close();

  return true;
}



//=============================================================================
//  Simulation::ConvertToCodeUnits
/// For any simulations loaded into memory via a snapshot file, all particle 
/// variables are converted into dimensionless code units here.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ConvertToCodeUnits(void)
{
  int i;                            ///< Particle counter
  int k;                            ///< Dimension counter

  debug2("[Simulation::ConvertToCodeUnits]");

  // Rescale all SPH particles
  //---------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    for (k=0; k<ndim; k++) sph->sphdata[i].r[k] /= simunits.r.inscale;
    for (k=0; k<ndim; k++) sph->sphdata[i].v[k] /= simunits.v.inscale;
    sph->sphdata[i].m /= simunits.m.inscale;
    sph->sphdata[i].h /= simunits.r.inscale;
    sph->sphdata[i].u /= simunits.u.inscale;
    sph->sphdata[i].rho /= simunits.rho.inscale;
    sph->sphdata[i].dudt /= simunits.dudt.inscale;
  }


  // Rescale all N-body particles
  //---------------------------------------------------------------------------
  for (i=0; i<nbody->Nstar; i++) {
    for (k=0; k<ndim; k++) nbody->stardata[i].r[k] /= simunits.r.inscale;
    for (k=0; k<ndim; k++) nbody->stardata[i].v[k] /= simunits.v.inscale;
    nbody->stardata[i].m /= simunits.m.inscale;
    nbody->stardata[i].h /= simunits.r.inscale;
    nbody->stardata[i].radius /= simunits.r.inscale;
  }

  t /= simunits.t.inscale;

  return;
}
