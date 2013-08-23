//=============================================================================
//  SimulationIO.hpp
//  ..
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics and Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G Rosotti
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
using namespace std;



//=============================================================================
//  SimulationBase::ReadSnapshotFile
//=============================================================================
bool SimulationBase::ReadSnapshotFile
(string filename,                   ///< [in] Name of input snapshot file
 string fileform)                   ///< [in] Format of input snapshot file
{
  debug2("[Simulation::ReadSnapshotFile]");

  cout << "Read file : " << filename << "   format : " << fileform << endl;

  if (fileform == "column")
    return ReadColumnSnapshotFile(filename);
  else if (fileform == "sf" || fileform == "seren_form")
    return ReadSerenFormSnapshotFile(filename);
  else {
    cout << "Unrecognised file format" << endl;
    return false;
  }
}



//=============================================================================
//  SimulationBase::WriteSnapshotFile
//=============================================================================
bool SimulationBase::WriteSnapshotFile
(string filename,                   ///< [in] Name of output snapshot file
 string fileform)                   ///< [in] Format of output snapshot file
{
  debug2("[Simulation::WriteSnapshotFile]");

  if (fileform == "column")
    return WriteColumnSnapshotFile(filename);
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
(string filename,                   ///< ..
 string format)                     ///< ..
{
  HeaderInfo info;
  ifstream infile;
  infile.open(filename.c_str());

  if (format == "column")
    ReadColumnHeaderFile(infile, info);
  else if (format == "sf" || format == "seren_form")
    ReadColumnHeaderFile(infile, info);
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

  // Open file and read header information
  infile >> info.Nsph;
  infile >> info.Nstar;
  infile >> info.ndim;
  infile >> info.t;

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
//=============================================================================
template <int ndim>
bool Simulation<ndim>::ReadColumnSnapshotFile(string filename)
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
  // --------------------------------------------------------------------------
  while (infile.good() && i < sph->Nsph) {
    if (ndim == 1) infile >> sph->sphdata[i].r[0] >> sph->sphdata[i].v[0] 
			  >> sph->sphdata[i].m >> sph->sphdata[i].h
			  >> sph->sphdata[i].rho >> sph->sphdata[i].u;
    else if (ndim == 2) infile >> sph->sphdata[i].r[0] >> sph->sphdata[i].r[1] 
			       >> sph->sphdata[i].v[0] >> sph->sphdata[i].v[1] 
			       >> sph->sphdata[i].m >> sph->sphdata[i].h
			       >> sph->sphdata[i].rho >> sph->sphdata[i].u;
    else if (ndim == 3) infile >> sph->sphdata[i].r[0] >> sph->sphdata[i].r[1] 
			       >> sph->sphdata[i].r[2] >> sph->sphdata[i].v[0] 
			       >> sph->sphdata[i].v[1] >> sph->sphdata[i].v[2] 
			       >> sph->sphdata[i].m >> sph->sphdata[i].h
			       >> sph->sphdata[i].rho >> sph->sphdata[i].u;
    i++;
  }

  nbody->Nstar = info.Nstar;
  nbody->AllocateMemory(nbody->Nstar);
  i = 0;

  // Read in data depending on dimensionality
  // --------------------------------------------------------------------------
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
//=============================================================================
template <int ndim>
bool Simulation<ndim>::WriteColumnSnapshotFile(string filename)
{
  int i;
  int ndimaux;
  int Npart;
  int Nstar;
  ofstream outfile;

  debug2("[Simulation::WriteColumnSnapshotFile]");

  cout << "Writing current data to snapshot file : " << filename << endl;

  // Open file and read header information
  outfile.open(filename.c_str());
  outfile << sph->Nsph << endl;
  outfile << nbody->Nstar << endl;
  outfile << ndim << endl;
  outfile << t << endl;

  // Write data for SPH particles
  // --------------------------------------------------------------------------
  for (i=0; i<sph->Nsph; i++) {
    if (ndim == 1) outfile << sph->sphdata[i].r[0] << "   " 
			   << sph->sphdata[i].v[0] << "   "
			   << sph->sphdata[i].m << "   " 
			   << sph->sphdata[i].h << "   "
			   << sph->sphdata[i].rho << "   " 
			   << sph->sphdata[i].u
			   << endl;
    else if (ndim == 2) outfile << sph->sphdata[i].r[0] << "   " 
				<< sph->sphdata[i].r[1] << "   " 
				<< sph->sphdata[i].v[0] << "   " 
				<< sph->sphdata[i].v[1] << "   "
				<< sph->sphdata[i].m << "   "
				<< sph->sphdata[i].h << "   "
				<< sph->sphdata[i].rho << "   "
				<< sph->sphdata[i].u
				<< endl;
    else if (ndim == 3) outfile << sph->sphdata[i].r[0] << "   "
				<< sph->sphdata[i].r[1] << "   " 
				<< sph->sphdata[i].r[2] << "   "
				<< sph->sphdata[i].v[0] << "   " 
				<< sph->sphdata[i].v[1] << "   " 
				<< sph->sphdata[i].v[2] << "   "
				<< sph->sphdata[i].m << "   "
				<< sph->sphdata[i].h << "   "
				<< sph->sphdata[i].rho << "   "
				<< sph->sphdata[i].u
				<< endl;
  }

  // Write data for SPH particles
  // --------------------------------------------------------------------------
  for (i=0; i<nbody->Nstar; i++) {
    if (ndim == 1) outfile << nbody->stardata[i].r[0] << "   " 
			   << nbody->stardata[i].v[0] << "   "
			   << nbody->stardata[i].m << "   " 
			   << nbody->stardata[i].h << "   "
			   << 0.0 << "   " 
			   << 0.0
			   << endl;
    else if (ndim == 2) outfile << nbody->stardata[i].r[0] << "   " 
				<< nbody->stardata[i].r[1] << "   " 
				<< nbody->stardata[i].v[0] << "   " 
				<< nbody->stardata[i].v[1] << "   "
				<< nbody->stardata[i].m << "   "
				<< nbody->stardata[i].h << "   "
				<< 0.0 << "   "
				<< 0.0
				<< endl;
    else if (ndim == 3) outfile << nbody->stardata[i].r[0] << "   "
				<< nbody->stardata[i].r[1] << "   " 
				<< nbody->stardata[i].r[2] << "   "
				<< nbody->stardata[i].v[0] << "   " 
				<< nbody->stardata[i].v[1] << "   " 
				<< nbody->stardata[i].v[2] << "   "
				<< nbody->stardata[i].m << "   "
				<< nbody->stardata[i].h << "   "
				<< 0.0 << "   "
				<< 0.0
				<< endl;
  }

  outfile.close();

  return true;
}



//=============================================================================
//  Simulation::ReadSerenFormHeaderFile
/// Function for reading the header file of a snapshot. Does not modify the
/// variables of the Simulation class, but rather returns information in a 
/// HeaderInfo struct.
//=============================================================================
template <int ndim>
void Simulation<ndim>::ReadSerenFormHeaderFile
(ifstream& infile,                  ///< ..
 HeaderInfo& info)                  ///< ..
{
  // Open file and read header information
  infile >> info.Nsph;
  infile >> info.Nstar;
  infile >> info.ndim;
  infile >> info.t;

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
//=============================================================================
template <int ndim>
bool Simulation<ndim>::ReadSerenFormSnapshotFile(string filename)
{
  int dim_check;               // Dimension check
  int dmdt_range_aux;          // Accretion history array size
  int i;                       // Aux. counter
  int ifirst;                  // ..
  int ilast;                   // ..
  int j;
  int k;
  //int j;                       // Aux. counter
  //int itemp;                   // Dummy integer variable
  int ndata;                   // No. of data arrays written
  int nunit;                   // No. of unit strings
  int p;                       // Particle counter
  int pr_check;                // Precision check
  int idata[50];               // Integer data
  int ilpdata[50];             // ..
  int typedata[50][5];         // ..
  FLOAT rdata[50];             // Real data array
  //FLOAT taux;                  // Aux. time variable
  FLOAT rtemp;                 // Dummy float variable
  DOUBLE ddata[50];            // ..
  string format_id;            // File format (for verification)
  string data_id[50];          // String ids of arrays written
  string unit_data[50];        // String ids of units written
  ifstream infile;             // ..
  FLOAT raux;                  // ..
  HeaderInfo info;             // ..
  int sink_data_length;        // ..
  string dummystring;          // ..
  bool booldummy;              // ..

  debug2("[Simulation::ReadColumnSnapshotFile]");

  cout << "Opening file : " << filename << endl;

  infile.open(filename.c_str());


  // Read information identifying format and precision of file.
  // Then check if each value corresponds to the current values.
  // --------------------------------------------------------------------------
  infile >> format_id;
  simparams->TrimWhiteSpace(format_id);

  cout << "Checking format : " << format_id << endl;
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

  // Output important info to screen
  cout << "SPH Particles  : " << sph->Nsph << endl;
  cout << "Star particles : " << nbody->Nstar << endl;

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
  //simunits.sigma.inunit = unit_data[6];
  simunits.press.inunit = unit_data[7];
  simunits.f.inunit = unit_data[8];
  simunits.E.inunit = unit_data[9];
  simunits.mom.inunit = unit_data[10];
  simunits.angmom.inunit = unit_data[11];
  simunits.angvel.inunit = unit_data[12];
  simunits.dmdt.inunit = unit_data[13];
  //simunits.L.inunit = unit_data[14];
  //simunits.kappa.inunit = unit_data[15];
  //simunits.B.inunit = unit_data[16];
  //simunits.Q.inunit = unit_data[17];
  //simunits.J.inunit = unit_data[18];

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
  // --------------------------------------------------------------------------
  for (j=0; j<ndata; j++) {
    cout << "Reading data for array : " << data_id[j] << endl;

    ifirst = typedata[i][1];
    ilast = typedata[i][2];

    // porig
    // ------------------------------------------------------------------------
    if (data_id[j] == "porig") {
      for (i=0; i<sph->Nsph; i++) infile >> sph->sphdata[i].iorig;
    }    

    // Positions
    // ------------------------------------------------------------------------
    else if (data_id[j] == "r") {
      if (ndim == 1)
	for (i=0; i<sph->Nsph; i++) infile >> sph->sphdata[i].r[0];
      else if (ndim == 2)
	for (i=0; i<sph->Nsph; i++) 
	  infile >> sph->sphdata[i].r[0] >> sph->sphdata[i].r[1];
      else if (ndim == 3)
	for (i=0; i<sph->Nsph; i++) 
	  infile >> sph->sphdata[i].r[0] >> sph->sphdata[i].r[1] 
		 >> sph->sphdata[i].r[2];
    }

    // Masses
    // ------------------------------------------------------------------------
    else if (data_id[j] == "m") {
      for (i=0; i<sph->Nsph; i++) infile >> sph->sphdata[i].m;
    }

    // Smoothing lengths
    // ------------------------------------------------------------------------
    else if (data_id[j] == "h") {
      for (i=0; i<sph->Nsph; i++) infile >> sph->sphdata[i].h;
    }

    // Velocities
    // ------------------------------------------------------------------------
    else if (data_id[j] == "v") {
      if (ndim == 1)
	for (i=0; i<sph->Nsph; i++) infile >> sph->sphdata[i].v[0];
      else if (ndim == 2)
	for (i=0; i<sph->Nsph; i++) 
	  infile >> sph->sphdata[i].v[0] >> sph->sphdata[i].v[1];
      else if (ndim == 3)
	for (i=0; i<sph->Nsph; i++) 
	  infile >> sph->sphdata[i].v[0] >> sph->sphdata[i].v[1] 
		 >> sph->sphdata[i].v[2];
    }

    // Other 1-D redundant information
    // ------------------------------------------------------------------------
    else if(data_id[j] == "temp") {
      for (i=0; i<sph->Nsph; i++) infile >> rtemp;
    }

    // Densities
    // ------------------------------------------------------------------------
    else if (data_id[j] == "rho") {
      for (i=0; i<sph->Nsph; i++) infile >> sph->sphdata[i].rho;
    }

    // Specific internal energies
    // ------------------------------------------------------------------------
    else if (data_id[j] == "u") {
      for (i=0; i<sph->Nsph; i++) infile >> sph->sphdata[i].u;
    }

    // Sinks/stars
    // ------------------------------------------------------------------------
    else if (data_id[j] == "sink_v1") {
      sink_data_length = 12 + 2*ndim + 2*dmdt_range_aux;
      int ii;
      FLOAT sdata[sink_data_length];
      for (ii=0; ii<6; ii++) infile >> idata[ii];
      if (nbody->Nstar > 0) {
        for (i=0; i<nbody->Nstar; i++) {
          for (ii=0; ii<2; ii++) infile >> booldummy;
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
    // ------------------------------------------------------------------------
    else if (typedata[i][0] >= 1) {
      int kk = typedata[i][0];
      for (i=ifirst-1; i<ilast; i++)
        for (k=0; k<kk; k++) infile >> dummystring;
    }

  }
  // --------------------------------------------------------------------------

  // Close file
  infile.close();
  
  return true;
}



//=============================================================================
//  Simulation::WriteSerenFormSnapshotFile
//=============================================================================
template <int ndim>
bool Simulation<ndim>::WriteSerenFormSnapshotFile(string filename)
{
  return true;
}
