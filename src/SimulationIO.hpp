//=============================================================================
//  SimulationIO.hpp
//=============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstring>
#include "Simulation.h"
#include "Parameters.h"
#include "Debug.h"
using namespace std;



//=============================================================================
//  SimulationBase::ReadSnapshotFile
//=============================================================================
bool SimulationBase::ReadSnapshotFile
(string filename,                   ///< [in] Name of input snapshot file
 string fileform)                   ///< [in] Format of input snapshot file
{
  debug1("[Simulation::ReadSnapshotFile]");

  if (fileform == "ascii")
    return ReadColumnSnapshotFile(filename);
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
  debug1("[Simulation::WriteSnapshotFile]");

  if (fileform == "column")
    return WriteColumnSnapshotFile(filename);
  else {
    cout << "Unrecognised file format" << endl;
    return false;
  }
}



//=============================================================================
//  Simulation::ReadColumnSnapshotFile
//=============================================================================
template <int ndim>
bool Simulation<ndim>::ReadColumnSnapshotFile(string filename)
{
  int i;
  int ndimaux;
  int Npart;
  int Nstar;
  ifstream infile;
  FLOAT raux;

  debug1("[Simulation::ReadColumnSnapshotFile]");

  // Open file and read header information
  infile.open(filename.c_str());
  infile >> Npart;
  infile >> Nstar;
  infile >> ndimaux;
  infile >> t;

  // Check dimensionality matches if using fixed dimensions
  if (ndimaux != ndim) {
    cout << "Incorrect no. of dimensions in file : " 
	 << ndimaux << "  [ndim : " << ndim << "]" << endl;
    return false;
  }

  sph->Nsph = Npart;
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

  nbody->Nstar = Nstar;
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

  debug1("[Simulation::WriteColumnSnapshotFile]");

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

