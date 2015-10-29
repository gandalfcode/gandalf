//=================================================================================================
//  CodeTiming.cpp
//  Contains class functions for recording CPU and wall-clock times for
//  marked sections of the code.
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


#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include "Exception.h"
#include "CodeTiming.h"
#include "Debug.h"
using namespace std;



//=================================================================================================
//  CodeTiming::CodeTiming
/// Constructor for CodeTiming class
//=================================================================================================
CodeTiming::CodeTiming()
{
  for (int i=0; i<Nlevelmax; i++) Nblock[i] = 0;
  level       = 0;
  Nlevel      = 0;
  ttot        = 0.0;
  ttot_wall   = 0.0;
  tstart      = clock();
  tstart_wall = WallClockTime();
#if defined MPI_PARALLEL
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  CPU_Master = rank==0 ? true : false;
#endif
}



//=================================================================================================
//  CodeTiming::~CodeTiming
/// CodeTiming destructor
//=================================================================================================
CodeTiming::~CodeTiming()
{
}



//=================================================================================================
//  CodeTiming::StartTimingSection
/// Start timing a block of code marked by the string 'newblock'.  If block is already on record
/// (e.g. from previous timestep), then timing is appended to previously recorded times.
//=================================================================================================
void CodeTiming::StartTimingSection
 (const string newblock)               ///< [in] String of new/existing timing block
{
  int iblock;                          // Integer id of existing timing block
#ifdef MPI_PARALLEL
  return;
#endif
  // If block string not in list, then create new entry to timing block
  if (blockmap[level].find(newblock) == blockmap[level].end()) {
    iblock                            = Nblock[level];
    blockmap[level][newblock]         = iblock;
    block[level][iblock].timing_level = level;
    block[level][iblock].block_name   = newblock;
    Nblock[level]++;
  }
  // Else, look-up existing timing block
  else {
    iblock = blockmap[level][newblock];
  }

  // Now record timein arrays
  block[level][iblock].tstart      = clock();
  block[level][iblock].tstart_wall = WallClockTime();

  // Add block to stack for book-keeping and error checking
  levelstack[level] = &(block[level][iblock]);
  level++;
  Nlevel = max(Nlevel, level);

  return;
}



//=================================================================================================
//  CodeTiming::EndTimingSection
/// Terminate timing a block of code signified by the string 's1' and record
/// time in main timing arrays.
//=================================================================================================
void CodeTiming::EndTimingSection
 (const string s1)                     ///< [in] String identifying block-end
{
  int iblock;                          // Integer i.d. of timing block in arrays
#ifdef MPI_PARALLEL
  return;
#endif
  // Check level is valid
  if (level <= 0) {
    ExceptionHandler::getIstance().raise("Error with timing levels");
  }
  level--;

  // If block not in list, then print error message and return
  if (blockmap[level].find(s1) == blockmap[level].end()) {
    string message = "Error : looking for incorrect timing block : " + s1;
    for (map<string, int>::iterator it = blockmap[level].begin();
         it != blockmap[level].end(); it++) {
      cout << it->first << " ";
    }
    cout << endl;
    ExceptionHandler::getIstance().raise(message);
  }
  // Else, look-up existing timing block
  else {
    iblock = blockmap[level][s1];
  }

  block[level][iblock].tend = clock();
  block[level][iblock].ttot +=
    (double) (block[level][iblock].tend - block[level][iblock].tstart) / (double) CLOCKS_PER_SEC;

  block[level][iblock].tend_wall = WallClockTime();
  block[level][iblock].ttot_wall +=
    (double) (block[level][iblock].tend_wall - block[level][iblock].tstart_wall);

  return;
}



//=================================================================================================
//  CodeTiming::ComputeTimingStatistics
/// Compute all timing statistics of the code and record statistics into an
/// external file 'run_id.timing'.
//=================================================================================================
void CodeTiming::ComputeTimingStatistics
 (string run_id)                       ///< [in] String i.d. of current simulation
{
  int iblock;                          // Timing block counter
  int l;                               // Timing block level
  DOUBLE tcount = 0.0;                 // Total time
  DOUBLE tcount_wall = 0.0;            // Total wall-clock time
  string filename;                     // Filename for writing statistics
  string fileend = "timing";           // Filename ending
  ofstream outfile;                    // Output file stream object

  filename  = run_id + "." + fileend;
  tend      = clock();
  tend_wall = WallClockTime();
  ttot      = (double) (tend - tstart) / (double) CLOCKS_PER_SEC;

#if defined MPI_PARALLEL
  // In this case, clock measures the cpu time on the local processor. We need to sum
  // the cpu time of ALL processors to get the cpu time of the simulation
  if (CPU_Master) {
    MPI_Reduce(MPI_IN_PLACE, &ttot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Reduce(&ttot, &ttot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  // Now we do the same for all the blocks. Saves the ttot of each single block
  // inside the temporary array ttot_temp (so that we call MPI_Reduce only once)
  int Nblocks=0;
  for (l=0; l<Nlevel; l++) {
    Nblocks += Nblock[l];
  }

  vector<double> ttot_temp(Nblocks);
  int counter=0;
  for (l=0; l< Nlevel; l++) {
    for (iblock=0; iblock<Nblock[l]; iblock++) {
      ttot_temp[counter]=block[l][iblock].ttot;
      counter++;
    }
  }

  if (CPU_Master) {
    MPI_Reduce(MPI_IN_PLACE, &ttot_temp[0], Nblocks, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Reduce(&ttot_temp[0], &ttot_temp[0], Nblocks, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  if (!CPU_Master) return;

  counter=0;
  for (l=0; l< Nlevel; l++) {
    for (iblock=0; iblock<Nblock[l]; iblock++) {
      block[l][iblock].ttot=ttot_temp[counter];
      counter++;
    }
  }
#endif
  ttot_wall = (double) (tend_wall - tstart_wall);

  outfile.open(filename.c_str());
  outfile << resetiosflags(ios::adjustfield);
  outfile << setiosflags(ios::left);
  std::string dashes(95,'-');
  outfile << dashes << endl;
  outfile << "Total simulation wall clock time : " << ttot_wall << endl;
  outfile << dashes << endl;

  // Output timing data on each hierarchical level
  //-----------------------------------------------------------------------------------------------
  for (l=0; l<Nlevel; l++) {
    tcount = 0.0;
    tcount_wall = 0.0;
    outfile << "Level : " << l << endl;
    outfile << "Block";
    outfile << std::string( 42, ' ' );
    outfile << "Wall time" << std::string(3,' ') << "%time";
    outfile << std::string(7,' ') << "CPU Time" << std::string(4,' ') << "%time"<<endl;
    outfile << dashes << endl;

    for (iblock=0; iblock<Nblock[l]; iblock++) {
      tcount      += block[l][iblock].ttot;
      tcount_wall += block[l][iblock].ttot_wall;
      block[l][iblock].tfraction      = block[l][iblock].ttot / ttot;
      block[l][iblock].tfraction_wall = block[l][iblock].ttot_wall / ttot_wall;
      outfile << setw(47) << block[l][iblock].block_name
              << setw(12) << block[l][iblock].ttot_wall
              << setw(12) << 100.0*block[l][iblock].tfraction_wall;
      outfile << setw(12) << block[l][iblock].ttot
              << setw(12) << 100.0*block[l][iblock].tfraction
              << endl;
    }
    outfile << setw(47) << "REMAINDER"
            << setw(12) << ttot_wall - tcount_wall
            << setw(12) << 100.0*(ttot_wall - tcount_wall)/ttot_wall;

    outfile << setw(12) << ttot - tcount
            << setw(12) << 100.0*(ttot - tcount)/ttot
            << endl;
    outfile << dashes << endl;
  }
  //-----------------------------------------------------------------------------------------------

  outfile << resetiosflags(ios::adjustfield);
  outfile.close();

  return;
}



//=================================================================================================
//  CodeTiming::WallClockTime
/// Returns current world clock time (from set point) in seconds
//=================================================================================================
double CodeTiming::WallClockTime(void)
{
  struct timeval tm;
  gettimeofday( &tm, NULL );
  return (double)tm.tv_sec + (double)tm.tv_usec / 1000000.0;
}
