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
return;
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
return;
  // Check level is valid
  if (level <= 0) {
    cout << "Error with timing levels" << endl;
    return;
  }
  level--;

  // If block not in list, then print error message and return
  if (blockmap[level].find(s1) == blockmap[level].end()) {
    cout << "Error : looking for incorrect timing block : " << s1 << endl;
    return;
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
 (string run_id)                       ///< String i.d. of current simulation
{
  int iblock;                          // Timing block counter
  int l;                               // Timing block level
  DOUBLE tcount = 0.0;                 // Total time
  DOUBLE tcount_wall = 0.0;            // Total wall-clock time
  string filename;                     // Filename for writing statistics
  string fileend = "timing";           // Filename ending
  ofstream outfile;                    // Output file stream object

  filename = run_id + "." + fileend;
  outfile.open(filename.c_str());

  tend      = clock();
  tend_wall = WallClockTime();
  ttot      = (double) (tend - tstart) / (double) CLOCKS_PER_SEC;
  ttot_wall = (double) (tend_wall - tstart_wall);

  outfile << resetiosflags(ios::adjustfield);
  outfile << setiosflags(ios::left);
  outfile << "----------------------------------------------------------------------------" << endl;
  outfile << "Total simulation wall clock time : " << ttot_wall << endl;
  outfile << "----------------------------------------------------------------------------" << endl;

  // Output timing data on each hierarchical level
  //-----------------------------------------------------------------------------------------------
  for (l=0; l<Nlevel; l++) {
    tcount = 0.0;
    tcount_wall = 0.0;
    outfile << "Level : " << l << endl;
    outfile << "Block                                          Wall time       %time" << endl;
    outfile << "----------------------------------------------------------------------------" << endl;

    for (iblock=0; iblock<Nblock[l]; iblock++) {
      tcount      += block[l][iblock].ttot;
      tcount_wall += block[l][iblock].ttot_wall;
      block[l][iblock].tfraction      = block[l][iblock].ttot / ttot;
      block[l][iblock].tfraction_wall = block[l][iblock].ttot_wall / ttot_wall;
      outfile << setw(47) << block[l][iblock].block_name
              << setw(16) << block[l][iblock].ttot_wall
              << setw(15) << 100.0*block[l][iblock].tfraction_wall << endl;
      /*outfile << setw(30) << block[l][iblock].block_name
              << setw(12) << block[l][iblock].ttot
              << setw(12) << 100.0*block[iblock][l].tfraction
              << setw(12) << block[l][iblock].ttot_wall
              << setw(12) << 100.0*block[l][iblock].tfraction_wall << endl;*/
    }
    outfile << setw(47) << "REMAINDER"
            << setw(16) << ttot_wall - tcount_wall
            << setw(15) << 100.0*(ttot_wall - tcount_wall)/ttot_wall << endl;

    /*outfile << setw(30) << "REMAINDER"
            << setw(12) << ttot - tcount
            << setw(12) << 100.0*(ttot - tcount)/ttot
            << setw(12) << ttot_wall - tcount_wall
            << setw(12) << 100.0*(ttot_wall - tcount_wall)/ttot_wall << endl;*/
    outfile << "----------------------------------------------------------------------------" << endl;
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
