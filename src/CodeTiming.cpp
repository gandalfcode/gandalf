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
  Nblock      = 0;
  Nlevelmax   = 0;
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
 (string newblock,                     ///< String of new/existing timing block
  int timing_level)                    ///< Timing level of block
{
  int iblock;                          // Integer id of existing timing block
return;
  // If block string not in list, then create new entry to timing block
  if (blockno.find(newblock) == blockno.end()) {
    iblock = Nblock;
    blockno[newblock] = iblock;
    block[iblock].timing_level = timing_level;
    block[iblock].block_name = newblock;
    Nblock++;
    Nlevelmax = max(Nlevelmax,timing_level);
  }
  // Else, look-up existing timing block
  else {
    iblock = blockno[newblock];
  }

  // Now record timein arrays
  block[iblock].tstart = clock();
  block[iblock].tstart_wall = WallClockTime();

  return;
}



//=================================================================================================
//  CodeTiming::EndTimingSection
/// Terminate timing a block of code signified by the string 's1' and record
/// time in main timing arrays.
//=================================================================================================
void CodeTiming::EndTimingSection
 (string s1)                           ///< String identifying block-end
{
  int iblock;                          // Integer i.d. of timing block in arrays
return;
  // If block not in list, then print error message and return
  if (blockno.find(s1) == blockno.end()) {
    cout << "Error : looking for incorrect timing block : " << s1 << endl;
    return;
  }
  // Else, look-up existing timing block
  else {
    iblock = blockno[s1];
  }

  block[iblock].tend = clock();
  block[iblock].ttot += (double) (block[iblock].tend - block[iblock].tstart) /
    (double) CLOCKS_PER_SEC;

  block[iblock].tend_wall = WallClockTime();
  block[iblock].ttot_wall += (double) (block[iblock].tend_wall - block[iblock].tstart_wall);

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
  int level;                           // Timing block level
  DOUBLE tcount = 0.0;                 // Total time
  DOUBLE tcount_wall = 0.0;            // Total wall-clock time
  string filename;                     // Filename for writing statistics
  string fileend = "timing";           // Filename ending
  ofstream outfile;                    // Output file stream object

  filename = run_id + "." + fileend;
  outfile.open(filename.c_str());

  tend = clock();
  tend_wall = WallClockTime();

  ttot = (double) (tend - tstart) / (double) CLOCKS_PER_SEC;
  ttot_wall = (double) (tend_wall - tstart_wall);

  outfile << resetiosflags(ios::adjustfield);
  outfile << setiosflags(ios::left);

  outfile << "----------------------------------------------------------------------------" << endl;
  outfile << "Total simulation time : " << ttot << "    " << ttot_wall << endl;
  outfile << "----------------------------------------------------------------------------" << endl;

  // Output timing data on each hierarchical level
  for (level=1; level<=Nlevelmax; level++) {
    tcount = 0.0;
    tcount_wall = 0.0;
    //outfile << "Level : " << level << endl;
    outfile << "Block                         Time        %time       Wall time   %time" << endl;
    outfile << "----------------------------------------------------------------------------" << endl;

    for (iblock=0; iblock<Nblock; iblock++) {
      if (block[iblock].timing_level != level) continue;
      tcount += block[iblock].ttot;
      tcount_wall += block[iblock].ttot_wall;
      block[iblock].tfraction = block[iblock].ttot / ttot;
      block[iblock].tfraction_wall = block[iblock].ttot_wall / ttot_wall;
      outfile << setw(30) << block[iblock].block_name
              << setw(12) << block[iblock].ttot
              << setw(12) << 100.0*block[iblock].tfraction
              << setw(12) << block[iblock].ttot_wall
              << setw(12) << 100.0*block[iblock].tfraction_wall << endl;
    }

    outfile << setw(30) << "REMAINDER"
            << setw(12) << ttot - tcount
            << setw(12) << 100.0*(ttot - tcount)/ttot
            << setw(12) << ttot_wall - tcount_wall
            << setw(12) << 100.0*(ttot_wall - tcount_wall)/ttot_wall << endl;
    outfile << "----------------------------------------------------------------------------" << endl;
  }

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
